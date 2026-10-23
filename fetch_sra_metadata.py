#!/usr/bin/env python3
"""
fetch_sra_metadata.py - build the pipeline samplesheet for CCLE/DepMap RNA-seq
runs held in SRA, joined to DepMap cell-line metadata and filtered by lineage.

Companion to fetch_gdc_metadata.py (same stdlib-only, no-dependency style), but
for the SRA entry point of main.nf (`--sra_samplesheet`).

Typical use (liver cell lines of the CCLE RNA-seq BioProject):

    python3 fetch_sra_metadata.py \\
        --bioproject PRJNA523380 \\
        --model-csv /data/CCLE/Model.csv \\
        --lineage Liver \\
        -o sra_samplesheet.csv

How the cell line is resolved
-----------------------------
SRA is queried for every run of the BioProject (esearch + efetch, XML). For
each run the cell-line name is taken from, in order:
  1. a SAMPLE_ATTRIBUTE named "cell line" / "cell_line" / "cell line name"
  2. the SAMPLE TITLE
  3. the SAMPLE alias
That name is then matched against DepMap's Model.csv on CCLEName,
StrippedCellLineName and ModelID (all normalised to letters+digits, uppercase,
and the "_LIVER"-style lineage suffix of CCLE names is tried both ways).
Only models whose OncotreeLineage equals --lineage are kept.

Why Model.csv rather than the SRA sample name alone: the lineage in a CCLE
name is a free-text suffix, while OncotreeLineage is the curated field the rest
of this stack already filters on (see build_ccle_reference.R).

Output columns
--------------
    sample,srr_run,model_id,ccle_name,stripped_cell_line_name,
    oncotree_lineage,oncotree_primary_disease,library_layout,library_strategy,
    total_spots,total_bases,biosample,sra_experiment,sra_sample,sra_gs_uri

  sample      = StrippedCellLineName, sanitised for use as a directory name
                (main.nf uses it as the output folder under --outdir).
  srr_run     = the SRA run accession the pipeline downloads.
  sra_gs_uri  = gs://sra-pub-run-odp/sra/<SRR>/<SRR>, SRA's public GCS mirror.
                NOTE this is the .sra archive, not FASTQ -- the pipeline still
                runs fasterq-dump on it. It is recorded because `prefetch`
                reads it natively when the job runs inside GCP, which avoids
                egress; see --gs-uris.

One run per cell line
---------------------
A cell line with several runs is normally several runs of ONE library, so
main.nf can merge them (`--merge_same_sample true`). Fused mode
(`--fused true`) needs exactly one run per sample; pass --one-run-per-line
here to keep only the largest run (by total_bases) for each cell line.

Exits 1 if any run could not be matched to a DepMap model, after writing the
good rows (same contract as fetch_gdc_metadata.py).
"""
import argparse
import csv
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
SRA_GS_PREFIX = "gs://sra-pub-run-odp/sra"
FETCH_BATCH = 200

# Model.csv column names have changed across DepMap releases; accept any.
MODEL_ID_COLS = ("ModelID", "DepMap_ID", "depmap_id")
CCLE_NAME_COLS = ("CCLEName", "CCLE_Name", "CCLE_name")
STRIPPED_COLS = ("StrippedCellLineName", "stripped_cell_line_name", "CellLineName")
LINEAGE_COLS = ("OncotreeLineage", "lineage", "Lineage")
DISEASE_COLS = ("OncotreePrimaryDisease", "primary_disease", "PrimaryDisease")

CELL_LINE_ATTRS = ("cell line", "cell_line", "cell line name", "cell-line", "cell_line_name")

OUT_COLUMNS = [
    "sample", "srr_run", "model_id", "ccle_name", "stripped_cell_line_name",
    "oncotree_lineage", "oncotree_primary_disease", "library_layout",
    "library_strategy", "total_spots", "total_bases", "biosample",
    "sra_experiment", "sra_sample", "sra_gs_uri",
]


# ---------------------------------------------------------------- utilities
def normalise(name):
    """Uppercase, letters+digits only. HEP G2 / Hep-G2 / HEPG2 all collapse."""
    return re.sub(r"[^A-Z0-9]", "", (name or "").upper())


def strip_lineage_suffix(name):
    """CCLE names are <LINE>_<TISSUE>; return the part before the first '_'."""
    return (name or "").split("_", 1)[0]


def sanitise_sample(name):
    """Make a string safe as a directory name (main.nf uses it as one)."""
    out = re.sub(r"[^A-Za-z0-9._-]+", "-", (name or "").strip())
    return out.strip("-._") or "UNKNOWN"


def first_col(row, names):
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return ""


def http_get(url, params, retries=3):
    query = urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(f"{url}?{query}", timeout=120) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            if attempt == retries - 1:
                raise SystemExit(f"NCBI request failed after {retries} tries: {exc}")
            time.sleep(2 ** attempt)


# ------------------------------------------------------------------ SRA side
def esearch(term, api_key=None):
    params = {"db": "sra", "term": term, "usehistory": "y", "retmax": 0}
    if api_key:
        params["api_key"] = api_key
    root = ET.fromstring(http_get(ESEARCH_URL, params))
    count = int(root.findtext("Count") or 0)
    return count, root.findtext("WebEnv"), root.findtext("QueryKey")


def efetch_runs(webenv, query_key, count, api_key=None, delay=0.34):
    """Yield EXPERIMENT_PACKAGE elements, paging through the history server."""
    for start in range(0, count, FETCH_BATCH):
        params = {
            "db": "sra", "WebEnv": webenv, "query_key": query_key,
            "retstart": start, "retmax": FETCH_BATCH, "retmode": "xml",
        }
        if api_key:
            params["api_key"] = api_key
        sys.stderr.write(f"  efetch {start + 1}-{min(start + FETCH_BATCH, count)} of {count}\n")

        # NCBI intermittently returns a truncated / malformed XML body (seen in
        # practice as "not well-formed (invalid token)" part-way through a
        # 2 MB batch). Without this retry a single bad batch kills the whole
        # run after several minutes of paging.
        root = None
        for attempt in range(3):
            try:
                root = ET.fromstring(http_get(EFETCH_URL, params))
                break
            except ET.ParseError as exc:
                if attempt == 2:
                    raise SystemExit(
                        f"efetch returned malformed XML for records "
                        f"{start + 1}-{start + FETCH_BATCH} after 3 tries: {exc}")
                sys.stderr.write(f"    malformed XML, retrying ({attempt + 1}/2)\n")
                time.sleep(2 ** attempt)

        for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            yield pkg
        time.sleep(delay)   # NCBI: 3 req/s without an API key, 10 with


def cell_line_from_sample(sample_el):
    """Best-effort cell line name for one SAMPLE element."""
    if sample_el is None:
        return ""
    for attr in sample_el.findall(".//SAMPLE_ATTRIBUTE"):
        tag = (attr.findtext("TAG") or "").strip().lower()
        if tag in CELL_LINE_ATTRS:
            value = (attr.findtext("VALUE") or "").strip()
            if value:
                return value
    for path in ("TITLE", ):
        value = (sample_el.findtext(path) or "").strip()
        if value:
            return value
    return (sample_el.get("alias") or "").strip()


def parse_package(pkg):
    """Flatten one EXPERIMENT_PACKAGE into one dict per RUN."""
    exp = pkg.find("EXPERIMENT")
    sample_el = pkg.find("SAMPLE")
    layout_el = pkg.find(".//LIBRARY_LAYOUT")
    layout = ""
    if layout_el is not None and len(layout_el):
        layout = layout_el[0].tag          # PAIRED or SINGLE
    strategy = pkg.findtext(".//LIBRARY_STRATEGY") or ""
    cell_line = cell_line_from_sample(sample_el)

    biosample = ""
    if sample_el is not None:
        for ext in sample_el.findall(".//EXTERNAL_ID"):
            if (ext.get("namespace") or "").lower() == "biosample":
                biosample = (ext.text or "").strip()
                break

    rows = []
    for run in pkg.findall(".//RUN"):
        acc = run.get("accession")
        if not acc:
            continue
        rows.append({
            "srr_run": acc,
            "cell_line_raw": cell_line,
            "library_layout": layout,
            "library_strategy": strategy,
            "total_spots": run.get("total_spots") or "",
            "total_bases": run.get("total_bases") or "",
            "biosample": biosample,
            "sra_experiment": (exp.get("accession") if exp is not None else "") or "",
            "sra_sample": (sample_el.get("accession") if sample_el is not None else "") or "",
        })
    return rows


# -------------------------------------------------------------- DepMap side
def load_models(path):
    """Return {normalised name -> model dict} covering several name forms."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        models = list(csv.DictReader(fh))
    if not models:
        raise SystemExit(f"No rows in {path}")

    sample_row = models[0]
    for label, cols in (("model id", MODEL_ID_COLS), ("lineage", LINEAGE_COLS)):
        if not any(c in sample_row for c in cols):
            raise SystemExit(
                f"{path} has no {label} column (looked for: {', '.join(cols)}). "
                "Is this a DepMap Model.csv?")

    index = {}
    for m in models:
        rec = {
            "model_id": first_col(m, MODEL_ID_COLS),
            "ccle_name": first_col(m, CCLE_NAME_COLS),
            "stripped": first_col(m, STRIPPED_COLS),
            "lineage": first_col(m, LINEAGE_COLS),
            "disease": first_col(m, DISEASE_COLS),
        }
        keys = {rec["model_id"], rec["ccle_name"], rec["stripped"],
                strip_lineage_suffix(rec["ccle_name"])}
        for k in keys:
            nk = normalise(k)
            # First writer wins; CCLEName/ModelID are unique, stripped names
            # can collide across models and we prefer the earlier (curated) row.
            if nk and nk not in index:
                index[nk] = rec
    return models, index


def match_model(cell_line_raw, index):
    for candidate in (cell_line_raw, strip_lineage_suffix(cell_line_raw)):
        rec = index.get(normalise(candidate))
        if rec:
            return rec
    return None


# -------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(
        description="Build the SRA samplesheet for CCLE RNA-seq, filtered by DepMap lineage.")
    p.add_argument("--bioproject", default="PRJNA523380",
                   help="SRA BioProject accession (default: PRJNA523380, CCLE RNA-seq)")
    p.add_argument("--model-csv",
                   help="DepMap Model.csv (ModelID/CCLEName/StrippedCellLineName/OncotreeLineage); "
                        "required unless --print-query")
    p.add_argument("--lineage", default="Liver",
                   help="OncotreeLineage to keep (default: Liver); 'all' disables the filter")
    p.add_argument("--strategy", default="RNA-Seq",
                   help="Keep only runs with this LIBRARY_STRATEGY (default: RNA-Seq; 'all' disables)")
    p.add_argument("--paired-only", action="store_true",
                   help="Drop runs whose LIBRARY_LAYOUT is not PAIRED (the pipeline needs paired reads)")
    p.add_argument("--one-run-per-line", action="store_true",
                   help="Keep only the largest run (by total_bases) per cell line; required for --fused")
    p.add_argument("--gs-uris", action="store_true",
                   help="Also print the gs:// mirror URIs to stderr as a summary")
    p.add_argument("--api-key", help="NCBI API key (raises the rate limit from 3/s to 10/s)")
    p.add_argument("--print-query", action="store_true",
                   help="Print the SRA search term and exit")
    p.add_argument("-o", "--output", default="sra_samplesheet.csv")
    args = p.parse_args()

    # Filter server-side where SRA supports it: PRJNA523380 is a mixed-assay
    # project (4550 packages, mostly WGS) and only ~1019 are RNA-Seq, so this
    # cuts the efetch paging by ~5x. The client-side checks below still run as
    # a safety net in case the index and the XML disagree.
    term = f"{args.bioproject}[BioProject]"
    if args.strategy.lower() != "all":
        term += f' AND "{args.strategy}"[Strategy]'
    if args.paired_only:
        term += ' AND "paired"[Layout]'
    if args.print_query:
        print(term)
        return 0
    if not args.model_csv:
        p.error("--model-csv is required (except with --print-query)")

    sys.stderr.write(f"Searching SRA: {term}\n")
    count, webenv, query_key = esearch(term, args.api_key)
    if count == 0:
        raise SystemExit(f"No SRA experiments found for {term}")
    sys.stderr.write(f"  {count} experiment package(s)\n")

    delay = 0.11 if args.api_key else 0.34
    runs = []
    for pkg in efetch_runs(webenv, query_key, count, args.api_key, delay):
        runs.extend(parse_package(pkg))
    sys.stderr.write(f"Parsed {len(runs)} run(s)\n")

    if args.strategy.lower() != "all":
        before = len(runs)
        runs = [r for r in runs if r["library_strategy"] == args.strategy]
        sys.stderr.write(f"  {len(runs)} of {before} are {args.strategy}\n")
    if args.paired_only:
        before = len(runs)
        runs = [r for r in runs if r["library_layout"] == "PAIRED"]
        sys.stderr.write(f"  {len(runs)} of {before} are PAIRED\n")

    sys.stderr.write(f"Loading DepMap models: {args.model_csv}\n")
    models, index = load_models(args.model_csv)
    sys.stderr.write(f"  {len(models)} models, {len(index)} lookup keys\n")

    rows, unmatched = [], []
    for r in runs:
        rec = match_model(r["cell_line_raw"], index)
        if rec is None:
            unmatched.append(r)
            continue
        if args.lineage.lower() != "all" and (rec["lineage"] or "").lower() != args.lineage.lower():
            continue
        name = rec["stripped"] or strip_lineage_suffix(rec["ccle_name"]) or rec["model_id"]
        rows.append({
            "sample": sanitise_sample(name),
            "srr_run": r["srr_run"],
            "model_id": rec["model_id"],
            "ccle_name": rec["ccle_name"],
            "stripped_cell_line_name": rec["stripped"],
            "oncotree_lineage": rec["lineage"],
            "oncotree_primary_disease": rec["disease"],
            "library_layout": r["library_layout"],
            "library_strategy": r["library_strategy"],
            "total_spots": r["total_spots"],
            "total_bases": r["total_bases"],
            "biosample": r["biosample"],
            "sra_experiment": r["sra_experiment"],
            "sra_sample": r["sra_sample"],
            "sra_gs_uri": f"{SRA_GS_PREFIX}/{r['srr_run']}/{r['srr_run']}",
        })

    if args.lineage.lower() != "all":
        sys.stderr.write(f"  {len(rows)} run(s) with OncotreeLineage == {args.lineage}\n")

    # Several runs per cell line are normally lanes of one library.
    per_line = Counter(r["sample"] for r in rows)
    multi = {k: v for k, v in per_line.items() if v > 1}
    if args.one_run_per_line and multi:
        best = {}
        for r in rows:
            prev = best.get(r["sample"])
            def size(x):
                try:
                    return int(x["total_bases"] or 0)
                except ValueError:
                    return 0
            if prev is None or size(r) > size(prev):
                best[r["sample"]] = r
        dropped = len(rows) - len(best)
        rows = sorted(best.values(), key=lambda x: x["sample"])
        sys.stderr.write(f"  --one-run-per-line: kept 1 run for each of {len(rows)} "
                         f"cell line(s), dropped {dropped}\n")
    else:
        rows.sort(key=lambda x: (x["sample"], x["srr_run"]))

    if not rows:
        raise SystemExit("No runs left after filtering -- check --lineage / --strategy.")

    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    cell_lines = len({r["sample"] for r in rows})
    sys.stderr.write(f"\nWrote {args.output}: {len(rows)} run(s) across {cell_lines} cell line(s)\n")

    if args.gs_uris:
        for r in rows:
            sys.stderr.write(f"  {r['sample']}\t{r['srr_run']}\t{r['sra_gs_uri']}\n")

    multi = {k: v for k, v in Counter(r["sample"] for r in rows).items() if v > 1}
    if multi:
        sys.stderr.write(
            f"\n{len(multi)} cell line(s) have more than one run: "
            f"{', '.join(f'{k} ({v})' for k, v in sorted(multi.items())[:10])}"
            f"{' ...' if len(multi) > 10 else ''}\n"
            "  Pass --merge_same_sample true to main.nf to merge them into one Salmon run,\n"
            "  or re-run this script with --one-run-per-line (required for --fused).\n")

    if unmatched:
        names = sorted({r['cell_line_raw'] or '<no name>' for r in unmatched})
        sys.stderr.write(
            f"\nWARNING: {len(unmatched)} run(s) across {len(names)} name(s) matched no DepMap model:\n"
            f"  {', '.join(names[:20])}{' ...' if len(names) > 20 else ''}\n"
            "  These are excluded. A missing liver line here would silently shrink the cohort.\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
