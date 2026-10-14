#!/usr/bin/env python3
"""
fetch_gdc_metadata.py - build the pipeline samplesheet from the GDC API.

Two ways to select files:

  1. A whole project (all STAR 2-Pass Genome RNA-Seq BAMs):
       python3 fetch_gdc_metadata.py --project TCGA-LIHC -o samplesheet.csv

  2. A list of GDC file IDs (text file with one ID per line, or a CSV/TSV with a
     column named gdc_file_id, file_id, "File ID" or id - e.g. a GDC manifest,
     a GDC cart sample sheet, or an old samplesheet):
       python3 fetch_gdc_metadata.py ids.txt -o samplesheet.csv

Add --print-query to just print the JSON request body (for use with curl) and exit.

Metadata queries are public (no token needed). Uses only the Python standard library.

Output columns:
  sample,gdc_file_id,case,tissue_type,sample_type,tumor_descriptor,project_id,file_name

  sample = the TCGA sample barcode, e.g. TCGA-BC-A110-01A (tumor) / TCGA-BC-A110-11A (normal)
"""
import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from collections import Counter

API_URL = "https://api.gdc.cancer.gov/files"
ID_BATCH_SIZE = 100
PAGE_SIZE = 500

FIELDS = [
    "file_id",
    "file_name",
    "data_type",
    "experimental_strategy",
    "analysis.workflow_type",
    "cases.case_id",
    "cases.submitter_id",
    "cases.project.project_id",
    "cases.samples.sample_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.samples.tissue_type",
    "cases.samples.tumor_descriptor",
]

ID_COLUMNS = ("gdc_file_id", "file_id", "File ID", "id")
PLAIN_HEADERS = ("id", "file_id", "gdc_file_id", "file id")

OUT_COLUMNS = [
    "sample", "gdc_file_id", "case", "tissue_type", "sample_type",
    "tumor_descriptor", "project_id", "file_name",
]


# --------------------------------------------------------------------------- queries
def in_filter(field, values):
    return {"op": "in", "content": {"field": field, "value": list(values)}}


def project_filters(project, workflow_type, strategy):
    return {"op": "and", "content": [
        in_filter("cases.project.project_id", [project]),
        in_filter("data_type", ["Aligned Reads"]),
        in_filter("data_format", ["BAM"]),
        in_filter("experimental_strategy", [strategy]),
        in_filter("analysis.workflow_type", [workflow_type]),
    ]}


def make_payload(filters, size, offset=0):
    return {
        "filters": filters,
        "fields": ",".join(FIELDS),
        "format": "JSON",
        "size": size,
        "from": offset,
    }


def query_gdc(payload, retries=3):
    request = urllib.request.Request(
        API_URL, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"}
    )
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.load(response)["data"]
        except urllib.error.HTTPError as err:
            sys.exit(f"ERROR: GDC API returned {err.code}: {err.read().decode(errors='replace')[:500]}")
        except urllib.error.URLError as err:
            if attempt == retries:
                sys.exit(f"ERROR: cannot reach the GDC API ({err})")
            time.sleep(2 * attempt)


def fetch_by_ids(ids):
    hits = []
    for i in range(0, len(ids), ID_BATCH_SIZE):
        batch = ids[i:i + ID_BATCH_SIZE]
        data = query_gdc(make_payload(in_filter("file_id", batch), len(batch)))
        hits.extend(data["hits"])
    return hits


def fetch_project(project, workflow_type, strategy, page_size=PAGE_SIZE):
    filters = project_filters(project, workflow_type, strategy)
    hits, offset = [], 0
    while True:
        data = query_gdc(make_payload(filters, page_size, offset))
        hits.extend(data["hits"])
        total = data["pagination"]["total"]
        offset += page_size
        if offset >= total or not data["hits"]:
            return hits, total


# --------------------------------------------------------------------------- input
def read_ids(path):
    with open(path, newline="") as fh:
        lines = [line for line in fh.read().splitlines() if line.strip()]
    if not lines:
        sys.exit(f"ERROR: {path} is empty")

    first = lines[0]
    delim = "\t" if "\t" in first else ("," if "," in first else None)

    if delim is None:  # plain list, one ID per line
        ids = [line.strip() for line in lines]
        if ids[0].lower() in PLAIN_HEADERS:
            ids = ids[1:]
        return ids

    reader = csv.DictReader(lines, delimiter=delim)
    col = next((c for c in ID_COLUMNS if c in (reader.fieldnames or [])), None)
    if col is None:
        sys.exit(f"ERROR: no file-ID column found in {path}. Expected one of: {', '.join(ID_COLUMNS)}")
    return [row[col].strip() for row in reader if (row.get(col) or "").strip()]


# --------------------------------------------------------------------------- parsing
def infer_tissue_from_barcode(barcode):
    """TCGA sample codes: 01-09 tumor, 10-19 normal (e.g. TCGA-BC-A110-01A / -11A)."""
    parts = barcode.split("-")
    if len(parts) >= 4 and parts[3][:2].isdigit():
        code = int(parts[3][:2])
        if 1 <= code <= 9:
            return "Tumor"
        if 10 <= code <= 19:
            return "Normal"
    return None


def build_rows(hits, requested):
    rows, problems, warnings = [], [], []
    found = set()

    for hit in hits:
        fid = hit.get("file_id")
        found.add(fid)

        cases = hit.get("cases") or []
        if len(cases) != 1:
            problems.append(f"{fid}: expected exactly 1 case, found {len(cases)}")
            continue
        case = cases[0]

        samples = case.get("samples") or []
        if len(samples) != 1:
            problems.append(f"{fid}: expected exactly 1 sample, found {len(samples)}")
            continue
        sample = samples[0]

        case_id = case.get("submitter_id") or ""
        barcode = sample.get("submitter_id") or f"{case_id}_{fid[:8]}"

        tissue = sample.get("tissue_type")
        if not tissue or tissue in ("Not Reported", "Unknown"):
            inferred = infer_tissue_from_barcode(barcode)
            if inferred:
                warnings.append(f"{fid}: tissue_type missing in GDC; inferred '{inferred}' from barcode {barcode}")
                tissue = inferred
            else:
                tissue = tissue or "Unknown"

        file_name = hit.get("file_name", "")
        if hit.get("data_type") != "Aligned Reads" or ".genomic." not in file_name:
            warnings.append(
                f"{fid}: '{file_name}' does not look like a genomic aligned-reads BAM "
                f"(data_type={hit.get('data_type')}); check that this is the file you want"
            )

        rows.append({
            "sample": barcode,
            "gdc_file_id": fid,
            "case": case_id,
            "tissue_type": tissue,
            "sample_type": sample.get("sample_type") or "",
            "tumor_descriptor": sample.get("tumor_descriptor") or "",
            "project_id": (case.get("project") or {}).get("project_id", ""),
            "file_name": file_name,
        })

    for fid in requested:
        if fid not in found:
            problems.append(f"{fid}: not returned by the GDC API (wrong ID?)")

    # Sample names must be unique (two BAMs of one sample vial would be merged otherwise)
    counts = Counter(r["sample"] for r in rows)
    for r in rows:
        if counts[r["sample"]] > 1:
            new_name = f"{r['sample']}_{r['gdc_file_id'][:8]}"
            warnings.append(f"{r['gdc_file_id']}: sample name {r['sample']} is shared by several files; renamed to {new_name}")
            r["sample"] = new_name

    return rows, problems, warnings


# --------------------------------------------------------------------------- filtering / ordering
def is_recurrence(row):
    """Recurrent tumors: tumor_descriptor 'Recurrence' / sample_type 'Recurrent Solid Tumor'."""
    text = f"{row['tumor_descriptor']} {row['sample_type']}".lower()
    return "recur" in text


def drop_recurrences(rows):
    kept = [r for r in rows if not is_recurrence(r)]
    dropped = [r for r in rows if is_recurrence(r)]
    return kept, dropped


def paired_cases(rows):
    """Cases that have at least one Tumor and one Normal sample in `rows`."""
    tissues = {}
    for r in rows:
        tissues.setdefault(r["case"], set()).add(r["tissue_type"])
    return {case for case, t in tissues.items() if {"Tumor", "Normal"} <= t}


# --------------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description="Build the pipeline samplesheet from the GDC API.")
    parser.add_argument("input", nargs="?", help="CSV/TSV with a file-ID column, or a text file with one ID per line")
    parser.add_argument("--project", help="GDC project ID, e.g. TCGA-LIHC (selects all matching BAMs)")
    parser.add_argument("--workflow-type", default="STAR 2-Pass Genome",
                        help="analysis.workflow_type for --project mode (default: %(default)s)")
    parser.add_argument("--experimental-strategy", default="RNA-Seq",
                        help="experimental_strategy for --project mode (default: %(default)s)")
    parser.add_argument("--include-recurrence", action="store_true",
                        help="keep recurrent-tumor samples (excluded by default)")
    parser.add_argument("--print-query", action="store_true",
                        help="print the JSON request body for --project mode and exit")
    parser.add_argument("-o", "--output", default="samplesheet_annotated.csv")
    args = parser.parse_args()

    if bool(args.input) == bool(args.project):
        parser.error("give exactly one of: an input file, or --project")

    if args.print_query:
        if not args.project:
            parser.error("--print-query needs --project")
        filters = project_filters(args.project, args.workflow_type, args.experimental_strategy)
        print(json.dumps(make_payload(filters, 1000), indent=2))
        return

    if args.project:
        hits, total = fetch_project(args.project, args.workflow_type, args.experimental_strategy)
        print(f"GDC reports {total} files matching {args.project} / {args.experimental_strategy} / "
              f"{args.workflow_type}", file=sys.stderr)
        requested = [h["file_id"] for h in hits]
    else:
        ids = read_ids(args.input)
        requested = list(dict.fromkeys(ids))
        if len(requested) != len(ids):
            print(f"WARNING: removed {len(ids) - len(requested)} duplicate file ID(s)", file=sys.stderr)
        hits = fetch_by_ids(requested)

    rows, problems, warnings = build_rows(hits, requested)

    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for p in problems:
        print(f"ERROR: {p}", file=sys.stderr)

    if not args.include_recurrence:
        rows, dropped = drop_recurrences(rows)
        for r in dropped:
            print(f"EXCLUDED (recurrence): {r['sample']}  {r['gdc_file_id']}  "
                  f"[{r['sample_type']} / {r['tumor_descriptor']}]", file=sys.stderr)
        print(f"Excluded {len(dropped)} recurrence sample(s)", file=sys.stderr)

    # Tumor tissue that is neither Primary nor a (now excluded) recurrence stays in, but is flagged
    other_tumors = [r for r in rows if r["tissue_type"] == "Tumor" and r["tumor_descriptor"] != "Primary"]
    if other_tumors:
        kinds = Counter(f"{r['sample_type']} / {r['tumor_descriptor']}" for r in other_tumors)
        print("NOTE: kept " + str(len(other_tumors)) + " non-primary tumor sample(s): "
              + ", ".join(f"{k} x{v}" for k, v in kinds.items()), file=sys.stderr)

    # Cases with both Tumor and Normal first; then by case, Normal before Tumor
    paired = paired_cases(rows)
    rows.sort(key=lambda r: (r["case"] not in paired, r["case"], r["tissue_type"], r["sample"]))
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    tissue_counts = Counter(r["tissue_type"] for r in rows)
    n_cases = len({r["case"] for r in rows})

    print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)
    print("Tissue types: " + ", ".join(f"{k}={v}" for k, v in sorted(tissue_counts.items())), file=sys.stderr)
    print(f"Cases: {n_cases} ({len(paired)} with both Tumor and Normal, listed first)", file=sys.stderr)

    if problems:
        sys.exit(f"ERROR: {len(problems)} file(s) skipped because of the problems above")


if __name__ == "__main__":
    main()
