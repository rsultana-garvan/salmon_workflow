#!/bin/bash
#PBS -N sra_download
#PBS -q copyq
#PBS -P vh48
#PBS -j oe
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l storage=scratch/vh48+gdata/vh48+gdata/io72
#PBS -l walltime=10:00:00
#PBS -l wd

# Pre-download SRA runs to FASTQ, then write the fastq samplesheet that
# main.nf's --fastq_samplesheet mode consumes.
#
# WHY THIS EXISTS (Gadi): only copyq has outbound network, but copyq is capped
# at 1 CPU, so it cannot also run Salmon. The normal queue has the CPUs but no
# network. So on Gadi the SRA path splits in two:
#
#   1. qsub download_sra_fastq.sh                      <- copyq, this script
#   2. qsub on normal:  nextflow run main.nf \
#          --fastq_samplesheet <OUTDIR>/fastq_samplesheet.csv ...
#
# Elsewhere (the GCP VM, a laptop) you do NOT need this: run main.nf with
# --sra_samplesheet and it does the download itself, or --fused true to do
# download + dump + Salmon in one task.
#
# Usage, anywhere:
#   bash download_sra_fastq.sh sra_samplesheet.csv /scratch/vh48/$USER/ccle_fastq
# or via PBS (args through -v, or edit the defaults below):
#   qsub -v SRA_SHEET=...,OUTDIR=... download_sra_fastq.sh
#
# Idempotent: a run whose R1+R2 already exist is skipped, so a job that hits
# the copyq walltime can simply be resubmitted and will carry on.
#
# Requires prefetch + fasterq-dump (sra-tools) on PATH. On Gadi load whatever
# module or conda env provides them, e.g.
#   module use /g/data/io72/apps/Modules/modulefiles && module load sratoolkit
# or activate a conda env that has sra-tools. The script checks and stops early.

set -euo pipefail

SRA_SHEET="${1:-${SRA_SHEET:-sra_samplesheet.csv}}"
OUTDIR="${2:-${OUTDIR:-fastq}}"
# copyq gives 1 CPU; elsewhere use what the job/host has.
THREADS="${THREADS:-${PBS_NCPUS:-4}}"

cd "${PBS_O_WORKDIR:-$PWD}"

command -v prefetch      >/dev/null || { echo "ERROR: prefetch not on PATH (load sra-tools)" >&2; exit 2; }
command -v fasterq-dump  >/dev/null || { echo "ERROR: fasterq-dump not on PATH (load sra-tools)" >&2; exit 2; }
[ -f "$SRA_SHEET" ] || { echo "ERROR: samplesheet not found: $SRA_SHEET" >&2; exit 2; }

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
SHEET_OUT="$OUTDIR/fastq_samplesheet.csv"
STATS="$OUTDIR/fastq_stats"
mkdir -p "$STATS"

# sra-tools refuses to start without a settings file and otherwise writes one
# into $HOME. Keep it beside the data so the job never depends on a writable
# home (and so two jobs cannot race on the same file).
export NCBI_SETTINGS="$OUTDIR/user-settings.mkfg"
if [ ! -f "$NCBI_SETTINGS" ]; then
    cat > "$NCBI_SETTINGS" <<'MKFG'
/LIBS/GUID = "5c0d1f66-0a0b-4e1a-9a5a-2d1f6b5a7c31"
/libs/cloud/report_instance_identity = "true"
MKFG
fi

echo "samplesheet : $SRA_SHEET"
echo "outdir      : $OUTDIR"
echo "threads     : $THREADS"
echo

# Column positions by name, so extra DepMap columns in the sheet do not matter.
read -r HEADER < "$SRA_SHEET"
col() { echo "$HEADER" | tr ',' '\n' | grep -nx "$1" | cut -d: -f1; }
SAMPLE_COL="$(col sample  || true)"
SRR_COL="$(col srr_run || true)"
[ -n "$SAMPLE_COL" ] && [ -n "$SRR_COL" ] || {
    echo "ERROR: $SRA_SHEET needs 'sample' and 'srr_run' columns; got: $HEADER" >&2; exit 2; }

printf 'sample,fastq_1,fastq_2\n' > "$SHEET_OUT"

n_done=0; n_skip=0; n_fail=0
while IFS=, read -r -a FIELDS; do
    SAMPLE="${FIELDS[$((SAMPLE_COL-1))]}"
    SRR="${FIELDS[$((SRR_COL-1))]}"
    [ -n "${SRR:-}" ] || continue

    R1="$OUTDIR/${SRR}_R1.fastq.gz"
    R2="$OUTDIR/${SRR}_R2.fastq.gz"

    if [ -s "$R1" ] && [ -s "$R2" ]; then
        echo "[skip] $SAMPLE $SRR (already downloaded)"
        n_skip=$((n_skip+1))
        printf '%s,%s,%s\n' "$SAMPLE" "$R1" "$R2" >> "$SHEET_OUT"
        continue
    fi

    echo "[get ] $SAMPLE $SRR"
    WORK="$OUTDIR/.tmp_$SRR"
    rm -rf "$WORK"; mkdir -p "$WORK"

    # --max-size u: prefetch's default 20 GB cap SKIPS larger runs with only a
    # warning, which would silently drop samples from the cohort.
    if ! prefetch --max-size u --output-directory "$WORK" "$SRR"; then
        echo "[FAIL] prefetch $SRR" >&2; n_fail=$((n_fail+1)); rm -rf "$WORK"; continue
    fi
    if ! fasterq-dump --split-3 --threads "$THREADS" \
            --temp "$WORK/tmp" --outdir "$WORK" "$WORK/$SRR/$SRR.sra"; then
        echo "[FAIL] fasterq-dump $SRR" >&2; n_fail=$((n_fail+1)); rm -rf "$WORK"; continue
    fi

    if [ ! -s "$WORK/${SRR}_1.fastq" ] || [ ! -s "$WORK/${SRR}_2.fastq" ]; then
        echo "[FAIL] $SRR produced no paired FASTQ (single-end?)" >&2
        n_fail=$((n_fail+1)); rm -rf "$WORK"; continue
    fi

    N_SINGLE=0
    [ -s "$WORK/${SRR}.fastq" ] && N_SINGLE=$(( $(wc -l < "$WORK/${SRR}.fastq") / 4 ))
    printf 'srr_run\tsingleton_reads\n%s\t%s\n' "$SRR" "$N_SINGLE" > "$STATS/${SRR}.fastq_stats.txt"

    if command -v pigz >/dev/null 2>&1; then
        pigz -p "$THREADS" -c "$WORK/${SRR}_1.fastq" > "$R1"
        pigz -p "$THREADS" -c "$WORK/${SRR}_2.fastq" > "$R2"
    else
        gzip -c "$WORK/${SRR}_1.fastq" > "$R1"
        gzip -c "$WORK/${SRR}_2.fastq" > "$R2"
    fi
    rm -rf "$WORK"

    printf '%s,%s,%s\n' "$SAMPLE" "$R1" "$R2" >> "$SHEET_OUT"
    n_done=$((n_done+1))
done < <(tail -n +2 "$SRA_SHEET")

echo
echo "downloaded : $n_done"
echo "skipped    : $n_skip (already present)"
echo "failed     : $n_fail"
echo "samplesheet: $SHEET_OUT"
echo
echo "Next, on the normal queue (no network needed):"
echo "  nextflow run main.nf --fastq_samplesheet $SHEET_OUT \\"
echo "      --salmon_index <index> --gene_map <gtf> --outdir results_ccle"
[ "$n_fail" -eq 0 ] || { echo "ERROR: $n_fail run(s) failed; see above." >&2; exit 1; }
