#!/bin/bash
#PBS -N sra_prefetch
#PBS -q copyq
#PBS -P vh48
#PBS -j oe
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l storage=gdata/vh48+gdata/io72
#PBS -l walltime=10:00:00
#PBS -l wd

# Stage 1 of 2 for Gadi: download the .sra files only (needs the network, so
# copyq). Stage 2, dump_sra_fastq.pbs, turns them into FASTQ on a many-CPU
# node of the normal queue. Splitting them keeps the slow single-threaded
# fasterq-dump + gzip work off the 1-CPU copyq.
#
# CCLE runs are ALIGNED (cSRA): the .sra holds differences against GRCh38, so
# prefetch also pulls ~72 reference sequences, and fasterq-dump needs them
# locally when it later runs with no network. Everything therefore stays under
# one directory ($OUTDIR/ncbi + $OUTDIR/sra), and both stages use the same
# NCBI_SETTINGS file. Do not delete $OUTDIR/ncbi between stages.
#
# Usage:
#   qsub -v SRA_SHEET=sra_samplesheet_liver.csv,OUTDIR=/g/data/vh48/AmaroqTx/ccle_fastq prefetch_sra.sh
#
# Idempotent: a run is skipped when its .done marker exists (written only
# after prefetch exited 0) or its final FASTQ pair already exists. Resubmit
# after a walltime kill.

set -euo pipefail

SRA_SHEET="${1:-${SRA_SHEET:-sra_samplesheet.csv}}"
OUTDIR="${2:-${OUTDIR:-fastq}}"

cd "${PBS_O_WORKDIR:-$PWD}"

if ! command -v prefetch >/dev/null && command -v module >/dev/null; then
    module use /g/data/io72/apps/Modules/modulefiles
    module load "${SRA_MODULE:-sra-tools/3.2.1}"
fi
command -v prefetch >/dev/null || { echo "ERROR: prefetch not on PATH (load sra-tools)" >&2; exit 2; }
[ -f "$SRA_SHEET" ] || { echo "ERROR: samplesheet not found: $SRA_SHEET" >&2; exit 2; }

mkdir -p "$OUTDIR/sra" "$OUTDIR/ncbi"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# Same file as dump_sra_fastq.pbs (it is only created if missing). The
# repository root keeps prefetch's reference downloads beside the data instead
# of in $HOME; UNVERIFIED that the key is honoured, so check the log for where
# the ncbi-acc: files land (find $OUTDIR -name 'NC_012920*').
export NCBI_SETTINGS="$OUTDIR/user-settings.mkfg"
if [ ! -f "$NCBI_SETTINGS" ]; then
    cat > "$NCBI_SETTINGS" <<MKFG
/LIBS/GUID = "5c0d1f66-0a0b-4e1a-9a5a-2d1f6b5a7c31"
/libs/cloud/report_instance_identity = "true"
/repository/user/main/public/root = "$OUTDIR/ncbi"
MKFG
fi

echo "samplesheet : $SRA_SHEET"
echo "outdir      : $OUTDIR"
echo

read -r HEADER < "$SRA_SHEET"
col() { echo "$HEADER" | tr ',' '\n' | grep -nx "$1" | cut -d: -f1; }
SRR_COL="$(col srr_run || true)"
[ -n "$SRR_COL" ] || { echo "ERROR: $SRA_SHEET needs a 'srr_run' column; got: $HEADER" >&2; exit 2; }

n_done=0; n_skip=0; n_fail=0
while IFS=, read -r -a FIELDS; do
    SRR="${FIELDS[$((SRR_COL-1))]}"
    [ -n "${SRR:-}" ] || continue

    if [ -e "$OUTDIR/sra/$SRR.done" ] || { [ -s "$OUTDIR/${SRR}_R1.fastq.gz" ] && [ -s "$OUTDIR/${SRR}_R2.fastq.gz" ]; }; then
        echo "[skip] $SRR"
        n_skip=$((n_skip+1)); continue
    fi

    echo "[get ] $SRR"
    rm -rf "$OUTDIR/sra/$SRR"   # a partial directory from a killed job
    # --max-size u: the default 20 GB cap would skip big runs with a warning.
    if prefetch --max-size u --output-directory "$OUTDIR/sra" "$SRR" \
            && [ -s "$OUTDIR/sra/$SRR/$SRR.sra" ]; then
        touch "$OUTDIR/sra/$SRR.done"
        n_done=$((n_done+1))
    else
        echo "[FAIL] prefetch $SRR" >&2
        rm -rf "$OUTDIR/sra/$SRR"; n_fail=$((n_fail+1))
    fi
done < <(tail -n +2 "$SRA_SHEET")

echo
echo "prefetched : $n_done"
echo "skipped    : $n_skip"
echo "failed     : $n_fail"
echo
echo "Next: qsub -v SRA_SHEET=$SRA_SHEET,OUTDIR=$OUTDIR dump_sra_fastq.pbs"
[ "$n_fail" -eq 0 ] || { echo "ERROR: $n_fail run(s) failed; resubmit to retry." >&2; exit 1; }
