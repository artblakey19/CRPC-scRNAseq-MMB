#!/usr/bin/env bash
# Build per-sample spliced/unspliced loom files from the 10x BAMs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
cd "${PROJECT_ROOT}"

OUT_DIR="${VELOCITY_LOOM_DIR:-Results/07_RNA_Velocity/loom_by_sample}"
THREADS="${THREADS:-8}"
SAMTOOLS_MEMORY="${SAMTOOLS_MEMORY:-4096}"
MAX_JOBS="${MAX_JOBS:-1}"
BAM_BASENAME="${BAM_BASENAME:-possorted_genome_bam.bam}"
GENES_GTF="${GENES_GTF:-${1:-}}"
REPEAT_MASK_GTF="${REPEAT_MASK_GTF:-${2:-}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-codex}"
mkdir -p "${MPLCONFIGDIR}"

if [[ -z "${GENES_GTF}" || ! -f "${GENES_GTF}" ]]; then
    cat >&2 <<EOF
Usage:
  GENES_GTF=/path/to/genes.gtf bash scripts/07_RNA_Velocity/01_make_velocity_looms.sh

Optional:
  REPEAT_MASK_GTF=/path/to/repeat_mask.gtf
  THREADS=16
  MAX_JOBS=2
  BAM_BASENAME=possorted_genome_bam.bam

The file Resources/hg38_gencode_v27.txt is not a GTF and cannot be used for
velocyto counting. Prefer the genes.gtf from the same Cell Ranger reference
used to align the BAMs.
EOF
    exit 2
fi

if [[ -n "${REPEAT_MASK_GTF}" && ! -f "${REPEAT_MASK_GTF}" ]]; then
    echo "Repeat-mask GTF does not exist: ${REPEAT_MASK_GTF}" >&2
    exit 2
fi

if command -v velocyto >/dev/null 2>&1; then
    VELOCYTO_BIN="$(command -v velocyto)"
elif [[ -x ".venv/bin/velocyto" ]]; then
    VELOCYTO_BIN=".venv/bin/velocyto"
else
    cat >&2 <<EOF
velocyto.py is not installed in this environment.

Install it into the active Python environment, then rerun this script, e.g.:
  .venv/bin/python -m pip install velocyto

scVelo is already installed here; this script only needs velocyto.py for the
BAM -> loom counting step.
EOF
    exit 127
fi

mkdir -p "${OUT_DIR}"

run_sample() {
    local sample="$1"
    bam="Raw_data/CRPC_bam/${sample}/${BAM_BASENAME}"
    barcodes="Raw_data/${sample}/filtered_feature_bc_matrix/barcodes.tsv.gz"

    if [[ ! -f "${bam}" ]]; then
        echo "Missing BAM for ${sample}: ${bam}" >&2
        exit 2
    fi
    if [[ ! -f "${barcodes}" ]]; then
        echo "Missing barcode file for ${sample}: ${barcodes}" >&2
        exit 2
    fi

    echo "[velocyto] ${sample}"
    cmd=(
        "${VELOCYTO_BIN}" run
        -b "${barcodes}"
        -o "${OUT_DIR}"
        -e "${sample}"
        -@ "${THREADS}"
        --samtools-memory "${SAMTOOLS_MEMORY}"
    )
    if [[ -n "${REPEAT_MASK_GTF}" ]]; then
        cmd+=(-m "${REPEAT_MASK_GTF}")
    fi
    cmd+=("${bam}" "${GENES_GTF}")

    "${cmd[@]}"
}

active_jobs=0
for sample in CRPC1 CRPC2 CRPC3; do
    log="${OUT_DIR}/${sample}.velocyto.log"
    echo "[submit] ${sample} -> ${log}"
    run_sample "${sample}" > "${log}" 2>&1 &
    active_jobs=$((active_jobs + 1))

    if (( active_jobs >= MAX_JOBS )); then
        wait -n
        active_jobs=$((active_jobs - 1))
    fi
done

while (( active_jobs > 0 )); do
    wait -n
    active_jobs=$((active_jobs - 1))
done

echo "Velocity loom files written under: ${OUT_DIR}"
