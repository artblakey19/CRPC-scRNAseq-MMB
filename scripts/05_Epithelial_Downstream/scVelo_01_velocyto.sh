#!/bin/bash
# scVelo step 1 — velocyto per-sample spliced/unspliced count quantification.
#
# Cellranger BAM 의 CB / UB tag 와 GENCODE GTF 를 사용해 spliced / unspliced /
# ambiguous count matrix 를 만들고 .loom 으로 저장. scVelo 가 입력으로 쓰는
# 표준 형식.
#
# Output (per sample):
#   Results/05_Epithelial_Downstream/scVelo/loom/{sample}.loom
#
# Prereqs:
#   - conda env "scvelo" with velocyto.py + samtools
#   - GENCODE v32 primary assembly GTF in ~/refs/velocyto/
#     (cellranger refdata-gex-GRCh38-2020-A 와 동일 버전 — barcode/gene id 호환)
#   - Cellranger BAM (Raw_data/CRPC_bam/{sample}/possorted_genome_bam.bam)
#   - Cellranger filtered barcodes (Raw_data/{sample}/filtered_feature_bc_matrix/barcodes.tsv.gz)

set -euo pipefail

# ============================== Config =======================================
PROJECT="/home/MMB/projects/Human CRPC scRNAseq"
REFS="$HOME/refs/velocyto"
SAMPLES=(CRPC1 CRPC2 CRPC3)
NTHREADS=30
SAMTOOLS_MEMORY=6000   # MB per samtools-sort thread → 30 × 6 GB = 180 GB peak (RAM 226 GB)
                       # 이전 실패의 진짜 원인은 메모리 한도 아니라 /tmp 가 tmpfs (RAM)
                       # 라서 samtools sort temp 파일이 RAM 을 또 잡아먹어 OOM 였음.
                       # 아래 TMPDIR override 로 sort temp 를 디스크로 빼서 해결.
CONDA_ENV="scvelo"

GTF="${REFS}/gencode.v32.primary_assembly.annotation.gtf"
# repeat mask (optional) — intronic 영역의 repeat-derived reads 마스킹용
RMSK_GTF="${REFS}/hg38_rmsk.gtf"

# ============================== Init ========================================
source ~/miniforge3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

if [ ! -f "$GTF" ] && [ -f "${GTF}.gz" ]; then
    echo "Decompressing GTF..."
    gunzip -k "${GTF}.gz"
fi
if [ ! -f "$GTF" ]; then echo "Missing GTF: $GTF"; exit 1; fi

cd "$PROJECT"
OUTDIR_BASE="Results/05_Epithelial_Downstream/scVelo/loom"
mkdir -p "$OUTDIR_BASE"

# /tmp 가 tmpfs (RAM 기반) 이므로 samtools sort temp 를 로컬 디스크로 강제.
# velocyto 내부의 samtools sort 호출도 TMPDIR 을 따라간다.
export TMPDIR="${PROJECT}/.tmp_sort"
mkdir -p "$TMPDIR"

RMSK_FLAG=""
if [ -f "$RMSK_GTF" ]; then RMSK_FLAG="--mask ${RMSK_GTF}"; fi

# ============================== Run =========================================
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="Raw_data/CRPC_bam/${SAMPLE}/possorted_genome_bam.bam"
    BARCODES="Raw_data/${SAMPLE}/filtered_feature_bc_matrix/barcodes.tsv.gz"
    LOOM_OUT="${OUTDIR_BASE}/${SAMPLE}.loom"

    if [ ! -f "$BAM" ]; then echo "Missing BAM: $BAM"; exit 1; fi
    if [ ! -f "$BARCODES" ]; then echo "Missing barcodes: $BARCODES"; exit 1; fi
    if [ -f "$LOOM_OUT" ]; then
        echo "Skip ${SAMPLE} — loom exists: ${LOOM_OUT}"
        continue
    fi

    # velocyto run 은 출력 디렉토리에 BAM 이름 기반 .loom 을 만듦. tmp 로 받은 뒤
    # 일관된 이름으로 옮겨준다.
    TMPOUT="${OUTDIR_BASE}/${SAMPLE}_tmp"
    mkdir -p "$TMPOUT"

    # decompress barcodes (velocyto 가 plain text 요구)
    BC_PLAIN="${TMPOUT}/barcodes.tsv"
    zcat "$BARCODES" > "$BC_PLAIN"

    echo ""
    echo "================================================================"
    echo "=== ${SAMPLE} velocyto start $(date) ==="
    echo "================================================================"
    velocyto run \
        -b "$BC_PLAIN" \
        -o "$TMPOUT" \
        -e "$SAMPLE" \
        ${RMSK_FLAG} \
        --samtools-threads "${NTHREADS}" \
        --samtools-memory  "${SAMTOOLS_MEMORY}" \
        "$BAM" \
        "$GTF" \
        2>&1 | tee "${OUTDIR_BASE}/${SAMPLE}_velocyto.log"

    mv "${TMPOUT}/${SAMPLE}.loom" "${LOOM_OUT}"
    rm -rf "$TMPOUT"
    echo "=== ${SAMPLE} velocyto done $(date) — ${LOOM_OUT} ==="
done

echo ""
echo "================================================================"
echo "=== ALL VELOCYTO LOOMS DONE $(date) ==="
echo "================================================================"
