#!/usr/bin/env bash
# Numbat step 1 — allele pileup (cellsnp-lite) + population phasing (Eagle2)
# Runs pkharchenkolab/numbat-rbase in Docker, one patient at a time.
# Output: Results/08_Numbat/<sample>/{pileup,phasing}/  +  <sample>_allele_counts.tsv.gz
set -euo pipefail

# 프로젝트 루트 = 이 스크립트(scripts/08_Numbat/) 기준 두 단계 위. 절대경로(docker --mount 용).
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${PROJ}/Results/08_Numbat"
IMAGE="pkharchenkolab/numbat-rbase:latest"
NCORES=30
SAMPLES=(CRPC1 CRPC2 CRPC3)

mkdir -p "${OUT}" "${PROJ}/logs"

for S in "${SAMPLES[@]}"; do
    mkdir -p "${OUT}/${S}"

    # cellsnp-lite gets a plain-text barcode list (CellRanger ships it gzipped)
    if [[ ! -s "${OUT}/${S}/barcodes.tsv" ]]; then
        zcat "${PROJ}/Raw_data/${S}/filtered_feature_bc_matrix/barcodes.tsv.gz" \
            > "${OUT}/${S}/barcodes.tsv"
    fi

    echo "[$(date '+%F %T')] === ${S}: pileup + phasing ==="
    # -e R_PROFILE_USER=/dev/null : 프로젝트 .Rprofile(renv/activate.R) 로드를 막아
    #   컨테이너 자체 R 라이브러리(optparse·numbat)를 그대로 쓰게 함
    docker run --rm \
        -u "$(id -u):$(id -g)" -e HOME=/tmp -e R_PROFILE_USER=/dev/null \
        --mount type=bind,source="${PROJ}",target=/mnt/project \
        -w /mnt/project \
        "${IMAGE}" \
        Rscript /numbat/inst/bin/pileup_and_phase.R \
            --label "${S}" \
            --samples "${S}" \
            --bams "/mnt/project/Raw_data/CRPC_bam/${S}/possorted_genome_bam.bam" \
            --barcodes "/mnt/project/Results/08_Numbat/${S}/barcodes.tsv" \
            --outdir "/mnt/project/Results/08_Numbat/${S}" \
            --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
            --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
            --paneldir /data/1000G_hg38 \
            --ncores ${NCORES} \
        2>&1 | tee "${PROJ}/logs/numbat_pileup_${S}.log"

    echo "[$(date '+%F %T')] === ${S}: done ==="
done