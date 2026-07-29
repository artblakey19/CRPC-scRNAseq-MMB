#!/usr/bin/env bash
# =============================================================================
# Numbat step 2b wrapper  —  CONTAINER / pkharchenkolab/numbat-rbase
#   02a 산출물(count_mat.rds, cell_annot_ref.csv) + step1 allele 로 run_numbat 실행.
#   query = Epithelial only, init_k 등 파라미터는 02b_run_numbat.R 안에 있다.
# 컨테이너엔 Seurat 이 없어 numbat 실행만 여기서 한다.
#
# ★ 파라미터(init_k / min_LLR / query) 튜닝 시 02a 재실행 없이 이 스크립트만 반복 실행하면 됨.
#   (02a 산출물은 파라미터와 무관 → Results/08_Numbat/<sample>/count_mat.rds 재사용)
# =============================================================================
set -euo pipefail

# 프로젝트 루트 = 이 스크립트(scripts/08_Numbat/) 기준 두 단계 위. 절대경로(docker --mount 용).
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="pkharchenkolab/numbat-rbase:latest"
cd "${PROJ}"
mkdir -p logs

echo "[$(date '+%F %T')] === 2b: run_numbat (container / numbat) ==="
# -e R_PROFILE_USER=/dev/null : 프로젝트 .Rprofile(renv) 로드 차단 → 컨테이너 numbat 사용
docker run --rm \
    -u "$(id -u):$(id -g)" -e HOME=/tmp -e R_PROFILE_USER=/dev/null \
    --mount type=bind,source="${PROJ}",target=/mnt/project \
    -w /mnt/project \
    "${IMAGE}" \
    Rscript /mnt/project/scripts/08_Numbat/02b_run_numbat.R \
    2>&1 | tee logs/numbat_run.log

echo "[$(date '+%F %T')] === 2b done -> Results/08_Numbat/<sample>/numbat/ ==="
