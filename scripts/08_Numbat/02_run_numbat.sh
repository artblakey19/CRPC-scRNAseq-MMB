#!/usr/bin/env bash
# =============================================================================
# Numbat step 2 orchestrator  —  2a(host/Seurat) → 2b(container/numbat) 순차 실행 (전체 재실행용)
#   2a: host renv 에서 Rscript 로 입력 준비 (numbat 파라미터와 무관)
#   2b: 02b_run_numbat.sh (docker 컨테이너에서 run_numbat)
# host 엔 numbat 이, 컨테이너엔 Seurat 이 없어 단계를 나눈다.
#
# ★ init_k / min_LLR / query 같은 파라미터만 조정할 땐 이 스크립트 대신
#   02b_run_numbat.sh 만 단독 실행하면 됨 (02a 산출물 재사용).
# =============================================================================
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJ}"
mkdir -p logs

# --- 2a : host (Seurat / renv) ---
echo "[$(date '+%F %T')] === 2a: export inputs (host / Seurat) ==="
Rscript scripts/08_Numbat/02a_prepare_numbat_inputs.R \
    2>&1 | tee logs/numbat_prepare.log

# --- 2b : container (numbat) ---
bash "${HERE}/02b_run_numbat.sh"

echo "[$(date '+%F %T')] === step2 done -> Results/08_Numbat/<sample>/numbat/ ==="
