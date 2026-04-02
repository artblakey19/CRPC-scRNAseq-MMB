#!/usr/bin/env bash
set -euo pipefail

SRC="/mnt/z/ZJ Sun Lab/incoming_COH_data/labs/ZJsun grp/Old Lab Members/WK"
DST="/mnt/c/Users/kcw78/Desktop/WK_Rscript"

mkdir -p "$DST"

# 디렉터리는 통과(*/), .R/.r 파일은 포함, 그 외는 모두 제외
rsync -av --include='*/' --include='*.R' --include='*.r' --exclude='*' "$SRC"/ "$DST"/
