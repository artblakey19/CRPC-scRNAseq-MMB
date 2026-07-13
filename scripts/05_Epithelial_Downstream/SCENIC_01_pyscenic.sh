#!/usr/bin/env bash
# pySCENIC 3-step pipeline on the exported epi loom.
# Run: bash SCENIC_01_pyscenic.sh
set -euo pipefail

PROJ="/home/MMB/projects/Human CRPC scRNAseq"
PY="$PROJ/.venv/bin/python3"
PYSCENIC="$PROJ/.venv/bin/pyscenic"

OUT="$PROJ/Results/05_Epithelial_Downstream/SCENIC"
RES="$PROJ/Resources/SCENIC"
LOOM_FULL="$OUT/epi_counts.loom"        # used for ctx + aucell
LOOM_GRN="$OUT/epi_counts_grn5k.loom"   # 5k-cell subsample for GRN (avoids dask msg-size limit)

TFS="$RES/hs_hgnc_curated_tfs.txt"
DB1="$RES/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
DB2="$RES/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOT="$RES/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

ADJ="$OUT/adjacencies.tsv"
REG="$OUT/regulons.csv"
AUC_LOOM="$OUT/epi_pyscenic_aucell.loom"

N_WORKERS=16   # keep below CPU count; large worker count amplifies serialization size

for f in "$LOOM_FULL" "$LOOM_GRN"; do
  [[ -f "$f" ]] || { echo "ERROR: $f missing" >&2; exit 1; }
done

echo "[1/3] GRN inference (GRNBoost2) on 5000-cell subsample"
"$PYSCENIC" grn \
  "$LOOM_GRN" "$TFS" \
  -o "$ADJ" \
  --num_workers "$N_WORKERS" \
  --seed 777

echo "[2/3] Context (motif enrichment + pruning) on full expression"
"$PYSCENIC" ctx \
  "$ADJ" "$DB1" "$DB2" \
  --annotations_fname "$ANNOT" \
  --expression_mtx_fname "$LOOM_FULL" \
  --output "$REG" \
  --num_workers "$N_WORKERS" \
  --mask_dropouts

echo "[3/3] AUCell scoring on full expression"
"$PYSCENIC" aucell \
  "$LOOM_FULL" "$REG" \
  --output "$AUC_LOOM" \
  --num_workers "$N_WORKERS"

echo "Done.  outputs:"
echo "  adjacencies: $ADJ"
echo "  regulons:    $REG"
echo "  AUC loom:    $AUC_LOOM"
