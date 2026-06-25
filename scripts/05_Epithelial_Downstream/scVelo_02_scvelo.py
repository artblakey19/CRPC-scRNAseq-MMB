"""scVelo step 2 — dynamical RNA velocity on Stage 3 filtered epithelial.

Inputs:
    - Results/05_Epithelial_Downstream/scVelo/loom/{CRPC1,CRPC2,CRPC3}.loom
        velocyto-derived spliced / unspliced count matrices per sample.
    - Results/05_Epithelial_Downstream/scVelo/epi_metadata_for_scvelo.csv
        Seurat metadata + UMAP coords (scVelo_00_export_seurat.R 가 생성).

Outputs:
    - Results/05_Epithelial_Downstream/scVelo/epi_scvelo.h5ad
    - Results/05_Epithelial_Downstream/scVelo/figures/*.png
        velocity_stream / velocity_arrows / latent_time / phase_portrait

Run from project root with scvelo conda env activated:
    conda activate scvelo
    python scripts/05_Epithelial_Downstream/scVelo_02_scvelo.py
"""

import os
import sys
import re
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad
import matplotlib as mpl

# ============================== Config =======================================
PROJECT = "/home/MMB/projects/Human CRPC scRNAseq"
SAMPLES = ["CRPC1", "CRPC2", "CRPC3"]
LOOM_DIR = f"{PROJECT}/Results/05_Epithelial_Downstream/scVelo/loom"
OUT_DIR = f"{PROJECT}/Results/05_Epithelial_Downstream/scVelo"
FIG_DIR = f"{OUT_DIR}/figures"
META_CSV = f"{OUT_DIR}/epi_metadata_for_scvelo.csv"
OUT_H5AD = f"{OUT_DIR}/epi_scvelo.h5ad"
# Re-derived 2026-06-25 after the deterministic (seeded) re-baseline of stages
# 02/04 — must match Epithelial_Annotation.R::cluster_to_label. Re-export the
# seurat object (scVelo_00) before re-running so cluster IDs are consistent.
CLUSTER_LABEL_MAP = {
    "0": "OE 1",
    "1": "OE 2",
    "2": "OE 3",
    "3": "Club-like",
    "4": "Hillock-like 1",
    "5": "Hillock-like 2",
    "6": "OE 4",
    "7": "BE 1",
    "8": "BE 2",
    "9": "ARPC",
    "10": "Ionocyte-like",
}
CLUSTER_LABEL_ORDER = [
    "ARPC", "Club-like", "Hillock-like 1", "Hillock-like 2", "BE 1", "BE 2",
    "OE 1", "OE 2", "OE 3", "OE 4", "Ionocyte-like",
]

os.makedirs(FIG_DIR, exist_ok=True)
scv.settings.figdir = FIG_DIR
sc.settings.figdir = FIG_DIR
scv.settings.set_figure_params("scvelo", dpi_save=200)
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"
mpl.rcParams["savefig.edgecolor"] = "white"
mpl.rcParams["savefig.transparent"] = False

# Velocyto loom barcode format: "{sample}:{barcode}x"  → Seurat 형식 "P{i}_{barcode}-1"
def loom_to_seurat_bc(loom_bcs: np.ndarray, sample_idx: int) -> np.ndarray:
    out = []
    for bc in loom_bcs:
        m = re.match(r"^[^:]+:([ACGT]+)x?$", bc)
        if not m:
            out.append(bc)
            continue
        out.append(f"P{sample_idx}_{m.group(1)}-1")
    return np.array(out)


# ============================== Load looms ===================================
print("=== Loading looms ===")
adatas = []
for i, s in enumerate(SAMPLES, start=1):
    loom_path = f"{LOOM_DIR}/{s}.loom"
    if not os.path.exists(loom_path):
        sys.exit(f"Missing loom: {loom_path} — run scVelo_01_velocyto.sh first.")
    a = ad.read_loom(loom_path, sparse=True, X_name="spliced")
    a.var_names_make_unique()
    a.obs_names = loom_to_seurat_bc(a.obs_names.to_numpy(), i)
    a.obs_names_make_unique()
    a.obs["orig.ident"] = s
    print(f"  {s}: {a.n_obs} cells × {a.n_vars} genes")
    adatas.append(a)

adata = ad.concat(adatas, join="outer", index_unique=None)
print(f"Merged: {adata.n_obs} cells × {adata.n_vars} genes")

# ============================== Join Seurat metadata =========================
print("=== Joining Seurat metadata + UMAP ===")
if not os.path.exists(META_CSV):
    sys.exit(f"Missing metadata CSV: {META_CSV} — run scVelo_00_export_seurat.R first.")

meta = pd.read_csv(META_CSV, index_col="cell")
common = adata.obs_names.intersection(meta.index)
n_loom, n_meta, n_common = adata.n_obs, len(meta), len(common)
print(f"  loom={n_loom}  seurat={n_meta}  overlap={n_common}")
if n_common < 0.5 * min(n_loom, n_meta):
    print("WARNING: <50% barcode overlap. Check loom_to_seurat_bc() mapping.")

adata = adata[common].copy()
for col in meta.columns:
    adata.obs[col] = meta.loc[adata.obs_names, col].values
adata.obs["epi_label"] = (
    adata.obs["seurat_clusters"].astype(str).map(CLUSTER_LABEL_MAP)
)
adata.obs["epi_label"] = pd.Categorical(
    adata.obs["epi_label"],
    categories=[x for x in CLUSTER_LABEL_ORDER if x in set(adata.obs["epi_label"])],
    ordered=True,
)

# Restore Seurat UMAP onto adata
adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)

# ============================== Preprocessing + dynamical model ==============
print("=== scVelo preprocess ===")
# scvelo 0.3.4 + scanpy 1.11.5: filter_and_normalize(..., n_top_genes=2000) 는
# 내부에서 normalize_per_cell(n_top_genes=...) 를 호출하는데 신버전 scanpy 의
# normalize_per_cell 시그니처에 n_top_genes 가 빠져 TypeError 발생.
# scvelo 권장대로 단계 분리.
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="seurat", subset=True
)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

print("=== Recover dynamics (slow) ===")
scv.tl.recover_dynamics(adata, n_jobs=20)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=20)

# ============================== Plots ========================================
print("=== Plots ===")
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="epi_label",
    save="velocity_stream_by_annotation.png", show=False
)
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="orig.ident",
    save="velocity_stream_by_patient.png", show=False
)
scv.pl.velocity_embedding(
    adata, basis="umap", color="epi_label", arrow_length=3,
    save="velocity_arrows_by_annotation.png", show=False
)

scv.tl.recover_latent_time(adata)
scv.pl.scatter(
    adata, color="latent_time", color_map="gnuplot",
    save="latent_time.png", show=False
)

# top driver genes
scv.tl.rank_velocity_genes(adata, groupby="seurat_clusters", min_corr=0.3)
df = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"]).head(30)
df.to_csv(f"{OUT_DIR}/top_velocity_genes_per_cluster.csv", index=False)

# ============================== Save ========================================
adata.write_h5ad(OUT_H5AD, compression="gzip")
print(f"=== Saved: {OUT_H5AD} ===")
