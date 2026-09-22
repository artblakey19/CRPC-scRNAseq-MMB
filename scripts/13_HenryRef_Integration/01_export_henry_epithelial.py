#!/usr/bin/env python3
# =============================================================================
# 13 Henry-ref integration — step 01 : Henry 2018(GSE120716) 정상 상피 전체 export
# -----------------------------------------------------------------------------
# CELLxGENE 처리본 h5ad(raw.X = 정수 UMI)에서 원 논문 라벨(obs$Population:
# BE/Club/Hillock/LE)을 그대로 들고 나와 R(02)에서 우리 상피와 통합.
#   - 09/07 과 달리 Hillock(CELLxGENE 에선 'epithelial cell of urethra')을 포함한다.
#   - BE 18k 가 kNN/label transfer 를 지배하지 않게 5000 으로 subsample(seed 42).
#   - NE(25개)는 너무 적어 제외.
# =============================================================================
import anndata as ad, numpy as np, scipy.io as sio, scipy.sparse as sp, gzip, os

H5   = "Raw_data/Henry_GSE120716/henry_epithelial.h5ad"
OUT  = "Raw_data/Henry_GSE120716/epi_matrix_all"
os.makedirs(OUT, exist_ok=True)
rng  = np.random.default_rng(42)
CAP  = {"BE": 5000, "Club": None, "Hillock": None, "LE": None}

a = ad.read_h5ad(H5)
counts = a.raw.X                              # cells x genes, 정수 UMI
genes  = a.raw.var                            # feature_name = symbol, index = Ensembl
assert sp.issparse(counts)

pop = a.obs["Population"].astype(str).values
idx_all = []
for lab, cap in CAP.items():
    ii = np.where(pop == lab)[0]
    if cap is not None and len(ii) > cap:
        ii = np.sort(rng.choice(ii, cap, replace=False))
    idx_all.append(ii)
    print(f"  {lab:8s} n={len(ii)}")
idx_all = np.sort(np.concatenate(idx_all))

M  = counts[idx_all, :].T.tocsc()            # genes x cells
bc = np.array([f"HENRY_{b}" for b in a.obs_names[idx_all]])   # 충돌 방지 접두사
print(f"matrix: {M.shape[0]} genes x {M.shape[1]} cells")

sio.mmwrite(os.path.join(OUT, "matrix.mtx"), M, field="integer")
os.system(f"gzip -f {os.path.join(OUT, 'matrix.mtx')}")
with gzip.open(os.path.join(OUT, "features.tsv.gz"), "wt") as f:
    for ens, sym in zip(genes.index, genes["feature_name"].astype(str)):
        f.write(f"{ens}\t{sym}\n")
with gzip.open(os.path.join(OUT, "barcodes.tsv.gz"), "wt") as f:
    f.write("\n".join(bc) + "\n")

meta = a.obs.iloc[idx_all][["Population", "donor_id", "Sample", "Region", "tissue"]].copy()
meta.insert(0, "cell", bc)
meta.to_csv(os.path.join(OUT, "cell_meta.csv"), index=False)
print("[01] done ->", OUT)
