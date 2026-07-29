#!/usr/bin/env python3
# =============================================================================
# 09 inferCNV — step 07 : Henry 2018(GSE120716) 정상 전립선 상피 reference 추출
# -----------------------------------------------------------------------------
# CELLxGENE 처리본 h5ad(raw.X = 정수 UMI) 에서 정상 상피(basal/luminal/secretory)를
# inferCNV 외부 reference 로 추출 → MTX 로 내보내 R(08)에서 우리 상피와 병합.
#   - basal 18k 이 평균을 지배하지 않게 3000 으로 subsample, luminal/secretory 전량.
#   - urethra(다른 조직)·neuroendocrine(25개) 제외.
#   - group 라벨(Henry_basal/luminal/secretory)로 inferCNV ref_group 3개 구성.
# 필드 표준(Kfoury 2021 CancerCell·Dong 2020 CommBiol 이 동일 Henry atlas 사용).
# =============================================================================
import anndata as ad, numpy as np, scipy.io as sio, scipy.sparse as sp, gzip, os

H5   = "Raw_data/Henry_GSE120716/henry_epithelial.h5ad"
OUT  = "Raw_data/Henry_GSE120716/ref_matrix"
os.makedirs(OUT, exist_ok=True)
rng  = np.random.default_rng(42)
KEEP = {"basal cell of prostate epithelium": ("Henry_basal", 3000),
        "luminal cell of prostate epithelium": ("Henry_luminal", None),
        "secretory cell": ("Henry_secretory", None)}

a = ad.read_h5ad(H5)
counts = a.raw.X                              # cells x genes, 정수 UMI
genes  = a.raw.var                            # feature_name = symbol, index = Ensembl
assert sp.issparse(counts)

# --- cell type 별 subsample -------------------------------------------------
idx_all, groups = [], []
ct = a.obs["cell_type"].astype(str).values
for lab, (grp, cap) in KEEP.items():
    ii = np.where(ct == lab)[0]
    if cap is not None and len(ii) > cap:
        ii = np.sort(rng.choice(ii, cap, replace=False))
    idx_all.append(ii); groups += [grp] * len(ii)
    print(f"  {lab:40s} -> {grp:16s} n={len(ii)}")
idx_all = np.concatenate(idx_all)
groups  = np.array(groups)
bc      = np.array([f"HENRY_{b}" for b in a.obs_names[idx_all]])   # 충돌 방지 접두사

M = counts[idx_all, :].T.tocsc()             # genes x cells
print(f"reference matrix: {M.shape[0]} genes x {M.shape[1]} cells")

# --- export (genes x cells MTX + features + barcodes + groups) ---------------
sio.mmwrite(os.path.join(OUT, "matrix.mtx"), M, field="integer")
with gzip.open(os.path.join(OUT, "features.tsv.gz"), "wt") as f:
    for ens, sym in zip(genes.index, genes["feature_name"].astype(str)):
        f.write(f"{ens}\t{sym}\n")
with gzip.open(os.path.join(OUT, "barcodes.tsv.gz"), "wt") as f:
    f.write("\n".join(bc) + "\n")
with open(os.path.join(OUT, "cell_groups.csv"), "w") as f:
    f.write("cell,group\n")
    for b, g in zip(bc, groups):
        f.write(f"{b},{g}\n")
os.system(f"gzip -f {os.path.join(OUT, 'matrix.mtx')}")
print("[07] done ->", OUT)
print("group counts:", {g: int((groups == g).sum()) for g in set(groups)})
