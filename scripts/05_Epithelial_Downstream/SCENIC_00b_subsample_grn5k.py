#!/usr/bin/env python3
"""SCENIC step 0b — create the 5000-cell GRN subsample loom.

pySCENIC GRNBoost2 (SCENIC_01 step 1) runs on a 5000-cell subsample to stay
under dask's serialized message-size limit. SCENIC_00 only writes the full
epi_counts.loom, so this script derives epi_counts_grn5k.loom from it.

Seeded (np.random.seed=777) for reproducibility across reseeded reruns.

Run: .venv/bin/python3 scripts/05_Epithelial_Downstream/SCENIC_00b_subsample_grn5k.py
"""
import os
import numpy as np
import loompy

PROJ = "/home/MMB/projects/Human CRPC scRNAseq"
OUT = f"{PROJ}/Results/05_Epithelial_Downstream/SCENIC"
SRC = f"{OUT}/epi_counts.loom"
DST = f"{OUT}/epi_counts_grn5k.loom"
N = 5000
SEED = 777

assert os.path.exists(SRC), f"missing {SRC} — run SCENIC_00_export_loom.R first"

np.random.seed(SEED)
with loompy.connect(SRC, mode="r") as ds:
    n_genes, n_cells = ds.shape
    k = min(N, n_cells)
    idx = np.sort(np.random.choice(n_cells, size=k, replace=False))
    mat = ds[:, idx]
    ra = {key: ds.ra[key] for key in ds.ra.keys()}
    ca = {key: ds.ca[key][idx] for key in ds.ca.keys()}
    print(f"source: {n_genes} genes x {n_cells} cells -> subsample {k} cells (seed {SEED})")

if os.path.exists(DST):
    os.remove(DST)
loompy.create(DST, mat, ra, ca)
print(f"wrote {DST}")
