# Cross-Species Integration with sysVI
# Mouse DNPC Model (DoubleTg/TripleTg) × Human CRPC Biopsy (CRPC1/2/3)
#
# Hypothesis:
#   CRPC1 (mostly normal/benign-PIN) → clusters with DoubleTg (PIN)
#   CRPC2, CRPC3 (advanced cancer)   → clusters with TripleTg (DNPC-like)
#
# Pipeline:
#   Step 1: sysVI integration (cross-species batch correction)
#   Step 2: KNN-based label transfer (Mouse EpiCellTypes → Human)

import os
import scanpy as sc
import scvi
from scvi.external import SysVI
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.io import mmread
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report

scvi.settings.seed = 42
sc.settings.figdir = "Results/Cross_Species"
os.makedirs("Results/Cross_Species", exist_ok=True)

# ============================================================
# Configuration
# ============================================================
SKIP_TRAINING = False          # Set True after first successful run
N_TOP_GENES = 2000             # HVGs for integration
N_LATENT = 30                  # Latent dimensions
MAX_EPOCHS = 200               # Training epochs
CYCLE_WEIGHT = 5               # Cycle-consistency loss weight (2-10 typical; up to 50)
KNN_K = 30                     # K for label transfer
CONFIDENCE_THRESHOLD = 0.7     # Minimum KNN voting fraction to accept label
RANDOM_SEEDS = [42, 0, 1]      # Multiple seeds for robustness (pick best)
BEST_SEED_IDX = 0              # Index into RANDOM_SEEDS to use (update after comparison)

fig_kw = dict(dpi=300, bbox_inches="tight")

# ============================================================
# Step 0: Load data exported from R (Seurat)
# ============================================================
print("=" * 60)
print("Loading data from R exports...")
print("=" * 60)

counts = mmread("Results/Cross_Species/scanvi_counts.mtx").T.tocsr()
genes = pd.read_csv("Results/Cross_Species/scanvi_genes.txt", header=None)[0].values
barcodes = pd.read_csv("Results/Cross_Species/scanvi_barcodes.txt", header=None)[0].values
metadata = pd.read_csv("Results/Cross_Species/scanvi_metadata.csv", index_col=0)

adata = sc.AnnData(X=counts, obs=metadata)
adata.var_names = genes
adata.obs_names = barcodes

print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Species: {adata.obs['species'].value_counts().to_dict()}")
print(f"Samples: {adata.obs['sample'].value_counts().to_dict()}")

# Verify mouse labels exist
mouse_mask = adata.obs["species"] == "Mouse"
human_mask = adata.obs["species"] == "Human"
print(f"Mouse EpiCellTypes: {adata.obs.loc[mouse_mask, 'scanvi_label'].nunique()} types")
print(f"  → {sorted(adata.obs.loc[mouse_mask, 'scanvi_label'].unique())}")

# ============================================================
# Step 1: Preprocessing for sysVI
# ============================================================
# sysVI requires: normalized + log-transformed data in .X
# HVGs selected per-system, then intersected

print("\n" + "=" * 60)
print("Preprocessing: per-species HVG selection...")
print("=" * 60)

sc.pp.filter_genes(adata, min_cells=10)

# Store raw counts
adata.layers["counts"] = adata.X.copy()

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Per-species HVG selection → intersection (as recommended by sysVI paper)
hvg_per_species = {}
for sp in ["Human", "Mouse"]:
    sp_adata = adata[adata.obs["species"] == sp].copy()
    sc.pp.highly_variable_genes(sp_adata, n_top_genes=N_TOP_GENES, flavor="seurat_v3",
                                 layer="counts")
    hvg_per_species[sp] = set(sp_adata.var_names[sp_adata.var["highly_variable"]])

shared_hvgs = hvg_per_species["Human"] & hvg_per_species["Mouse"]
print(f"HVGs: Human={len(hvg_per_species['Human'])}, Mouse={len(hvg_per_species['Mouse'])}")
print(f"Shared HVGs (intersection): {len(shared_hvgs)}")

# If intersection is too small, fall back to union-based selection
if len(shared_hvgs) < 1000:
    print("⚠ Intersection < 1000 — using batch-aware HVG selection instead")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, batch_key="species")
    adata = adata[:, adata.var["highly_variable"]].copy()
else:
    adata = adata[:, list(shared_hvgs)].copy()

print(f"Final gene set: {adata.n_vars} genes")

# ============================================================
# Step 2: sysVI Integration
# ============================================================
# batch_key = "species" (the substantial batch effect)
# categorical_covariate_keys = ["sample"] (weaker batch: individual samples)

if not SKIP_TRAINING:
    print("\n" + "=" * 60)
    print("Training sysVI model...")
    print("=" * 60)

    # Encode species as integer for sysVI
    adata.obs["system"] = (adata.obs["species"] == "Human").astype(int)
    # system: 0 = Mouse, 1 = Human

    # Setup anndata
    SysVI.setup_anndata(
        adata=adata,
        batch_key="system",
        categorical_covariate_keys=["sample"],
    )

    # Use best seed
    scvi.settings.seed = RANDOM_SEEDS[BEST_SEED_IDX]

    # Initialize model
    model = SysVI(adata=adata)

    # Train with cycle-consistency
    model.train(
        max_epochs=MAX_EPOCHS,
        check_val_every_n_epoch=1,
        plan_kwargs={"z_distance_cycle_weight": CYCLE_WEIGHT},
    )

    # --- Plot training losses ---
    epochs_detail = MAX_EPOCHS // 2
    losses = ["reconstruction_loss_train", "kl_local_train", "cycle_loss_train"]

    fig, axs = plt.subplots(2, len(losses), figsize=(len(losses) * 4, 5))
    for ax_i, l_train in enumerate(losses):
        l_val = l_train.replace("_train", "_validation")
        l_name = l_train.replace("_train", "")

        l_train_vals = model.trainer.logger.history[l_train].copy()
        l_train_vals.index = l_train_vals.index + 1

        try:
            l_val_vals = model.trainer.logger.history[l_val].copy()
            l_val_vals.index = l_val_vals.index + 1
            has_val = True
        except KeyError:
            has_val = False

        # Full history
        axs[0, ax_i].plot(l_train_vals.index, l_train_vals.values.ravel(),
                          c="tab:blue", alpha=1, label="train")
        if has_val:
            axs[0, ax_i].plot(l_val_vals.index, l_val_vals.values.ravel(),
                              c="tab:orange", alpha=0.5, label="val")
        axs[0, ax_i].set_title(l_name)
        axs[0, ax_i].legend(fontsize=8)

        # Zoom-in (last half)
        axs[1, ax_i].plot(l_train_vals.index[epochs_detail:],
                          l_train_vals.values.ravel()[epochs_detail:],
                          c="tab:blue")
        if has_val:
            axs[1, ax_i].plot(l_val_vals.index[epochs_detail:],
                              l_val_vals.values.ravel()[epochs_detail:],
                              c="tab:orange", alpha=0.5)
        axs[1, ax_i].set_xlabel("Epoch")

    fig.suptitle(f"sysVI Training (seed={RANDOM_SEEDS[BEST_SEED_IDX]}, "
                 f"cycle_weight={CYCLE_WEIGHT})", fontsize=12)
    fig.tight_layout()
    plt.savefig("Results/Cross_Species/sysVI_training_loss.png", **fig_kw)
    plt.close()
    print("✓ Training loss plot saved")

    # --- Get latent representation ---
    adata.obsm["X_sysVI"] = model.get_latent_representation(adata=adata)

    # Save
    adata.write("Results/Cross_Species/combined_sysVI.h5ad")
    print("✓ Model trained and saved")

else:
    print("\n" + "=" * 60)
    print("Loading pre-trained sysVI results...")
    print("=" * 60)
    adata = sc.read("Results/Cross_Species/combined_sysVI.h5ad")
    mouse_mask = adata.obs["species"] == "Mouse"
    human_mask = adata.obs["species"] == "Human"

# ============================================================
# Step 3: UMAP on sysVI latent space
# ============================================================
print("\n" + "=" * 60)
print("Computing UMAP on sysVI latent space...")
print("=" * 60)

sc.pp.neighbors(adata, use_rep="X_sysVI", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# ============================================================
# Step 4: KNN Label Transfer (Mouse → Human)
# ============================================================
print("\n" + "=" * 60)
print("KNN Label Transfer: Mouse EpiCellTypes → Human cells...")
print("=" * 60)

# Mouse cells = training set (labeled), Human cells = prediction set
mouse_labels = adata.obs.loc[mouse_mask, "scanvi_label"].values
mouse_latent = adata[mouse_mask].obsm["X_sysVI"]
human_latent = adata[human_mask].obsm["X_sysVI"]

# Weighted KNN classifier
knn = KNeighborsClassifier(n_neighbors=KNN_K, weights="distance", metric="euclidean")
knn.fit(mouse_latent, mouse_labels)

# Predict labels
human_predictions = knn.predict(human_latent)

# Prediction probabilities (voting fractions)
human_proba = knn.predict_proba(human_latent)
human_confidence = human_proba.max(axis=1)

# Store results
adata.obs["knn_prediction"] = "reference"
adata.obs.loc[human_mask, "knn_prediction"] = human_predictions

adata.obs["knn_confidence"] = 1.0
adata.obs.loc[human_mask, "knn_confidence"] = human_confidence

# Apply confidence threshold
adata.obs["celltype_final"] = adata.obs["knn_prediction"].copy()
low_conf = adata.obs["knn_confidence"] < CONFIDENCE_THRESHOLD
adata.obs.loc[human_mask & low_conf, "celltype_final"] = "Uncertain"

# For mouse cells, keep original EpiCellTypes
adata.obs.loc[mouse_mask, "celltype_final"] = adata.obs.loc[mouse_mask, "scanvi_label"]

print(f"\nLabel Transfer Results (K={KNN_K}, threshold={CONFIDENCE_THRESHOLD}):")
print(f"  Total human cells: {human_mask.sum()}")
print(f"  High confidence (≥{CONFIDENCE_THRESHOLD}): {(~low_conf[human_mask]).sum()}")
print(f"  Low confidence (Uncertain): {(low_conf[human_mask]).sum()}")

# ============================================================
# Step 5: Visualizations
# ============================================================
print("\n" + "=" * 60)
print("Generating visualizations...")
print("=" * 60)

species_cols = {"Human": "#E64B35", "Mouse": "#4DBBD5"}
sample_cols = {
    "CRPC1": "#E64B35", "CRPC2": "#F39B7F", "CRPC3": "#B09C85",
    "DoubleTg": "#4DBBD5", "TripleTg": "#3C5488",
}

# --- 5a. Overview: Species / Cell Type / Sample ---
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
sc.pl.umap(adata, color="species", palette=species_cols,
           ax=axes[0], title="Species", show=False)
sc.pl.umap(adata, color="celltype_final",
           ax=axes[1], title="Cell Type (sysVI + KNN)", show=False)
sc.pl.umap(adata, color="sample", palette=sample_cols,
           ax=axes[2], title="Sample", show=False)
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_sysVI_overview.png", **fig_kw)
plt.close()
print("✓ Overview UMAP saved")

# --- 5b. Species side-by-side with cell types ---
fig, axes = plt.subplots(1, 2, figsize=(20, 8))
for ax_i, (sp, title) in enumerate([
    ("Mouse", "Mouse EpiCellTypes (Reference)"),
    ("Human", "Human EpiCellTypes (KNN Predicted)")
]):
    mask = adata.obs["species"] == sp
    # Background: all cells grey
    sc.pl.umap(adata, ax=axes[ax_i], show=False, na_color="lightgray", size=3)
    color_key = "scanvi_label" if sp == "Mouse" else "knn_prediction"
    sc.pl.umap(adata[mask], color=color_key, ax=axes[ax_i],
               title=title, size=10, show=False)
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_sysVI_celltype_by_species.png", **fig_kw)
plt.close()
print("✓ Species side-by-side UMAP saved")

# --- 5c. Split by sample (highlight each on grey background) ---
samples = ["CRPC1", "CRPC2", "CRPC3", "DoubleTg", "TripleTg"]
fig, axes = plt.subplots(1, 5, figsize=(30, 6))
for ax, s in zip(axes, samples):
    mask = adata.obs["sample"] == s
    sc.pl.umap(adata, ax=ax, show=False, na_color="lightgray", size=3)
    sc.pl.umap(adata[mask], color="celltype_final", ax=ax,
               title=s, legend_loc="none", size=10, show=False)
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_sysVI_split_by_sample.png", **fig_kw)
plt.close()
print("✓ Split-by-sample UMAP saved")

# --- 5d. Human prediction confidence ---
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata[human_mask], color="knn_confidence", ax=axes[0],
           title="Human: KNN Confidence", vmin=0, vmax=1, cmap="RdYlGn", show=False)
sc.pl.umap(adata[human_mask], color="celltype_final", ax=axes[1],
           title="Human: Predicted EpiCellTypes", show=False)
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_sysVI_human_confidence.png", **fig_kw)
plt.close()
print("✓ Human confidence UMAP saved")

# --- 5e. Prediction distribution per human sample (bar plot) ---
epi_types = sorted(adata.obs.loc[mouse_mask, "scanvi_label"].unique())
fig, axes = plt.subplots(1, 3, figsize=(24, 6))
for ax, patient in zip(axes, ["CRPC1", "CRPC2", "CRPC3"]):
    mask = (adata.obs["sample"] == patient) & human_mask
    ct_counts = adata.obs.loc[mask, "knn_prediction"].value_counts()
    ct_counts = ct_counts.reindex(epi_types, fill_value=0)
    ct_counts.plot.bar(ax=ax, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.set_title(f"{patient}: Predicted EpiCellType Distribution")
    ax.set_ylabel("Number of cells")
    ax.tick_params(axis="x", rotation=45)
plt.tight_layout()
plt.savefig("Results/Cross_Species/sysVI_prediction_barplot.png", **fig_kw)
plt.close()
print("✓ Prediction bar plot saved")

# ============================================================
# Step 6: Quantitative Summary
# ============================================================
print("\n" + "=" * 60)
print("=== KNN Label Transfer Summary (Human cells) ===")
print("=" * 60)

human_data = adata.obs[human_mask]
print(f"\nTotal human cells: {len(human_data)}")
print(f"High confidence (≥{CONFIDENCE_THRESHOLD}): {(human_data['knn_confidence'] >= CONFIDENCE_THRESHOLD).sum()}")
print(f"Low confidence (Uncertain): {(human_data['knn_confidence'] < CONFIDENCE_THRESHOLD).sum()}")

print(f"\nMean confidence per predicted type:")
conf_summary = (human_data
                .groupby("knn_prediction")["knn_confidence"]
                .agg(["mean", "median", "count"])
                .sort_values("count", ascending=False)
                .round(3))
print(conf_summary)

# Prediction composition table (sample × celltype)
pred_comp = pd.crosstab(
    adata.obs.loc[human_mask, "sample"],
    adata.obs.loc[human_mask, "knn_prediction"]
)
pred_comp.to_csv("Results/Cross_Species/sysVI_knn_prediction_by_sample.csv")

print(f"\n=== Prediction Composition by Sample ===")
print(pred_comp)

# Proportion table (normalized)
pred_prop = pred_comp.div(pred_comp.sum(axis=1), axis=0).round(3)
pred_prop.to_csv("Results/Cross_Species/sysVI_knn_prediction_proportion.csv")
print(f"\n=== Prediction Proportions by Sample ===")
print(pred_prop)

# ============================================================
# Step 7: Hypothesis Test — Co-clustering Analysis
# ============================================================
print("\n" + "=" * 60)
print("=== Co-clustering Analysis (Hypothesis Test) ===")
print("=" * 60)

# Leiden clustering on sysVI latent space
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_sysVI")

# Cross-tab: cluster × sample
cluster_sample = pd.crosstab(
    adata.obs["leiden_sysVI"],
    adata.obs["sample"],
    normalize="index"  # row-normalized → fraction of each sample per cluster
).round(3)
cluster_sample.to_csv("Results/Cross_Species/sysVI_cluster_sample_composition.csv")

print("\nCluster composition (fraction of each sample per cluster):")
print(cluster_sample)

# Identify clusters enriched for specific sample combinations
print("\n--- Hypothesis evaluation ---")
for cl in cluster_sample.index:
    row = cluster_sample.loc[cl]
    doubleTg = row.get("DoubleTg", 0)
    tripleTg = row.get("TripleTg", 0)
    crpc1 = row.get("CRPC1", 0)
    crpc2 = row.get("CRPC2", 0)
    crpc3 = row.get("CRPC3", 0)

    # Flag clusters where hypothesis samples co-occur
    if crpc1 > 0.1 and doubleTg > 0.1:
        print(f"  Cluster {cl}: CRPC1={crpc1:.1%} + DoubleTg={doubleTg:.1%} ← PIN-like co-clustering")
    if (crpc2 > 0.1 or crpc3 > 0.1) and tripleTg > 0.1:
        print(f"  Cluster {cl}: CRPC2={crpc2:.1%}, CRPC3={crpc3:.1%} + TripleTg={tripleTg:.1%} ← DNPC-like co-clustering")

# Cluster UMAP
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color="leiden_sysVI", ax=axes[0],
           title="sysVI Leiden Clusters", legend_loc="on data", show=False)
sc.pl.umap(adata, color="sample", palette=sample_cols, ax=axes[1],
           title="Samples", show=False)
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_sysVI_leiden_clusters.png", **fig_kw)
plt.close()
print("✓ Leiden cluster UMAP saved")

# Save final object
adata.write("Results/Cross_Species/combined_sysVI.h5ad")
print("\n✓ Final h5ad saved: Results/Cross_Species/combined_sysVI.h5ad")
print("Done!")