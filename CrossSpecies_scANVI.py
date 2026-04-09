# Cross-Species Integration with scANVI
# Semi-supervised: Mouse BE/LE labels guide human cell mapping

import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sc.settings.figdir = "Results/Cross_Species"
scvi.settings.seed = 42

# ============================================================
# SKIP_TRAINING: Set to True to load saved h5ad and skip model training.
#                Set to False when you modify training parameters above.
# ============================================================
SKIP_TRAINING = True

if not SKIP_TRAINING:
    # --- Data loading & preprocessing ---
    from scipy.io import mmread

    counts = mmread("Results/Cross_Species/scanvi_counts.mtx").T.tocsr()
    genes = pd.read_csv("Results/Cross_Species/scanvi_genes.txt", header=None)[0].values
    barcodes = pd.read_csv("Results/Cross_Species/scanvi_barcodes.txt", header=None)[0].values
    metadata = pd.read_csv("Results/Cross_Species/scanvi_metadata.csv", index_col=0)

    adata = sc.AnnData(X=counts, obs=metadata)
    adata.var_names = genes
    adata.obs_names = barcodes

    # Basic filtering
    sc.pp.filter_genes(adata, min_cells=10)

    # Use all mouse EpiCellTypes as scANVI labels, Human cells remain "Unknown"
    adata.obs["scanvi_label"] = adata.obs["scanvi_label"].astype(str)

    # HVG selection on raw counts
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, batch_key="species", subset=True
    )
    # Restore raw counts for scVI
    adata.X = adata.layers["counts"].copy()

    # --- Model training ---
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="species",
        categorical_covariate_keys=["sample"],
    )

    vae = scvi.model.SCVI(
        adata,
        n_latent=30,
        n_layers=3,
        gene_likelihood="nb",
    )
    vae.train(max_epochs=200, early_stopping=True)

    # Train scANVI with mouse labels
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        unlabeled_category="Unknown",
        labels_key="scanvi_label",
    )
    scanvi_model.train(max_epochs=200, n_samples_per_label=200)

    # Get latent representation
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

    # Predict labels and get prediction probabilities
    adata.obs["scanvi_prediction"] = scanvi_model.predict()
    pred_probs = scanvi_model.predict(soft=True)
    adata.obs["scanvi_confidence"] = pred_probs.max(axis=1).values

    # Mark low-confidence predictions
    confidence_threshold = 0.5
    adata.obs["celltype_final"] = adata.obs["scanvi_prediction"].copy()
    low_conf = adata.obs["scanvi_confidence"] < confidence_threshold
    adata.obs.loc[low_conf, "celltype_final"] = "Uncertain"

    # For mouse cells, keep original labels
    mouse_mask = adata.obs["species"] == "Mouse"
    human_mask = adata.obs["species"] == "Human"
    adata.obs.loc[mouse_mask, "celltype_final"] = adata.obs.loc[mouse_mask, "scanvi_label"]

    # Save trained model and data
    adata.write("Results/Cross_Species/combined_scANVI.h5ad")

else:
    # --- Load previously saved results ---
    adata = sc.read("Results/Cross_Species/combined_scANVI.h5ad")
    mouse_mask = adata.obs["species"] == "Mouse"
    human_mask = adata.obs["species"] == "Human"
    confidence_threshold = 0.5
    low_conf = adata.obs["scanvi_confidence"] < confidence_threshold

# --- UMAP directly from scANVI latent space (no Harmony) ---
sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30, key_added="scANVI_nn")
sc.tl.umap(adata, min_dist=0.2, spread=0.8, neighbors_key="scANVI_nn")
adata.obsm["X_umap_scANVI"] = adata.obsm["X_umap"].copy()

# --- Harmony on scANVI latent space (visualization only — predictions already finalized) ---
import harmonypy as hm

harmony_out = hm.run_harmony(adata.obsm["X_scANVI"], adata.obs, "species")
adata.obsm["X_scANVI_harmony"] = harmony_out.Z_corr

# UMAP on harmony-corrected latent space
sc.pp.neighbors(adata, use_rep="X_scANVI_harmony", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.2, spread=0.8)
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()

# Visualizations ----
fig_kw = dict(dpi=300, bbox_inches="tight")

# --- scANVI-direct UMAP (no Harmony) ---
adata.obsm["X_umap"] = adata.obsm["X_umap_scANVI"]

fig, axes = plt.subplots(1, 3, figsize=(24, 7))
sc.pl.umap(adata, color="species", ax=axes[0], title="Species", show=False)
sc.pl.umap(adata, color="celltype_final", ax=axes[1], title="Cell Type (scANVI)", show=False)
sc.pl.umap(adata, color="sample", ax=axes[2], title="Sample", show=False)
fig.suptitle("scANVI Latent Space (no Harmony)", fontsize=14, y=1.02)
plt.savefig("Results/Cross_Species/UMAP_scANVI_noHarmony_overview.png", **fig_kw)
plt.close()

# --- Harmony-corrected UMAP ---
adata.obsm["X_umap"] = adata.obsm["X_umap_harmony"]

# Comparison: scANVI vs scANVI+Harmony side by side
fig, axes = plt.subplots(2, 3, figsize=(24, 14))
for col, (color, label) in enumerate(
    [("species", "Species"), ("celltype_final", "Cell Type"), ("sample", "Sample")]
):
    adata.obsm["X_umap"] = adata.obsm["X_umap_scANVI"]
    sc.pl.umap(adata, color=color, ax=axes[0, col], title=f"{label} (scANVI only)", show=False)
    adata.obsm["X_umap"] = adata.obsm["X_umap_harmony"]
    sc.pl.umap(adata, color=color, ax=axes[1, col], title=f"{label} (scANVI + Harmony)", show=False)
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_scANVI_vs_Harmony_comparison.png", **fig_kw)
plt.close()

# Overview (Harmony): species, predicted cell types, samples
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
sc.pl.umap(adata, color="species", ax=axes[0], title="Species", show=False)
sc.pl.umap(adata, color="celltype_final", ax=axes[1], title="Cell Type (scANVI)", show=False)
sc.pl.umap(adata, color="sample", ax=axes[2], title="Sample", show=False)
plt.savefig("Results/Cross_Species/UMAP_scANVI_overview.png", **fig_kw)
plt.close()

# Split by sample
samples = ["CRPC1", "CRPC2", "CRPC3", "DoubleTg", "TripleTg"]
fig, axes = plt.subplots(1, 5, figsize=(28, 6))
for ax, s in zip(axes, samples):
    mask = adata.obs["sample"] == s
    sc.pl.umap(adata, ax=ax, show=False, na_color="lightgray", size=5)
    sc.pl.umap(
        adata[mask],
        color="celltype_final",
        ax=ax,
        title=s,
        legend_loc="none",
        size=10,
        show=False,
    )
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_scANVI_split_by_sample.png", **fig_kw)
plt.close()

# Prediction confidence UMAP
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata[human_mask], color="scanvi_confidence", ax=axes[0],
           title="Human: Prediction Confidence", vmin=0, vmax=1, cmap="RdYlGn", show=False)
sc.pl.umap(adata[human_mask], color="celltype_final", ax=axes[1],
           title="Human: Predicted EpiCellTypes", show=False)
plt.savefig("Results/Cross_Species/UMAP_scANVI_human_confidence.png", **fig_kw)
plt.close()

# Mouse reference vs Human predicted (side by side)
fig, axes = plt.subplots(1, 2, figsize=(20, 8))
sc.pl.umap(adata[mouse_mask], color="scanvi_label", ax=axes[0],
           title="Mouse EpiCellTypes (Reference)", show=False)
sc.pl.umap(adata[human_mask], color="scanvi_prediction", ax=axes[1],
           title="Human EpiCellTypes (Predicted)", show=False)
plt.savefig("Results/Cross_Species/UMAP_scANVI_celltype_prediction.png", **fig_kw)
plt.close()

# Prediction distribution per human sample (bar plot)
human_types = sorted(adata.obs.loc[mouse_mask, "scanvi_label"].unique())
fig, axes = plt.subplots(1, 3, figsize=(24, 6))
for ax, p in zip(axes, ["CRPC1", "CRPC2", "CRPC3"]):
    mask = (adata.obs["sample"] == p) & human_mask
    ct_counts = adata.obs.loc[mask, "scanvi_prediction"].value_counts()
    ct_counts = ct_counts.reindex(human_types, fill_value=0)
    ct_counts.plot.bar(ax=ax)
    ax.set_title(f"{p}: Predicted EpiCellType")
    ax.set_ylabel("Number of cells")
    ax.tick_params(axis="x", rotation=45)
plt.tight_layout()
plt.savefig("Results/Cross_Species/scANVI_prediction_barplot.png", **fig_kw)
plt.close()

# Confidence summary
human_data = adata.obs[human_mask]
print("\n=== Prediction Confidence Summary (Human cells) ===")
print(f"Total human cells: {len(human_data)}")
print(f"High confidence (>={confidence_threshold}): {(~low_conf[human_mask]).sum()}")
print(f"Low confidence (<{confidence_threshold}): {low_conf[human_mask].sum()}")
print(f"\nMean confidence per predicted type:")
print(human_data.groupby("scanvi_prediction")["scanvi_confidence"].agg(["mean", "count"])
      .sort_values("count", ascending=False).round(3))

# Prediction composition table
pred_comp = pd.crosstab(adata.obs.loc[human_mask, "sample"],
                        adata.obs.loc[human_mask, "scanvi_prediction"])
pred_comp.to_csv("Results/Cross_Species/scANVI_prediction_by_sample.csv")

print("\n=== scANVI Prediction by Sample ===")
print(pred_comp)

