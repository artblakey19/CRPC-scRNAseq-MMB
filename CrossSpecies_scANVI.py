# Cross-Species Integration with scANVI
# Semi-supervised: Mouse BE/LE labels guide human cell mapping

import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sc.settings.figdir = "Results/Cross_Species"
scvi.settings.seed = 42

# Load data exported from R (mtx + metadata)
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

# HVG selection on raw counts
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000, batch_key="species", subset=True
)
# Restore raw counts for scVI
adata.X = adata.layers["counts"].copy()

# Setup scVI model (base for scANVI)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="species",
    categorical_covariate_keys=["sample"],
)

vae = scvi.model.SCVI(
    adata,
    n_latent=30,
    n_layers=2,
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
scanvi_model.train(max_epochs=50, n_samples_per_label=100)

# Get latent representation
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

# Predict labels for human cells
adata.obs["scanvi_prediction"] = scanvi_model.predict()

# UMAP and clustering on scANVI latent space
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scanvi")

# Save trained model and data
adata.write("Results/Cross_Species/combined_scANVI.h5ad")

# Visualizations ----
fig_kw = dict(dpi=300, bbox_inches="tight")

# Overview: species, clusters, samples
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
sc.pl.umap(adata, color="species", ax=axes[0], title="scANVI: Species", show=False)
sc.pl.umap(
    adata,
    color="leiden_scanvi",
    ax=axes[1],
    title="scANVI: Clusters",
    legend_loc="on data",
    show=False,
)
sc.pl.umap(adata, color="sample", ax=axes[2], title="scANVI: Samples", show=False)
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
        color="leiden_scanvi",
        ax=ax,
        title=s,
        legend_loc="none",
        size=10,
        show=False,
    )
plt.tight_layout()
plt.savefig("Results/Cross_Species/UMAP_scANVI_split_by_sample.png", **fig_kw)
plt.close()

# scANVI predictions: mouse reference + human predicted labels
fig, axes = plt.subplots(1, 2, figsize=(20, 8))
mouse_mask = adata.obs["species"] == "Mouse"
human_mask = adata.obs["species"] == "Human"

sc.pl.umap(
    adata[mouse_mask],
    color="scanvi_label",
    ax=axes[0],
    title="Mouse EpiCellTypes (Reference)",
    show=False,
)
sc.pl.umap(
    adata[human_mask],
    color="scanvi_prediction",
    ax=axes[1],
    title="Human Predicted EpiCellTypes (scANVI)",
    show=False,
)
plt.savefig("Results/Cross_Species/UMAP_scANVI_celltype_prediction.png", **fig_kw)
plt.close()

# Prediction distribution per sample (stacked bar)
fig, axes = plt.subplots(1, 3, figsize=(24, 6))
for ax, p in zip(axes, ["CRPC1", "CRPC2", "CRPC3"]):
    mask = adata.obs["sample"] == p
    counts = adata.obs.loc[mask, "scanvi_prediction"].value_counts()
    counts = counts.reindex(
        [c for c in adata.obs["scanvi_label"].cat.categories if c != "Unknown"],
        fill_value=0,
    )
    counts.plot.bar(ax=ax)
    ax.set_title(f"{p}: Predicted EpiCellType")
    ax.set_ylabel("Number of cells")
    ax.tick_params(axis="x", rotation=45)
plt.tight_layout()
plt.savefig("Results/Cross_Species/scANVI_prediction_barplot.png", **fig_kw)
plt.close()

# Cluster composition table
comp = pd.crosstab(adata.obs["leiden_scanvi"], adata.obs["sample"])
comp.to_csv("Results/Cross_Species/scANVI_cluster_composition.csv")

pred_comp = pd.crosstab(adata.obs["sample"], adata.obs["scanvi_prediction"])
pred_comp.to_csv("Results/Cross_Species/scANVI_prediction_by_sample.csv")

print("\n=== scANVI Prediction Summary ===")
print(pred_comp)
print("\n=== Cluster Composition ===")
print(comp)
