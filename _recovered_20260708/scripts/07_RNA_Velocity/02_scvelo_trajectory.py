#!/usr/bin/env python3
"""Run scVelo trajectory analysis on epithelial CRPC cells.

Inputs:
  - Per-sample loom files from velocyto.py.
  - Seurat-exported metadata with annotation and UMAP coordinates.

Outputs are written to Results/07_RNA_Velocity/scvelo by default.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv


DEFAULT_SAMPLE_MAP = {
    "CRPC1": "P1",
    "CRPC2": "P2",
    "CRPC3": "P3",
    "P1": "P1",
    "P2": "P2",
    "P3": "P3",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--loom-dir",
        default="Results/07_RNA_Velocity/loom_by_sample",
        help="Directory containing velocyto loom files.",
    )
    parser.add_argument(
        "--metadata",
        default="Results/07_RNA_Velocity/seurat_velocity_metadata.csv",
        help="CSV from 00_export_seurat_velocity_inputs.R.",
    )
    parser.add_argument(
        "--out-dir",
        default="Results/07_RNA_Velocity/scvelo",
        help="Output directory.",
    )
    parser.add_argument(
        "--mode",
        choices=("deterministic", "stochastic", "dynamical"),
        default="deterministic",
        help="Velocity mode. Use dynamical for latent_time, but it is slower.",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=3000,
        help="Highly variable genes to keep for velocity moments.",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of PCs for neighbors/moments.",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Neighbor count for moments.",
    )
    parser.add_argument(
        "--min-shared-counts",
        type=int,
        default=20,
        help="scVelo gene filter threshold.",
    )
    parser.add_argument(
        "--sample-map",
        nargs="*",
        default=[],
        help="Optional sample mappings like CRPC1=P1 CRPC2=P2 CRPC3=P3.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="Figure DPI.",
    )
    return parser.parse_args()


def build_sample_map(items: list[str]) -> dict[str, str]:
    sample_map = dict(DEFAULT_SAMPLE_MAP)
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --sample-map item: {item!r}")
        left, right = item.split("=", 1)
        sample_map[left] = right
    return sample_map


def discover_looms(loom_dir: Path) -> list[Path]:
    looms = sorted(loom_dir.glob("*.loom"))
    if not looms:
        raise FileNotFoundError(f"No .loom files found in {loom_dir}")
    return looms


def sample_from_path(path: Path, sample_map: dict[str, str]) -> str | None:
    stem = path.stem
    for key in sorted(sample_map, key=len, reverse=True):
        if key in stem:
            return sample_map[key]
    return None


def read_loom(path: Path, sample_map: dict[str, str]) -> ad.AnnData:
    adata = sc.read_loom(str(path), sparse=True, cleanup=False, X_name="spliced")
    fallback_sample = sample_from_path(path, sample_map)
    if fallback_sample is None:
        fallback_sample = path.stem
    adata.obs["velocity_file"] = path.name
    adata.obs["velocity_sample"] = fallback_sample
    adata.var_names_make_unique()
    return adata


def normalize_barcode(barcode: str) -> str:
    bc = barcode.strip()
    bc = re.sub(r"^[A-Za-z0-9_.-]+:", "", bc)
    if "_" in bc:
        maybe_sample, rest = bc.split("_", 1)
        if maybe_sample in DEFAULT_SAMPLE_MAP.values() or maybe_sample in DEFAULT_SAMPLE_MAP:
            bc = rest
    bc = bc.replace("-1x", "-1")
    if bc.endswith("x") and "-" not in bc:
        bc = bc[:-1] + "-1"
    if not re.search(r"-\d+$", bc):
        bc = bc + "-1"
    return bc


def sample_from_cell_id(cell_id: str, fallback: str, sample_map: dict[str, str]) -> str:
    raw = str(cell_id)
    if ":" in raw:
        prefix = raw.split(":", 1)[0]
        if prefix in sample_map:
            return sample_map[prefix]
    if "_" in raw:
        prefix = raw.split("_", 1)[0]
        if prefix in sample_map:
            return sample_map[prefix]
    return sample_map.get(fallback, fallback)


def attach_seurat_ids(adata: ad.AnnData, sample_map: dict[str, str]) -> ad.AnnData:
    original_ids = adata.obs_names.astype(str).to_numpy()
    fallback_samples = adata.obs["velocity_sample"].astype(str).to_numpy()
    samples = [
        sample_from_cell_id(cell_id, fallback, sample_map)
        for cell_id, fallback in zip(original_ids, fallback_samples, strict=False)
    ]
    barcodes = [normalize_barcode(cell_id) for cell_id in original_ids]
    seurat_cells = [f"{sample}_{barcode}" for sample, barcode in zip(samples, barcodes, strict=False)]

    adata.obs["velocity_cell_id"] = original_ids
    adata.obs["sample_prefix"] = samples
    adata.obs["barcode"] = barcodes
    adata.obs["seurat_cell"] = seurat_cells
    adata.obs_names = pd.Index(seurat_cells, name="cell")
    return adata


def load_velocity_data(loom_dir: Path, sample_map: dict[str, str]) -> ad.AnnData:
    looms = discover_looms(loom_dir)
    adatas = [read_loom(path, sample_map) for path in looms]
    if len(adatas) == 1:
        adata = adatas[0]
    else:
        adata = ad.concat(
            adatas,
            join="outer",
            merge="same",
            label="loom_batch",
            keys=[path.stem for path in looms],
            index_unique=None,
        )
    adata = attach_seurat_ids(adata, sample_map)
    adata.var_names_make_unique()
    return adata


def attach_metadata(adata: ad.AnnData, metadata_csv: Path, out_dir: Path) -> ad.AnnData:
    meta = pd.read_csv(metadata_csv)
    if "cell" not in meta.columns:
        raise ValueError(f"Metadata is missing required column 'cell': {metadata_csv}")
    meta = meta.drop_duplicates("cell").set_index("cell")

    matched = adata.obs_names[adata.obs_names.isin(meta.index)]
    report = {
        "velocity_cells": int(adata.n_obs),
        "seurat_metadata_cells": int(meta.shape[0]),
        "matched_cells": int(len(matched)),
        "match_fraction_velocity": float(len(matched) / max(adata.n_obs, 1)),
        "match_fraction_seurat": float(len(matched) / max(meta.shape[0], 1)),
    }
    (out_dir / "barcode_match_report.json").write_text(json.dumps(report, indent=2) + "\n")
    if len(matched) == 0:
        raise ValueError(
            "No cells matched between velocyto loom IDs and Seurat metadata. "
            "Check sample prefixes and barcode suffixes."
        )

    adata = adata[matched].copy()
    meta = meta.loc[adata.obs_names]
    for col in meta.columns:
        adata.obs[col] = meta[col].to_numpy()

    adata.obsm["X_umap"] = meta[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
    for col in ("annotation", "sample", "orig.ident", "seurat_cluster", "sample_prefix"):
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")
    return adata


def write_obs_summaries(adata: ad.AnnData, out_dir: Path) -> None:
    obs_cols = [
        col
        for col in (
            "sample",
            "orig.ident",
            "seurat_cluster",
            "annotation",
            "velocity_pseudotime",
            "latent_time",
            "velocity_length",
            "velocity_confidence",
        )
        if col in adata.obs.columns
    ]
    adata.obs[obs_cols].to_csv(out_dir / "cell_velocity_metrics.csv")

    numeric_cols = [
        col
        for col in (
            "velocity_pseudotime",
            "latent_time",
            "velocity_length",
            "velocity_confidence",
        )
        if col in adata.obs.columns
    ]
    if numeric_cols and "annotation" in adata.obs.columns:
        summary = (
            adata.obs.groupby("annotation", observed=True)[numeric_cols]
            .agg(["count", "mean", "median", "std"])
            .sort_index()
        )
        summary.to_csv(out_dir / "annotation_velocity_summary.csv")

    if "sample" in adata.obs.columns and "annotation" in adata.obs.columns:
        comp = pd.crosstab(adata.obs["annotation"], adata.obs["sample"])
        comp.to_csv(out_dir / "annotation_by_sample_cell_counts.csv")


def save_velocity_gene_tables(adata: ad.AnnData, out_dir: Path) -> None:
    if "rank_velocity_genes" not in adata.uns:
        return
    result = adata.uns["rank_velocity_genes"]
    groups = result.get("names", {}).dtype.names
    if not groups:
        return
    rows = []
    for group in groups:
        names = result["names"][group]
        scores = result["scores"][group] if "scores" in result else [np.nan] * len(names)
        for rank, (gene, score) in enumerate(zip(names, scores, strict=False), start=1):
            rows.append({"annotation": group, "rank": rank, "gene": gene, "score": score})
    pd.DataFrame(rows).to_csv(out_dir / "rank_velocity_genes_by_annotation.csv", index=False)


def plot_outputs(adata: ad.AnnData, out_dir: Path, dpi: int, mode: str) -> None:
    basis = "umap"

    def save_current(name: str) -> None:
        plt.savefig(out_dir / name, dpi=dpi, bbox_inches="tight")
        plt.close("all")

    color = "annotation" if "annotation" in adata.obs.columns else None
    if color is not None:
        scv.pl.velocity_embedding_stream(
            adata,
            basis=basis,
            color=color,
            legend_loc="right margin",
            dpi=dpi,
            show=False,
        )
        save_current("velocity_embedding_stream_by_annotation.png")
        scv.pl.velocity_embedding_grid(
            adata,
            basis=basis,
            color=color,
            legend_loc="right margin",
            dpi=dpi,
            show=False,
        )
        save_current("velocity_embedding_grid_by_annotation.png")

    for metric in ("velocity_pseudotime", "latent_time", "velocity_confidence", "velocity_length"):
        if metric in adata.obs.columns:
            scv.pl.scatter(
                adata,
                basis=basis,
                color=metric,
                color_map="viridis",
                dpi=dpi,
                show=False,
            )
            save_current(f"{metric}.png")

    if "root_cells" in adata.obs.columns or "end_points" in adata.obs.columns:
        terminal_colors = [c for c in ("root_cells", "end_points") if c in adata.obs.columns]
        scv.pl.scatter(
            adata,
            basis=basis,
            color=terminal_colors,
            dpi=dpi,
            show=False,
        )
        save_current("terminal_states.png")

    if mode == "dynamical":
        top_genes = list(adata.var_names[:6])
        if "fit_likelihood" in adata.var.columns:
            top_genes = (
                adata.var["fit_likelihood"]
                .sort_values(ascending=False)
                .head(6)
                .index.tolist()
            )
        if top_genes:
            scv.pl.velocity(
                adata,
                var_names=top_genes,
                color=color,
                ncols=2,
                dpi=dpi,
                show=False,
            )
            save_current("top_dynamic_genes.png")


def main() -> None:
    args = parse_args()
    loom_dir = Path(args.loom_dir)
    metadata_csv = Path(args.metadata)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scv.settings.verbosity = 3
    scv.settings.presenter_view = False
    sc.set_figure_params(frameon=False, dpi=args.dpi)

    sample_map = build_sample_map(args.sample_map)
    adata = load_velocity_data(loom_dir, sample_map)
    adata = attach_metadata(adata, metadata_csv, out_dir)

    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=args.min_shared_counts,
    )
    if args.n_top_genes > 0 and adata.n_vars > args.n_top_genes:
        try:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=args.n_top_genes,
                flavor="seurat",
            )
            keep_hvg = adata.var["highly_variable"].to_numpy()
            if keep_hvg.sum() > 0:
                adata = adata[:, keep_hvg].copy()
        except Exception as exc:
            (out_dir / "highly_variable_genes_warning.txt").write_text(
                str(exc) + "\n"
            )
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=8)
        scv.tl.velocity(adata, mode="dynamical")
    elif args.mode == "deterministic":
        scv.tl.velocity(adata, mode="deterministic")
    else:
        scv.tl.velocity(adata, mode="stochastic")

    scv.tl.velocity_graph(
        adata,
        n_jobs=1,
        backend="threading",
        show_progress_bar=False,
    )
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_confidence(adata)

    try:
        scv.tl.terminal_states(adata)
    except Exception as exc:  # terminal state estimation can fail on sparse graphs
        (out_dir / "terminal_states_warning.txt").write_text(str(exc) + "\n")

    if args.mode == "dynamical":
        scv.tl.latent_time(adata)

    if "annotation" in adata.obs.columns:
        try:
            scv.tl.rank_velocity_genes(adata, groupby="annotation", min_corr=0.3)
        except Exception as exc:
            (out_dir / "rank_velocity_genes_warning.txt").write_text(str(exc) + "\n")

    write_obs_summaries(adata, out_dir)
    save_velocity_gene_tables(adata, out_dir)
    plot_outputs(adata, out_dir, args.dpi, args.mode)

    adata.write_h5ad(out_dir / "crpc_epithelial_scvelo.h5ad")
    print(f"Wrote scVelo outputs to {out_dir}")


if __name__ == "__main__":
    main()
