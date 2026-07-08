#!/usr/bin/env python3
"""Format stress-regression reclustering results as a confusion matrix."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUT_DIR = Path(
    "Results/06_QC_Artifact_Diagnostics/"
    "02b_RegressOut_Denisenko17_BaselineControl"
)
IN_CSV = OUT_DIR / "crosstab_regress_vs_annotation.csv"


def read_counts(path: Path) -> tuple[list[str], list[str], list[list[int]]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        cols = header[1:]
        rows: list[str] = []
        counts: list[list[int]] = []
        for line in reader:
            rows.append(line[0])
            counts.append([int(float(x)) for x in line[1:]])
    return rows, cols, counts


def write_matrix(path: Path, rows: list[str], cols: list[str], matrix: list[list[object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + cols)
        for row_name, values in zip(rows, matrix):
            writer.writerow([row_name] + values)


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    rows, cols, counts = read_counts(IN_CSV)

    row_sums = [sum(row) for row in counts]
    col_sums = [sum(counts[i][j] for i in range(len(rows))) for j in range(len(cols))]

    row_pct = [
        [100.0 * value / row_sums[i] if row_sums[i] else 0.0 for value in row]
        for i, row in enumerate(counts)
    ]
    col_pct = [
        [100.0 * counts[i][j] / col_sums[j] if col_sums[j] else 0.0 for j in range(len(cols))]
        for i in range(len(rows))
    ]

    write_matrix(OUT_DIR / "confusion_original_annotation_vs_regress_counts.csv",
                 rows, cols, counts)
    write_matrix(OUT_DIR / "confusion_original_annotation_vs_regress_row_percent.csv",
                 rows, cols, [[f"{x:.1f}" for x in row] for row in row_pct])
    readable = [
        [f"{counts[i][j]} ({row_pct[i][j]:.1f}%)" for j in range(len(cols))]
        for i in range(len(rows))
    ]
    write_matrix(OUT_DIR / "confusion_original_annotation_vs_regress_readable.csv",
                 rows, cols, readable)

    composition_rows: list[dict[str, object]] = []
    col_labels: list[str] = []
    for j, cluster in enumerate(cols):
        sorted_parts = sorted(
            [(counts[i][j], rows[i]) for i in range(len(rows))],
            reverse=True,
        )
        total = col_sums[j]
        top = []
        for n, annotation in sorted_parts[:3]:
            frac = n / total if total else 0.0
            top.append((annotation, n, frac))
        composition_rows.append({
            "regress_cluster": cluster,
            "n_cells": total,
            "top1_annotation": top[0][0],
            "top1_n": top[0][1],
            "top1_fraction": f"{top[0][2]:.3f}",
            "top2_annotation": top[1][0],
            "top2_n": top[1][1],
            "top2_fraction": f"{top[1][2]:.3f}",
            "top3_annotation": top[2][0],
            "top3_n": top[2][1],
            "top3_fraction": f"{top[2][2]:.3f}",
        })
        if top[1][2] >= 0.20:
            ann_label = f"{top[0][0]} + {top[1][0]}"
        else:
            ann_label = top[0][0]
        col_labels.append(f"post {cluster}\n{ann_label}")

    with (OUT_DIR / "regress_cluster_composition.csv").open("w", newline="") as handle:
        fieldnames = list(composition_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(composition_rows)

    long_rows: list[dict[str, object]] = []
    for i, row_name in enumerate(rows):
        for j, cluster in enumerate(cols):
            n = counts[i][j]
            if n == 0:
                continue
            long_rows.append({
                "original_annotation": row_name,
                "regress_cluster": cluster,
                "n_cells": n,
                "pct_of_original_annotation": f"{row_pct[i][j]:.1f}",
                "pct_of_regress_cluster": f"{col_pct[i][j]:.1f}",
            })
    with (OUT_DIR / "confusion_original_annotation_vs_regress_long.csv").open(
        "w", newline=""
    ) as handle:
        fieldnames = [
            "original_annotation",
            "regress_cluster",
            "n_cells",
            "pct_of_original_annotation",
            "pct_of_regress_cluster",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(long_rows)

    fig, ax = plt.subplots(figsize=(13, 8.5))
    image = ax.imshow(row_pct, cmap="mako_r" if "mako_r" in plt.colormaps() else "viridis",
                      vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=9)
    ax.set_xlabel("Cluster after Denisenko17 regression")
    ax.set_ylabel("Original annotation / cluster before regression")
    ax.set_title("Before vs After Stress Regression Confusion Matrix")

    for i in range(len(rows)):
        for j in range(len(cols)):
            if counts[i][j] == 0:
                continue
            text_color = "white" if row_pct[i][j] >= 45 else "black"
            ax.text(
                j,
                i,
                f"{counts[i][j]}\n{row_pct[i][j]:.1f}%",
                ha="center",
                va="center",
                fontsize=7,
                color=text_color,
            )

    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("% of original annotation")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Confusion_original_annotation_vs_regress_heatmap.png",
                dpi=220)
    plt.close(fig)

    print(f"Wrote formatted confusion matrix outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
