# RNA Velocity Trajectory Workflow

This workflow adds RNA velocity to the existing filtered epithelial analysis.
It keeps the Seurat UMAP and curated epithelial annotation from:

`Results/05_Epithelial_Downstream/epi_annotated.rds`

## 1. Export Seurat metadata

```bash
Rscript scripts/07_RNA_Velocity/00_export_seurat_velocity_inputs.R
```

This writes:

- `Results/07_RNA_Velocity/seurat_velocity_metadata.csv`
- `Results/07_RNA_Velocity/seurat_umap.csv`
- `Results/07_RNA_Velocity/seurat_annotation_summary.csv`

## 2. Build velocity loom files

`velocyto.py` is required for BAM to loom counting. Pass a real GTF file,
preferably `genes.gtf` from the same Cell Ranger reference used for alignment.

If regular pip installation fails during build isolation, install the build
helpers first and then install velocyto without build isolation:

```bash
.venv/bin/python -m pip install wheel cython
.venv/bin/python -m pip install --no-build-isolation velocyto
```

```bash
GENES_GTF=/path/to/refdata-gex-GRCh38-*/genes/genes.gtf \
THREADS=16 \
bash scripts/07_RNA_Velocity/01_make_velocity_looms.sh
```

Sample-level parallel run, usually safer than running all three BAMs at once:

```bash
GENES_GTF=/path/to/refdata-gex-GRCh38-*/genes/genes.gtf \
THREADS=4 \
SAMTOOLS_MEMORY=1024 \
MAX_JOBS=2 \
bash scripts/07_RNA_Velocity/01_make_velocity_looms.sh
```

Optional repeat-mask GTF:

```bash
GENES_GTF=/path/to/genes.gtf \
REPEAT_MASK_GTF=/path/to/repeat_mask.gtf \
bash scripts/07_RNA_Velocity/01_make_velocity_looms.sh
```

The local `Resources/hg38_gencode_v27.txt` is not a GTF and is not sufficient
for this step.

## 3. Run scVelo trajectory analysis

Fast deterministic model, recommended for this environment:

```bash
.venv/bin/python scripts/07_RNA_Velocity/02_scvelo_trajectory.py \
  --mode deterministic \
  --n-top-genes 0
```

Stochastic model:

```bash
.venv/bin/python scripts/07_RNA_Velocity/02_scvelo_trajectory.py \
  --mode stochastic
```

Dynamical model with latent time:

```bash
.venv/bin/python scripts/07_RNA_Velocity/02_scvelo_trajectory.py \
  --mode dynamical
```

Main outputs:

- `Results/07_RNA_Velocity/scvelo/velocity_embedding_stream_by_annotation.png`
- `Results/07_RNA_Velocity/scvelo/velocity_pseudotime.png`
- `Results/07_RNA_Velocity/scvelo/latent_time.png` when `--mode dynamical`
- `Results/07_RNA_Velocity/scvelo/cell_velocity_metrics.csv`
- `Results/07_RNA_Velocity/scvelo/annotation_velocity_summary.csv`
- `Results/07_RNA_Velocity/scvelo/crpc_epithelial_scvelo.h5ad`
