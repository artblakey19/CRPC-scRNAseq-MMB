# Numbat CNV inference — plan / TODO

Planned once BAM files are obtained. Goal: run haplotype-aware CNV inference
(Numbat), compare against copyKAT, and fold it into the copyKAT-vs-monocle3
malignisation comparison.

## Motivation

- copyKAT (expression-only) per-patient subclone clustering has **weak
  silhouette (~0.1 at every K)** — can't tell whether the CNV structure is
  genuine discrete subclones, a continuum (ongoing CIN), or just noise.
- All expression-only tools (copyKAT, inferCNV, SCEVAN) share this limit; the
  raw data only has 10x count matrices.
- **Numbat is haplotype-aware** — adds allele-imbalance signal (orthogonal,
  much less noisy than expression magnitude) and reconstructs a clonal
  phylogeny. It decisively separates "noise" from "real clonal structure".

## Prerequisites — obtain BAMs

- [ ] Request from the sequencing facility: cellranger
      `possorted_genome_bam.bam` + `.bai` per sample (CRPC1, CRPC2, CRPC3);
      ideally the whole `outs/` folder. Same 10x libraries as the
      `filtered_feature_bc_matrix` already in `Raw data/` → barcodes match.
- [ ] Phasing reference files:
  - 1000G SNP VCF (for cellsnp-lite pileup)
  - Eagle2 phasing reference panel
  - Numbat ships download instructions / a prebuilt 1000G reference.

## Pipeline

1. [ ] Install: R package `numbat` + `cellsnp-lite` + `Eagle2`
       (Numbat also has a Docker/conda image bundling these).
2. [ ] Per sample — `pileup_and_phase.R`: BAM + barcodes → allele counts at
       phased SNP sites.
3. [ ] Per sample — `run_numbat()`: expression matrix (from
       `filtered_feature_bc_matrix`) + allele dataframe → per-cell CNV,
       aneuploid/diploid call, subclones, clonal phylogeny tree.
   - Reference / diploid anchor: **cluster 8**, same as copyKAT & monocle3,
     for consistency.
4. [ ] Run **per-patient** (cohort = 3 independent late-stage DNPC patients;
       copyKAT and monocle3 are both per-patient).

## Compare against copyKAT

- [ ] aneuploid/diploid concordance, per cell.
- [ ] subclone concordance — do copyKAT and Numbat agree on the number of
      forms and on membership? (alluvial / contingency table)
- [ ] Does Numbat's clone separation beat copyKAT's ~0.1 silhouette?
      → this is the noise-vs-continuum verdict.
- [ ] Numbat clonal phylogeny vs copyKAT subclone CNA tree.

## Fold into the malignisation comparison

Extend `Epithelial_copyKAT_vs_Monocle3_compare.R` (or a Numbat variant):

- [ ] Dual-progression scatter with **Numbat** as the genomic axis
      (allele-based CNA distance / clone position) vs monocle3 pseudotime.
- [ ] Numbat subclone / phylogeny ↔ monocle3 pseudotime + trajectory branches.
- [ ] **Decision point:** if Numbat resolves crisp clones → the "2-3 discrete
      CNV forms" framing holds; if it also shows a continuum → keep the
      continuous-axis framing (continuous CNA progression vs pseudotime).

## Open questions

- Use Numbat's own subclone calls, or re-cluster for consistency with the
  copyKAT script's silhouette-based form definition?
- inferCNV subcluster mode is a free cross-check (count matrix only) but is a
  lateral move, not a fundamental upgrade — optional, low priority.

## Status

- BAMs: **not yet obtained** (user expects to be able to get them).
- Until BAMs arrive: copyKAT "forms" are treated as **exploratory /
  heatmap-validated**, not hard clonal claims.

## Related

- `Epithelial_copyKAT_perpatient.R` — current copyKAT per-patient pipeline
- `Epithelial_copyKAT_vs_Monocle3_compare.R` — malignisation comparison
- `Epithelial_Monocle3_persample.R` — per-patient pseudotime
