# Single-Cell RNA-seq Analysis Pipeline for A-319 Treatment in SLE

## Overview

This repository contains the R scripts used for processing and analyzing single-cell RNA sequencing (scRNA-seq) data from peripheral blood mononuclear cells (PBMCs) of patients with systemic lupus erythematosus (SLE) treated with A-319, a CD3xCD19 T cell engager, as described in the accompanying manuscript. The analysis pipeline covers the full workflow from raw data ingestion through cell type annotation, interferon (IFN) gene signature scoring, clonotype analysis, and integration with external CD19 CAR-T cell therapy data.

PBMC samples were collected at three timepoints: Baseline (Day 0), Day 28 (intra-treatment), and Month 3 (post-treatment) and processed using the 10X Genomics Chromium Single Cell Immune Profiling platform (5-prime v2). Sequencing data were demultiplexed and aligned to GRCh38 using Cell Ranger v7.2.0. Downstream analysis was performed in R using Seurat v5.

## Scripts

The scripts are numbered in the order they should be executed:

---

### 1_mkseurats.R - Create and Annotate Seurat Objects

**Purpose:** This script ingests the combined Cell Ranger output (a pre-assembled Seurat object containing all patients and timepoints) and prepares it for downstream analysis.

**Key steps:**
- Loads the combined Seurat object generated from Cell Ranger count matrices
- Calculates per-cell mitochondrial gene percentage (PercentageFeatureSet)
- Performs log-normalization (NormalizeData) and scaling (ScaleData)
- Adds clinical metadata to each cell, including patient ID, cohort assignment (Cohorts 1-3), dose level (0.3, 0.6, or 1.2 ug/kg), disease diagnosis, age, and disease duration
- Corrects timepoint labels (e.g., M13 to M3)
- Saves the annotated Seurat object as well as subsetted objects by patient, cohort, and timepoint for sample-level analyses
- Performs cell cycle scoring using Regev lab cell cycle gene sets (S-phase and G2M-phase genes)

**Primary output:** `all.patients.w.mdata.rds` - the master Seurat object with clinical metadata for all downstream scripts.

---

### 2_global_umap.R - Global UMAP Construction and Cell Type Annotation

**Purpose:** This script constructs a unified UMAP embedding of all cells across all patients and timepoints, integrates cell type annotations from separate B cell, T cell, and myeloid analyses, and assigns final cell type labels.

**Key steps:**
- Loads the master Seurat object and previously annotated B cell, T cell, and myeloid sub-objects
- Reads BCR and TCR filtered contig annotation files (from 10X VDJ) and flags cells with detectable B-cell or T-cell receptor sequences
- Runs the standard Seurat dimensionality reduction pipeline: FindVariableFeatures (3,000 genes), ScaleData, RunPCA (25 PCs), RunUMAP (20 PCs), FindNeighbors, FindClusters (multiple resolutions: 0.5, 1, 2)
- Integrates cell type annotations from separate B cell, T cell, and myeloid sub-analyses into the global object
- Applies quality control filters: mitochondrial RNA < 10%, UMI count > 500
- Manually annotates unannotated clusters based on marker gene expression (e.g., NK cells, cycling T cells, HSPCs, plasma cells) and differential expression (FindMarkers)
- Imputes remaining unannotated cells using K-nearest neighbors (k=50) on UMAP coordinates
- Assigns each cell to a major lineage trajectory (B, T, or Myeloid)
- Computes cell type proportions by timepoint
- Saves final subsetted objects for B cells, T cells, and myeloid cells

**Primary outputs:** `all.patients.w.mdata.with.anno1.updated.rds`, `b.final.rds`, `t.final.rds`, `m.final.rds`

---

### 3_global_umap_IFN.R - Interferon Signature Scoring on the Global UMAP

**Purpose:** This script computes interferon (IFN) gene signature activity scores across all cells in the global UMAP using AUCell, and visualizes IFN signature changes across timepoints and cell types for both A-319 and CAR-T treated patients.

**Key steps:**
- Defines multiple published IFN gene signatures (Banchereau, M1.2, Feng, Yao, Higgs, Landolt, and others) as well as plasmablast gene signatures
- Builds gene expression rankings and calculates AUCell scores using AUCell_buildRankings and AUCell_calcAUC with aucMaxRank = 25% of genes
- Stores AUC scores as a new assay (AUC_IFN) in the Seurat object
- Updates cell annotations (consolidating annotations, resolving NKT, gamma-delta T cells)
- Generates IFN score FeaturePlots on the global UMAP split by timepoint (Baseline vs. Month 3)
- Summarizes mean IFN scores per cell type and timepoint as heatmaps (pheatmap) and line plots
- Integrates CAR-T patient IFN scores for comparative visualization between A-319 and CAR-T treated cohorts
- Performs paired t-tests comparing IFN scores between Baseline and Month 3

**Primary output:** `Lupus_Reference_UMAP_Fixed.with.AUC.IFN.rds`

---

### 4_project_v5.R - Reference Dataset Projection

**Purpose:** This script projects each patient's cells onto two external reference atlases using Seurat's label transfer framework to obtain reference-based cell type annotations.

**Key steps:**
- **Bone Marrow Reference Projection:** Projects each patient's sample onto a healthy pediatric bone marrow/thymus reference dataset using FindTransferAnchors (20 PCs) and MapQuery to assign reference cell type labels (cell.type.short, cell.type.binned, trajectory)
- **Tonsil Reference Projection:** Projects each patient's sample onto the Human Cell Atlas tonsil reference atlas using FindTransferAnchors and TransferData to assign tonsil-derived annotations (annotation_level_1, annotation_figure_1, annotation_20230508)
- For each projection, sample-level Seurat objects are normalized (NormalizeData), variable features are identified (FindVariableFeatures) and scaled (ScaleData), and PCA is computed (RunPCA) before transfer anchor identification
- Saves projected metadata and UMAP coordinates for each patient

**Primary outputs:** Per-patient projected Seurat objects and metadata files (*_BM_proj.rds, *_tonsil.predictions.rds)

---

### 5_reanalyze_T.R - T Cell Compartment Re-analysis and Fine Annotation

**Purpose:** This script performs iterative sub-clustering of the T/NK cell compartment to achieve fine-grained cell type annotations across naive T cells, activated T cells, cycling T cells, NK cells, NKT cells, and gamma-delta T cells.

**Key steps:**
- Loads the T cell subset from the global UMAP (~131K cells)
- Re-runs the dimensionality reduction pipeline (FindVariableFeatures, ScaleData, RunPCA, RunUMAP) and clustering
- Removes contaminating non-T cell types (RBC, platelets, mast cells, neutrophils, megakaryocytes, HSPCs, B cell doublets) based on cluster identity and annotation
- **Naive T cell sub-analysis:** Subsets naive/central memory T cells, re-clusters, and annotates populations including nCD4, nCD8, CD4 Th1 CCR7+, CD4 Th2 CCR7-low, Treg Fr-like, Treg Eff-like, IFI+ CD4, and quiescent T cells using marker gene expression and differential expression analysis
- **Activated T cell sub-analysis:**
  - Cycling T/NK cells: Subsets and re-clusters proliferating cells, annotating Cycling CD4, Cycling CD8, and Cycling NK populations
  - NK and gamma-delta T cells: Subsets and re-clusters to identify NCAM1+/FCGR3A+ NK subsets, NKT cells, and gamma-delta T cells
  - Activated non-cycling T cells: Subsets and annotates GZMB+ aCD8, GZMK+ aCD8, and aCD4 populations
- Merges all T cell sub-compartments back into a single re-annotated T cell object
- Creates summarized annotations (anno1.updated) that group fine annotations into broader categories (e.g., nCD4, nCD8, CD4 Th, Cycling, aCD8, aCD4, NK, NKT, gd, Treg)

**Primary output:** `t.final.112044.cells.updated.anno.w.clonotype.rds`

---

### 6_analyze_BCR.R - BCR and TCR Clonotype Overlap Analysis

**Purpose:** This script analyzes B-cell receptor (BCR) and T-cell receptor (TCR) clonotype data from 10X VDJ sequencing to assess clonal overlap across treatment timepoints.

**Key steps:**
- Imports BCR and TCR filtered contig annotation files for all patients and timepoints
- Extracts patient and timepoint metadata from file paths
- For BCR clonotypes: groups by patient and CDR3 amino acid sequence to calculate clonotype frequency, proportion, and number of timepoints represented (Baseline, Day 28, Month 3)
- For TCR clonotypes: performs the same clonotype overlap analysis
- Calculates summary statistics: number of unique clonotypes, maximum/mean clonal proportions, and number of cells per overlap category
- Generates visualizations showing clonotype overlap by patient and across timepoints (stacked bar plots of clonotype sharing)
- Extracts chain information (heavy/light chain for BCR; alpha/beta chain for TCR) from CDR3 sequences
- Saves processed clonotype data for integration with cell-level metadata

**Primary outputs:** `b.clonotypes.rds`, `t.clonotypes.rds`

---

### 7_read_cart_analyze_cart.R - CAR-T Patient Data Processing and Integration

**Purpose:** This script reads, processes, and integrates external CD19 CAR-T cell therapy scRNA-seq data (from Wilheim et al., GEO: GSE263931) with the A-319 reference dataset to enable cross-treatment comparisons.

**Key steps:**
- Reads 10X raw data files (matrix, barcodes, features) for multiple CAR-T patients (pre- and post-treatment)
- Creates and merges Seurat objects for each patient-timepoint
- Standard preprocessing: NormalizeData, FindVariableFeatures, ScaleData, RunPCA, RunUMAP
- Calculates mitochondrial percentage and filters cells (mitochondrial < 10%)
- Updates the A-319 reference annotations by integrating the refined B cell and T cell annotations from scripts 5 and 6
- Projects CAR-T patient cells onto the A-319 reference using FindTransferAnchors (20 PCs) and TransferData to assign predicted cell types (trajectory, anno1.updated, anno1)
- Computes IFN signature scores (AUCell) on both the CAR-T data and the A-319 B cell subset
- Performs differential IFN expression analysis in plasma cells (FindMarkers on the AUC_IFN assay) between pre- and post-treatment timepoints for both CAR-T and A-319 cohorts
- Computes B-cell subtype proportions over time for both CAR-T and A-319 cohorts and generates comparative boxplots with statistical tests (Wilcoxon rank-sum, paired t-tests)
- Bins B cell subtypes into broader categories (Naive, Memory, Plasma) for simplified comparison

**Primary outputs:** `CART.patient.scRNA.with.predictions.rds`, `lupus.predictions.mt.lt.10.with.AUC.IFN.rds`

---

### 8_CART_IFN_do_CART_score.R - IFN Signature Scoring for CAR-T T Cells

**Purpose:** This script computes IFN activity scores specifically in the T cell compartment of CAR-T treated patients and compares IFN levels between pre- and post-treatment timepoints.

**Key steps:**
- Subsets CAR-T cells to the T cell compartment based on predicted annotations matching the A-319 T cell reference labels
- Reclassifies IFI+ CD4 cells as nCD4 for consistency with the A-319 annotation scheme
- Builds gene expression rankings and computes AUCell IFN scores using the same molecular signatures as Script 3
- Generates violin plots of IFN scores by predicted cell type and timepoint
- Extracts AUC scores into a long-format data frame for statistical analysis
- Computes mean IFN scores per cell type and timepoint, and generates trend line plots showing IFN changes pre- vs. post-CAR-T
- Performs paired t-tests comparing mean IFN scores across cell types between pre- and post-treatment
- Computes per-patient mean IFN scores and visualizes patient-level IFN trends

**Primary output:** `CART.t.only.with.AUC.IFN.rds`

---

## Dependencies

- **R** (v4.0.5 or later)
- **Seurat** v5
- **Azimuth** / **SeuratData**
- **tidyverse**
- **AUCell** (v1.12.0)
- **pheatmap**
- **patchwork**
- **ggpubr** (for stat_compare_means)
- **ggrepel** (for geom_text_repel)
- **FNN** (for get.knn in K-nearest neighbor imputation)
- **RColorBrewer**

## Data Availability

- Processed Seurat objects are available as an interactive Cell x Gene object at: https://cellxgene.cziscience.com/collections/6c147b03-08e4-4bb7-afb7-2bcc0a814b65
- Raw sequencing data are deposited under accession number [TBD]
- CAR-T comparison data (Wilheim et al.) were obtained from GEO accession GSE263931
