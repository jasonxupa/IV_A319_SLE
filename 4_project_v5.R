# project_v5.R
# Purpose: Project SLE PBMC scRNA-seq data onto two reference atlases
#          (bone marrow and tonsil) using Seurat v5 label transfer to obtain
#          cell-type predictions for each patient.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(magrittr)

# ---- Set working directory ---------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/02_Project_to_tonsil")

# ==============================================================================
# SECTION 1: Project each patient onto the bone marrow (BM) reference
# ==============================================================================

# ---- Load the BM reference object with pre-computed UMAP for projection ------
ref.traj.umap4projection <- readRDS(
  "~/Library/CloudStorage/Box-Box/Jason_Files/scRNA/14_project_using_seurat_v4/2025-02-23_ProjectionFiles_6Samples_25PC_931VEG_PASKMG/ref.traj.umap4projection.rds"
)

# ---- List per-patient Seurat RDS files ---------------------------------------
patient.files <- list.files(
  "/Users/jasonxu/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/01_mkseurats/by_patient/",
  full.names = TRUE
)

# ---- Loop: Normalize, PCA, find anchors, and map to BM reference ------------
for (x in 1:length(patient.files)) {
  print(x)
  patient_seurat <- readRDS(patient.files[x])
  
  # Standard pre-processing
  
  patient_seurat <- patient_seurat %>%
    NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
    RunPCA() %>% RunUMAP(dims = 1:25)
  
  # Find transfer anchors to BM reference
  patient.anchors.bm <- FindTransferAnchors(
    reference            = ref.traj.umap4projection,
    query                = patient_seurat,
    k.filter             = NA,
    reference.reduction  = "pca",
    reference.neighbors  = "pca.annoy.neighbors",
    dims                 = 1:20
  )
  
  # Map query onto reference UMAP and transfer cell-type labels
  patient_seurat <- MapQuery(
    anchorset           = patient.anchors.bm,
    query               = patient_seurat,
    reference           = ref.traj.umap4projection,
    refdata             = list(cell.type.short  = "cell.type.short",
                               cell.type.binned = "cell.type.binned",
                               trajectory       = "trajectory"),
    reference.reduction = "pca",
    reduction.model     = "umap4projection"
  )
  
  filename <- patient_seurat@meta.data$patient %>% unique() %>% print()
  saveRDS(patient_seurat,            file = paste0(filename, "_BM_proj.rds"))
  saveRDS(patient_seurat@meta.data,  paste0(filename, "_metadata_BM_proj.rds"))
}

# ---- Save reference UMAP coordinates into metadata for each projection -------
projection.files <- list.files(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/02_Project_to_tonsil/",
  pattern = "[0-9]_BM_proj", full.names = TRUE
)

for (x in 1:length(projection.files)) {
  print(x)
  patient_seurat <- readRDS(projection.files[x])
  
  UMAP.coords <- patient_seurat@reductions$ref.umap@cell.embeddings %>%
    as_tibble(rownames = "cell.name")
  
  patient_seurat@meta.data <- patient_seurat@meta.data %>% left_join(UMAP.coords)
  rownames(patient_seurat@meta.data) <- patient_seurat$cell.name
  
  filename <- patient_seurat@meta.data$patient %>% unique() %>% print()
  saveRDS(patient_seurat,           file = paste0(filename, "_BM_proj_w_ref_coords_mdata.rds"))
  saveRDS(patient_seurat@meta.data, paste0(filename, "_metadata_BM_proj_w_ref_coords.rds"))
}

# ==============================================================================
# SECTION 2: Project each patient onto the tonsil reference atlas
# ==============================================================================

# ---- Load tonsil reference ---------------------------------------------------
rm(ref.traj.umap4projection)
tonsil <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/02_Project_to_tonsil/scRNA-seq/20230911_tonsil_atlas_rna_seurat_obj.rds"
)
tonsil <- tonsil %>% FindVariableFeatures()

# ---- Loop: Find anchors and transfer tonsil cell-type labels -----------------
for (x in 1:length(patient.files)) {
  print(x)
  patient_seurat <- readRDS(patient.files[x])
  
  patient_seurat <- patient_seurat %>%
    NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
    RunPCA()
  
  patient.anchors.tonsil <- FindTransferAnchors(
    reference           = tonsil,
    query               = patient_seurat,
    k.filter            = NA,
    reference.reduction = "pca",
    dims                = 1:20
  )
  
  tonsil.predictions <- TransferData(
    anchorset = patient.anchors.tonsil,
    query     = patient_seurat,
    reference = tonsil,
    refdata   = list(annotation_level_1  = "annotation_level_1",
                     annotation_figure_1 = "annotation_figure_1",
                     annotation_20230508 = "annotation_20230508")
  )
  
  filename <- patient_seurat@meta.data$patient %>% unique() %>% print()
  saveRDS(tonsil.predictions, file = paste0(filename, "_tonsil.predictions.rds"))
}

# ---- Save UMAP coordinates for tonsil projections ----------------------------
projection.files.tonsil <- list.files(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/02_Project_to_tonsil/",
  pattern = "[0-9]_tonsil", full.names = TRUE
)

for (x in 1:length(projection.files.tonsil)) {
  print(x)
  patient_seurat <- readRDS(projection.files.tonsil[x])
  
  UMAP.coords <- patient_seurat@reductions$ref.umap@cell.embeddings %>%
    as_tibble(rownames = "cell.name")
  
  patient_seurat@meta.data <- patient_seurat@meta.data %>% left_join(UMAP.coords)
  rownames(patient_seurat@meta.data) <- patient_seurat$cell.name
  
  filename <- patient_seurat@meta.data$patient %>% unique() %>% print()
  saveRDS(patient_seurat,           file = paste0(filename, "_BM_proj_w_ref_coords_mdata.rds"))
  saveRDS(patient_seurat@meta.data, paste0(filename, "_metadata_BM_proj_w_ref_coords.rds"))
}