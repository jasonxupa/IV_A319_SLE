# mkseurats.R
# Purpose: Create Seurat objects from combined scRNA-seq data for SLE patients
#          treated with A-319. Adds patient metadata, computes QC metrics,
#          and saves per-patient/cohort/timepoint subsets.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)

# ---- Set working directory ---------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/01_mkseurats")

# ---- Read in the pre-combined Seurat object ----------------------------------
combined.seurat <- readRDS("/Users/jasonxu/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/Ptdata/combined_seurat02.rds")

# ---- Copy CellRanger web summary files to a single directory -----------------
# Consolidates web_summary.html files from each sample into one folder
# for convenient QC review.

source_dir <- "/Users/jasonxu/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/Ptdata/Ptdata"
target_dir <- file.path(source_dir, "web_summary_combined")

if (!dir.exists(target_dir)) {
  dir.create(target_dir)
}

web_files <- list.files(source_dir, recursive = TRUE,
                        pattern = "web_summary.html", full.names = TRUE)

for (file in web_files) {
  prefix        <- basename(dirname(file) %>% gsub(pattern = "/outs", replacement = ""))
  new_file_name <- paste0(prefix, "_web_summary.html")
  new_file_path <- file.path(target_dir, new_file_name)
  file.copy(file, new_file_path)
}

cat("Files copied and renamed successfully to:", target_dir, "\\n")

# ---- Calculate mitochondrial gene percentage (QC metric) ---------------------
mito.features <- grep(pattern = "^MT-",
                      x = rownames(x = combined.seurat@assays$RNA), value = TRUE)
combined.seurat[["percent.mito"]] <- PercentageFeatureSet(combined.seurat,
                                                          features = mito.features)

hist(combined.seurat$percent.mito)
summary(combined.seurat$percent.mito)
ggplot(combined.seurat@meta.data, aes(y = sample, x = percent.mito)) +
  geom_violin() + theme_bw()

table(combined.seurat$percent.mito > 10)

# ---- Normalize and scale data ------------------------------------------------
combined.seurat <- combined.seurat %>% NormalizeData() %>% ScaleData()

# ---- Add patient-level clinical metadata -------------------------------------
# Metadata table linking sample IDs to patient demographics and dosing cohorts.
mdata.to.add <- tibble(
  patient.id.nejm  = 1:11,
  patient          = c("S01","S02","S03","S04","S05","S06","S08",
                       "S09","S14","S16","S17"),
  cohort           = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L),
  dose             = c("0.3 ug/kg","0.3 ug/kg","0.3 ug/kg","0.3 ug/kg",
                       "0.3 ug/kg","0.3 ug/kg","0.6 ug/kg","0.6 ug/kg",
                       "1.2 ug/kg","1.2 ug/kg","1.2 ug/kg"),
  disease          = c("SLE/LN/SSc","SLE/LN/RA","SLE/LN","SLE/LN","SLE/LN",
                       "SLE/LN/NPSLE","SLE/LN/sSS","SLE/LN","SLE/LN/sSS",
                       "SLE/LN/sSS","SLE/MDA5+DM/APS/sSS/ILD/PAH"),
  age              = c(38L, 41L, 21L, 30L, 45L, 26L, 20L, 49L, 38L, 38L, 54L),
  disease.duration = c(10, 9, 2.5, 10, 6, 0.5, 7, 18, 6.5, 6.5, 26)
)

combined.seurat@meta.data$cell.name <- rownames(combined.seurat@meta.data)

combined.seurat@meta.data <- combined.seurat@meta.data %>%
  mutate(patient   = str_split_fixed(sample, pattern = "_", n = 2)[, 1],
         timepoint = str_split_fixed(sample, pattern = "_", n = 2)[, 2]) %>%
  left_join(mdata.to.add)

rownames(combined.seurat@meta.data) <- combined.seurat@meta.data$cell.name

# ---- Join RNA layers and fix timepoint label ---------------------------------
combined.seurat <- JoinLayers(combined.seurat, assay = "RNA")
combined.seurat@meta.data$timepoint[combined.seurat@meta.data$timepoint == "M13"] <- "M3"

# ---- Save the full combined Seurat object ------------------------------------
saveRDS(combined.seurat, "all.patients.w.mdata.rds")
saveRDS(as.matrix(GetAssayData(combined.seurat, layer = "counts")),
        file = "counts.mtx.rds")

# ---- Split and save Seurat objects by patient --------------------------------
dir.create("by_patient", showWarnings = FALSE)
metadata        <- combined.seurat@meta.data
unique_patients <- unique(metadata$patient)

for (patient_id in unique_patients) {
  patient_seurat    <- subset(combined.seurat, subset = patient == patient_id)
  cell.num          <- ncol(patient_seurat)
  patient_file_name <- paste0("by_patient/", patient_id, "_", cell.num, "_cells.rds")
  saveRDS(patient_seurat, file = patient_file_name)
  cat("Saved:", patient_file_name, "\\n")
}

# ---- Split and save Seurat objects by cohort ---------------------------------
dir.create("by_cohort", showWarnings = FALSE)
unique_cohorts <- unique(combined.seurat$cohort)

for (cohort_id in unique_cohorts) {
  cohort_seurat    <- subset(combined.seurat, subset = cohort == cohort_id)
  cell.num         <- ncol(cohort_seurat)
  cohort_file_name <- paste0("by_cohort/cohort_", cohort_id, "_", cell.num, "cells_.rds")
  saveRDS(cohort_seurat, file = cohort_file_name)
  cat("Saved:", cohort_file_name, "\\n")
}

# ---- Split and save Seurat objects by timepoint ------------------------------
dir.create("by_timepoint", showWarnings = FALSE)
unique_timepoints <- unique(metadata$timepoint)

for (timepoint_id in unique_timepoints) {
  timepoint_seurat    <- subset(combined.seurat, subset = timepoint == timepoint_id)
  cell.num            <- ncol(timepoint_seurat)
  timepoint_file_name <- paste0("by_timepoint/", timepoint_id, "_", cell.num, "_cells.rds")
  saveRDS(timepoint_seurat, file = timepoint_file_name)
  cat("Saved:", timepoint_file_name, "\\n")
}

# ---- Cell cycle scoring using Regev lab gene sets ----------------------------
cycle3    <- readLines("/Users/jasonxu/Library/CloudStorage/Box-Box/ETP_SingleCell_Project/GCB_Courses/536_kim/traf6/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
s.genes   <- cycle3[1:43]
g2m.genes <- cycle3[44:97]

combined.seurat <- CellCycleScoring(object       = combined.seurat,
                                    s.features   = s.genes,
                                    g2m.features = g2m.genes,
                                    set.ident    = FALSE)