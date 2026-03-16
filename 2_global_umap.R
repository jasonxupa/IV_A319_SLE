# global_umap.R
# Purpose: Construct a global UMAP embedding for all SLE PBMC cells across
#          patients and timepoints. Integrate cell-type annotations from
#          separate T-cell, B-cell, and myeloid analyses, assign labels to
#          unannotated cells via KNN imputation, and compute cell-type
#          proportions by timepoint.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(FNN)           # for get.knn()
library(RColorBrewer)

# ---- Set working directory ---------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/11_make_global_umap")

# ---- Load the merged Seurat object with patient metadata ---------------------
all.patients.w.mdata <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/01_mkseurats/all.patients.w.mdata.rds"
)

# ---- Load lineage-specific annotated objects ---------------------------------
b.only.jx.anno.with.AUCell <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/05_Analyze_B/b.only.jx.anno.with.AUCell.rds"
)
b.in.m.final.570.cells <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/07_analyze_myeloid/b.in.m.final.570.cells.rds"
)

# ==============================================================================
# SECTION 1: Import BCR and TCR filtered contig annotations
# ==============================================================================
# Used to flag cells with productive VDJ rearrangements.

# ---- BCR contig annotations --------------------------------------------------
b.clonotype.files <- list.files(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/09_BCR_TCR/TCR and BCR",
  recursive = TRUE, all.files = TRUE, full.names = TRUE
) %>%
  grep(pattern = "vdj_b", value = TRUE) %>%
  grep(pattern = "filtered_contig_annotations.csv", value = TRUE)

b.clonotype.list <- list()
for (file in b.clonotype.files) {
  b.clonotype.list[[file]] <- read_csv(file)
  b.clonotype.list[[file]]$file.name <- file
}
b.clonotypes <- do.call(rbind, b.clonotype.list)

b.clonotypes <- b.clonotypes %>%
  mutate(
    sample.timepoint = str_split_fixed(file.name, pattern = "/", n = 15)[, 13],
    patient          = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 1],
    timepoint        = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 2],
    cell.name        = paste0(patient, "_", timepoint, "_", barcode)
  )

# ---- TCR contig annotations --------------------------------------------------
t.clonotype.files <- list.files(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/09_BCR_TCR/TCR and BCR",
  recursive = TRUE, all.files = TRUE, full.names = TRUE
) %>%
  grep(pattern = "vdj_t", value = TRUE) %>%
  grep(pattern = "filtered_contig_annotations.csv", value = TRUE)

t.clonotype.list <- list()
for (file in t.clonotype.files) {
  t.clonotype.list[[file]] <- read_csv(file)
  t.clonotype.list[[file]]$file.name <- file
}
t.clonotypes <- do.call(rbind, t.clonotype.list)

t.clonotypes <- t.clonotypes %>%
  mutate(
    sample.timepoint = str_split_fixed(file.name, pattern = "/", n = 15)[, 13],
    patient          = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 1],
    timepoint        = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 2],
    cell.name        = paste0(patient, "_", timepoint, "_", barcode)
  )

# ==============================================================================
# SECTION 2: Dimensionality reduction and clustering (initial)
# ==============================================================================

all.patients.w.mdata <- all.patients.w.mdata %>%
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)

all.patients.w.mdata <- all.patients.w.mdata %>%
  FindNeighbors() %>%
  FindClusters(resolution = c(0.5, 1, 2))

# ==============================================================================
# SECTION 3: Consolidate cell-type annotations from sub-analyses
# ==============================================================================
# Merge annotations from independently analyzed T-cell, myeloid, and B-cell
# objects into a single annotation table.

all.t.combined.with.oxphos.AUC <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/06_analyze_T/all.t.combined.with.oxphos.AUC.rds"
)
m.only.metadata <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/07_analyze_myeloid/m.only.metadata.rds"
)
mono.only.metadata <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/07_analyze_myeloid/mono.only.metadata.rds"
)
neut.only.metadata <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/07_analyze_myeloid/neut.only.metadata.rds"
)

annos <- rbind(
  neut.only.metadata %>%
    select(cell.name, anno1 = cell.type.jx.1.m, anno2 = cell.type.jx.1.m.binned,
           anno3 = cell.type.jx.0.5.neut.binned, anno4 = cell.type.jx.0.5.neut),
  m.only.metadata %>%
    filter(!cell.name %in% c(neut.only.metadata$cell.name, mono.only.metadata$cell.name)) %>%
    select(cell.name, anno1 = cell.type.jx.1.m, anno2 = cell.type.jx.1.m.binned,
           anno3 = cell.type.jx.1.m.binned, anno4 = cell.type.jx.1.m.binned),
  mono.only.metadata %>%
    select(cell.name, anno1 = cell.type.jx.1.m, anno2 = cell.type.jx.1.m.binned,
           anno3 = cell.type.jx.0.25.mono.binned2, anno4 = cell.type.jx.0.25.mono.binned),
  all.t.combined.with.oxphos.AUC@meta.data %>%
    select(cell.name, anno1 = cell.type.jx.1.binned2, anno2 = cell.type.jx.1.binned,
           anno3 = cell.type.jx.1, anno4 = cell.type.jx.1),
  b.only.jx.anno.with.AUCell@meta.data %>%
    select(cell.name, anno1 = cell.type.jx.0.5, anno2 = cell.type.jx.0.5,
           anno3 = cell.type.jx.0.5, anno4 = cell.type.jx.0.5)
)

all.patients.w.mdata@meta.data <- all.patients.w.mdata@meta.data %>%
  left_join(annos)
rownames(all.patients.w.mdata@meta.data) <- all.patients.w.mdata@meta.data$cell.name

# ---- QC filtering: remove high-mito and low-count cells ---------------------
all.patients.w.mdata <- subset(all.patients.w.mdata, percent.mito < 10)
all.patients.w.mdata <- subset(all.patients.w.mdata, nCount_RNA > 500)

# ---- Flag cells with productive BCR / TCR rearrangements --------------------
all.patients.w.mdata$b.clone <- all.patients.w.mdata$cell.name %in% b.clonotypes$cell.name
all.patients.w.mdata$t.clone <- all.patients.w.mdata$cell.name %in% t.clonotypes$cell.name

# ==============================================================================
# SECTION 4: Identify and annotate unresolved clusters
# ==============================================================================

Idents(all.patients.w.mdata) <- all.patients.w.mdata$RNA_snn_res.1

c19.markers <- FindMarkers(all.patients.w.mdata, ident.1 = "19",
                           group.by = "RNA_snn_res.1", max.cells.per.ident = 1000)
c13.markers <- FindMarkers(all.patients.w.mdata, ident.1 = "13",
                           group.by = "RNA_snn_res.1", max.cells.per.ident = 1000)
c27.markers <- FindMarkers(all.patients.w.mdata, ident.1 = "27",
                           group.by = "RNA_snn_res.1", max.cells.per.ident = 1000)

# Assign labels to unannotated clusters based on marker gene expression
all.patients.w.mdata$anno1[all.patients.w.mdata$RNA_snn_res.1 %in% c(13, 19)] <- "T/M doublet"
all.patients.w.mdata$anno1[all.patients.w.mdata$RNA_snn_res.1 %in% c(9, 14)]  <- "NK"
all.patients.w.mdata$anno1[all.patients.w.mdata$RNA_snn_res.1 %in% c(21, 17)] <- "cycling T"
all.patients.w.mdata$anno1[all.patients.w.mdata$RNA_snn_res.1 %in% c(15)]     <- "NBC/A-NBC/csNBC"
all.patients.w.mdata$anno1[all.patients.w.mdata$RNA_snn_res.1 %in% c(25)]     <- "Plasma"
all.patients.w.mdata$anno1[all.patients.w.mdata$RNA_snn_res.1 %in% c(27)]     <- "HSPC"

# ==============================================================================
# SECTION 5: Optimized UMAP with multiple feature/dimension settings
# ==============================================================================
# Test different numbers of variable features and PCA dimensions to find the
# best visual separation of cell types.

# 1000 variable features
all.patients.w.mdata <- all.patients.w.mdata %>%
  FindVariableFeatures(nfeatures = 1000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap20_1k") %>%
  RunUMAP(dims = 1:15, reduction.name = "umap15_1k") %>%
  RunUMAP(dims = 1:10, reduction.name = "umap10_1k")

# 3000 variable features (selected for final visualization)
all.patients.w.mdata <- all.patients.w.mdata %>%
  FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap20_3k") %>%
  RunUMAP(dims = 1:15, reduction.name = "umap15_3k") %>%
  RunUMAP(dims = 1:10, reduction.name = "umap10_3k")

# ==============================================================================
# SECTION 6: KNN imputation for remaining NA annotations
# ==============================================================================
# Use k=50 nearest neighbors in UMAP space to assign cell-type labels to cells
# that were not captured by any of the lineage-specific analyses.

umap_coords <- all.patients.w.mdata@reductions$umap20_3k@cell.embeddings
annotations <- all.patients.w.mdata$anno1
na_cells    <- which(is.na(annotations))
k           <- 50
neighbors   <- get.knn(umap_coords, k = k)

for (i in na_cells) {
  neighbor_indices     <- neighbors$nn.index[i, ]
  neighbor_annotations <- annotations[neighbor_indices]
  most_common          <- names(sort(table(neighbor_annotations), decreasing = TRUE))[1]
  annotations[i]       <- most_common
}

all.patients.w.mdata$anno.1.updated <- annotations

# ---- Save annotated global UMAP PDF -----------------------------------------
pdf("global.umap.20.3k.pdf", width = 10, height = 8)
DimPlot(all.patients.w.mdata, group.by = "anno.1.updated", label = TRUE,
        reduction = "umap20_3k")
DimPlot(all.patients.w.mdata, group.by = "anno.1.updated", label = TRUE,
        reduction = "umap20_3k") + NoLegend()
dev.off()

# ==============================================================================
# SECTION 7: Assign lineage trajectory and compute cell-type proportions
# ==============================================================================

# ---- Create ordered factor for cell-type annotation -------------------------
all.patients.w.mdata$anno.1.updated.factor <- factor(
  all.patients.w.mdata$anno.1.updated,
  levels = c("Plasma", "A-NBC & csNBC", "NBC/A-NBC/csNBC",
             "T/M doublet", "aCD8", "CD4.CMT", "Treg", "nCD4", "nCD8",
             "NK", "cycling T",
             m.final$anno.1.updated %>% table() %>% sort(decreasing = TRUE) %>% names())
)

# ---- Assign broad lineage trajectory (B, T, or Myeloid) --------------------
all.patients.w.mdata@meta.data <- all.patients.w.mdata@meta.data %>%
  mutate(trajectory = case_when(
    anno.1.updated %in% c("Plasma", "A-NBC & csNBC", "NBC/A-NBC/csNBC") ~ "B",
    anno.1.updated %in% c("NK", "aCD8", "CD4.CMT", "Treg", "nCD4", "nCD8",
                          "T/M doublet", "cycling T") ~ "T",
    TRUE ~ "M"
  ))

rownames(all.patients.w.mdata@meta.data) <- all.patients.w.mdata$cell.name

# ---- Compute cell-type proportions by timepoint -----------------------------
prop.traj <- all.patients.w.mdata@meta.data %>%
  group_by(trajectory, timepoint) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(
    all.patients.w.mdata@meta.data %>%
      group_by(timepoint) %>% summarize(total = n(), .groups = "drop")
  ) %>%
  mutate(prop = n / total)

prop.traj2 <- all.patients.w.mdata@meta.data %>%
  group_by(anno.1.updated.factor, timepoint) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(
    all.patients.w.mdata@meta.data %>%
      group_by(timepoint) %>% summarize(total = n(), .groups = "drop")
  ) %>%
  mutate(prop = n / total)

# ---- Define color palette for cell types ------------------------------------
colors <- c(
  RColorBrewer::brewer.pal(3, "Blues"),
  colorRampPalette(RColorBrewer::brewer.pal(8, "Greens"))(10)[3:10],
  colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel1"))(14) %>% rev()
)

# ---- Plot trajectory proportions by timepoint --------------------------------
pdf("traj.by.timepoints.pdf")
ggplot(all.patients.w.mdata@meta.data, aes(x = timepoint, fill = trajectory)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip()
ggplot(prop.traj, aes(x = timepoint, y = prop, fill = trajectory)) +
  geom_col() +
  geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) +
  theme_bw() + coord_flip()
ggplot(prop.traj2, aes(x = timepoint, y = prop, fill = anno.1.updated.factor)) +
  geom_col() +
  geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) +
  theme_bw() + coord_flip() + scale_fill_manual(values = colors)
dev.off()

# ---- Final color-coded global UMAP ------------------------------------------
pdf("global.umap.20.3k.nolegend.pdf", width = 10, height = 10)
DimPlot(all.patients.w.mdata, group.by = "anno.1.updated.factor", label = TRUE,
        reduction = "umap20_3k", cols = colors) + NoLegend()
DimPlot(all.patients.w.mdata, group.by = "anno.1.updated.factor", label = TRUE,
        reduction = "umap20_3k", cols = colors)
dev.off()

# ==============================================================================
# SECTION 8: Subset and save lineage-specific objects
# ==============================================================================

b.final <- subset(all.patients.w.mdata,
                  b.clone == TRUE | anno.1.updated %in% c("Plasma", "A-NBC & csNBC", "NBC/A-NBC/csNBC"))
t.final <- subset(all.patients.w.mdata,
                  t.clone == TRUE | anno.1.updated %in% c("NK", "aCD8", "CD4.CMT", "Treg", "nCD4",
                                                          "nCD8", "T/M doublet", "cycling T"))
m.final <- subset(all.patients.w.mdata,
                  !anno.1.updated %in% c("NK", "aCD8", "CD4.CMT", "Treg", "nCD4", "nCD8",
                                         "T/M doublet", "cycling T", "Plasma", "A-NBC & csNBC",
                                         "NBC/A-NBC/csNBC"))

saveRDS(b.final, "b.final.6044.cells.rds")
saveRDS(t.final, "t.final.131227.cells.rds")
saveRDS(m.final, "m.final.168202.cells.rds")
saveRDS(all.patients.w.mdata, "all.patients.w.mdata.with.anno1.updated.rds")