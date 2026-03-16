# reanalyze_T.R
# Purpose: Sub-cluster the T/NK cell compartment from the global UMAP to
#          remove doublets/contaminating cells, identify naive T, activated T,
#          cycling T, NK, NKT, and other T subsets, and assign refined
#          cell-type annotations at multiple resolution levels.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# ---- Set working directory ---------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/13_reanalyze_t_add_clono")

# ---- Load T/NK lineage subset from global UMAP analysis ---------------------
t.final.131227.cells <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/11_make_global_umap/t.final.131227.cells.rds"
)

# ==============================================================================
# SECTION 1: Initial dimensionality reduction and clustering
# ==============================================================================

t.final.131227.cells <- t.final.131227.cells %>%
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)

t.final.131227.cells <- FindNeighbors(t.final.131227.cells,
                                      reduction = "pca", dims = 1:25)
t.final.131227.cells <- t.final.131227.cells %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 2)

# ---- Diagnostic plots to identify doublet/contaminating clusters -------------
pdf("T_doublets_17014_cells.pdf")
DimPlot(t.final.131227.cells, group.by = "RNA_snn_res.1", label = TRUE) + NoLegend()
DimPlot(t.final.131227.cells, group.by = "anno.1.updated", label = TRUE) + NoLegend()
DimPlot(t.final.131227.cells, group.by = "RNA_snn_res.0.5", label = TRUE) +
  DimPlot(t.final.131227.cells, group.by = "anno.1.updated", label = TRUE) + NoLegend()
DimPlot(t.final.131227.cells, group.by = "t.clone", label = TRUE)
dev.off()

# ==============================================================================
# SECTION 2: Remove doublets and non-T/NK contaminants
# ==============================================================================
# Keep clusters consistent with T/NK identity; exclude cells annotated as
# non-T/NK lineages (PB, RBC, Platelet, Mast, Neut, etc.).

t.final.no.doublet <- subset(
  t.final.131227.cells,
  RNA_snn_res.1 %in% c(1, 24, 27, 3, 16, 9, 2, 28, 25, 8, 19, 22,
                       0, 13, 4, 5, 15, 23, 10, 6, 14, 20, 26, 12, 21) &
    !anno.1.updated %in% c("PB (Cycling)", "RBC", "Platelet", "Mast", "Neut",
                           "Megakaryocyte", "Plasma", "T/M doublet", "HSPC",
                           "NBC/A-NBC/csNBC") &
    !anno2 %in% c("PB (Cycling)", "RBC", "Platelet", "Mast", "Neut",
                  "Megakaryocyte", "Plasma", "T/M doublet", "HSPC",
                  "NBC/A-NBC/csNBC")
)

cells.removed <- subset(
  t.final.131227.cells,
  !RNA_snn_res.1 %in% c(1, 24, 27, 3, 16, 9, 2, 28, 25, 8, 19, 22,
                        0, 13, 4, 5, 15, 23, 10, 6, 14, 20, 26, 12, 21) &
    anno.1.updated %in% c("PB (Cycling)", "RBC", "Platelet", "Mast", "Neut",
                          "Megakaryocyte", "Plasma", "T/M doublet", "HSPC",
                          "NBC/A-NBC/csNBC") &
    anno2 %in% c("PB (Cycling)", "RBC", "Platelet", "Mast", "Neut",
                 "Megakaryocyte", "Plasma", "T/M doublet", "HSPC",
                 "NBC/A-NBC/csNBC")
)

saveRDS(cells.removed, "cells.removed.17k.t.traj.rds")
saveRDS(t.final.no.doublet, "t.final.no.doublet.114k.rds")

# ---- Re-cluster the cleaned T/NK object -------------------------------------
t.final.no.doublet <- t.final.no.doublet %>%
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)
t.final.no.doublet <- t.final.no.doublet %>%
  FindClusters(resolution = c(0.5, 1))

# ---- Find markers to identify remaining contaminants ------------------------
Idents(t.final.no.doublet) <- t.final.no.doublet$RNA_snn_res.0.5
markers <- FindAllMarkers(t.final.no.doublet)

# Remove platelet-contaminated cluster (cluster 14 at res 0.5)
t.final.no.doublet.no.platelet <- subset(t.final.no.doublet,
                                         !RNA_snn_res.0.5 %in% c(14))

t.final.no.doublet.no.platelet <- t.final.no.doublet.no.platelet %>%
  FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:25)
t.final.no.doublet.no.platelet <- t.final.no.doublet.no.platelet %>%
  FindClusters(resolution = c(0.5, 1))

# ==============================================================================
# SECTION 3: Annotate the naive T-cell compartment
# ==============================================================================
# Isolate naive T cells (nCD4, nCD8, Tregs, Th subsets) from the broader
# T/NK object based on cluster identity and prior annotations.

naive.T <- subset(
  t.final.no.doublet.no.platelet,
  (RNA_snn_res.0.5 %in% c(14, 3, 5, 4, 17, 7) &
     !anno2 %in% c("CD16-CD56+ NK", "CD16+CD56- NK", "GNLY+ EM T",
                   "IFN+ CD8", "ZEB2+ EM.CD8", "GNLY-low aCD8")) |
    anno2 %in% c("Treg Fr-like", "Treg Eff-like", "CD4.CM.CCR7.low",
                 "CD4.CM.CCR7-int", "nCD8", "nCD4")
)

naive.T <- naive.T %>%
  FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:25)
naive.T <- naive.T %>%
  FindNeighbors() %>% FindClusters(resolution = c(0.5, 1))

# ---- Identify markers for each naive T cluster ------------------------------
Idents(naive.T) <- naive.T$RNA_snn_res.1
markers.n.t <- FindAllMarkers(naive.T)

# ---- Assign anno1 cell-type labels (coarse) ---------------------------------
naive.T@meta.data <- naive.T@meta.data %>%
  mutate(anno1 = case_when(
    RNA_snn_res.1 %in% c(2)                          ~ "CD4 Th2 CCR7-low",
    RNA_snn_res.1 %in% c(1, 3, 15, 13)               ~ "nCD4",
    RNA_snn_res.1 %in% c(14, 9, 5, 7, 17, 10, 12)    ~ "nCD8",
    RNA_snn_res.1 %in% c(0, 6)                        ~ "CD4 Th1 CCR7+",
    RNA_snn_res.1 == 4                                 ~ "Treg Fr-like",
    RNA_snn_res.1 == 8                                 ~ "Treg Eff-like",
    RNA_snn_res.1 == 11                                ~ "Mito-high FOXP1+ qT",
    RNA_snn_res.1 == 18                                ~ "IFI+ CD4",
    RNA_snn_res.1 == 16                                ~ "nCD4"
  ))

# ---- Assign anno2 cell-type labels (fine) ------------------------------------
naive.T@meta.data <- naive.T@meta.data %>%
  mutate(anno2 = case_when(
    RNA_snn_res.1 == 0  ~ "TH1.RB.high",
    RNA_snn_res.1 == 1  ~ "nCD4",
    RNA_snn_res.1 == 2  ~ "TH2.GATA3.hi",
    RNA_snn_res.1 == 3  ~ "nCD4",
    RNA_snn_res.1 == 4  ~ "Treg Fr-like",
    RNA_snn_res.1 == 5  ~ "nCD8",
    RNA_snn_res.1 == 6  ~ "TH1.IRF1.PD1+AP+.ag.presenting.CD4",
    RNA_snn_res.1 == 7  ~ "nCD8",
    RNA_snn_res.1 == 8  ~ "Treg-Eff-like",
    RNA_snn_res.1 == 9  ~ "nCD8",
    RNA_snn_res.1 == 10 ~ "nCD8",
    RNA_snn_res.1 == 11 ~ "FOXP1-hi.IFI.low qT",
    RNA_snn_res.1 == 12 ~ "nCD8",
    RNA_snn_res.1 == 13 ~ "nCD4",
    RNA_snn_res.1 == 14 ~ "nCD8",
    RNA_snn_res.1 == 15 ~ "nCD4",
    RNA_snn_res.1 == 16 ~ "nCD4",
    RNA_snn_res.1 == 17 ~ "nCD8",
    RNA_snn_res.1 == 18 ~ "IFI+ CD4"
  ))

# ---- Save naive T diagnostic plots ------------------------------------------
pdf("NaiveT_compartment.pdf")
DimPlot(naive.T, group.by = "anno2", label = TRUE) +
  DimPlot(naive.T, group.by = "patient", label = TRUE) + NoLegend() +
  DimPlot(naive.T, group.by = "timepoint", label = TRUE) + NoLegend() +
  DimPlot(naive.T, group.by = "anno1", label = TRUE) + NoLegend()
dev.off()

saveRDS(naive.T, "naive.T.rds")

# ==============================================================================
# SECTION 4: Annotate the activated/effector T-cell compartment
# ==============================================================================

activated.T <- subset(t.final.no.doublet.no.platelet,
                      !cell.name %in% naive.T$cell.name)

activated.T <- activated.T %>%
  FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:25)
activated.T <- activated.T %>%
  FindNeighbors() %>% FindClusters(resolution = c(0.5, 1))

Idents(activated.T) <- activated.T$RNA_snn_res.0.5
markers.activated.t <- FindAllMarkers(activated.T)

# ---- 4a: Cycling T/NK subset ------------------------------------------------
activated.T.cycling <- subset(activated.T, RNA_snn_res.0.5 %in% c(7, 13))
activated.T.cycling <- activated.T.cycling %>%
  FindVariableFeatures(nfeatures = 1000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:10)
activated.T.cycling <- activated.T.cycling %>%
  FindNeighbors() %>% FindClusters(resolution = c(0.1, 0.25, 0.5, 1))

Idents(activated.T.cycling) <- activated.T.cycling$RNA_snn_res.0.25
markers.cycling.t <- FindAllMarkers(activated.T.cycling)

# Annotate cycling T cells
activated.T.cycling@meta.data <- activated.T.cycling@meta.data %>%
  mutate(
    anno1 = case_when(
      RNA_snn_res.0.25 %in% c(0, 1, 2, 4, 5, 7) ~ "Cycling.T",
      RNA_snn_res.0.25 %in% c(3, 6)              ~ "Cycling.NK"
    ),
    anno2 = case_when(
      RNA_snn_res.0.25 %in% c(0, 1, 2, 5, 7) ~ "Cycling.CD8",
      RNA_snn_res.0.25 == 3                    ~ "Cycling.NK",
      RNA_snn_res.0.25 == 4                    ~ "Cycling.CD4",
      RNA_snn_res.0.25 == 6                    ~ "Cycling.NK"
    )
  )

saveRDS(activated.T.cycling, "activated.T.cycling.rds")

# ---- 4b: NK / NKT / gamma-delta T subset ------------------------------------
activated.T.NK.gd <- subset(activated.T,
                            RNA_snn_res.0.5 %in% c(3, 8, 10, 5, 12, 16))
activated.T.NK.gd <- activated.T.NK.gd %>%
  FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:10)
activated.T.NK.gd <- activated.T.NK.gd %>%
  FindNeighbors() %>% FindClusters(resolution = c(0.1, 0.25, 0.5, 1))

# Add BM/tonsil projection annotations for context
mdata.bm.tonsil.combined <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/03_changes_over_time/mdata.bm.tonsil.combined.rds"
)
to.add <- mdata.bm.tonsil.combined %>%
  select(cell.name, predicted.cell.type.short, predicted.annotation_20230508)
activated.T.NK.gd@meta.data <- activated.T.NK.gd@meta.data %>% left_join(to.add)
rownames(activated.T.NK.gd@meta.data) <- activated.T.NK.gd@meta.data$cell.name

Idents(activated.T.NK.gd) <- activated.T.NK.gd$RNA_snn_res.0.25
markers.activated.T.NK.gd <- FindAllMarkers(activated.T.NK.gd)

# Annotate NK / NKT / gamma-delta subsets
activated.T.NK.gd@meta.data <- activated.T.NK.gd@meta.data %>%
  mutate(
    anno1 = case_when(
      RNA_snn_res.0.25 %in% c(0, 1)    ~ "NCAM1-low FCGR3A+ NK",
      RNA_snn_res.0.25 %in% c(2, 3, 6) ~ "NKT",
      RNA_snn_res.0.25 %in% c(4, 7)    ~ "NCAM1+ FCGR3A-low NK",
      RNA_snn_res.0.25 == 5             ~ "gd"
    ),
    anno2 = case_when(
      RNA_snn_res.0.25 == 0             ~ "FCGR3A+ NK",
      RNA_snn_res.0.25 == 1             ~ "HLA-hi NCAM1-low FCGR3A+ NK",
      RNA_snn_res.0.25 %in% c(2, 3, 6) ~ "NKT",
      RNA_snn_res.0.25 %in% c(4, 7)    ~ "NCAM1+ FCGR3A-low NK",
      RNA_snn_res.0.25 == 5             ~ "gd"
    )
  )

saveRDS(activated.T.NK.gd, "activated.T.NK.gd.rds")

# ---- 4c: Activated CD8 / CD4 effectors (non-cycling, non-NK) ----------------
activated.T.nonNK.noncycling <- subset(
  activated.T,
  !cell.name %in% c(activated.T.cycling$cell.name, activated.T.NK.gd$cell.name)
)
activated.T.nonNK.noncycling <- activated.T.nonNK.noncycling %>%
  FindVariableFeatures(nfeatures = 1000) %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:10)
activated.T.nonNK.noncycling <- activated.T.nonNK.noncycling %>%
  FindNeighbors() %>% FindClusters(resolution = c(0.1, 0.25, 0.5, 1))

Idents(activated.T.nonNK.noncycling) <- activated.T.nonNK.noncycling$RNA_snn_res.0.1
markers.activated.T.nonNK.noncycling <- FindAllMarkers(activated.T.nonNK.noncycling)

# Annotate activated effector T cells
activated.T.nonNK.noncycling@meta.data <- activated.T.nonNK.noncycling@meta.data %>%
  mutate(
    anno1 = case_when(
      RNA_snn_res.0.25 %in% c(1, 2, 3, 4, 0) ~ "aCD8",
      RNA_snn_res.0.25 %in% c(5, 6)           ~ "aCD4"
    ),
    anno2 = case_when(
      RNA_snn_res.0.25 == 3             ~ "GZMK aCD8",
      RNA_snn_res.0.25 %in% c(0, 1, 2, 4) ~ "GZMB aCD8",
      RNA_snn_res.0.25 %in% c(5, 6)    ~ "aCD4"
    )
  )

pdf("activated.T.nonNK.noncycling_featurePlot_w_anno.pdf")
FeaturePlot(activated.T.nonNK.noncycling,
            c("CD4", "CD8A", "CD3E", "NCAM1", "HIST1H4C", "KLRC2",
              "FCGR3A", "FCGR3A", "TOP2A", "GNLY", "GZMB", "GZMK",
              "CXCR6", "TRGV2", "TRGC2", "CD2")) +
  DimPlot(activated.T.nonNK.noncycling, group.by = "anno2", label = TRUE) +
  DimPlot(activated.T.nonNK.noncycling, group.by = "anno1", label = TRUE) +
  DimPlot(activated.T.nonNK.noncycling, group.by = "patient", label = TRUE) + NoLegend() +
  DimPlot(activated.T.nonNK.noncycling, group.by = "t.clone", label = TRUE) + NoLegend()
dev.off()

saveRDS(activated.T.nonNK.noncycling, "activated.T.nonNK.noncycling.rds")

# ==============================================================================
# SECTION 5: Merge all T/NK subsets and create final annotated object
# ==============================================================================

t.final.112044.cells.updated.anno <- merge(
  naive.T,
  list(activated.T.nonNK.noncycling, activated.T.NK.gd, activated.T.cycling)
)
t.final.112044.cells.updated.anno <- JoinLayers(t.final.112044.cells.updated.anno)

t.final.112044.cells.updated.anno <- t.final.112044.cells.updated.anno %>%
  FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>%
  RunPCA(npcs = 30) %>% RunUMAP(dims = 1:30)

# ---- Create a simplified (updated) annotation collapsing subsets -------------
t.final.112044.cells.updated.anno@meta.data <- t.final.112044.cells.updated.anno@meta.data %>%
  mutate(anno1.updated = case_when(
    anno1 %in% c("nCD4", "IFI+ CD4")                         ~ "nCD4",
    anno1 == "nCD8"                                            ~ "nCD8",
    anno1 %in% c("CD4 Th1 CCR7+", "CD4 Th2 CCR7-low")       ~ "CD4 Th",
    anno1 %in% c("Cycling.NK", "Cycling.T")                   ~ "Cycling",
    anno1 == "Mito-high FOXP1+ qT"                            ~ "Activating",
    anno1 == "aCD8"                                            ~ "aCD8",
    anno1 == "aCD4"                                            ~ "aCD4",
    anno1 == "gd"                                              ~ "gd",
    anno1 %in% c("NCAM1+ FCGR3A-low NK", "NCAM1-low FCGR3A+ NK") ~ "NK",
    anno1 == "NKT"                                             ~ "NKT",
    anno1 %in% c("Treg Fr-like", "Treg Eff-like")             ~ "Treg"
  ))

# ---- Save final annotated T/NK object ---------------------------------------
saveRDS(t.final.112044.cells.updated.anno,
        "t.final.112044.cells.updated.anno.w.clonotype.rds")