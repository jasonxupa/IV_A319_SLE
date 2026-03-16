# read_cart_analyze_cart.R
# Purpose: Import CAR-T patient scRNA-seq data (10x scRNA), create Seurat
#          objects, map cells to the SLE/A-319 reference for label transfer,
#          compute IFN gene signature scores via AUCell, and compare B-cell
#          subset proportions and IFN activity pre- vs. post-treatment.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(AUCell)
library(ggpubr)

# ---- Set working directory ---------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/14_analyze_CART")

# ==============================================================================
# SECTION 1: Read in CAR-T 10x count matrices and create Seurat objects
# ==============================================================================

# ---- Parse raw data file listing ---------------------------------------------
raw.data.files <- tibble(filename = list.files("RAW_DATA")) %>%
  mutate(
    patient           = str_split_fixed(filename, n = 4, pattern = "_")[, 2],
    timepoint         = str_split_fixed(filename, n = 4, pattern = "_")[, 3],
    filetype          = str_split_fixed(filename, n = 4, pattern = "_")[, 4],
    patient.timepoint = paste0(patient, ".", timepoint)
  )

# ---- Read count matrices per patient-timepoint -------------------------------
patient.scrna.list <- list()

for (pt.timepoint in unique(raw.data.files$patient.timepoint)) {
  print(pt.timepoint)
  to.read <- raw.data.files %>%
    filter(patient.timepoint == pt.timepoint) %>% arrange(filetype) %>%
    mutate(filename.full = paste0(
      "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/14_analyze_CART/RAW_DATA/",
      filename))

  patient.scrna.list[[pt.timepoint]] <- ReadMtx(
    mtx      = to.read$filename.full[4],
    cells    = to.read$filename.full[1],
    features = to.read$filename.full[2]
  )
}

# ---- Create individual Seurat objects and add sample metadata ----------------
patient.scrna.objects <- lapply(patient.scrna.list, CreateSeuratObject)

for (x in 1:length(patient.scrna.objects)) {
  object <- patient.scrna.objects[x]
  object[[1]]@meta.data$pt.timepoint <- names(patient.scrna.objects[x])
  object[[1]]@meta.data$og.cell.name <- colnames(patient.scrna.objects[x][[1]])

  colnames(patient.scrna.objects[x][[1]]) <-
    paste0(names(patient.scrna.objects[x]), "_", object[[1]]@meta.data$og.cell.name)

  object[[1]]@meta.data$cell.name <- colnames(patient.scrna.objects[x][[1]])
  patient.scrna.objects[x] <- object[[1]]
}

# ---- Merge all patient Seurat objects ----------------------------------------
patient.scrna.merge <- merge(patient.scrna.objects$pat1.post,
                              patient.scrna.objects[-1])
patient.scrna.merge <- JoinLayers(patient.scrna.merge)

# ==============================================================================
# SECTION 2: Pre-processing and QC
# ==============================================================================

patient.scrna.merge <- patient.scrna.merge %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)
patient.scrna.merge <- FindNeighbors(patient.scrna.merge,
                                      reduction = "pca", dims = 1:25)

# ---- Calculate mitochondrial percentage --------------------------------------
patient.scrna.merge$percent.mito <- PercentageFeatureSet(
  object = patient.scrna.merge, pattern = "^MT"
)

# ---- Feature plots for lineage marker QC ------------------------------------
pdf("CART_feature_Plot_raw.pdf", height = 20, width = 20)
FeaturePlot(patient.scrna.merge,
            features = c("CD3E", "CD19", "CD20", "LYZ", "NCAM1", "CD4",
                          "CD8A", "PF4", "CD14", "HBB", "TNFRSF17",
                          "percent.mito"))
dev.off()

saveRDS(patient.scrna.merge, "CART.patient.scrna.merge.unfiltered.105037.rds")

# ---- Filter cells with >10% mitochondrial reads -----------------------------
patient.scrna.merge.mito.lt.10 <- subset(patient.scrna.merge, percent.mito <= 10)

# ==============================================================================
# SECTION 3: Reference mapping - transfer A-319 cell-type labels to CAR-T data
# ==============================================================================

# ---- Load the A-319 reference with updated annotations ----------------------
all.patients.w.mdata.with.anno1.updated <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/11_make_global_umap/all.patients.w.mdata.with.anno1.updated.rds"
)

# ---- Load reanalyzed B-cell and T-cell objects for annotation updates --------
b.final.2.no.doublet.4829.cells.w.clonotype <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/12_reanalyze_b_add_clono/b.final.2.no.doublet.4829.cells.w.clonotype.rds"
)
t.final.112044.cells.updated.anno.with.clonotypes <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/13_reanalyze_t_add_clono/t.final.112044.cells.updated.anno.with.clonotypes.rds"
)

# ---- Reconcile annotations between original global and B/T re-analyses ------
og.annos <- all.patients.w.mdata.with.anno1.updated@meta.data %>%
  select(cell.name, patient, timepoint,
         "anno1.og" = anno1, "anno.1.updated.og" = anno.1.updated)

tb.updated <- rbind(
  b.final.2.no.doublet.4829.cells.w.clonotype@meta.data %>%
    select(cell.name, patient, timepoint,
           "anno1.bt" = anno1, "anno.1.updated.bt" = anno.1.updated),
  t.final.112044.cells.updated.anno.with.clonotypes@meta.data %>%
    select(cell.name, patient, timepoint,
           "anno1.bt" = anno1, "anno.1.updated.bt" = anno.1.updated)
)

annos.tb.update <- og.annos %>% left_join(tb.updated) %>% as_tibble()

# Prioritize B/T re-analysis annotations where available
annos.tb.update <- annos.tb.update %>%
  mutate(
    anno1 = case_when(
      anno1.og != anno1.bt & !is.na(anno1.bt)           ~ anno1.bt,
      anno1.og == anno1.bt & !is.na(anno1.bt)           ~ anno1.bt,
      is.na(anno1.bt)                                     ~ anno1.og,
      cell.name %in% duplicated.cells$cell.name           ~ anno1.og
    ),
    anno1.updated = case_when(
      anno.1.updated.og != anno.1.updated.bt & !is.na(anno.1.updated.bt) ~ anno.1.updated.bt,
      anno.1.updated.og == anno.1.updated.bt & !is.na(anno.1.updated.bt) ~ anno.1.updated.bt,
      is.na(anno.1.updated.bt)                                             ~ anno.1.updated.og,
      cell.name %in% duplicated.cells$cell.name                           ~ anno.1.updated.og
    )
  )

# Handle duplicated cell names (cells present in both B and T re-analyses)
duplicated.cells  <- annos.tb.update[duplicated(annos.tb.update$cell.name), ]
duplicated.cells2 <- annos.tb.update %>% filter(cell.name %in% duplicated.cells$cell.name)
annos.tb.update.non.dup <- annos.tb.update %>%
  filter(!cell.name %in% duplicated.cells$cell.name)
annos.tb.update.final   <- rbind(annos.tb.update.non.dup,
                                  duplicated.cells2[c(2, 4, 6, 8), ])

# ---- Update the reference with reconciled annotations -----------------------
all.patients.w.mdata.with.anno1.updated@meta.data <-
  all.patients.w.mdata.with.anno1.updated@meta.data %>%
  select(-anno1, -anno.1.updated) %>%
  left_join(annos.tb.update.final %>%
              select(cell.name, patient, timepoint, anno1, anno1.updated))

# ---- Assign lineage trajectory -----------------------------------------------
all.patients.w.mdata.with.anno1.updated@meta.data <-
  all.patients.w.mdata.with.anno1.updated@meta.data %>%
  mutate(trajectory = case_when(
    anno1.updated %in% c("Plasma", "A-NBC & csNBC", "NBC/A-NBC/csNBC") ~ "B",
    anno1.updated %in% c("NK", "aCD8", "CD4.CMT", "Treg", "nCD4", "nCD8",
                          "T/M doublet", "cycling T") ~ "T",
    TRUE ~ "M"
  ))

rownames(all.patients.w.mdata.with.anno1.updated@meta.data) <-
  all.patients.w.mdata.with.anno1.updated$cell.name

# Minor annotation corrections
all.patients.w.mdata.with.anno1.updated@meta.data$anno1.updated[
  all.patients.w.mdata.with.anno1.updated@meta.data$anno1.updated == "A-NBC & csNBC"
] <- "Plasma"
all.patients.w.mdata.with.anno1.updated@meta.data$anno1.updated[
  all.patients.w.mdata.with.anno1.updated@meta.data$anno1 == "gd"
] <- "gd"
all.patients.w.mdata.with.anno1.updated@meta.data$anno1.updated[
  all.patients.w.mdata.with.anno1.updated@meta.data$anno1 == "NKT"
] <- "NKT"

saveRDS(all.patients.w.mdata.with.anno1.updated,
        "all.patients.w.mdata.with.anno1.updated.11_9.rds")

# ---- Find transfer anchors and project CAR-T data ---------------------------
patient.anchors <- FindTransferAnchors(
  reference           = all.patients.w.mdata.with.anno1.updated,
  query               = patient.scrna.merge,
  k.filter            = NA,
  reference.reduction = "pca",
  dims                = 1:20
)

lupus.predictions <- TransferData(
  anchorset = patient.anchors,
  query     = patient.scrna.merge,
  reference = all.patients.w.mdata.with.anno1.updated,
  refdata   = list(trajectory    = "trajectory",
                   anno1.updated = "anno1.updated",
                   anno1         = "anno1")
)

# ---- Filter to mitochondrial <10% and visualize predictions -----------------
lupus.predictions.mt.lt.10 <- subset(lupus.predictions, percent.mito <= 10)

pdf("CART_mapped_to_lupus_ref.pdf")
DimPlot(lupus.predictions, group.by = "predicted.anno1", label = TRUE)
DimPlot(lupus.predictions, group.by = "predicted.anno1.updated", label = TRUE)
DimPlot(lupus.predictions, group.by = "predicted.trajectory", label = TRUE)
FeaturePlot(lupus.predictions,
            c("predicted.anno1.score", "predicted.anno1.updated.score",
              "predicted.trajectory.score"), label = TRUE)
dev.off()

pdf("CART_mapped_to_lupus_ref_mito_lt10.pdf")
DimPlot(lupus.predictions.mt.lt.10, group.by = "predicted.anno1",
        label = TRUE, raster = TRUE) + NoLegend()
DimPlot(lupus.predictions.mt.lt.10, group.by = "predicted.anno1.updated",
        label = TRUE, raster = TRUE) + NoLegend()
DimPlot(lupus.predictions.mt.lt.10, group.by = "predicted.trajectory",
        label = TRUE, raster = TRUE) + NoLegend()
FeaturePlot(lupus.predictions.mt.lt.10,
            c("predicted.anno1.score", "predicted.anno1.updated.score",
              "predicted.trajectory.score"),
            raster = TRUE, label = TRUE) + NoLegend()
dev.off()

saveRDS(lupus.predictions,
        "CART.patient.scRNA.with.predictions.105037.rds")
saveRDS(lupus.predictions.mt.lt.10,
        "CART.patient.scRNA.mt.10.filter.with.predictions.95k.rds")

# ==============================================================================
# SECTION 4: IFN gene signature scoring with AUCell
# ==============================================================================

# ---- Define published IFN and related molecular signatures -------------------
molecular_signatures <- list(
  Bennett = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "DNAPT6", "EIF2AK2",
              "FCGR1A", "HERC6"),
  Baechler = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "EIF2AK2", "EPSTI1",
               "FCGR1A", "FLJ20035", "HERC5", "HERC6", "IFI6", "IFI27",
               "IFI35", "IFI44"),
  Kirou   = c("IFI6", "IFI27", "IFI35"),
  Feng    = c("IFI6", "IFI27", "IFI35", "IFI44", "IFI44L"),
  Nikpour = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "DNAPT6", "EIF2AK2",
              "EPSTI1", "FCGR1A", "FLJ20035", "HERC5", "HERC6", "IFI6",
              "IFI27", "IFI35", "IFI44", "IFI44L", "IFIH1", "ISG15", "LAMP3",
              "LAP3", "LGALS3BP", "LGP2", "LOC129607", "LY6E", "MX1", "MX2",
              "OAS1", "OAS2", "OAS3", "OASL", "PML", "PLSCR1", "RNASE2",
              "RSAD2", "RTP4", "SERPING1", "SIGLEC1", "SP110", "STAT1",
              "STAT2", "TAP1", "UBE2L6", "USP18", "XAF1"),
  Landolt = c("IFI6", "IFI27", "IFI35", "IFI44", "IFI44L"),
  Petri   = c("IFI6", "IFI27", "IFI35"),
  Yao     = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "DNAPT6", "EIF2AK2",
              "EPSTI1", "FCGR1A", "FLJ20035", "HERC5", "HERC6", "IFI6",
              "IFI27", "IFI35", "IFI44", "IFI44L", "IFIH1", "ISG15", "LAMP3",
              "LAP3", "LGALS3BP", "LGP2"),
  Higgs   = c("IFI6", "IFI27", "IFI35", "IFI44", "RSAD2"),
  banchereau.ifn = c("ABCA1", "BATF2", "BTN3A1", "CARD16", "CARD17", "CCL8",
                     "CEACAM1", "CMPK2", "DHRS9", "DYNLT1", "GALM", "GBP1",
                     "GBP1", "GBP2", "HELZ2", "IFI44", "IFIH1", "IFIT3",
                     "IFITM1", "IRF7", "ISG20", "MX1", "NBN", "OAS1", "OAS1",
                     "OAS3", "PARP9", "SP110", "TMEM140", "TNFAIP6", "TNFSF10",
                     "TRIM21", "TRIM22", "TRIM5", "TRIM5", "UBE2L6", "ZBP1",
                     "ZNF684"),
  banchereau.pc = c("CRELD2", "DNAJB11", "FAM98A", "MYDGF", "NME1", "RFC4",
                    "SDF2L1", "SRM", "STT3A"),
  M1.2 = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "OAS1",
           "OAS2", "OAS3", "MX1", "MX2", "RSAD2", "STAT1", "STAT2", "STAT3",
           "STAT4", "STAT5A", "STAT5B", "IRF1", "IRF7", "IRF9", "IFNAR1",
           "IFNAR2", "IFNGR1", "IFNGR2", "IL6", "IL10", "IL12A", "IL12B",
           "IL18", "TNF", "TNFRSF1A", "TNFRSF1B", "CXCL10", "CXCL11",
           "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
           "CCL23", "CCL24", "CCL26", "CCL27", "CCL28", "CXCL9", "CXCL12",
           "CXCL13", "CXCL14", "CXCL15", "CXCL16", "CXCL17", "CXCL18",
           "CXCL19", "CXCL20", "CXCL21", "CXCL22", "CXCL23", "CXCL24",
           "CXCL25", "CXCL26", "CXCL27", "CXCL28", "CXCL29"),
  M3.4 = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "OAS1",
           "OAS2", "OAS3", "MX1", "MX2", "RSAD2", "STAT1", "STAT2", "STAT3",
           "STAT4", "STAT5A", "STAT5B", "IRF1", "IRF7", "IRF9", "IFNAR1",
           "IFNAR2", "IFNGR1", "IFNGR2", "IL6", "IL10", "IL12A", "IL12B",
           "IL18", "TNF", "TNFRSF1A", "TNFRSF1B", "CXCL10", "CXCL11",
           "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
           "CCL23", "CCL24", "CCL26", "CCL27", "CCL28", "CXCL9", "CXCL12",
           "CXCL13", "CXCL14", "CXCL15", "CXCL16", "CXCL17", "CXCL18",
           "CXCL19", "CXCL20", "CXCL21", "CXCL22", "CXCL23", "CXCL24",
           "CXCL25", "CXCL26", "CXCL27", "CXCL28", "CXCL29"),
  M5.12 = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "OAS1",
            "OAS2", "OAS3", "MX1", "MX2", "RSAD2", "STAT1", "STAT2", "STAT3",
            "STAT4", "STAT5A", "STAT5B", "IRF1", "IRF7", "IRF9", "IFNAR1",
            "IFNAR2", "IFNGR1", "IFNGR2", "IL6", "IL10", "IL12A", "IL12B",
            "IL18", "TNF", "TNFRSF1A", "TNFRSF1B", "CXCL10", "CXCL11",
            "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
            "CCL23", "CCL24", "CCL26", "CCL27", "CCL28", "CXCL9", "CXCL12",
            "CXCL13", "CXCL14", "CXCL15", "CXCL16", "CXCL17", "CXCL18",
            "CXCL19", "CXCL20", "CXCL21", "CXCL22", "CXCL23", "CXCL24",
            "CXCL25", "CXCL26", "CXCL27", "CXCL28", "CXCL29")
)

names(molecular_signatures) <- paste0(names(molecular_signatures), ".",
                                       lapply(molecular_signatures, length))

# ---- Compute AUCell rankings and scores on CAR-T data -----------------------
expr.matrix        <- lupus.predictions.mt.lt.10@assays$RNA$counts
cells_rankings_vst <- AUCell_buildRankings(as.matrix(expr.matrix))
cells_AUC.IFN      <- AUCell_calcAUC(molecular_signatures, cells_rankings_vst,
                                       aucMaxRank = nrow(cells_rankings_vst) * 0.25,
                                       verbose = TRUE)

lupus.predictions.mt.lt.10[["AUC_IFN"]] <- CreateAssayObject(
  cells_AUC.IFN@assays@data$AUC
)

saveRDS(lupus.predictions.mt.lt.10, "lupus.predictions.mt.lt.10.with.AUC.IFN.rds")

# ---- Define signatures to focus on for downstream analysis -------------------
mole.sig.keep <- c("banchereau.ifn.38", "M1.2.67", "Feng.5",
                    "Yao.22", "Higgs.5", "Landolt.5")

# ---- Parse patient and timepoint from sample labels --------------------------
lupus.predictions.mt.lt.10$timepoint <- str_split_fixed(
  lupus.predictions.mt.lt.10$pt.timepoint, pattern = "\\\\.", n = 2)[, 2]
lupus.predictions.mt.lt.10$patient <- str_split_fixed(
  lupus.predictions.mt.lt.10$pt.timepoint, pattern = "\\\\.", n = 2)[, 1]
lupus.predictions.mt.lt.10$timepoint <- factor(
  lupus.predictions.mt.lt.10$timepoint, levels = c("pre", "post"))

# ==============================================================================
# SECTION 5: Visualize IFN scores in CAR-T plasma cells (pre vs. post)
# ==============================================================================

pdf("CART_PlasmaCell_IFN_scoring.pdf")
VlnPlot(subset(lupus.predictions.mt.lt.10,
               predicted.anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)")),
        group.by = "timepoint",
        features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint") +
  ggtitle("IFN scores in plasmablasts, CAR-T patients, n=7")

VlnPlot(subset(lupus.predictions.mt.lt.10,
               predicted.anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)")),
        group.by = "patient",
        features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint") +
  ggtitle("IFN scores in plasmablasts by patient, CAR-T patients")

VlnPlot(subset(lupus.predictions.mt.lt.10,
               predicted.anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)")),
        group.by = "timepoint",
        features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "patient") +
  ggtitle("IFN scores in plasmablasts split by patient")
dev.off()

# ---- Differential IFN score testing (pre vs. post) in CAR-T plasma cells ----
plasma.cells <- subset(lupus.predictions.mt.lt.10,
                        predicted.anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)"))

ifn.markers <- FindMarkers(plasma.cells, ident.1 = "pre", ident.2 = "post",
                             group.by = "timepoint", assay = "AUC_IFN",
                             logfc.threshold = 0, pseudocount.use = 0)
ifn.markers <- ifn.markers[mole.sig.keep, ]
ifn.markers %>% mutate(fc = 2^avg_log2FC)

# ==============================================================================
# SECTION 6: IFN scoring in B-cell compartment (CAR-T)
# ==============================================================================

b.cell.compartment <- subset(
  lupus.predictions.mt.lt.10,
  predicted.anno1 %in% c("NBC", "NBC/Activated NBC", "Activated NBC",
                           "csNBC", "IgA+/IgG+ PC", "PB (Cycling)")
)
b.cell.compartment$predicted.anno1 <- factor(
  b.cell.compartment$predicted.anno1,
  levels = c("NBC", "NBC/Activated NBC", "Activated NBC",
             "csNBC", "IgA+/IgG+ PC", "PB (Cycling)")
)

pdf("IFN.expression.b.cells.CART.pdf")
VlnPlot(b.cell.compartment, group.by = "predicted.anno1",
        features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint") +
  ggtitle("IFN scores in B-cell subsets, CAR-T patients")
dev.off()

# ==============================================================================
# SECTION 7: IFN scoring in A-319 plasma cells for comparison
# ==============================================================================

b.final.2.no.doublet.4829.cells.w.clonotype <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/12_reanalyze_b_add_clono/b.final.2.no.doublet.4829.cells.w.clonotype.rds"
)

cells_rankings_vst <- AUCell_buildRankings(
  as.matrix(b.final.2.no.doublet.4829.cells.w.clonotype@assays$RNA$counts)
)
cells_AUC.IFN <- AUCell_calcAUC(molecular_signatures, cells_rankings_vst,
                                  aucMaxRank = nrow(cells_rankings_vst) * 0.25,
                                  verbose = TRUE)
b.final.2.no.doublet.4829.cells.w.clonotype[["AUC_IFN"]] <- CreateAssayObject(
  cells_AUC.IFN@assays@data$AUC
)

plasma.cells.a319 <- subset(b.final.2.no.doublet.4829.cells.w.clonotype,
                              anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)"))

ifn.markers2 <- FindMarkers(plasma.cells.a319, ident.1 = "Base", ident.2 = "M3",
                              group.by = "timepoint", assay = "AUC_IFN",
                              logfc.threshold = 0, pseudocount.use = 0)
ifn.markers2 <- ifn.markers2[mole.sig.keep, ]
ifn.markers2 %>% mutate(fc = 2^avg_log2FC)

# ==============================================================================
# SECTION 8: Read CAR-T VDJ (BCR) filtered contig annotations
# ==============================================================================

clones <- raw.data.files %>% filter(filetype == "filtered_contig_annotations.csv")

clones.list <- list()
for (x in 1:nrow(clones)) {
  print(x)
  clones.read <- clones[x, ] %>%
    mutate(filename.full = paste0(
      "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/14_analyze_CART/RAW_DATA/",
      filename))
  clones.list[[clones.read$patient.timepoint]] <- read_csv(clones.read$filename.full)
}

clones.list.bound <- do.call(rbind, clones.list)

# ==============================================================================
# SECTION 9: Compare B-cell subset proportions (A-319 vs. CAR-T, pre vs. post)
# ==============================================================================

# ---- B-cell proportion summary: CAR-T (detailed subtypes) -------------------
b.summ.stats.binned <- b.cell.compartment@meta.data %>%
  group_by(patient, timepoint, predicted.anno1, .drop = FALSE) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(
    b.cell.compartment@meta.data %>%
      group_by(patient, timepoint, .drop = FALSE) %>%
      summarize(n.total = n(), .groups = "drop")
  ) %>%
  mutate(prop = n / n.total)

# ---- B-cell proportion summary: A-319 (detailed subtypes) -------------------
b.summ.stats.binned2 <- b.final.2.no.doublet.4829.cells.w.clonotype@meta.data %>%
  group_by(patient, timepoint, anno1.factor, .drop = FALSE) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(
    b.final.2.no.doublet.4829.cells.w.clonotype@meta.data %>%
      group_by(patient, timepoint, .drop = FALSE) %>%
      summarize(n.total = n(), .groups = "drop")
  ) %>%
  mutate(prop = n / n.total)

# ---- Plot B-cell trajectory changes -----------------------------------------
pdf("B.traj.changes.319.vs.CAR.pdf")
ggplot(b.summ.stats.binned2 %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) + geom_point() +
  geom_line(aes(group = patient)) +
  facet_wrap(~anno1.factor, nrow = 1) +
  stat_compare_means(comparisons = list(c("Base", "M3")), paired = FALSE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("A319")

ggplot(b.summ.stats.binned, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~predicted.anno1, nrow = 1) +
  stat_compare_means(comparisons = list(c("pre", "post")), paired = FALSE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("CART")

ggplot(b.summ.stats.binned, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_line(aes(group = patient)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~predicted.anno1, nrow = 1) +
  stat_compare_means(comparisons = list(c("pre", "post")), paired = FALSE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("CART")
dev.off()

# ==============================================================================
# SECTION 10: Binned B-cell proportions (Naive / Memory / Plasma)
# ==============================================================================

# ---- A-319 binned categories ------------------------------------------------
b.final.2.no.doublet.4829.cells.w.clonotype@meta.data <-
  b.final.2.no.doublet.4829.cells.w.clonotype@meta.data %>%
  mutate(anno1.binned = recode_factor(
    anno1,
    "NBC" = "Naive", "NBC/Activated NBC" = "Naive",
    "Activated NBC" = "Memory", "csNBC" = "Memory",
    "IgA+/IgG+ PC" = "Plasma", "PB (Cycling)" = "Plasma"
  ))

b.summ.stats.binned3 <- b.final.2.no.doublet.4829.cells.w.clonotype@meta.data %>%
  group_by(patient, timepoint, anno1.binned, .drop = FALSE) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(
    b.final.2.no.doublet.4829.cells.w.clonotype@meta.data %>%
      group_by(patient, timepoint, .drop = FALSE) %>%
      summarize(n.total = n(), .groups = "drop")
  ) %>%
  mutate(prop = n / n.total)

# ---- CAR-T binned categories ------------------------------------------------
b.cell.compartment@meta.data <- b.cell.compartment@meta.data %>%
  mutate(anno1.binned = recode_factor(
    predicted.anno1,
    "NBC" = "Naive", "NBC/Activated NBC" = "Naive",
    "Activated NBC" = "Memory", "csNBC" = "Memory",
    "IgA+/IgG+ PC" = "Plasma", "PB (Cycling)" = "Plasma"
  ))

b.summ.stats.binned4 <- b.cell.compartment@meta.data %>%
  group_by(patient, timepoint, anno1.binned, .drop = FALSE) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(
    b.cell.compartment@meta.data %>%
      group_by(patient, timepoint, .drop = FALSE) %>%
      summarize(n.total = n(), .groups = "drop")
  ) %>%
  mutate(prop = n / n.total)

# ---- Boxplots of binned B-cell proportions -----------------------------------
pdf("B.cell.traj.change.over.time.binned.pdf")
ggplot(b.summ.stats.binned3 %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~anno1.binned, nrow = 1) +
  stat_compare_means(comparisons = list(c("Base", "M3")), paired = TRUE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("A319")

ggplot(b.summ.stats.binned4,
       aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~anno1.binned, nrow = 1) +
  stat_compare_means(comparisons = list(c("pre", "post")), paired = TRUE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("CART")

ggplot(b.summ.stats.binned3 %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_line(aes(group = patient)) + geom_point() +
  facet_wrap(~anno1.binned, nrow = 1) +
  stat_compare_means(comparisons = list(c("Base", "M3")), paired = TRUE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("A319")

ggplot(b.summ.stats.binned4,
       aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_line(aes(group = patient)) + geom_point() +
  facet_wrap(~anno1.binned, nrow = 1) +
  stat_compare_means(comparisons = list(c("pre", "post")), paired = TRUE) +
  theme_bw() + scale_fill_manual(values = blue.colors) + ggtitle("CART")
dev.off()

# ---- Save final objects ------------------------------------------------------
saveRDS(b.cell.compartment, "b.cell.compartment.CART.rds")
saveRDS(plasma.cells, "plasma.cells.CART.rds")

save.image("11_5_all_CART_data_duplicated_image_with_all_A319_data.rdata")