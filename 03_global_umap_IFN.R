# global_umap_IFN.R
# Purpose: Compute IFN gene signature scores (AUCell) on the global UMAP
#          object for A-319 treated patients, visualize IFN activity across
#          cell types and timepoints, and generate heatmaps and trend-line
#          plots comparing IFN reduction in A-319 vs. CAR-T cohorts.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(AUCell)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# ==============================================================================
# SECTION 1: Load the global reference UMAP with fixed coordinates
# ==============================================================================

Lupus_Reference_UMAP_Fixed <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/14_analyze_CART/Lupus_Reference_UMAP_Fixed.rds"
)

# ==============================================================================
# SECTION 2: Define published IFN and related molecular signatures
# ==============================================================================

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

mole.sig.keep <- c("banchereau.ifn.38", "M1.2.67", "Feng.5",
                    "Yao.22", "Higgs.5", "Landolt.5")

# ==============================================================================
# SECTION 3: Compute AUCell IFN scores on the A-319 reference
# ==============================================================================

cells_rankings_vst <- AUCell_buildRankings(
  Lupus_Reference_UMAP_Fixed@assays$RNA$counts
)
cells_AUC.IFN <- AUCell_calcAUC(molecular_signatures, cells_rankings_vst,
                                  aucMaxRank = nrow(cells_rankings_vst) * 0.25,
                                  verbose = TRUE)

Lupus_Reference_UMAP_Fixed[["AUC_IFN"]] <- CreateAssayObject(
  cells_AUC.IFN@assays@data$AUC
)

# ==============================================================================
# SECTION 4: Visualize IFN scores on the global UMAP (A-319 data)
# ==============================================================================

pdf("IFN_score_global_umap.pdf", width = 8, height = 6.5)
FeaturePlot(subset(Lupus_Reference_UMAP_Fixed, timepoint == "Base"),
            features = c("banchereau.ifn.38", "M1.2.67", "Yao.22"),
            reduction = "umap20_3k_fix", order = TRUE, raster = TRUE,
            raster.dpi = c(500, 500), split.by = "timepoint")
FeaturePlot(subset(Lupus_Reference_UMAP_Fixed, timepoint == "M3"),
            features = c("banchereau.ifn.38", "M1.2.67", "Yao.22"),
            reduction = "umap20_3k_fix", order = TRUE, raster = TRUE,
            raster.dpi = c(500, 500))
dev.off()

pdf("IFN_score_global_umap_2.pdf", width = 12, height = 5)
FeaturePlot(subset(Lupus_Reference_UMAP_Fixed, timepoint == "Base"),
            features = c("banchereau.ifn.38", "M1.2.67"),
            reduction = "umap20_3k_fix", order = FALSE, raster = TRUE,
            raster.dpi = c(500, 500))
FeaturePlot(subset(Lupus_Reference_UMAP_Fixed, timepoint == "D28"),
            features = c("banchereau.ifn.38", "M1.2.67"),
            reduction = "umap20_3k_fix", order = FALSE, raster = TRUE,
            raster.dpi = c(500, 500))
FeaturePlot(subset(Lupus_Reference_UMAP_Fixed, timepoint == "M3"),
            features = c("banchereau.ifn.38", "M1.2.67"),
            reduction = "umap20_3k_fix", order = FALSE, raster = TRUE,
            raster.dpi = c(500, 500))
dev.off()

# ==============================================================================
# SECTION 5: Update cell-type annotations on the reference
# ==============================================================================

# Merge reconciled B/T annotations back into the reference metadata
annos.tb.update.final <- rbind(annos.tb.update.non.dup,
                                duplicated.cells2[c(2, 4, 6, 8), ])

Lupus_Reference_UMAP_Fixed$cell.name <- colnames(Lupus_Reference_UMAP_Fixed@assays$RNA)
Lupus_Reference_UMAP_Fixed@meta.data <- Lupus_Reference_UMAP_Fixed@meta.data %>%
  select(-anno1, -anno1.updated) %>%
  left_join(annos.tb.update.final %>% select(cell.name, anno1, anno1.updated))

# ---- Assign lineage trajectory -----------------------------------------------
Lupus_Reference_UMAP_Fixed@meta.data <- Lupus_Reference_UMAP_Fixed@meta.data %>%
  mutate(trajectory = case_when(
    anno1.updated %in% c("Plasma", "A-NBC & csNBC", "NBC/A-NBC/csNBC") ~ "B",
    anno1.updated %in% c("NK", "aCD8", "CD4.CMT", "Treg", "nCD4", "nCD8",
                          "T/M doublet", "cycling T") ~ "T",
    TRUE ~ "M"
  ))

# ---- Minor annotation corrections -------------------------------------------
Lupus_Reference_UMAP_Fixed@meta.data$anno1.updated[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1.updated == "A-NBC & csNBC"
] <- "Plasma"
Lupus_Reference_UMAP_Fixed@meta.data$anno1.updated[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1 == "gd"
] <- "gd"
Lupus_Reference_UMAP_Fixed@meta.data$anno1.updated[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1 == "NKT"
] <- "NKT"

# ---- Rename specific anno1 labels for consistency ----------------------------
Lupus_Reference_UMAP_Fixed@meta.data$anno1[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1 == "Plasma"
] <- "IgA+/IgG+ PC"
Lupus_Reference_UMAP_Fixed@meta.data$anno1[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1 == "IFI+ CD4"
] <- "nCD4"
Lupus_Reference_UMAP_Fixed@meta.data$anno1[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1 == "Treg"
] <- "Treg Fr-like"
Lupus_Reference_UMAP_Fixed@meta.data$anno1[
  Lupus_Reference_UMAP_Fixed@meta.data$anno1 == "CD4.CMT"
] <- "CD4 Th1 CCR7+"

rownames(Lupus_Reference_UMAP_Fixed@meta.data) <-
  Lupus_Reference_UMAP_Fixed$cell.name

# ---- Create ordered factor for cell-type annotation -------------------------
Lupus_Reference_UMAP_Fixed$anno.1.updated.factor <- factor(
  Lupus_Reference_UMAP_Fixed$anno1.updated,
  levels = c("Plasma", "NBC/A-NBC/csNBC",
             "T/M doublet", "aCD8", "CD4.CMT", "Treg", "nCD4", "nCD8",
             "NK", "cycling T", "NKT", "gd",
             m.final$anno.1.updated %>% table() %>% sort(decreasing = TRUE) %>% names())
)

saveRDS(Lupus_Reference_UMAP_Fixed,
        "Lupus_Reference_UMAP_Fixed.with.AUC.IFN.rds")

# ==============================================================================
# SECTION 6: Heatmap of mean IFN scores by cell type and timepoint (A-319)
# ==============================================================================

to.heatmap <- cbind(Lupus_Reference_UMAP_Fixed@meta.data,
                     t(Lupus_Reference_UMAP_Fixed@assays$AUC_IFN@counts))

# ---- Mean Banchereau IFN score per cell type and timepoint -------------------
matrix <- to.heatmap %>%
  select(banchereau.ifn.38, cell.name, anno.1.updated.factor, timepoint) %>%
  group_by(anno.1.updated.factor, timepoint) %>%
  summarize(mean.ifn = mean(banchereau.ifn.38, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(values_from = mean.ifn, names_from = timepoint)

mtx <- matrix[, 2:4] %>% as.matrix()
rownames(mtx) <- matrix$anno.1.updated.factor

# ---- Heatmap: clustered and unclustered views --------------------------------
pheatmap::pheatmap(mtx, cluster_cols = FALSE, cluster_rows = FALSE, scale = "row")
pheatmap::pheatmap(mtx, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
pheatmap::pheatmap(mtx, cluster_cols = FALSE, cluster_rows = TRUE, cutree_rows = 4)
pheatmap::pheatmap(mtx[, c(1, 3)], cluster_cols = FALSE, cluster_rows = FALSE,
                   scale = "row")

# ---- Heatmap using fine (anno1) annotations ----------------------------------
matrix2 <- to.heatmap %>%
  select(banchereau.ifn.38, cell.name, anno1, timepoint) %>%
  group_by(anno1, timepoint) %>%
  summarize(mean.ifn = mean(banchereau.ifn.38, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(values_from = mean.ifn, names_from = timepoint)

mtx2 <- matrix2[, 2:4] %>% as.matrix()
rownames(mtx2) <- matrix2$anno1

pheatmap::pheatmap(mtx2, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")

# ==============================================================================
# SECTION 7: Trend-line plots of IFN score changes (A-319)
# ==============================================================================

# ---- Define color palette ----------------------------------------------------
colors <- c(
  RColorBrewer::brewer.pal(3, "Blues")[c(1, 3)],
  colorRampPalette(RColorBrewer::brewer.pal(8, "Greens"))(10)[3:10],
  colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel1"))(14) %>% rev()
)

# ---- Mean IFN per cell type (A-319) -----------------------------------------
to.ggplot <- to.heatmap %>%
  select(banchereau.ifn.38, cell.name, anno.1.updated.factor, timepoint) %>%
  group_by(anno.1.updated.factor, timepoint) %>%
  summarize(mean.ifn = mean(banchereau.ifn.38, na.rm = TRUE),
            n = n(), .groups = "drop")

pdf("IFN.reduction.all.319.pdf")
ggplot(to.ggplot %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor)) +
  geom_point(size = 4) +
  geom_line(aes(group = anno.1.updated.factor)) +
  theme_bw() + scale_color_manual(values = colors) +
  geom_text(aes(label = anno.1.updated.factor)) +
  stat_compare_means(comparisons = list(c("Base", "M3")),
                     method = "t.test", paired = TRUE)

ggplot(to.ggplot %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor,
           group = anno.1.updated.factor)) +
  geom_point(size = 4) + geom_line() +
  theme_bw() + scale_color_manual(values = colors)
dev.off()

# ==============================================================================
# SECTION 8: Trend-line plots of IFN score changes (CAR-T)
# ==============================================================================

lupus.predictions.mt.lt.10.with.AUC.IFN <- readRDS(
  "~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/14_analyze_CART/lupus.predictions.mt.lt.10.with.AUC.IFN.rds"
)

lupus.predictions.mt.lt.10.with.AUC.IFN@meta.data <- cbind(
  lupus.predictions.mt.lt.10.with.AUC.IFN@meta.data,
  lupus.predictions.mt.lt.10.with.AUC.IFN@assays$AUC_IFN@counts %>% t()
)

lupus.predictions.mt.lt.10.with.AUC.IFN@meta.data$anno.1.updated.factor <- factor(
  lupus.predictions.mt.lt.10.with.AUC.IFN@meta.data$predicted.anno1.updated,
  levels = c("Plasma", "NBC/A-NBC/csNBC", "T/M doublet", "aCD8", "CD4.CMT",
             "Treg", "nCD4", "nCD8", "NK", "cycling T", "NKT", "gd",
             m.final$anno.1.updated %>% table() %>% sort(decreasing = TRUE) %>% names())
)

lupus.predictions.mt.lt.10.with.AUC.IFN$timepoint <- str_split_fixed(
  lupus.predictions.mt.lt.10.with.AUC.IFN$pt.timepoint, pattern = "\\\\.", n = 2)[, 2]
lupus.predictions.mt.lt.10.with.AUC.IFN$timepoint <- factor(
  lupus.predictions.mt.lt.10.with.AUC.IFN$timepoint, levels = c("pre", "post"))

# ---- Mean IFN per cell type (CAR-T) -----------------------------------------
to.ggplot.cart <- lupus.predictions.mt.lt.10.with.AUC.IFN@meta.data %>%
  select(banchereau.ifn.38, cell.name, anno.1.updated.factor, timepoint) %>%
  group_by(anno.1.updated.factor, timepoint) %>%
  summarize(mean.ifn = mean(banchereau.ifn.38, na.rm = TRUE),
            n = n(), .groups = "drop")

pdf("IFN.reduction.all.CART.pdf")
ggplot(to.ggplot.cart,
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor,
           group = anno.1.updated.factor)) +
  geom_point(size = 4) + geom_line() +
  theme_bw() + scale_color_manual(values = colors) +
  geom_text(aes(label = anno.1.updated.factor))

ggplot(to.ggplot.cart,
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor,
           group = anno.1.updated.factor)) +
  geom_point(size = 4) + geom_line() +
  theme_bw() + scale_color_manual(values = colors) +
  geom_text_repel(aes(label = anno.1.updated.factor))
dev.off()

# ==============================================================================
# SECTION 9: Combined A-319 vs. CAR-T IFN reduction comparison
# ==============================================================================

pdf("IFN.reduction.all.CART.a319.revised.colors.pdf")

# A-319: Base vs. M3
ggplot(to.ggplot %>%
         filter(timepoint != "D28", !anno.1.updated.factor %in% c("NKT", "gd")),
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor,
           group = anno.1.updated.factor)) +
  geom_point(size = 4) + geom_line() +
  theme_bw() + scale_color_manual(values = colors) +
  geom_text(aes(label = anno.1.updated.factor)) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), limits = c(0.1, 0.45)) +
  stat_compare_means(comparisons = list(c("Base", "M3")),
                     method = "t.test", paired = TRUE)

# CAR-T: pre vs. post
ggplot(to.ggplot.cart %>% filter(!anno.1.updated.factor %in% c("NKT", "gd")),
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor)) +
  geom_point(size = 4) +
  geom_line(aes(group = anno.1.updated.factor)) +
  theme_bw() + scale_color_manual(values = colors) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), limits = c(0.1, 0.45)) +
  stat_compare_means(comparisons = list(c("pre", "post")),
                     method = "t.test", paired = TRUE)

dev.off()

# ==============================================================================
# SECTION 10: Additional IFN module score comparisons (M1.2 module)
# ==============================================================================

to.ggplot2 <- to.heatmap %>%
  select(M1.2.67, anno.1.updated.factor, timepoint) %>%
  group_by(anno.1.updated.factor, timepoint) %>%
  summarize(mean.ifn = mean(M1.2.67, na.rm = TRUE),
            n = n(), .groups = "drop")

to.ggplot.cart2 <- lupus.predictions.mt.lt.10.with.AUC.IFN@meta.data %>%
  select(M1.2.67, cell.name, anno.1.updated.factor, timepoint) %>%
  group_by(anno.1.updated.factor, timepoint) %>%
  summarize(mean.ifn = mean(M1.2.67, na.rm = TRUE),
            n = n(), .groups = "drop")

pdf("IFN.reduction.all.CART.a319.M1.2.67.pdf")
ggplot(to.ggplot2 %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor,
           group = anno.1.updated.factor)) +
  geom_point(size = 4) + geom_line() +
  theme_bw() + scale_color_manual(values = colors) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), limits = c(0.1, 0.45))

ggplot(to.ggplot.cart2,
       aes(x = timepoint, y = mean.ifn, color = anno.1.updated.factor,
           group = anno.1.updated.factor)) +
  geom_point(size = 4) + geom_line() +
  theme_bw() + scale_color_manual(values = colors) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), limits = c(0.1, 0.45))
dev.off()