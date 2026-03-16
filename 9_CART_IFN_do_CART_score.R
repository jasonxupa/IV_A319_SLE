# CART_IFN_do_CART_score.R
# Purpose: Compute IFN gene signature AUCell scores specifically for the
#          T-cell compartment in CAR-T treated patients, compare pre- vs.
#          post-treatment scores across cell types and patients, and generate
#          violin/trend-line plots.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(AUCell)
library(ggpubr)
library(RColorBrewer)

# ==============================================================================
# SECTION 1: Subset CAR-T data to T cells only
# ==============================================================================
# Keep only cells whose predicted annotation matches a T-cell label present
# in the A-319 T-cell reference.

CART.t.only <- subset(
  CART.patient.scRNA.mt.10.filter.with.predictions.95k,
  predicted.anno1 %in% c(t.final.112044.cells.updated.anno.with.clonotypes$anno1)
)

# ---- Reclassify IFI+ CD4 cells as nCD4 for consistency ----------------------
CART.t.only@meta.data$predicted.anno1[
  CART.t.only$predicted.anno1 == "IFI+ CD4"
] <- "nCD4"
CART.t.only@meta.data$predicted.anno1.updated[
  CART.t.only$predicted.anno1.updated == "IFI+ CD4"
] <- "nCD4"

# ---- Set cell-type factor levels to match the A-319 T-cell reference ---------
CART.t.only$predicted.anno1.factor <- factor(
  CART.t.only$predicted.anno1,
  levels = levels(t.final.112044.cells.updated.anno.with.clonotypes$anno1.factor)
)

# ==============================================================================
# SECTION 2: Compute AUCell IFN scores for CAR-T T cells
# ==============================================================================

cells_rankings_vst <- AUCell_buildRankings(CART.t.only@assays$RNA$counts)

cells_AUC.IFN <- AUCell_calcAUC(molecular_signatures, cells_rankings_vst,
                                  aucMaxRank = nrow(cells_rankings_vst) * 0.25,
                                  verbose = TRUE)

CART.t.only[["AUC_IFN"]] <- CreateAssayObject(cells_AUC.IFN@assays@data$AUC)

saveRDS(CART.t.only, "CART.t.only.with.AUC.IFN.rds")

# ==============================================================================
# SECTION 3: Violin plots of IFN scores by cell type and timepoint
# ==============================================================================

CART.t.only$timepoint <- factor(CART.t.only$timepoint, levels = c("pre", "post"))

pdf("CART_IFN_expression_AUC.pdf")

VlnPlot(CART.t.only, group.by = "predicted.anno1.factor",
        features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint")

VlnPlot(subset(CART.t.only,
               predicted.anno1.updated.factor %in% c("Cycling", "Treg")),
        group.by = "predicted.anno1.factor", features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint")

VlnPlot(CART.t.only, group.by = "predicted.anno1.updated.factor",
        features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint")

VlnPlot(CART.t.only, group.by = "patient", features = mole.sig.keep,
        stack = TRUE, flip = TRUE, fill.by = "ident", split.by = "timepoint")

dev.off()

# ==============================================================================
# SECTION 4: Prepare AUC scores for statistical comparison
# ==============================================================================

# ---- Merge AUC scores into metadata -----------------------------------------
CART.t.only@meta.data <- cbind(CART.t.only@meta.data,
                                t(CART.t.only@assays$AUC_IFN@counts))

# ---- Parse patient and timepoint from sample labels --------------------------
CART.t.only$timepoint <- str_split_fixed(CART.t.only$pt.timepoint,
                                          pattern = "\\\\.", n = 2)[, 2]
CART.t.only$patient   <- str_split_fixed(CART.t.only$pt.timepoint,
                                          pattern = "\\\\.", n = 2)[, 1]

# ---- Pivot scores to long format for ggplot ----------------------------------
scores.auc <- CART.t.only@meta.data %>%
  dplyr::select(predicted.anno1.factor, timepoint, patient,
                `Bennett.8`:`M5.12.67`) %>%
  pivot_longer(cols = `Bennett.8`:`M5.12.67`,
               names_to = "score.name", values_to = "score.value")

scores.auc <- scores.auc %>% filter(score.name %in% mole.sig.keep)
scores.auc$timepoint <- factor(scores.auc$timepoint, levels = c("pre", "post"))

# ==============================================================================
# SECTION 5: Statistical comparison of IFN scores (violin plots with p-values)
# ==============================================================================

pdf("CART_scores.pval.t.pdf", width = 20, height = 20)

ggplot(scores.auc, aes(x = timepoint, y = score.value)) +
  facet_grid(predicted.anno1.factor ~ score.name) +
  geom_violin() +
  stat_compare_means(comparisons = list(c("Base", "D28"),
                                         c("D28", "M3"),
                                         c("Base", "M3"))) +
  theme_bw()

ggplot(scores.auc %>% filter(timepoint != "D28"),
       aes(x = timepoint, y = score.value)) +
  facet_grid(score.name ~ predicted.anno1.factor) +
  geom_violin() +
  stat_compare_means(comparisons = list(c("Base", "M3")), label.y = 1) +
  theme_bw() + expand_limits(y = 1.3)

ggplot(scores.auc %>%
         filter(predicted.anno1.factor %in% c("cycling T", "Treg Eff-like",
                                                "Treg Fr-like")),
       aes(x = timepoint, y = score.value)) +
  facet_grid(score.name ~ dose) +
  geom_violin() +
  stat_compare_means(comparisons = list(c("Base", "D28"),
                                         c("D28", "M3"),
                                         c("Base", "M3")),
                     label.y = 1) +
  theme_bw() + expand_limits(y = 1.5) + ggtitle("Cycling T + Treg")

dev.off()

# ==============================================================================
# SECTION 6: Trend-line plots of mean IFN scores by cell type
# ==============================================================================

mean.scores <- scores.auc %>%
  group_by(score.name, predicted.anno1.factor, timepoint) %>%
  summarize(mean.score = mean(score.value), .groups = "drop")

colors <- c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
mean.scores$timepoint <- factor(mean.scores$timepoint, levels = c("pre", "post"))

pdf("CART_IFN.t.trend.line.mean.auc.pdf")
ggplot(mean.scores,
       aes(y = mean.score, x = timepoint,
           group = predicted.anno1.factor, color = predicted.anno1.factor)) +
  geom_point() + geom_line() + facet_wrap(~score.name) +
  scale_color_manual(values = colors) + theme_bw() +
  stat_compare_means(comparisons = list(c("pre", "post")), paired = TRUE)

ggplot(mean.scores %>% filter(timepoint != "D28"),
       aes(y = mean.score, x = timepoint,
           group = predicted.anno1.factor, color = predicted.anno1.factor)) +
  geom_point() + geom_line() + facet_wrap(~score.name) +
  scale_color_manual(values = colors) + theme_bw() +
  stat_compare_means(comparisons = list(c("pre", "post")), paired = TRUE)
dev.off()

# ==============================================================================
# SECTION 7: Per-patient trend-line plots
# ==============================================================================

mean.scores.pt <- scores.auc %>%
  group_by(score.name, patient, timepoint) %>%
  summarize(mean.score = mean(score.value), .groups = "drop")

pdf("IFN.t.trend.line.mean.auc.pt.pdf")
ggplot(mean.scores.pt %>% filter(timepoint != "D28"),
       aes(y = mean.score, x = timepoint, group = patient, color = patient)) +
  geom_point() + geom_line() + facet_wrap(~score.name) +
  theme_bw() +
  stat_compare_means(comparisons = list(c("Base", "M3")), paired = TRUE) +
  ggtitle("IFN change in T-cells by patient")
dev.off()