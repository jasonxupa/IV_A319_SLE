setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/33_recalc_T_cell_proportions")


total.cells = Lupus_Reference_UMAP_Fixed.with.AUC.IFN@meta.data %>% group_by(patient.id, timepoint) %>% summarize(total = length(patient.id))
t.final.112044.cells.updated.anno.with.clonotypes@meta.data %>% group_by(patient.id, timepoint) %>% summarize(total = length(patient.id),
                                                                                                              cycling.T = sum(anno1.updated == "Cycling.T"),
                                                                                                              cycling.NK = sum(anno1.updated == "Cycling.NK"))

t.final.112044.cells.updated.anno.with.clonotypes@meta.data$anno1 %>% table()

t.final.112044.cells.updated.anno.with.clonotypes@meta.data$anno2 %>% table()



prop.changes = t.final.112044.cells.updated.anno.with.clonotypes@meta.data %>% group_by(patient.id, timepoint, dose, anno2) %>% summarize(n = length(patient.id)) %>% left_join(total.cells) %>%
  mutate(prop = n/total)



pdf("CyclingT.Cycling.NK.prop.use.graph4.pdf")

ggplot(prop.changes, aes(x = timepoint, y = prop)) + geom_boxplot() + stat_compare_means() + facet_wrap(~anno2, scales = "free")
ggplot(prop.changes %>% filter(anno2 %in% c("Cycling.CD4", "Cycling.CD8", "Cycling.NK")), 
       aes(x = timepoint, y = prop*100)) + geom_boxplot(outlier.size = -1) + geom_point(aes(shape = dose)) + stat_compare_means() + facet_wrap(~anno2, scales = "free") + theme_bw()
ggplot(prop.changes %>% filter(anno2 %in% c("Cycling.CD4", "Cycling.CD8", "Cycling.NK")), 
       aes(x = timepoint, y = prop*100)) + geom_boxplot(outlier.size = -1) + geom_point(aes(shape = dose), position = position_jitter(width = 0.1)) + stat_compare_means() + facet_wrap(~anno2, scales = "free") + theme_bw()
ggplot(prop.changes %>% filter(anno2 %in% c("Cycling.CD4", "Cycling.CD8", "Cycling.NK")), 
       aes(x = timepoint, y = prop*100, fill = timepoint)) + geom_boxplot(outlier.size = -1) + geom_point(aes(shape = dose), position = position_jitter(width = 0.1)) + stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) + facet_wrap(~anno2, scales = "free") + theme_bw()


ggplot(prop.changes %>% filter(anno2 %in% c("Cycling.CD4", "Cycling.CD8", "Cycling.NK")), 
       aes(x = timepoint, y = prop*100)) + geom_boxplot(outlier.size = -1) + geom_line(aes(group = patient.id)) + geom_point(aes(shape = dose)) + stat_compare_means() + facet_wrap(~anno2, scales = "free") + theme_bw()
dev.off()


#save the cycling object for upload to cellx gene
Cycling.T.NK = subset(t.final.112044.cells.updated.anno.with.clonotypes, anno2 %in% c("Cycling.CD4", "Cycling.CD8", "Cycling.NK"))

prop.changes2 = t.final.112044.cells.updated.anno.with.clonotypes@meta.data %>% group_by(patient.id, timepoint, dose, anno1) %>% summarize(n = length(patient.id)) %>% left_join(total.cells) %>%
  mutate(prop = n/total)

ggplot(prop.changes2, aes(x = timepoint, y = prop)) + geom_boxplot() + stat_compare_means() + facet_wrap(~anno1, scales = "free")




#global T-cell GEXP changes
# For CD4.CMT

saveRDS(Cycling.T.NK, file = "2_14_Cycling.T.NK.4860.cells.4_upload.rds")

D28.vs.base.Cycling.all = FindMarkers(Cycling.T.NK, ident.1 = "D28", ident.2 = "Base", group.by = "timepoint")
D28.vs.base.Cycling.all = D28.vs.base.Cycling.all %>% as_tibble(rownames = "gene")
D28.vs.base.Cycling.all = D28.vs.base.Cycling.all %>% filter(p_val_adj < 0.05)
pdf("Cycling.Base.vs.D28.pdf")
ggplot(D28.vs.base.Cycling.all, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = D28.vs.base.Cycling.all %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = D28.vs.base.Cycling.all %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = D28.vs.base.Cycling.all %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("CD4.CMT Base vs D28") +
  geom_text_repel(aes(label = gene), data = D28.vs.base.Cycling.all %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()

VlnPlot(Cycling.T.NK, features = c("LYZ", "CD8A", "GZMB", "PRF1"), group.by = "timepoint", stack = T, flip = T)



t.final.112044.cells.updated.anno.with.clonotypes
saveRDS(t.final.112044.cells.updated.anno.with.clonotypes, "t.final.112044.cells.updated.anno.with.clonotypes.rds")
D28.vs.base.T.all = FindMarkers(t.final.112044.cells.updated.anno.with.clonotypes, ident.1 = "D28", ident.2 = "Base", group.by = "timepoint")
D28.vs.base.T.all = D28.vs.base.T.all %>% as_tibble(rownames = "gene")
D28.vs.base.T.all = D28.vs.base.T.all %>% filter(p_val_adj < 0.05)
pdf("D28.vs.base.T.all.pdf", height = 20, width = 20)
ggplot(D28.vs.base.T.all, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = D28.vs.base.T.all %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = D28.vs.base.T.all %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = D28.vs.base.T.all %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("All T Base vs D28") +
  geom_text_repel(aes(label = gene), data = D28.vs.base.T.all %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()

VlnPlot(t.final.112044.cells.updated.anno.with.clonotypes, features = c("LYZ", "CD8A", "GZMB", "PRF1"), group.by = "timepoint", stack = T, flip = T)


pdf("CD3e.expression.pdf")
VlnPlot(subset(t.final.112044.cells.updated.anno.with.clonotypes, anno1 %in% c("aCD8", "aCD4", "NKT", "Cycling.T", "Cycling.NK")), 
        features = c("CD3E", "CD8A", "GZMB", "PRF1"), group.by = "anno1", 
        split.by = "timepoint", stack = T, flip = T)

VlnPlot(subset(t.final.112044.cells.updated.anno.with.clonotypes, anno1 %in% c("NKT", "Cycling.T", "Cycling.NK")), 
        features = c("CD3E", "CD8A", "GZMB", "PRF1"), group.by = "anno1", 
        split.by = "timepoint", stack = T, flip = T)


VlnPlot(subset(t.final.112044.cells.updated.anno.with.clonotypes, anno2 %in% c("GZMB aCD8", "GZMK aCD8", "aCD4", "NKT", "Cycling.CD4", "Cycling.CD8", "Cycling.NK", "NK")), 
        features = c("CD3E", "CD8A", "GZMB", "PRF1"), group.by = "anno2", 
        split.by = "timepoint", stack = T, flip = T)
dev.off()


#activation score
cycling.t.nk = subset(t.final.112044.cells.updated.anno.with.clonotypes, anno1 %in% c("aCD8", "aCD4", "NKT", "Cycling.T", "Cycling.NK"))

GOBP_GRANZYME_MEDIATED_PROGRAMMED_CELL_DEATH_SIGNALING_PATHWAY.v2026.1.Hs=  c("UBE4B","GSDME","GZMA","GZMB","GZMK","LAMP1","NKG7","PRF1","SRGN","GSDMB","BNIP3")
jx.activation = c("CD8A", "GZMB", "PRF1", "GZMK", "FASLG", "GZMH", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA", "CD160")
activation = c("GZMB", "PRF1", "GZMK", "FASLG", "GZMH", "GZMA", "GZMM")
checkpoints = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA", "CD160")

activation.signatures = list(jx.activation, GOBP_GRANZYME_MEDIATED_PROGRAMMED_CELL_DEATH_SIGNALING_PATHWAY.v2026.1.Hs)
names(activation.signatures) = c("jx.activation", "GO:0140507")



library(AUCell)
cells_rankings_vst <- AUCell_buildRankings(cycling.t.nk@assays$RNA$counts)

cells_AUC.activation<- AUCell_calcAUC(activation.signatures, cells_rankings_vst, 
                                aucMaxRank=nrow(cells_rankings_vst)*0.25, verbose = T)


cycling.t.nk[["Activation"]] <- CreateAssayObject(cells_AUC.activation@assays@data$AUC)

cycling.t.nk@assays$Activation

pdf("activation.sig.pdf")
cycling.t.nk$anno1 = factor(cycling.t.nk$anno1, levels = c("Cycling.T", "aCD4", "aCD8", "Cycling.NK", "NKT"))

VlnPlot(cycling.t.nk, 
        features = c("CD3E", "jx.activation", "GO:0140507"), group.by = "anno1", 
        split.by = "timepoint", stack = T, flip = T)
VlnPlot(cycling.t.nk, 
        features = c("CD3E", "jx.activation"), group.by = "anno1", 
        split.by = "timepoint", stack = T, flip = T)
dev.off()

#compare remnant plasma B at D28 to cohort 1 plasma Baseline

remnant.b = subset(b.final.2.no.doublet.4829.cells.w.clonotype, timepoint == "D28")
saveRDS(remnant.b, "remnant.b.rds")
remnant.b$patient.id %>% unique()
remnant.b$patient.id %>% table()
remnant.b$anno1 %>% table()


remnant.b.with.baseline.all = subset(b.final.2.no.doublet.4829.cells.w.clonotype, timepoint %in% c("Base", "D28") & anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)"))
saveRDS(remnant.b.with.baseline.all, "remnant.b.with.baseline.all.rds")


remnant.b.with.baseline = subset(b.final.2.no.doublet.4829.cells.w.clonotype, timepoint %in% c("Base", "D28") & 
                                   patient.id %in% c(2, 3, 6) & anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)"))

remnant.b.with.baseline.cohort1 = subset(b.final.2.no.doublet.4829.cells.w.clonotype, timepoint %in% c("Base", "D28") & 
                                   patient.id %in% c(1:6) & anno1 %in% c("IgA+/IgG+ PC", "PB (Cycling)"))

saveRDS(remnant.b.with.baseline, "remnant.b.with.baseline.rds")
saveRDS(remnant.b.with.baseline.cohort1, "remnant.b.with.baseline.cohort1.rds")


#three sets of DEG, whichever is better
D28.vs.base.plasma.pt.236 = FindMarkers(remnant.b.with.baseline.cohort1, ident.1 = "D28", ident.2 = "Base", group.by = "timepoint")
D28.vs.base.plasma.pt.236 = D28.vs.base.plasma.pt.236 %>% as_tibble(rownames = "gene")
D28.vs.base.plasma.pt.236 = D28.vs.base.plasma.pt.236 %>% filter(p_val_adj < 0.25)

pdf("remnant.b.changes.pdf")
ggplot(D28.vs.base.plasma.pt.236, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = D28.vs.base.plasma.pt.236 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = D28.vs.base.plasma.pt.236 %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = D28.vs.base.plasma.pt.236 %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Vs P 236 Plasma") +
  geom_text_repel(aes(label = gene), data = D28.vs.base.plasma.pt.236 %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

VlnPlot(remnant.b.with.baseline.cohort1, features = c("CD19", "CD38", "TNFRSF17", "CD27"), group.by = "timepoint", stack = T, flip = T)
VlnPlot(remnant.b.with.baseline.cohort1, features = c("CD19", "CD38", "TNFRSF17", "CD27"), group.by = "timepoint", split.by = "anno1", stack = T, flip = T)
VlnPlot(remnant.b.with.baseline.cohort1, features = c("CD19", "CD38", "TNFRSF17", "CD27"), group.by = "anno1", split.by = "timepoint", stack = T, flip = T)
dev.off()

D28.vs.base.plasma.all.cohort.1 = FindMarkers(remnant.b.with.baseline.cohort1, ident.1 = "D28", ident.2 = "Base", group.by = "timepoint")

D28.vs.base.plasma.all.pts = FindMarkers(remnant.b.with.baseline.all, ident.1 = "D28", ident.2 = "Base", group.by = "timepoint")
D28.vs.base.plasma.all.pts = D28.vs.base.plasma.all.pts %>% as_tibble(rownames = "gene")
D28.vs.base.plasma.all.pts = D28.vs.base.plasma.all.pts %>% filter(p_val_adj < 0.25)
D28.vs.base.plasma.all.pts %>% filter(gene == "CD19")

pdf("remnant.b.vsall.pts.pdf")
ggplot(D28.vs.base.plasma.all.pts, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = D28.vs.base.plasma.all.pts %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = D28.vs.base.plasma.all.pts %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = D28.vs.base.plasma.all.pts %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Vs all Plasma") +
  geom_text_repel(aes(label = gene), data = D28.vs.base.plasma.all.pts %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()



save.image("2_15_all_image.rdata")
