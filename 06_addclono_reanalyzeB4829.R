
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/12_reanalyze_b_add_clono")

b.final = subset(all.patients.w.mdata, b.clone == T | anno.1.updated %in% c("Plasma", "A-NBC & csNBC", "NBC/A-NBC/csNBC"))
b.final = b.final %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20) 
b.final = FindNeighbors(b.final, reduction = 'pca', dims = 1:25)

b.final = b.final %>% FindClusters(resolution = 1)
b.final = b.final %>% FindClusters(resolution = 2)

DimPlot(b.final, group.by = "RNA_snn_res.1", label = T) + NoLegend()
DimPlot(b.final, group.by = "anno.1.updated", label = T) + NoLegend()



b.final.2.no.doublet = subset(b.final, RNA_snn_res.1 %in% c(10, 2, 0, 5, 7, 4, 1, 3))

DimPlot(b.final.2.no.doublet, group.by = "anno2", label = T) + NoLegend()
b.final.2.no.doublet



b.final.2.no.doublet = b.final.2.no.doublet %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20) 
b.final.2.no.doublet = b.final.2.no.doublet %>% FindClusters(resolution = c(0.5, 1))

DimPlot(b.final.2.no.doublet, group.by = "RNA_snn_res.0.5", label = T) + NoLegend() +
DimPlot(b.final.2.no.doublet, group.by = "anno2", label = T) + NoLegend()



b.final.2.no.doublet@meta.data  = 
  b.final.2.no.doublet@meta.data %>% mutate(anno1 = 
                                 case_when(RNA_snn_res.0.5 == 0 ~ "NBC/Activated NBC",
                                           RNA_snn_res.0.5 == 1 ~ "Activated NBC",
                                           RNA_snn_res.0.5 == 3 ~ "NBC",
                                           RNA_snn_res.0.5 == 2 ~ "IgA+/IgG+ PC",
                                           RNA_snn_res.0.5 == 4 ~ "csNBC",
                                           RNA_snn_res.0.5 == 5 ~ "PB (Cycling)"))


DimPlot(b.final.2.no.doublet, group.by = "anno1", label = T) + NoLegend()
b.final.2.no.doublet$anno1.factor = factor(b.final.2.no.doublet$anno1, levels = c("NBC", "NBC/Activated NBC", "Activated NBC", "csNBC", "IgA+/IgG+ PC", "PB (Cycling)"))

DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", label = T, cols = custom_palette) + NoLegend()


b.summ.stats2 = b.only.jx.anno@meta.data %>% group_by(timepoint, anno1.factor, .drop = F) %>% summarize(n = length(patient)) %>% 
  left_join(b.only.jx.anno@meta.data %>% group_by(timepoint) %>% summarize(n.total = length(patient))) %>% mutate(prop = n/n.total)


custom_palette <- c(
  "NBC" = "#1F77B4",                # Blue for NBC
  "NBC/Activated NBC" = "#3B9AE1",  # Blend of Blue and Orange for NBC/Activated NBC
  "Activated NBC" = "#FFBB78",      # Orange for Activated NBC
  "csNBC" = "#FF7F0E",              # Light Orange for csNBC
  "IgA+/IgG+ PC" = "#2CA02C",       # Green for IgA+/IgG+ PC
  "PB (Cycling)" = "#A6D8B8"        # Light Green for PB (Cycling)
)


pdf("B.updated.umap.pdf")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", label = T,
        cells.highlight = b.final.2.no.doublet$cell.name[!b.final.2.no.doublet$cell.name %in% b.only.jx.anno.with.AUCell$cell.name] )

DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", label = T, cols = custom_palette)

dev.off()



#saveRDS(b.final.2.no.doublet, "b.final.2.no.doublet.4829.cells.rds")



b.final.2.no.doublet <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/12_reanalyze_b_add_clono/b.final.2.no.doublet.4829.cells.rds")

b.clonotype.files = list.files("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/09_BCR_TCR/TCR and BCR", recursive = T, all.files = T, full.names = T) %>% grep(pattern = "vdj_b", value = T) %>% grep(pattern = "filtered_contig_annotations.csv", value = T)
#import all these files

b.clonotype.list = list()

for(file in b.clonotype.files){
  
  b.clonotype.list[[file]] = read_csv(file)
  b.clonotype.list[[file]]$file.name = file
  
}

b.clonotypes = do.call(rbind, b.clonotype.list)

lapply(X = b.clonotype.list, nrow) %>% unlist() %>% sum()
 
b.clonotypes$file.name %>% str_split_fixed(pattern = "/", n=15)
b.clonotypes = b.clonotypes %>% mutate(sample.timepoint = str_split_fixed(file.name, pattern = "/", n=15)[,13],
                                       patient = str_split_fixed(sample.timepoint, pattern = " ", n=2)[,1],
                                       timepoint =  str_split_fixed(sample.timepoint, pattern = " ", n=2)[,2])


b.clonotypes = b.clonotypes %>% mutate(cell.name = paste0(patient, "_", timepoint, "_", barcode))

table(b.clonotypes$cell.name %>% duplicated())

b.clonotypes.filter = b.clonotypes %>% filter(cell.name %in% b.final.2.no.doublet$cell.name)


b.clonotypes.filter[b.clonotypes.filter$cell.name %>% duplicated(),] %>% select(cell.name)

b.clonotypes.filter %>% filter(cell.name == "S01_Base_AAGACCTGTCGTCTTC-1")
clonotype.ids = b.clonotypes.filter %>% select(cell.name, raw_clonotype_id) %>% unique()
clonotype.ids = b.clonotypes.filter %>% select(cell.name, raw_clonotype_id) %>% na.omit() %>% unique()

clonotype.ids[clonotype.ids$cell.name %>% duplicated(),]



b.final.2.no.doublet@meta.data = b.final.2.no.doublet@meta.data %>% left_join(clonotype.ids)

rownames(b.final.2.no.doublet@meta.data) = b.final.2.no.doublet@meta.data$cell.name


b.clonotypes.aa <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/09_BCR_TCR/b.clonotypes.rds")

b.clonotypes.aa = b.clonotypes.aa %>% rename(raw_clonotype_id = clonotype_id)%>% select(-sample)

b.final.2.no.doublet@meta.data  = b.final.2.no.doublet@meta.data %>% left_join(b.clonotypes.aa )
rownames(b.final.2.no.doublet@meta.data) = b.final.2.no.doublet@meta.data$cell.name

saveRDS(b.final.2.no.doublet, "b.final.2.no.doublet.4829.cells.w.clonotype.rds")

b.final.2.no.doublet@meta.data %>% filter(patient == "S08", chain.1 == "IGH:CARGQYDSGSASEYW")


pdf("top10.clonotypes.pdf")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", 
        cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$proportion > 0.01 & b.final.2.no.doublet$timepoint == "Base" ], 
        raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle(">1% at base")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 5 & b.final.2.no.doublet$timepoint == "Base" ], raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 5 clonotypes at base")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 25 & b.final.2.no.doublet$timepoint == "Base" ], raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 10 clonotypes at base")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 25 & b.final.2.no.doublet$timepoint == "Base"], raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 25 clonotypes at base")

DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", 
        cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$proportion > 0.01 & b.final.2.no.doublet$timepoint == "M3" ], 
        raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle(">1% at M3")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 5 & b.final.2.no.doublet$timepoint == "M3" ], 
        raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 5 clonotypes at M3")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 10 & b.final.2.no.doublet$timepoint == "M3" ], 
        raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 10 clonotypes at M3")


DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 25 & b.final.2.no.doublet$timepoint == "M3" ], 
        raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 25 clonotypes at M3")

DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$chain.1 == "IGH:CARGQYDSGSASEYW"], raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("S08 IGH:CARGQYDSGSASEYW")

dev.off()


pdf("top25clones.pdf")
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cells.highlight = b.final.2.no.doublet$cell.name[b.final.2.no.doublet$clonotype_num <= 25 & b.final.2.no.doublet$timepoint == "Base"], raster = T, raster.dpi = c(800,800), pt.size = 5) + ggtitle("top 25 clonotypes at base")
dev.off()


pdf("top.clonotype.prop.pdf")
b.final.2.no.doublet@meta.data %>% filter(timepoint == "Base") %>% 
  dplyr::group_by(proportion > 0.01, anno1.factor) %>% summarize(n=length(patient)) %>% na.omit() %>% 
  left_join(b.final.2.no.doublet@meta.data  %>% filter(timepoint == "Base") %>%  dplyr::group_by(proportion > 0.01) %>% summarize(total=length(patient)) %>% na.omit()) %>% mutate(prop = n/total) %>%
  ggplot(aes(x = `proportion > 0.01`, y = prop, fill = anno1.factor)) + 
  geom_col() + geom_text(aes(label = paste0(n, "\n (", prop %>% round(2), ")")), position = position_stack(vjust = 0.5)) + theme_bw() + scale_fill_manual(values = custom_palette) + ggtitle(">1% @ base")


b.final.2.no.doublet@meta.data %>% filter(timepoint == "M3") %>% 
  dplyr::group_by(proportion > 0.01, anno1.factor) %>% summarize(n=length(patient)) %>% na.omit() %>% 
  left_join(b.final.2.no.doublet@meta.data  %>% filter(timepoint == "M3") %>%  dplyr::group_by(proportion > 0.01) %>% summarize(total=length(patient)) %>% na.omit()) %>% mutate(prop = n/total) %>%
  ggplot(aes(x = `proportion > 0.01`, y = prop, fill = anno1.factor)) + 
  geom_col() + geom_text(aes(label = paste0(n, "\n (", prop %>% round(2), ")")), position = position_stack(vjust = 0.5)) + theme_bw() + scale_fill_manual(values = custom_palette) + ggtitle(">1% @ M3")



b.final.2.no.doublet@meta.data  %>% filter(timepoint == "Base") %>% dplyr::group_by(clonotype_num <= 25, anno1.factor) %>% summarize(n=length(patient)) %>% na.omit() %>% 
  left_join(b.final.2.no.doublet@meta.data  %>% filter(timepoint == "Base")  %>% dplyr::group_by(clonotype_num <= 25) %>% summarize(total=length(patient)) %>% na.omit()) %>% mutate(prop = n/total) %>%
  ggplot(aes(x = `clonotype_num <= 25`, y = prop, fill = anno1.factor)) + geom_col() + geom_text(aes(label = paste0(n, "\n (", prop %>% round(2), ")")), position = position_stack(vjust = 0.5)) + theme_bw() + 
  scale_fill_manual(values = custom_palette) + ggtitle("top 25 clones @ base")

b.final.2.no.doublet@meta.data %>% filter(timepoint == "Base") %>% dplyr::group_by(clonotype_num <= 10, anno1.factor) %>% summarize(n=length(patient)) %>% na.omit() %>% 
  left_join(b.final.2.no.doublet@meta.data %>% filter(timepoint == "Base")  %>% dplyr::group_by(clonotype_num <= 10) %>% summarize(total=length(patient)) %>% na.omit()) %>% mutate(prop = n/total) %>%
  ggplot(aes(x = `clonotype_num <= 10`, y = prop, fill = anno1.factor)) + geom_col() + geom_text(aes(label = paste0(n, "\n (", prop %>% round(2), ")")), position = position_stack(vjust = 0.5)) + 
  theme_bw() + scale_fill_manual(values = custom_palette)   + ggtitle("top 10 clones @ base")


b.final.2.no.doublet@meta.data  %>% filter(timepoint == "M3") %>% dplyr::group_by(clonotype_num <= 25, anno1.factor) %>% summarize(n=length(patient)) %>% na.omit() %>% 
  left_join(b.final.2.no.doublet@meta.data %>% filter(timepoint == "M3")  %>% dplyr::group_by(clonotype_num <= 25) %>% summarize(total=length(patient)) %>% na.omit()) %>% mutate(prop = n/total) %>%
  ggplot(aes(x = `clonotype_num <= 25`, y = prop, fill = anno1.factor)) + geom_col() + geom_text(aes(label = paste0(n, "\n (", prop %>% round(2), ")")), position = position_stack(vjust = 0.5)) + 
  theme_bw() + scale_fill_manual(values = custom_palette) +  ggtitle("top 25 clones @ m3")

b.final.2.no.doublet@meta.data %>% filter(timepoint == "M3") %>% dplyr::group_by(clonotype_num <= 10, anno1.factor) %>% summarize(n=length(patient)) %>% na.omit() %>% 
  left_join(b.final.2.no.doublet@meta.data  %>% filter(timepoint == "M3") %>% dplyr::group_by(clonotype_num <= 10) %>% summarize(total=length(patient)) %>% na.omit()) %>% mutate(prop = n/total) %>%
  ggplot(aes(x = `clonotype_num <= 10`, y = prop, fill = anno1.factor)) + geom_col() + geom_text(aes(label = paste0(n, "\n (", prop %>% round(2), ")")), position = position_stack(vjust = 0.5)) + 
  theme_bw() + scale_fill_manual(values = custom_palette) +
  ggtitle("top 10 clones @ M3")

ggplot(b.final.2.no.doublet@meta.data %>% filter(proportion > 0.1), aes(x = patient, fill = anno1.factor)) + geom_bar() + theme_bw() + ggtitle("prop > 0.1") + scale_fill_manual(values = custom_palette) + facet_wrap(~timepoint)
ggplot(b.final.2.no.doublet@meta.data %>% filter(proportion > 0.05), aes(x = patient, fill = anno1.factor)) + geom_bar() + theme_bw() + ggtitle("prop > 0.05")+ scale_fill_manual(values = custom_palette) 
ggplot(b.final.2.no.doublet@meta.data %>% filter(clonotype_num <= 5), aes(x = patient, fill = anno1.factor)) + geom_bar() + theme_bw() + ggtitle("clonotype_num <= 5")+ scale_fill_manual(values = custom_palette) +facet_wrap(~timepoint)
ggplot(b.final.2.no.doublet@meta.data %>% filter(clonotype_num <= 10), aes(x = patient, fill = anno1.factor)) + geom_bar() + theme_bw() + ggtitle("clonotype_num <= 10")+ scale_fill_manual(values = custom_palette) + facet_wrap(~timepoint)
ggplot(b.final.2.no.doublet@meta.data %>% filter(clonotype_num <= 25), aes(x = patient, fill = anno1.factor)) + geom_bar() + theme_bw() + ggtitle("clonotype_num <= 25")+ scale_fill_manual(values = custom_palette) + facet_wrap(~timepoint)


dev.off()











b.final.2.no.doublet@meta.data %>% filter(proportion > 0.1)

c6.genes$gene %>% head(50) %>% cat()



b.final.2.no.doublet$anno1.factor = factor(b.final.2.no.doublet$anno1, levels = c("NBC", "NBC/Activated NBC", "Activated NBC", "csNBC", "IgA+/IgG+ PC", "PB (Cycling)"))
b.final.2.no.doublet$timepoint = factor(b.final.2.no.doublet$timepoint, levels = c("Base", "D28", "M3"))
b.final.2.no.doublet$patient = factor(b.final.2.no.doublet$patient)

b.final.2.no.doublet$patient %>% table()

b.summ.stats = b.final.2.no.doublet@meta.data %>% dplyr::group_by(patient, timepoint, anno1.factor, .drop = F) %>% dplyr::summarize(n = length(patient), .drop = F) %>% 
  left_join(b.final.2.no.doublet@meta.data %>% dplyr::group_by(patient, timepoint, .drop = F) %>% dplyr::summarize(n.total = length(patient), .drop = F)) %>% mutate(prop = n/n.total)

b.summ.stats$prop[b.summ.stats$n==0] = 0

library(ggpubr)





b.summ.stats2 = b.final.2.no.doublet@meta.data %>% dplyr::group_by(timepoint, anno1.factor, .drop = F) %>% summarize(n = length(patient)) %>% 
  left_join(b.final.2.no.doublet@meta.data %>% group_by(timepoint) %>% summarize(n.total = length(patient))) %>% mutate(prop = n/n.total)
b.summ.stats2$anno1.factor = factor(b.summ.stats2$anno1.factor, levels = c("NBC", "NBC/Activated NBC", "Activated NBC", "csNBC", "IgA+/IgG+ PC", "PB (Cycling)"))


custom_palette <- c(
  "NBC" = "#1F77B4",                # Blue for NBC
  "NBC/Activated NBC" = "#3B9AE1",  # Blend of Blue and Orange for NBC/Activated NBC
  "Activated NBC" = "#FFBB78",      # Orange for Activated NBC
  "csNBC" = "#FF7F0E",              # Light Orange for csNBC
  "IgA+/IgG+ PC" = "#2CA02C",       # Green for IgA+/IgG+ PC
  "PB (Cycling)" = "#A6D8B8"        # Light Green for PB (Cycling)
)


pdf("b.cell.changes.over.tx2.pdf")
ggplot(b.summ.stats, aes(x = timepoint, y = prop)) + geom_boxplot() + geom_point(position = position_jitter(width = 0.1)) + facet_wrap(~anno1.factor) + stat_compare_means() + theme_bw()
ggplot(b.summ.stats, aes(x = timepoint, y = prop)) + geom_boxplot(outlier.size = -1) + geom_point(position = position_jitter(width = 0.1)) + facet_wrap(~anno1.factor) + stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) + theme_bw()
ggplot(b.summ.stats2, aes(x = timepoint, y = prop, fill = anno1.factor)) + geom_col() + theme_bw() + scale_fill_manual(values = custom_palette) 
dev.off()


pdf("b.cell.umap.pdf", height = 8, width = 8)
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", cols = custom_palette, label = T, raster = T, raster.dpi = c(800,800), pt.size = 2)
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", split.by = "timepoint", cols = custom_palette, label = T, raster = T, raster.dpi = c(800,800),pt.size = 2) + theme_classic() + theme(legend.position = "bottom")
dev.off()



pdf("b.cell.umap2.pdf", height = 8, width = 15)
DimPlot(b.final.2.no.doublet, group.by = "anno1.factor", split.by = "timepoint", cols = custom_palette, label = T, raster = T, raster.dpi = c(800,800),pt.size = 2) + theme_classic() + theme(legend.position = "bottom")
dev.off()

pdf("b.cell.changes.over.tx3.pdf")
ggplot(b.summ.stats2, aes(x = timepoint, y = prop, fill = anno1.factor)) + geom_col() + theme_bw() + scale_fill_manual(values = custom_palette) + geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) + coord_flip()
dev.off()


# Create a color palette with 3 shades of blue
blue.colors <- brewer.pal(3, "PuBuGn") %>% rev()


pdf("b.cell.changes.over.tx2.narrow.pdf")
ggplot(b.summ.stats, aes(x = timepoint, y = prop, fill = timepoint)) + geom_boxplot(outlier.size = -1) + geom_point(position = position_jitter(width = 0.1)) + facet_wrap(~anno1.factor, nrow = 1) + stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) + theme_bw() +  scale_fill_manual(values = blue.colors)
dev.off()



#b.final.2.no.doublet@meta.data$anno1.factor = factor(b.final.2.no.doublet@meta.data$anno1.factor, levels = c("NBC", "NBC/Activated NBC", "Activated NBC", "csNBC", "IgA+/IgG+ PC", "PB (Cycling)"))


b.final.2.no.doublet@meta.data <- b.final.2.no.doublet@meta.data %>%
  mutate(anno1.factor.binned = recode_factor(
    anno1.factor,
    "NBC" = "NBC/A-NBC Int",
    "NBC/Activated NBC" = "NBC/A-NBC Int",
    "Activated NBC" = "A-NBC & csNBC",
    "csNBC" = "A-NBC & csNBC",
    "IgA+/IgG+ PC" = "Plasma",
    "PB (Cycling)" = "Plasma"
  ))

b.summ.stats.binned <- b.final.2.no.doublet@meta.data %>%
  group_by(patient, timepoint, anno1.factor.binned, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(b.final.2.no.doublet@meta.data %>%
              group_by(patient, timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

# Create the plot
ggplot(b.summ.stats.binned, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~anno1.factor.binned, nrow = 1) +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) +
  theme_bw() +
  scale_fill_manual(values = blue.colors)

ggplot(b.summ.stats, aes(x = timepoint, y = prop, fill = timepoint)) + geom_boxplot(outlier.size = -1) + geom_point(position = position_jitter(width = 0.1)) + facet_wrap(~anno1.factor, nrow = 1) + stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) + theme_bw() +  scale_fill_manual(values = blue.colors)

pdf("CD19.and.other.expression.pdf")
VlnPlot(b.final.2.no.doublet, group.by = "anno1.factor", features = c("CD19"), fill.by = "ident", cols = custom_palette, pt.size = 0) + NoLegend()
VlnPlot(b.final.2.no.doublet, group.by = "anno1.factor", features = c("CD19", "IGHD", "IGHM", "IGHA1", "IGHG1", "JCHAIN",  "CD69", 'FCRL5',  "CD27", "TNFRSF13B", "TNFRSF17", "MS4A1"), stack = T, flip = T, fill.by = "ident", cols = custom_palette) + NoLegend()
dev.off()

cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "IgA+/IgG+ PC", group.by = "anno1.factor", features = "CD19")
cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "PB (Cycling)", group.by = "anno1.factor", features = "CD19")

cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "csNBC", group.by = "anno1.factor", features = "CD19")


cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "IgA+/IgG+ PC", ident.2 = "NBC", group.by = "anno1.factor", features = "CD19") %>% print()
cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "IgA+/IgG+ PC", ident.2 = "NBC", group.by = "anno1.factor", features = "CD19") %>% print()

cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "PB (Cycling)",  ident.2 = "NBC", group.by = "anno1.factor", features = "CD19") %>% print()

cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "csNBC",  ident.2 = "NBC",group.by = "anno1.factor", features = "CD19")

b.final.2.no.doublet$anno1.factor %>% table()

FindMarkers(b.final.2.no.doublet, ident.1 = "NBC/Activated NBC",  ident.2 = "NBC",group.by = "anno1.factor", features = "CD19")
FindMarkers(b.final.2.no.doublet, ident.1 = "NBC/Activated NBC",  ident.2 = "NBC",group.by = "anno1.factor", features = "CD19", logfc.threshold = 0)
FindMarkers(b.final.2.no.doublet, ident.1 = "Activated NBC",  ident.2 = "NBC",group.by = "anno1.factor", features = "CD19", logfc.threshold = 0)

cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "csNBC", group.by = "anno1.factor", features = "CD19")
cd19.change = FindMarkers(b.final.2.no.doublet, ident.1 = "csNBC", group.by = "anno1.factor", features = "CD19")

