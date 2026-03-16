myeloid.in.t_6271_cells <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/06_analyze_T/myeloid.in.t_6271_cells.rds")
m.cells.in.b_364 <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/05_Analyze_B/m.cells.in.b_364.rds")
Myeloid_177754_cells <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/04_split_combined_seurat_by_traj/Myeloid_177754_cells.rds")
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/07_analyze_myeloid")


#Myeloid_177754_cells = merge(Myeloid_177754_cells, myeloid.in.t_6271_cells)
#Myeloid_177754_cells = merge(Myeloid_177754_cells, m.cells.in.b_364)


Myeloid_177754_cells = JoinLayers(Myeloid_177754_cells)
Myeloid_177754_cells$predicted.annotation_20230508 %>% table() %>% sort()
Myeloid_177754_cells$predicted.cell.type.short %>% table() %>% sort()


non.neut.mono = subset(Myeloid_177754_cells, !predicted.annotation_20230508 %in% c("Monocytes", "Neutrophils"))


non.neut.mono = non.neut.mono %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)
non.neut.mono = non.neut.mono %>% FindNeighbors() %>% FindClusters(resolution = c(0.5, 1, 2))

DimPlot(non.neut.mono, group.by = "predicted.annotation_20230508", label = T) + NoLegend()
DimPlot(non.neut.mono, group.by = "RNA_snn_res.1", label = T)
DimPlot(non.neut.mono, group.by = "RNA_snn_res.0.5", label = T) + NoLegend()

non.neut.mono$predicted.annotation_20230508



DimPlot(non.neut.mono, group.by = "RNA_snn_res.1", label = T) + NoLegend()
FeaturePlot(non.neut.mono, c("CD4", "CD8A", "CD19"), label = T) + NoLegend()
FeaturePlot(non.neut.mono, c("CD4", "CD8A", "IGHM"), label = T) + NoLegend()

non.neut.mono@meta.data %>%
  group_by(RNA_snn_res.1, predicted.annotation_20230508) %>%
  summarize(count = n(), .groups = 'drop') %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Create the cell.types string
  summarize(cell.types = paste(predicted.annotation_20230508, count, sep = " (", collapse = "), "), ")") %>%
  ungroup() %>% print(n=30)



myeloid_tonsil <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/02_Project_to_tonsil/scRNA-seq/20230911_myeloid_seurat_obj.rds")

PDC_tonsil <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/02_Project_to_tonsil/scRNA-seq/20230911_PDC_seurat_obj.rds")
#myeloid_tonsil = merge(myeloid_tonsil , PDC_tonsil)

myeloid_tonsil = myeloid_tonsil %>% FindVariableFeatures()
myeloid_tonsil = myeloid_tonsil %>% ScaleData() %>% RunPCA(npcs = 25)
saveRDS(myeloid_tonsil, "myeloid.tonsil.joined.pdc.rds")

patient.anchors.tonsil.m <- FindTransferAnchors(
  reference = myeloid_tonsil,
  query = non.neut.mono,
  k.filter = NA,
  reference.reduction = "pca", 
  #reference.neighbors = "pca.annoy.neighbors", 
  dims = 1:20
)

tonsil.predictions.m  <- TransferData(
  anchorset = patient.anchors.tonsil.m, 
  query = non.neut.mono,
  reference = myeloid_tonsil, 
  refdata = list(annotation_figure_1_M.only = "annotation_figure_1",
                 annotation_20230508_M.only = "annotation_20230508"))


tonsil.predictions.m@meta.data

saveRDS(tonsil.predictions.m, "tonsil.predictions.m.withPDC.rds")




tonsil.predictions.m@meta.data %>%
  group_by(RNA_snn_res.1, predicted.annotation_20230508_M.only) %>%
  summarize(count = n(), .groups = 'drop') %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string
  summarize(cell.types = paste(predicted.annotation_20230508_M.only, count, sep = " (", collapse = "), "), ")") %>%
  ungroup() %>% 
  print(n = 30)


tonsil.predictions.m@meta.data %>%
  group_by(RNA_snn_res.1, predicted.annotation_20230508) %>%
  summarize(count = n(), .groups = 'drop') %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string
  summarize(cell.types = paste(predicted.annotation_20230508, count, sep = " (", collapse = "), "), ")") %>%
  ungroup() %>% 
  print(n = 30)


tonsil.predictions.m@meta.data %>%
  group_by(RNA_snn_res.1, predicted.cell.type.short) %>%
  summarize(count = n(), .groups = 'drop') %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string
  summarize(cell.types = paste(predicted.cell.type.short, count, sep = " (", collapse = "), "), ")") %>%
  ungroup() %>% 
  print(n = 30)

tonsil.predictions.m@meta.data %>%
  group_by(RNA_snn_res.1, predicted.cell.type.short) %>%
  summarize(
    count = n(),
    mean_score = mean(predicted.cell.type.short.score, na.rm = TRUE),  # Calculate mean score
    .groups = 'drop'
  ) %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string with count and mean score included
  summarize(
    cell.types = paste(predicted.cell.type.short, count, 
                       round(mean_score, 2),  # Round mean score to 2 decimal places
                       sep = " (", collapse = "), "), ")"
  ) %>%
  ungroup() %>% 
  print(n = 30)



tonsil.predictions.m@meta.data %>%
  group_by(RNA_snn_res.1, predicted.annotation_20230508_M.only) %>%
  summarize(
    count = n(),
    mean_score = mean(predicted.annotation_20230508_M.only.score, na.rm = TRUE),  # Calculate mean score
    .groups = 'drop'
  ) %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string with count and mean score included
  summarize(
    cell.types = paste(predicted.annotation_20230508_M.only, count, 
                       round(mean_score, 2),  # Round mean score to 2 decimal places
                       sep = " (", collapse = "), "), ")"
  ) %>%
  ungroup() %>% 
  print(n = 30)



tonsil.predictions.m@meta.data %>%
  group_by(RNA_snn_res.1, predicted.annotation_20230508) %>%
  summarize(
    count = n(),
    mean_score = mean(predicted.annotation_20230508.score, na.rm = TRUE),  # Calculate mean score
    .groups = 'drop'
  ) %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string with count and mean score included
  summarize(
    cell.types = paste(predicted.annotation_20230508, count, 
                       round(mean_score, 2),  # Round mean score to 2 decimal places
                       sep = " (", collapse = "), "), ")"
  ) %>%
  ungroup() %>% 
  print(n = 30)





tonsil.predictions.m@meta.data <- 
  tonsil.predictions.m@meta.data %>% mutate(cell.type.jx.1 = 
                                       case_when(
                                         RNA_snn_res.1 == 0 ~ "Neut",
                                         RNA_snn_res.1 == 1 ~ "Neut",
                                         RNA_snn_res.1 == 2 ~ "EffectorT",
                                         RNA_snn_res.1 == 3 ~ "Platelet",
                                         RNA_snn_res.1 == 4 ~ "Mono",
                                         RNA_snn_res.1 == 5 ~ "Platelet",
                                         RNA_snn_res.1 == 6 ~ "Mono",
                                         RNA_snn_res.1 == 7 ~ "DC2",
                                         RNA_snn_res.1 == 8 ~ "Macrophage",
                                         RNA_snn_res.1 == 9 ~ "Neutrophil/Cycling",
                                         RNA_snn_res.1 == 10 ~ "EffectorT",
                                         RNA_snn_res.1 == 11 ~ "MMP8+ OLR1+ Neut",
                                         RNA_snn_res.1 == 12 ~ "NaiveT",
                                         RNA_snn_res.1 == 13 ~ "Platelet",
                                         RNA_snn_res.1 == 14 ~ "Platelet",
                                         RNA_snn_res.1 == 15 ~ "Macrophage",
                                         RNA_snn_res.1 == 16 ~ "EffectorT",
                                         RNA_snn_res.1 == 17 ~ "PDC",
                                         RNA_snn_res.1 == 18 & predicted.cell.type.short != "Mono" ~ "EffectorT",
                                         RNA_snn_res.1 == 18 & predicted.cell.type.short == "Mono" ~ "Mono",
                                         RNA_snn_res.1 == 19 ~ "RBC",
                                         RNA_snn_res.1 == 20 ~ "MEP",
                                         RNA_snn_res.1 == 21 ~ "Platelet",
                                         RNA_snn_res.1 == 22 ~ "Treg",
                                         RNA_snn_res.1 == 23 ~ "MEP",
                                         RNA_snn_res.1 == 24 ~ "B"
                                       ))
DimPlot(tonsil.predictions.m, group.by = "RNA_snn_res.1", label = T) + NoLegend()
DimPlot(tonsil.predictions.m, group.by = "cell.type.jx.1", label = T) + NoLegend()
saveRDS(tonsil.predictions.m, "tonsil.predictions.m.with.jx.anno.rds")

# Combined genes vector
all_genes <- c(
  # DC1 precursor
  "CLEC9A", "IDO1", "CADM1", "CST3", "SNX3", "C1orf54", "HLA-DPB1", "IRF8", "CCND1", "HLA-DPA1",
  
  # DC1 mature
  "CLEC9A", "XCR1",
  
  # DC2
  "CD1C", "CLEC10A", "FCER1A", "HLA-DQB1", "HLA-DRB1", "HLA-DQA1", "HLA-DPB1",
  
  # DC3
  "CD1C", "CLEC10A", "FCER1A", "S100A8", "S100A9", "VCAN", "LYZ", "ANXA1",
  
  # DC4
  "FCGR3A", "FTL", "SIGLEC10", "SERPINA1", "LST1", "AIF1",
  
  # DC5
  "SIGLEC6", "AXL", "CDH1", "CD22", "LILRA4",
  
  # IL7R DC
  "IL7R", "NR4A3", "CLDN1", "NRARP", "CLEC10A", "AXL", "IL4I1",
  
  # aDC1
  "SLCO5A1", "CCR7", "DUSP5", "LAD1", "TREML1", "LAMP3",
  
  # aDC2
  "PRDM16", "CCNA1", "SLC16A2", "S100A2", "AIRE", "CCL22",
  
  # aDC3
  "CCL19", "MARCKSL1", "FSCN1", "ACTG1", "SYNPO", "TMSB10",
  
  # M1 Macrophages
  "CD68", "CCL4", "IL1B", "CCL3", "TNF", "IL6", "IL1A", "PTX3",
  
  # Monocytes
  "CD14", "VCAN", "S100A9", "S100A8", "FCN1", "S100A4", "S100A12",
  
  # Mast
  "KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2", "LTC4S", "HPGD", "HPGDS", "LMO4",
  
  # Neutrophils
  "CXCL8", "G0S2", "PI3", "CCL3L1", "FCGR3B", "PTGS2", "IFITM2",
  
  # Cycling
  "CCNB2", "PCNA", "MKI67", "TOP2A",
  
  # MMP Slan-like
  "MMP9", "MMP12", "MMP14", "HLA-DRA", "HLA-DRB", "CCR1", "FCGR2A", "ADAMDEC1", "C1QC", "NECTIN2", "FUCA1",
  
  # C1Q Slan-like
  "C1QC", "APOE", "MMP9", "FUCA1", "APOC1", "ADAMDEC1", "HLA-DQA1", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
  
  # SELENOP Slan-like
  "SELENOP", "FOLR2", "PTGDS", "FUCA1", "APOE", "APOC1", "CCL18",
  
  # ITGAX Slan-like
  "FUCA1", "NR1H3", "APOE", "ITGAX", "ZEB2", "IL6R", "IL18", "MAF", "MAFB", "MAFG"
)

pdf("marker.gene.exp.myeloid.pdf", height = 15, width = 15)
VlnPlot(tonsil.predictions.m, features = all_genes, stack = T, flip = T, fill.by = "ident", group.by = "RNA_snn_res.1") + NoLegend()
dev.off()


genes_by_cell_type <- list(
  DC1_precursor = c("CLEC9A", "IDO1", "CADM1", "CST3", "SNX3", "C1orf54", "HLA-DPB1", "IRF8", "CCND1", "HLA-DPA1"),
  
  DC1_mature = c("CLEC9A", "XCR1"),
  
  DC2 = c("CD1C", "CLEC10A", "FCER1A", "HLA-DQB1", "HLA-DRB1", "HLA-DQA1", "HLA-DPB1"),
  
  DC3 = c("CD1C", "CLEC10A", "FCER1A", "S100A8", "S100A9", "VCAN", "LYZ", "ANXA1"),
  
  DC4 = c("FCGR3A", "FTL", "SIGLEC10", "SERPINA1", "LST1", "AIF1"),
  
  DC5 = c("SIGLEC6", "AXL", "CDH1", "CD22", "LILRA4"),
  
  IL7R_DC = c("IL7R", "NR4A3", "CLDN1", "NRARP", "CLEC10A", "AXL", "IL4I1"),
  
  aDC1 = c("SLCO5A1", "CCR7", "DUSP5", "LAD1", "TREML1", "LAMP3"),
  
  aDC2 = c("PRDM16", "CCNA1", "SLC16A2", "S100A2", "AIRE", "CCL22"),
  
  aDC3 = c("CCL19", "MARCKSL1", "FSCN1", "ACTG1", "SYNPO", "TMSB10"),
  
  M1_Macrophages = c("CD68", "CCL4", "IL1B", "CCL3", "TNF", "IL6", "IL1A", "PTX3"),
  
  Monocytes = c("CD14", "VCAN", "S100A9", "S100A8", "FCN1", "S100A4", "S100A12"),
  
  Mast = c("KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2", "LTC4S", "HPGD", "HPGDS", "LMO4"),
  
  Neutrophils = c("CXCL8", "G0S2", "PI3", "CCL3L1", "FCGR3B", "PTGS2", "IFITM2"),
  
  Cycling = c("CCNB2", "PCNA", "MKI67", "TOP2A"),
  
  MMP_Slan_like = c("MMP9", "MMP12", "MMP14", "HLA-DRA", "HLA-DRB", "CCR1", "FCGR2A", "ADAMDEC1", "C1QC", "NECTIN2", "FUCA1"),
  
  C1Q_Slan_like = c("C1QC", "APOE", "MMP9", "FUCA1", "APOC1", "ADAMDEC1", "HLA-DQA1", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"),
  
  SELENOP_Slan_like = c("SELENOP", "FOLR2", "PTGDS", "FUCA1", "APOE", "APOC1", "CCL18"),
  
  ITGAX_Slan_like = c("FUCA1", "NR1H3", "APOE", "ITGAX", "ZEB2", "IL6R", "IL18", "MAF", "MAFB", "MAFG")
)

  FeaturePlot(tonsil.predictions.m, features = c("MMP8", "OLR1", "CD24", "DEFA3", "MKI67")) + NoLegend()
  FeaturePlot(tonsil.predictions.m, features = c("KIT", "RUNX1", "GATA1", "CD34", "SPINK2", "HOXA7", "IRF7")) + NoLegend()
  
  
c11.markers = FindMarkers(tonsil.predictions.m, ident.1 = "11", group.by = "RNA_snn_res.1")
c11.markers %>% filter(avg_log2FC > 0)


c10.markers = FindMarkers(tonsil.predictions.m, ident.1 = "10", group.by = "RNA_snn_res.1")
c10.markers %>% filter(avg_log2FC > 0)
c10.markers %>% filter(avg_log2FC < 0)




c3.markers = FindMarkers(tonsil.predictions.m, ident.1 = "3", group.by = "RNA_snn_res.1")
c3.markers %>% filter(avg_log2FC > 0)
c10.markers %>% filter(avg_log2FC < 0)


c21.markers = FindMarkers(tonsil.predictions.m, ident.1 = "21", group.by = "RNA_snn_res.1")
c21.markers %>% filter(avg_log2FC > 0)
c10.markers %>% filter(avg_log2FC < 0)

c24.markers = FindMarkers(tonsil.predictions.m, ident.1 = "24", group.by = "RNA_snn_res.1")
c24.markers %>% filter(avg_log2FC > 0)


c22.markers = FindMarkers(tonsil.predictions.m, ident.1 = "22", group.by = "RNA_snn_res.1")
c22.markers %>% filter(avg_log2FC > 0)


c19.markers = FindMarkers(tonsil.predictions.m, ident.1 = "19", group.by = "RNA_snn_res.1")
c19.markers %>% filter(avg_log2FC > 0)

c24.markers = FindMarkers(tonsil.predictions.m, ident.1 = "20", group.by = "RNA_snn_res.1")
c24.markers %>% filter(avg_log2FC > 0)


c24.markers = FindMarkers(tonsil.predictions.m, ident.1 = "8", group.by = "RNA_snn_res.1")
c24.markers %>% filter(avg_log2FC > 0)
c24.markers = FindMarkers(tonsil.predictions.m, ident.1 = "15", group.by = "RNA_snn_res.1")
c24.markers %>% filter(avg_log2FC > 0)

t.b.in.m = subset(tonsil.predictions.m, cell.type.jx.1 %in% c("B", "Treg", "NaiveT", "EffectorT"))
saveRDS(t.b.in.m, "t.b.in.m.rds")

t.in.m = subset(tonsil.predictions.m, cell.type.jx.1 %in% c("Treg", "NaiveT", "EffectorT"))
saveRDS(t.in.m, "t.only.in.m.rds")

b.in.m = subset(tonsil.predictions.m, cell.type.jx.1 %in% c("B"))
saveRDS(b.in.m, "b.only.in.m.rds")

#myeloid all
Myeloid.all = merge(subset(Myeloid_177754_cells, predicted.annotation_20230508 %in% c("Monocytes", "Neutrophils")), tonsil.predictions.m)
Myeloid.all@meta.data

Myeloid.all = Myeloid.all %>% JoinLayers()

Myeloid.all = Myeloid.all %>% FindVariableFeatures() %>% ScaleData()

Myeloid.all = Myeloid.all %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)
Myeloid.all = Myeloid.all %>% FindNeighbors() %>% FindClusters(resolution = c(0.5, 1, 2))

DimPlot(Myeloid.all, group.by = "cell.type.jx.1", label = T)

DimPlot(Myeloid.all, group.by = "predicted.cell.type.short", label = T)
DimPlot(Myeloid.all, group.by = "predicted.annotation_20230508", label = T) + NoLegend()
DimPlot(Myeloid.all, group.by = "RNA_snn_res.1", label = T) + NoLegend()


Myeloid.all@meta.data %>%
  group_by(RNA_snn_res.1, predicted.annotation_20230508) %>%
  summarize(
    count = n(),
    mean_score = mean(predicted.annotation_20230508.score, na.rm = TRUE),  # Calculate mean score
    .groups = 'drop'
  ) %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string with count and mean score included
  summarize(
    cell.types = paste(predicted.annotation_20230508, count, 
                       round(mean_score, 2),  # Round mean score to 2 decimal places
                       sep = " (", collapse = "), "), ")"
  ) %>%
  ungroup() %>% 
  print(n = 30)



Myeloid.all@meta.data %>%
  group_by(RNA_snn_res.1, predicted.cell.type.short) %>%
  summarize(
    count = n(),
    mean_score = mean(predicted.annotation_20230508.score, na.rm = TRUE),  # Calculate mean score
    .groups = 'drop'
  ) %>%
  group_by(RNA_snn_res.1) %>%
  # Sort by count within each cluster
  arrange(desc(count)) %>%
  # Select the top three cell types for each cluster
  slice_head(n = 3) %>%
  # Create the cell.types string with count and mean score included
  summarize(
    cell.types = paste(predicted.cell.type.short, count, 
                       round(mean_score, 2),  # Round mean score to 2 decimal places
                       sep = " (", collapse = "), "), ")"
  ) %>%
  ungroup() %>% 
  print(n = 30)




Myeloid.all@meta.data <- 
  Myeloid.all@meta.data %>% mutate(cell.type.jx.1.m = 
                                              case_when(
                                                RNA_snn_res.1 == 0 ~ "Neut",
                                                RNA_snn_res.1 == 1 ~ "Neut",
                                                RNA_snn_res.1 == 2 ~ "Monocyte",
                                                RNA_snn_res.1 == 3 ~ "Monocyte",
                                                RNA_snn_res.1 == 4 ~ "Monocyte",
                                                RNA_snn_res.1 == 5 ~ "Monocyte",
                                                RNA_snn_res.1 == 6 ~ "Monocyte",
                                                RNA_snn_res.1 == 7 ~ "Neut",
                                                RNA_snn_res.1 == 8 ~ "SIGLEC10+/TCF7L2+ Monocyte",
                                                RNA_snn_res.1 == 9 ~ "Monocyte",
                                                RNA_snn_res.1 == 10 ~ "MMP8+ CEACAM8+ Neut",
                                                RNA_snn_res.1 == 11 ~ "EffectorT",
                                                RNA_snn_res.1 == 12 ~ "Platelet",
                                                RNA_snn_res.1 == 13 ~ "EffectorT",
                                                RNA_snn_res.1 == 14 ~ "EffectorT",
                                                RNA_snn_res.1 == 15 ~ "Megakaryocyte",
                                                RNA_snn_res.1 == 16 ~ "NaiveT",
                                                RNA_snn_res.1 == 17 ~ "NaiveT",
                                                RNA_snn_res.1 == 18 ~ "DC2",
                                                RNA_snn_res.1 == 19 ~ "Neut",
                                                RNA_snn_res.1 == 20 ~ "ZEB2-hi Monocyte",
                                                RNA_snn_res.1 == 21 ~ "Platelet",
                                                RNA_snn_res.1 == 22 ~ "Platelet",
                                                RNA_snn_res.1 == 23 ~ "B",
                                                RNA_snn_res.1 == 24 ~ "RBC",
                                                RNA_snn_res.1 == 25 ~ "PDC",
                                                RNA_snn_res.1 == 26 ~ "Granulocyte",
                                                RNA_snn_res.1 == 27 ~ "Treg",
                                                RNA_snn_res.1 == 28 ~ "HSPC/LMPP",
                                                RNA_snn_res.1 == 29 ~ "Mast",
                                                
                                              ))


c28.markers = FindMarkers(Myeloid.all, ident.1 = "28", group.by = "RNA_snn_res.1", max.cells.per.ident = 5000)
c28.markers %>% filter(avg_log2FC > 0)
c10.markers %>% filter(avg_log2FC < 0)


c29.markers = FindMarkers(Myeloid.all, ident.1 = "10", group.by = "RNA_snn_res.1",  max.cells.per.ident = 5000)
c29.markers %>% filter(avg_log2FC > 0)

DimPlot(Myeloid.all, group.by = "RNA_snn_res.1", label = T) + NoLegend()
DimPlot(Myeloid.all, group.by = "cell.type.jx.1.m", label = T) + NoLegend()

FeaturePlot(Myeloid.all, c("MMP8", "ZEB2", "DEFA3"))



Myeloid.all@meta.data <- 
  Myeloid.all@meta.data %>% mutate(cell.type.jx.1.m.binned = 
                                     case_when(
                                       RNA_snn_res.1 == 0 ~ "Neut",
                                       RNA_snn_res.1 == 1 ~ "Neut",
                                       RNA_snn_res.1 == 2 ~ "Monocyte",
                                       RNA_snn_res.1 == 3 ~ "Monocyte",
                                       RNA_snn_res.1 == 4 ~ "Monocyte",
                                       RNA_snn_res.1 == 5 ~ "Monocyte",
                                       RNA_snn_res.1 == 6 ~ "Monocyte",
                                       RNA_snn_res.1 == 7 ~ "Neut",
                                       RNA_snn_res.1 == 8 ~ "Monocyte",
                                       RNA_snn_res.1 == 9 ~ "Monocyte",
                                       RNA_snn_res.1 == 10 ~ "MMP8+ CEACAM8+ Neut",
                                       RNA_snn_res.1 == 11 ~ "EffectorT",
                                       RNA_snn_res.1 == 12 ~ "Platelet",
                                       RNA_snn_res.1 == 13 ~ "EffectorT",
                                       RNA_snn_res.1 == 14 ~ "EffectorT",
                                       RNA_snn_res.1 == 15 ~ "Megakaryocyte",
                                       RNA_snn_res.1 == 16 ~ "NaiveT",
                                       RNA_snn_res.1 == 17 ~ "NaiveT",
                                       RNA_snn_res.1 == 18 ~ "DC2",
                                       RNA_snn_res.1 == 19 ~ "Neut",
                                       RNA_snn_res.1 == 20 ~ "Monocyte",
                                       RNA_snn_res.1 == 21 ~ "Platelet",
                                       RNA_snn_res.1 == 22 ~ "Platelet",
                                       RNA_snn_res.1 == 23 ~ "B",
                                       RNA_snn_res.1 == 24 ~ "RBC",
                                       RNA_snn_res.1 == 25 ~ "PDC",
                                       RNA_snn_res.1 == 26 ~ "Granulocyte",
                                       RNA_snn_res.1 == 27 ~ "Treg",
                                       RNA_snn_res.1 == 28 ~ "HSPC/LMPP",
                                       RNA_snn_res.1 == 29 ~ "Mast",
                                       
                                     ))


saveRDS(Myeloid.all, "Myeloid.all_184389_cells.rds")



t.in.m.final = subset(Myeloid.all, cell.type.jx.1.m %in% c("HSPC/LMPP", "Treg", "NaiveT", "EffectorT"))

b.in.m.final = subset(Myeloid.all, cell.type.jx.1.m %in% c("B"))
saveRDS(t.in.m.final, "t.in.m.final.16661.cells.rds")
saveRDS(b.in.m.final, "b.in.m.final.570.cells.rds")

m.only = subset(Myeloid.all, !cell.type.jx.1.m %in% c("HSPC/LMPP", "Treg", "NaiveT", "EffectorT", "B"))
m.only$cell.type.jx.1.m.binned %>% table()

m.only = m.only %>% FindVariableFeatures() %>% ScaleData()
m.only = m.only %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)
m.only = m.only %>% FindNeighbors() %>% FindClusters(resolution = c(0.5, 1, 2))

DimPlot(m.only, group.by = "cell.type.jx.1.m", label = T) + NoLegend()

saveRDS(m.only, "m.only.final.rds")

m.only$cell.type.jx.1.m.binned %>% table() 
m.only$cell.type.jx.1.m %>% table() %>% length()

#m.only$cell.type.jx.1.m[m.only$RNA_snn_res.1 == 21] = "RBC"

#m.only$cell.type.jx.1.m.factor = factor(m.only$cell.type.jx.1.m, levels = c("Monocyte", "SIGLEC10+/TCF7L2+ Monocyte", "ZEB2-hi Monocyte", "Neut", "MMP8+ CEACAM8+ Neut", "DC2", "PDC" ,"Mast", "Megakaryocyte", "Granulocyte", "Platelet", "RBC"))
m.only$cell.type.jx.1.m.binned.factor = factor(m.only$cell.type.jx.1.m.binned, levels = c("Monocyte", "Neut", "MMP8+ CEACAM8+ Neut", "DC2", "PDC" ,"Mast", "Megakaryocyte", "Granulocyte", "Platelet", "RBC")) 

saveRDS(m.only, "m.only.final.rds")

                              
DimPlot(m.only, group.by = "RNA_snn_res.1", label = T)



myeloid.colors = colorRampPalette(RColorBrewer::brewer.pal(n = 12, name = "Paired"))(13)
myeloid.colors.binned = colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = "Accent"))(10)

pdf("myeloid.umap.2.pdf", height = 6, width = 8)
DimPlot(m.only, group.by = "cell.type.jx.1.m.factor", label = T, cols = myeloid.colors) + ggtitle("m.only. 167158 cells")
DimPlot(m.only, group.by = "cell.type.jx.1.m.binned", label = T, cols = myeloid.colors.binned)
dev.off()




m.summ.stats.res.1.m.binned.2 <- m.only@meta.data %>%
  group_by(patient, timepoint, cell.type.jx.1.m.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(m.only@meta.data %>%
              group_by(patient, timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

# Create plots for myeloid cells
pdf("m.summ.stats.m.2.pdf", height = 5, width = 15)
ggplot(m.summ.stats.res.1.m.binned.2, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~cell.type.jx.1.m.factor, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) +
  theme_bw()   + ggtitle("with unpaired p-values")

ggplot(m.summ.stats.res.1.m.binned.2, aes(x = timepoint, y = prop, fill = timepoint, group = patient)) +
  geom_boxplot(outlier.size = -1, aes(x = timepoint, y = prop, fill = timepoint), inherit.aes = F) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~cell.type.jx.1.m.factor, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3")), paired = T) +
  theme_bw() + ggtitle("with paired p-values")

dev.off()


#de DE for the following cell types

m.only$cell.type.jx.1.m.factor %>% table()

# For Treg
siglec.mono = FindMarkers(m.only, ident.1 = "SIGLEC10+/TCF7L2+ Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)

VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c("CD1A", "CD14", "LYZ", "ITGAX", "SIGLEC10", "HES4", "TCF7L2", "FCGR3A", "FCGR3B", "VCAN", "S100A8", "S100A9", "LST1", "CD14", "MS4A7", "IFITM3", "IFITM1", "CX3CR1", "SAT1", "NCF1", "CD52", "KLF2", "CD99", "CD79B", "B2M", "SELL", "C1QA", "LYN", "FCER1G", "IFNG"), stack = T, flip = T, fill.by = "ident")
VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c("CD14", "FCGR3A", "CLEC4C", "ITGAX", "IFI30", "S100A8", "S100A9", "CD33", "HLA-DRA", "CD74","CD36", "CCR2", "IL12A", "CX3CR1", "NRP1", "NRP2", "CD1C", "CD3E", "ITGAM", "CD40", "ADGRE1"), stack = T, flip = T, fill.by = "ident")

siglec.mono = siglec.mono %>% as_tibble(rownames = "gene")
siglec.mono = siglec.mono %>% filter(p_val_adj < 0.05)
ggplot(siglec.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = siglec.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = siglec.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("CD16 Mono vs all") +
  #geom_text_repel(aes(label = gene), data = siglec.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-10, 10, by = 1))  # Adjust the range as needed



siglec.mono2 = FindMarkers(m.only, ident.1 = "SIGLEC10+/TCF7L2+ Monocyte", ident.2 = "Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
siglec.mono2 = siglec.mono2 %>% as_tibble(rownames = "gene")
siglec.mono2 = siglec.mono2 %>% filter(p_val_adj < 0.05)
ggplot(siglec.mono2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = siglec.mono2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = siglec.mono2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("CD16+ Mono vs CD14+ Mono") +
  #geom_text_repel(aes(label = gene), data = siglec.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-10, 10, by = 1))  # Adjust the range as needed


mast = FindMarkers(m.only, ident.1 = "Mast", ident.2 = "Neut", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
mast = mast %>% as_tibble(rownames = "gene")
mast = mast %>% filter(p_val_adj < 0.05)
ggplot(mast, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mast %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mast %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Mast vs Neut") +
  #geom_text_repel(aes(label = gene), data = siglec.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed



mast2 = FindMarkers(m.only, ident.1 = "Mast", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
mast2 = mast2 %>% as_tibble(rownames = "gene")
mast2 = mast2 %>% filter(p_val_adj < 0.05)
ggplot(mast2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mast2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mast2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-30), max.overlaps = 50) +
  ggtitle("Mast vs Neut") +
  #geom_text_repel(aes(label = gene), data = siglec.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed




# 1. ZEB2-hi Monocyte
zeb2_mono = FindMarkers(m.only, ident.1 = "ZEB2-hi Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
zeb2_mono = zeb2_mono %>% as_tibble(rownames = "gene")
zeb2_mono = zeb2_mono %>% filter(p_val_adj < 0.05)
ggplot(zeb2_mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = zeb2_mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = zeb2_mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-75), max.overlaps = 50) +
  ggtitle("ZEB2-hi Monocyte") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# 1. ZEB2-hi Monocyte vs Mono
zeb2_mono2 = FindMarkers(m.only, ident.1 = "ZEB2-hi Monocyte", ident.2 = "Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
zeb2_mono2 = zeb2_mono2 %>% as_tibble(rownames = "gene")
zeb2_mono2 = zeb2_mono2 %>% filter(p_val_adj < 0.05)
ggplot(zeb2_mono2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = zeb2_mono2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = zeb2_mono2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("ZEB2-hi Monocyte vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed



# 2. MMP8+ CEACAM8+ Neut
mmp8_neut = FindMarkers(m.only, ident.1 = "MMP8+ CEACAM8+ Neut", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
mmp8_neut = mmp8_neut %>% as_tibble(rownames = "gene")
mmp8_neut = mmp8_neut %>% filter(p_val_adj < 0.05)
ggplot(mmp8_neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mmp8_neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mmp8_neut %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("MMP8+ CEACAM8+ Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed


# 2. MMP8+ CEACAM8+ Neut vs Neut
mmp8_neut2 = FindMarkers(m.only, ident.1 = "MMP8+ CEACAM8+ Neut", ident.2 = "Neut", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
mmp8_neut2 = mmp8_neut2 %>% as_tibble(rownames = "gene")
mmp8_neut2 = mmp8_neut2 %>% filter(p_val_adj < 0.05)
ggplot(mmp8_neut2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mmp8_neut2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mmp8_neut2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-75), max.overlaps = 50) +
  ggtitle("MMP8+ CEACAM8+ Neut vs Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed


# 3. DC2
dc2 = FindMarkers(m.only, ident.1 = "DC2", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
dc2 = dc2 %>% as_tibble(rownames = "gene")
dc2 = dc2 %>% filter(p_val_adj < 0.05)
ggplot(dc2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = dc2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = dc2 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("DC2") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# 3. DC2 vs Mono
dc2 = FindMarkers(m.only, ident.1 = "DC2", ident.2 = "Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
dc2 = dc2 %>% as_tibble(rownames = "gene")
dc2 = dc2 %>% filter(p_val_adj < 0.05)
ggplot(dc2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = dc2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = dc2 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("DC2 vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# 4. PDC
pdc = FindMarkers(m.only, ident.1 = "PDC", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
pdc = pdc %>% as_tibble(rownames = "gene")
pdc = pdc %>% filter(p_val_adj < 0.05)
ggplot(pdc, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = pdc %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = pdc %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("PDC") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed


# 4. PDC vs Mono
pdc2 = FindMarkers(m.only, ident.1 = "PDC", ident.2 = "DC2", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
pdc2 = pdc2 %>% as_tibble(rownames = "gene")
pdc2 = pdc2 %>% filter(p_val_adj < 0.05)
ggplot(pdc2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = pdc2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = pdc2 %>% filter(abs(avg_log2FC) > 3 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("PDC vs DC2") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed


# 5. Megakaryocyte
megakaryocyte = FindMarkers(m.only, ident.1 = "Megakaryocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
megakaryocyte = megakaryocyte %>% as_tibble(rownames = "gene")
megakaryocyte = megakaryocyte %>% filter(p_val_adj < 0.05)
ggplot(megakaryocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = megakaryocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = megakaryocyte %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Megakaryocyte") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# 5. Megakaryocyte vs PLT
megakaryocyte2 = FindMarkers(m.only, ident.1 = "Megakaryocyte", ident.2 = "Platelet", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
megakaryocyte2 = megakaryocyte2 %>% as_tibble(rownames = "gene")
megakaryocyte2 = megakaryocyte2 %>% filter(p_val_adj < 0.05)

megakaryocyte2 <- megakaryocyte2 %>%
  mutate(is_ribosomal = grepl("^RP", gene))

# Create the plot
pdf("volcano.meg.vs.plt.pdf", height = 15, width = 15)
ggplot(megakaryocyte2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(aes(color = is_ribosomal), alpha = 1) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +  # Set colors for ribosomal and non-ribosomal genes
  theme_bw() + 
  geom_text_repel(aes(label = gene, color = is_ribosomal), 
                  data = megakaryocyte2 %>% filter(abs(avg_log2FC) > 3 & p_val_adj < 1e-100), 
                  max.overlaps = 50) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +  # Set colors for text as well
  ggtitle("Megakaryocyte vs Platelet") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1)) +  # Adjust the range as needed
  labs(color = "Ribosomal Gene")  # Add a legend title

ggplot(megakaryocyte2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(aes(color = is_ribosomal), alpha = 1) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +  # Set colors for ribosomal and non-ribosomal genes
  theme_bw() + 
  geom_text_repel(aes(label = gene), 
                  data = megakaryocyte2 %>% filter(abs(avg_log2FC) > 3 & p_val_adj < 1e-100 & !is_ribosomal), 
                  max.overlaps = 50) +  # Exclude ribosomal genes from labeling
  ggtitle("Megakaryocyte vs Platelet") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1)) +  # Adjust the range as needed
  labs(color = "Ribosomal Gene")  # Add a legend title
dev.off()


# 5. Megakaryocyte vs Mono
megakaryocyte3 = FindMarkers(m.only, ident.1 = "Megakaryocyte", ident.2 = "Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
megakaryocyte3 = megakaryocyte3 %>% as_tibble(rownames = "gene")
megakaryocyte3 = megakaryocyte3 %>% filter(p_val_adj < 0.05)
pdf("volcano.meg.vs.mono.pdf")
ggplot(megakaryocyte3, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = megakaryocyte3 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = megakaryocyte3 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Megakaryocyte vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed
dev.off()

# 6. Granulocyte
granulocyte = FindMarkers(m.only, ident.1 = "Granulocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
granulocyte = granulocyte %>% as_tibble(rownames = "gene")
granulocyte = granulocyte %>% filter(p_val_adj < 0.05)
ggplot(granulocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = granulocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = granulocyte %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-10), max.overlaps = 50) +
  ggtitle("Granulocyte") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed


# 6. Granulocyte vs Mono
granulocyte2 = FindMarkers(m.only, ident.1 = "Granulocyte", ident.2 = "Monocyte", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
granulocyte2 = granulocyte2 %>% as_tibble(rownames = "gene")
granulocyte2 = granulocyte2 %>% filter(p_val_adj < 0.05)
pdf("volcano.gran.vs.mono.pdf")
ggplot(granulocyte2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = granulocyte2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = granulocyte2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Granulocyte vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed
dev.off()

# 6. Granulocyte vs Neut
granulocyte3 = FindMarkers(m.only, ident.1 = "Granulocyte", ident.2 = "Neut", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
granulocyte3 = granulocyte3 %>% as_tibble(rownames = "gene")
granulocyte3 = granulocyte3 %>% filter(p_val_adj < 0.05)
granulocyte3 <- granulocyte3 %>%
  mutate(is_ribosomal = grepl("^RP", gene))

# Create the plot
pdf("volcano.gran.vs.neut.pdf")
ggplot(granulocyte3, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(aes(color = is_ribosomal), alpha = 1) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +  # Set colors for ribosomal and non-ribosomal genes
  theme_bw() + 
  geom_text_repel(aes(label = gene), 
                  data = granulocyte3 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100 & !is_ribosomal | avg_log2FC < -2 & p_val_adj < 1e-50), 
                  max.overlaps = 50) +  # Exclude ribosomal genes from labeling
  ggtitle("Granulocyte vs Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1)) +  # Adjust the range as needed
  labs(color = "Ribosomal Gene")  # Add a legend title
dev.off()

# 7. Platelet
platelet = FindMarkers(m.only, ident.1 = "Platelet", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
platelet = platelet %>% as_tibble(rownames = "gene")
platelet = platelet %>% filter(p_val_adj < 0.05)
ggplot(platelet, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = platelet %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = platelet %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Platelet") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# 8. RBC
rbc = FindMarkers(m.only, ident.1 = "RBC", group.by = "cell.type.jx.1.m.factor", max.cells.per.ident = 1000)
rbc = rbc %>% as_tibble(rownames = "gene")
rbc = rbc %>% filter(p_val_adj < 0.05)
ggplot(rbc, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = rbc %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = rbc %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("RBC") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed



# Start the PDF device----

pdf("all.deg.volcano.plots.pdf", height = 20, width = 20)

# Plot for SIGLEC10+/TCF7L2+ Monocyte vs all
ggplot(siglec.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = siglec.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = siglec.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("CD16 Mono vs all") +
  scale_x_continuous(breaks = seq(-10, 10, by = 1))  # Adjust the range as needed

# Plot for SIGLEC10+/TCF7L2+ Monocyte vs Monocyte
ggplot(siglec.mono2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = siglec.mono2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = siglec.mono2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("CD16+ Mono vs CD14+ Mono") +
  scale_x_continuous(breaks = seq(-10, 10, by = 1))  # Adjust the range as needed

# Plot for Mast vs Neut
ggplot(mast, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mast %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mast %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Mast vs Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for Mast vs all
ggplot(mast2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mast2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mast2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-30), max.overlaps = 50) +
  ggtitle("Mast vs Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for ZEB2-hi Monocyte
ggplot(zeb2_mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = zeb2_mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = zeb2_mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-75), max.overlaps = 50) +
  ggtitle("ZEB2-hi Monocyte") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for ZEB2-hi Monocyte vs Mono
ggplot(zeb2_mono2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = zeb2_mono2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = zeb2_mono2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("ZEB2-hi Monocyte vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for MMP8+ CEACAM8+ Neut
ggplot(mmp8_neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mmp8_neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mmp8_neut %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("MMP8+ CEACAM8+ Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for MMP8+ CEACAM8+ Neut vs Neut
ggplot(mmp8_neut2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = mmp8_neut2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = mmp8_neut2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-75), max.overlaps = 50) +
  ggtitle("MMP8+ CEACAM8+ Neut vs Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for DC2
ggplot(dc2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = dc2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = dc2 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("DC2") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for DC2 vs Mono
ggplot(dc2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = dc2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = dc2 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("DC2 vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for PDC
ggplot(pdc, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = pdc %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = pdc %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("PDC") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for PDC vs DC2
ggplot(pdc2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = pdc2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = pdc2 %>% filter(abs(avg_log2FC) > 3 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("PDC vs DC2") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for Megakaryocyte
ggplot(megakaryocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = megakaryocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = megakaryocyte %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Megakaryocyte") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for Megakaryocyte vs Platelet
ggplot(megakaryocyte2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(aes(color = is_ribosomal), alpha = 1) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +  # Set colors for ribosomal and non-ribosomal genes
  theme_bw() + 
  geom_text_repel(aes(label = gene), 
                  data = megakaryocyte2 %>% filter(abs(avg_log2FC) > 3 & p_val_adj < 1e-100 & !is_ribosomal), 
                  max.overlaps = 50) +  # Exclude ribosomal genes from labeling
  ggtitle("Megakaryocyte vs Platelet") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1)) +  # Adjust the range as needed
  labs(color = "Ribosomal Gene")  # Add a legend title

# Plot for Megakaryocyte vs Monocyte
ggplot(megakaryocyte3, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = megakaryocyte3 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = megakaryocyte3 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Megakaryocyte vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for Granulocyte
ggplot(granulocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = granulocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = granulocyte %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-10), max.overlaps = 50) +
  ggtitle("Granulocyte") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for Granulocyte vs Monocyte
ggplot(granulocyte2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = granulocyte2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = granulocyte2 %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Granulocyte vs Mono") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for Granulocyte vs Neut
ggplot(granulocyte3, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(aes(color = is_ribosomal), alpha = 1) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +  # Set colors for ribosomal and non-ribosomal genes
  theme_bw() + 
  geom_text_repel(aes(label = gene), 
                  data = granulocyte3 %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100 & !is_ribosomal | avg_log2FC < -2 & p_val_adj < 1e-50), 
                  max.overlaps = 50) +  # Exclude ribosomal genes from labeling
  ggtitle("Granulocyte vs Neut") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1)) +  # Adjust the range as needed
  labs(color = "Ribosomal Gene")  # Add a legend title

# Plot for Platelet
ggplot(platelet, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = platelet %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = platelet %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Platelet") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Plot for RBC
ggplot(rbc, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = rbc %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = rbc %>% filter(abs(avg_log2FC) > 2 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("RBC") +
  scale_x_continuous(breaks = seq(-10, 20, by = 1))  # Adjust the range as needed

# Close the PDF device
dev.off()

#make vlnplot----
VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c("ZBTB16", "DDIT4", "IRS2", "FLT3", "CXCR4", "CEBPD", "CD74", "", "HBB", "HBA2"), split.by = "timepoint", stack = T, flip = T, fill.by = "ident")




#DE Base vs M3
VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c("ZBTB16", "DDIT4", "IRS2", "FLT3", "CXCR4", "CEBPD", "CD74", "B2M", "CTSS", "S100A9"), split.by = "timepoint", stack = T, flip = T, fill.by = "ident")

# For Monocyte
Idents(m.only) = m.only$cell.type.jx.1.m.binned.factor
m.only$cell.type.jx.1.m.binned.factor %>% table() %>% names()
isg.consensus = molecular_signatures %>% unlist() %>% unique()

m3.vs.base.Monocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Monocyte", max.cells.per.ident = 5000)
m3.vs.base.Monocyte <- m3.vs.base.Monocyte %>% as_tibble(rownames = "gene")
m3.vs.base.Monocyte <- m3.vs.base.Monocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Monocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Monocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Monocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-25), max.overlaps = 50) +
  ggtitle("Monocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed



# For Neut
m3.vs.base.Neut <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Neut", max.cells.per.ident = 5000)
m3.vs.base.Neut <- m3.vs.base.Neut %>% as_tibble(rownames = "gene")
m3.vs.base.Neut <- m3.vs.base.Neut %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Neut %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-25), max.overlaps = 50) +
  ggtitle("Neut Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For MMP8+ CEACAM8+ Neut
m3.vs.base.MMP8_CEACAM8_Neut <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "MMP8+ CEACAM8+ Neut", max.cells.per.ident = 5000)
m3.vs.base.MMP8_CEACAM8_Neut <- m3.vs.base.MMP8_CEACAM8_Neut %>% as_tibble(rownames = "gene")
m3.vs.base.MMP8_CEACAM8_Neut <- m3.vs.base.MMP8_CEACAM8_Neut %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.MMP8_CEACAM8_Neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.MMP8_CEACAM8_Neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.MMP8_CEACAM8_Neut %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5), max.overlaps = 50) +
  ggtitle("MMP8+ CEACAM8+ Neut Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For DC2
m3.vs.base.DC2 <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "DC2", max.cells.per.ident = 5000)
m3.vs.base.DC2 <- m3.vs.base.DC2 %>% as_tibble(rownames = "gene")
m3.vs.base.DC2 <- m3.vs.base.DC2 %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.DC2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.DC2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.DC2 %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-10), max.overlaps = 50) +
  ggtitle("DC2 Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For PDC
m3.vs.base.PDC <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "PDC", max.cells.per.ident = 5000)
m3.vs.base.PDC <- m3.vs.base.PDC %>% as_tibble(rownames = "gene")
m3.vs.base.PDC <- m3.vs.base.PDC %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.PDC, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.PDC %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.PDC %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-01), max.overlaps = 50) +
  ggtitle("PDC Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Mast
m3.vs.base.Mast <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Mast", max.cells.per.ident = 5000)
m3.vs.base.Mast <- m3.vs.base.Mast %>% as_tibble(rownames = "gene")
m3.vs.base.Mast <- m3.vs.base.Mast %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Mast, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Mast %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Mast %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-25), max.overlaps = 50) +
  ggtitle("Mast Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Megakaryocyte
m3.vs.base.Megakaryocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Megakaryocyte", max.cells.per.ident = 5000)
m3.vs.base.Megakaryocyte <- m3.vs.base.Megakaryocyte %>% as_tibble(rownames = "gene")
m3.vs.base.Megakaryocyte <- m3.vs.base.Megakaryocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Megakaryocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Megakaryocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Megakaryocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5), max.overlaps = 50) +
  ggtitle("Megakaryocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Granulocyte
m3.vs.base.Granulocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Granulocyte", max.cells.per.ident = 5000)
m3.vs.base.Granulocyte <- m3.vs.base.Granulocyte %>% as_tibble(rownames = "gene")
m3.vs.base.Granulocyte <- m3.vs.base.Granulocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Granulocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Granulocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Granulocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-1), max.overlaps = 50) +
  ggtitle("Granulocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Platelet
m3.vs.base.Platelet <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Platelet", max.cells.per.ident = 5000)
m3.vs.base.Platelet <- m3.vs.base.Platelet %>% as_tibble(rownames = "gene")
m3.vs.base.Platelet <- m3.vs.base.Platelet %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Platelet, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Platelet %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Platelet %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5), max.overlaps = 50) +
  ggtitle("Platelet Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For RBC
m3.vs.base.RBC <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "RBC", max.cells.per.ident = 5000)
m3.vs.base.RBC <- m3.vs.base.RBC %>% as_tibble(rownames = "gene")
m3.vs.base.RBC <- m3.vs.base.RBC %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.RBC, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.RBC %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.RBC %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-1), max.overlaps = 50) +
  ggtitle("RBC Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed



m.only
m.only@meta.data
# Identify ribosomal and mitochondrial genes
ribosomal_genes <- grep("^RPL|^RPS", rownames(m.only), value = F)
mitochondrial_genes <- grep("^MT-", rownames(m.only), value = F)
anti.sense =  grep("^HLA-", rownames(m.only), value = T) %>% grep(pattern = "AS", value = F)
HLA.genes <- grep("^HLA-", rownames(m.only), value = F)[-anti.sense]

anti.sense.dr =  grep("^HLA-D", rownames(m.only), value = T) %>% grep(pattern = "AS", value = F)
HLA.D.genes <- grep("^HLA-D", rownames(m.only), value = F)[-anti.sense.dr]

# Sum the expression values for ribosomal genes
m.only@meta.data$RB.genes <- Matrix::colSums(m.only@assays$RNA@layers$data[ribosomal_genes, ])

# Sum the expression values for mitochondrial genes
m.only@meta.data$MT.genes <- Matrix::colSums(m.only@assays$RNA@layers$data[mitochondrial_genes, ])

m.only@meta.data$HLA.genes <- Matrix::colSums(m.only@assays$RNA@layers$data[HLA.genes, ])
m.only@meta.data$HLA.D.genes <- Matrix::colSums(m.only@assays$RNA@layers$data[HLA.D.genes, ])



pdf("marker.gene.expression.pdf", height = 12, width = 10)
VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c('nCount_RNA', "FCGR3B", "CEACAM8", "LTF", "S100A12", "IL3RA", "HDC", "CLC"), stack = T, flip = T, fill.by = "ident")

VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c('nCount_RNA', 'CXCL8', "IFI27", "IL1R2", "FCGR3A", "FLT3", "CXCR4", "CEBPD", "CD74", "B2M", "CTSS", "S100A9", "IRF1", "TLR7"), stack = T, flip = T, fill.by = "ident")
VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", c("CD1A", "CD14", "LYZ", "ITGAX", "SIGLEC10", "HES4", "TCF7L2", "ZEB2", "FCGR3A",  "FCGR3B", "VCAN", 
                                                        "S100A8", "S100A9", "LST1", "MMP8", "CEACAM8", "DEFA3", "NAMPT", "SOD2",  "MS4A7", "IFITM3", "IFITM1", "CX3CR1", "SAT1", 
                                                        "NCF1", "CD52", "HLA-DQA1", "CD74", "FCER1A", "CD1C", "CLEC10A", "IGKC", "JCHAIN", "GZMB", "IRF8", 
                                                        "CD99", "CD79B", "B2M", "SELL", "C1QA", "LYN", "FCER1G", "GATA2", "HDC", "IL3RA", 
                                                        "HBA2", "PF4", "NRGN", "PPBP", "RB.genes", "MT.genes", "nCount_RNA"), stack = T, flip = T, fill.by = "ident", cols = myeloid.colors)
dev.off()




pdf("UMAP.over.time.pdf", height = 8, width = 24)
DimPlot(m.only, group.by = "cell.type.jx.1.m.factor", label = T, cols = myeloid.colors, split.by = "timepoint")

dev.off()


# Summarize statistics for the binned cell types
m.summ.stats.res.1 <- m.only@meta.data %>%
  group_by(patient, timepoint, cell.type.jx.1.m.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(m.only@meta.data %>%
              group_by(timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)


pdf("bar.changes.over.time.pdf")
ggplot(m.summ.stats.res.1, aes(x = timepoint, y = prop, fill = cell.type.jx.1.m.factor)) + geom_col() + theme_bw() + scale_fill_manual(values = myeloid.colors)
dev.off()

m.summ.stats.res.2 <- m.only@meta.data %>%
  group_by(timepoint, cell.type.jx.1.m.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(m.only@meta.data %>%
              group_by(timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)


pdf("bar.changes.over.time.pdf")
ggplot(m.summ.stats.res.2, aes(x = timepoint, y = prop, fill = cell.type.jx.1.m.factor)) + geom_col() + theme_bw() + scale_fill_manual(values = myeloid.colors) + geom_text(aes(label = round(prop, 2)),position = position_stack(vjust = 0.5))
dev.off()


m.only$cell.type.jx.1.m.binned.factor = m.only$cell.type.jx.1.m.binned
m.only$cell.type.jx.1.m.binned.factor[m.only$cell.type.jx.1.m.factor == "SIGLEC10+/TCF7L2+ Monocyte"] = "CD16+ Monocyte"

m.only$cell.type.jx.1.m.binned.factor = factor(m.only$cell.type.jx.1.m.binned.factor, levels =  c("Monocyte", "CD16+ Monocyte", "Neut","MMP8+ CEACAM8+ Neut", "DC2",   "PDC",               
                                                  "Mast", "Megakaryocyte", "Granulocyte", "Platelet",  "RBC"))

m.only$cell.type.jx.1.m.binned.factor %>% table()



Idents(m.only) = m.only$cell.type.jx.1.m.binned.factor
Idents(m.only) %>% table()
m.only$cell.type.jx.1.m %>% table()

# For CD16 Monocyte
m3.vs.base.CD16_Monocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "CD16+ Monocyte")
m3.vs.base.CD16_Monocyte <- m3.vs.base.CD16_Monocyte %>% as_tibble(rownames = "gene")
m3.vs.base.CD16_Monocyte <- m3.vs.base.CD16_Monocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.CD16_Monocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.CD16_Monocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.CD16_Monocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-10), max.overlaps = 50) +
  ggtitle("CD16 Monocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For ZEB2-hi Monocyte
m3.vs.base.ZEB2_hi_Monocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "ZEB2-hi Monocyte")
m3.vs.base.ZEB2_hi_Monocyte <- m3.vs.base.ZEB2_hi_Monocyte %>% as_tibble(rownames = "gene")
m3.vs.base.ZEB2_hi_Monocyte <- m3.vs.base.ZEB2_hi_Monocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.ZEB2_hi_Monocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.ZEB2_hi_Monocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.ZEB2_hi_Monocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.ZEB2_hi_Monocyte %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("ZEB2-hi Monocyte Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.ZEB2_hi_Monocyte %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Monocyte
m3.vs.base.Monocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Monocyte", max.cells.per.ident = 10000)
m3.vs.base.Monocyte <- m3.vs.base.Monocyte %>% as_tibble(rownames = "gene")
m3.vs.base.Monocyte <- m3.vs.base.Monocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Monocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Monocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Monocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-25), max.overlaps = 50) +
  ggtitle("Monocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

pdf("m3.vs.base.Mono.binned.pdf")
ggplot(m3.vs.base.Monocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Monocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Monocyte %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100 | p_val_adj < 1e-200), max.overlaps = 50) +
  ggtitle("Monocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('CEBPD', 'DUSP1', "FOS", "S100A9", "NAMPT", "FKBP5", "CCND3", "SOCS3", "EGR1", "MT.genes", "HLA.D.genes", "B2M", "CD74", "ITM2B"), stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Monocyte", "CD16+ Monocyte"))

VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c("MT.genes", "HLA.D.genes", "B2M", "CD74", "ITM2B"), stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Monocyte", "CD16+ Monocyte"))
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('CEBPD', 'DUSP1', "FOS", "S100A9", "NAMPT", "FKBP5", "CCND3", "SOCS3", "EGR1"), stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Monocyte", "CD16+ Monocyte"))
dev.off()


# For Neut
m3.vs.base.Neut <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Neut", max.cells.per.ident = 5000)
m3.vs.base.Neut <- m3.vs.base.Neut %>% as_tibble(rownames = "gene")
m3.vs.base.Neut <- m3.vs.base.Neut %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Neut %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-50 | p_val_adj < 1e-50, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.Neut %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Neut Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.Neut %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

pdf("m3.vs.base.Neut.binned.pdf")
ggplot(m3.vs.base.Neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Neut %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Neut Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('IL1R2', 'FKBP5', "FCGR3A", "TLR2", "CEBPD", "CXCR4", "CEBPD", "CCND3"), stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Neut"))
dev.off()



# For MMP8+ CEACAM8+ Neut
m3.vs.base.MMP8_CEACAM8_Neut <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "MMP8+ CEACAM8+ Neut")
m3.vs.base.MMP8_CEACAM8_Neut <- m3.vs.base.MMP8_CEACAM8_Neut %>% as_tibble(rownames = "gene")
m3.vs.base.MMP8_CEACAM8_Neut <- m3.vs.base.MMP8_CEACAM8_Neut %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.MMP8_CEACAM8_Neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.MMP8_CEACAM8_Neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.MMP8_CEACAM8_Neut %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.MMP8_CEACAM8_Neut %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("MMP8+ CEACAM8+ Neut Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.MMP8_CEACAM8_Neut %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For DC2
m3.vs.base.DC2 <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "DC2")
m3.vs.base.DC2 <- m3.vs.base.DC2 %>% as_tibble(rownames = "gene")
m3.vs.base.DC2 <- m3.vs.base.DC2 %>% filter(p_val_adj < 0.05)

pdf("m3.vs.base.DC2.binned.pdf")
ggplot(m3.vs.base.DC2, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.DC2 %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.DC2 %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  ggtitle("DC2 Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()

pdf("m3.vs.base.DC2.binned.expression.pdf")
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('CXCR4', 'DDIT4', "IL1R2", "CEBPD", "CEBPB", "MT.genes"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("DC2"))

dev.off()


# For PDC
m.only$cell.type.jx.1.m.binned.factor %>% table() %>% names()
m.only@meta.data %>% group_by(cell.type.jx.1.m.binned.factor, timepoint) %>% summarize(n=length(timepoint)) %>% print(n=30)

m3.vs.base.PDC <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "PDC")
m3.vs.base.PDC <- m3.vs.base.PDC %>% as_tibble(rownames = "gene")
m3.vs.base.PDC <- m3.vs.base.PDC %>% filter(p_val_adj < 0.8)

pdf("m3.vs.base.PDC.binned.pdf")
ggplot(m3.vs.base.PDC, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.PDC %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.PDC %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-2), max.overlaps = 50) +
  geom_point(data = m3.vs.base.PDC %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("PDC Base vs M3, base = 92, m3 = 174, d28 = 191") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.PDC %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()



m3.vs.d28.PDC <- FindMarkers(m.only, ident.1 = "D28", ident.2 = "M3", group.by = "timepoint", subset.ident = "PDC")
m3.vs.d28.PDC <- m3.vs.d28.PDC %>% as_tibble(rownames = "gene")
m3.vs.d28.PDC <- m3.vs.d28.PDC %>% filter(p_val_adj < 0.8)

pdf("m3.vs.D28.PDC.binned.pdf")
ggplot(m3.vs.d28.PDC, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.d28.PDC %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.d28.PDC %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05), max.overlaps = 50) +
  ggtitle("PDC D28 vs M3, base = 92, m3 = 174, d28 = 191") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()


d28.vs.base.PDC <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "D28", group.by = "timepoint", subset.ident = "PDC")
d28.vs.base.PDC <- d28.vs.base.PDC %>% as_tibble(rownames = "gene")
d28.vs.base.PDC <- d28.vs.base.PDC %>% filter(p_val_adj < 0.8)

pdf("d28.vs.base.PDC.binned.pdf")
ggplot(d28.vs.base.PDC, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = d28.vs.base.PDC %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = d28.vs.base.PDC %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-2), max.overlaps = 50) +
  ggtitle("PDC base vs D28, base = 92, m3 = 174, d28 = 191") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()

pdf("PDC.changes.over.time.pdf")
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('CXCR4', 'RASD1', "DDIT4", "IRF4", "FCRLA"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("PDC")) + ggtitle("PDC changes; cell # base = 92, m3 = 174, d28 = 191")
dev.off()

pdf("DC2.changes.over.time.pdf")
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('CXCR4', 'DDIT4', "IL1R2", "CEBPD", "CEBPB", "MT.genes", "NRF1", "AREG"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("DC2")) + ggtitle("DC changes base 488, d28 901, m3 884")
dev.off()


FeaturePlot(m.only, features = c("CSF3R", "NAMPT", "PF4"))


#do these

# List of cell types to iterate over
cell_types <- c("Monocyte", "CD16+ Monocyte", "Neut", "MMP8+ CEACAM8+ Neut", 
                "DC2")

Idents(m.only) %>% table()

# Open a PDF device to save all plots in one file
plot.list <- list()

# Loop through each cell type
for (cell_type in cell_types) {
  
  # Find markers for M3 vs Base
  m3.vs.base <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", 
                            group.by = "timepoint", subset.ident = cell_type)
  m3.vs.base <- m3.vs.base %>% as_tibble(rownames = "gene") %>% filter(p_val_adj < 0.05)
  
  # Identify top 75 genes based on adjusted p-value
  top_genes_m3_base <- m3.vs.base %>% arrange(p_val_adj) %>% head(75)
  
  # Determine cutoffs for p-value and log2FC
  p_val_cutoff_m3_base <- max(top_genes_m3_base$p_val_adj)
  log2FC_cutoff_m3_base <- quantile(abs(top_genes_m3_base$avg_log2FC), 0.75)  # 75th percentile
  
  # Create plot for M3 vs Base and save to list
  plot_m3_base <- ggplot(m3.vs.base, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
    geom_point(color = "gray", alpha = 0.5) + 
    geom_point(data = m3.vs.base %>% filter(p_val_adj < p_val_cutoff_m3_base & abs(avg_log2FC) > log2FC_cutoff_m3_base)) + 
    theme_bw() + 
    geom_text_repel(aes(label = gene), data = m3.vs.base %>% filter(abs(avg_log2FC) > log2FC_cutoff_m3_base & p_val_adj < p_val_cutoff_m3_base), max.overlaps = 50) +
    ggtitle(paste("Base vs M3 for", cell_type)) +
    scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
  
  plot.list[[paste0("m3_vs_base_", cell_type)]] <- plot_m3_base
  
  # Find markers for M3 vs D28
  m3.vs.d28 <- FindMarkers(m.only, ident.1 = "D28", ident.2 = "M3", 
                           group.by = "timepoint", subset.ident = cell_type)
  m3.vs.d28 <- m3.vs.d28 %>% as_tibble(rownames = "gene") %>% filter(p_val_adj < 0.05)
  
  # Identify top 75 genes for M3 vs D28
  top_genes_m3_d28 <- m3.vs.d28 %>% arrange(p_val_adj) %>% head(75)
  
  # Determine cutoffs for p-value and log2FC
  p_val_cutoff_m3_d28 <- max(top_genes_m3_d28$p_val_adj)
  log2FC_cutoff_m3_d28 <- quantile(abs(top_genes_m3_d28$avg_log2FC), 0.05)  # 75th percentile
  
  # Create plot for M3 vs D28 and save to list
  plot_m3_d28 <- ggplot(m3.vs.d28, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
    geom_point(color = "gray", alpha = 0.5) + 
    geom_point(data = m3.vs.d28 %>% filter(p_val_adj < p_val_cutoff_m3_d28 & abs(avg_log2FC) > log2FC_cutoff_m3_d28)) + 
    theme_bw() + 
    geom_text_repel(aes(label = gene), data = m3.vs.d28 %>% filter(abs(avg_log2FC) > log2FC_cutoff_m3_d28 & p_val_adj < p_val_cutoff_m3_d28), max.overlaps = 50) +
    ggtitle(paste("D28 vs M3 for", cell_type)) +
    scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
  
  plot.list[[paste0("m3_vs_d28_", cell_type)]] <- plot_m3_d28
  
  # Find markers for D28 vs Base
  d28.vs.base <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "D28", 
                             group.by = "timepoint", subset.ident = cell_type)
  d28.vs.base <- d28.vs.base %>% as_tibble(rownames = "gene") %>% filter(p_val_adj < 0.8)
  
  # Identify top 75 genes for D28 vs Base
  top_genes_d28_base <- d28.vs.base %>% arrange(p_val_adj) %>% head(75)
  
  # Determine cutoffs for p-value and log2FC
  p_val_cutoff_d28_base <- max(top_genes_d28_base$p_val_adj)
  log2FC_cutoff_d28_base <- quantile(abs(top_genes_d28_base$avg_log2FC), 0.75)  # 75th percentile
  
  # Create plot for D28 vs Base and save to list
  plot_d28_base <- ggplot(d28.vs.base, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
    geom_point(color = "gray", alpha = 0.5) + 
    geom_point(data = d28.vs.base %>% filter(p_val_adj < p_val_cutoff_d28_base & abs(avg_log2FC) > log2FC_cutoff_d28_base)) + 
    theme_bw() + 
    geom_text_repel(aes(label = gene), data = d28.vs.base %>% filter(abs(avg_log2FC) > log2FC_cutoff_d28_base & p_val_adj < p_val_cutoff_d28_base), max.overlaps = 50) +
    ggtitle(paste("Base vs D28 for", cell_type)) +
    scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
  
  plot.list[[paste0("d28_vs_base_", cell_type)]] <- plot_d28_base
}

# Optionally, save all plots to a single PDF file
pdf("all_cell_types_plots.pdf")
for (plot_name in names(plot.list)) {
  print(plot.list[[plot_name]])
}
dev.off()


# For Mast
m3.vs.base.Mast <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Mast")
m3.vs.base.Mast <- m3.vs.base.Mast %>% as_tibble(rownames = "gene")
m3.vs.base.Mast <- m3.vs.base.Mast %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Mast, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Mast %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Mast %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.Mast %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Mast Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.Mast %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Megakaryocyte
m3.vs.base.Megakaryocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Megakaryocyte")
m3.vs.base.Megakaryocyte <- m3.vs.base.Megakaryocyte %>% as_tibble(rownames = "gene")
m3.vs.base.Megakaryocyte <- m3.vs.base.Megakaryocyte %>% filter(p_val_adj < 0.05)

pdf("M3.vs.Base.Megakaryocyte.w.violin.pdf")
ggplot(m3.vs.base.Megakaryocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Megakaryocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Megakaryocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.Megakaryocyte %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Megakaryocyte Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.Megakaryocyte %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

ggplot(m3.vs.base.Megakaryocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Megakaryocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Megakaryocyte %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5), max.overlaps = 50) +
  ggtitle("Megakaryocyte Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('ZBTB16', "IFI27", "SIGLEC1", "MX1", "IRF7", "LGALS1", "PF4", "GP1BB", "nCount_RNA", "MT.genes", "RB.genes"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Megakaryocyte", "Platelet")) + ggtitle("Megakaryocyte changes")


VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c('ZBTB16', "IFI27", "SIGLEC1", "MX1", "IRF7", "LGALS1"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Megakaryocyte")) + ggtitle("Megakaryocyte changes")

VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", c("IFI27", "PF4", "GP1BB", "MT.genes"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = c("Platelet")) + ggtitle("Megakaryocyte changes")

dev.off()



# For Granulocyte
m3.vs.base.Granulocyte <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Granulocyte")
m3.vs.base.Granulocyte <- m3.vs.base.Granulocyte %>% as_tibble(rownames = "gene")
m3.vs.base.Granulocyte <- m3.vs.base.Granulocyte %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Granulocyte, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Granulocyte %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene)) +
  geom_point(data = m3.vs.base.Granulocyte %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Granulocyte Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.Granulocyte %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For Platelet
m3.vs.base.Platelet <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Platelet")
m3.vs.base.Platelet <- m3.vs.base.Platelet %>% as_tibble(rownames = "gene")
m3.vs.base.Platelet <- m3.vs.base.Platelet %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.Platelet, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Platelet %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Platelet %>% filter(abs(avg_log2FC) > 0.25 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.Platelet %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("Platelet Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.Platelet %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

# For RBC
m3.vs.base.RBC <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "RBC")
m3.vs.base.RBC <- m3.vs.base.RBC %>% as_tibble(rownames = "gene")
m3.vs.base.RBC <- m3.vs.base.RBC %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.base.RBC, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.RBC %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.RBC %>% filter(abs(avg_log2FC) > 0.25 & p_val_adj < 1e-2 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.RBC %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("RBC Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.RBC %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed



save.image("Sep23_myeloid_image.with.m.only.rb.mt.scored.rdata")


#score by IFN



#rank by AUCell, score IFN-----


molecular_signatures <- list(
  Bennett = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "DNAPT6", "EIF2AK2", "FCGR1A", "HERC6"),
  
  Baechler = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "EIF2AK2", "EPSTI1", "FCGR1A", 
               "FLJ20035", "HERC5", "HERC6", "IFI6", "IFI27", "IFI35", "IFI44"),
  
  Kirou = c("IFI6", "IFI27", "IFI35"),
  
  Feng = c("IFI6", "IFI27", "IFI35", "IFI44", "IFI44L"),
  
  Nikpour = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "DNAPT6", "EIF2AK2", "EPSTI1", 
              "FCGR1A", "FLJ20035", "HERC5", "HERC6", "IFI6", "IFI27", "IFI35", 
              "IFI44", "IFI44L", "IFIH1", "ISG15", "LAMP3", "LAP3", "LGALS3BP", 
              "LGP2", "LOC129607", "LY6E", "MX1", "MX2", "OAS1", "OAS2", "OAS3", 
              "OASL", "PML", "PLSCR1", "RNASE2", "RSAD2", "RTP4", "SERPING1", 
              "SIGLEC1", "SP110", "STAT1", "STAT2", "TAP1", "UBE2L6", "USP18", 
              "XAF1"),
  
  Landolt = c("IFI6", "IFI27", "IFI35", "IFI44", "IFI44L"),
  
  Petri = c("IFI6", "IFI27", "IFI35"),
  
  Yao = c("APOBEC1-like", "BST2", "C7orf6", "CIC", "DNAPT6", "EIF2AK2", "EPSTI1", 
          "FCGR1A", "FLJ20035", "HERC5", "HERC6", "IFI6", "IFI27", "IFI35", 
          "IFI44", "IFI44L", "IFIH1", "ISG15", "LAMP3", "LAP3", "LGALS3BP", 
          "LGP2"),
  
  Higgs = c("IFI6", "IFI27", "IFI35", "IFI44", "RSAD2"),
  
  banchereau.ifn = c("ABCA1", "BATF2", "BTN3A1", "CARD16", "CARD17", "CCL8", "CEACAM1", "CMPK2", "DHRS9", "DYNLT1", "GALM", "GBP1", "GBP1", "GBP2", "HELZ2", "IFI44", "IFIH1", "IFIT3", 
                     "IFITM1", "IRF7", "ISG20", "MX1", "NBN", "OAS1", "OAS1", "OAS3", "PARP9", "SP110", "TMEM140", "TNFAIP6", "TNFSF10", "TRIM21", "TRIM22", "TRIM5", "TRIM5", "UBE2L6", "ZBP1", "ZNF684"),
  
  banchereau.pc = c("CRELD2", "DNAJB11", "FAM98A", "MYDGF", "NME1", "RFC4", "SDF2L1", "SRM", "STT3A"),
  
  M1.2 = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "OAS1", "OAS2", 
           "OAS3", "MX1", "MX2", "RSAD2", "STAT1", "STAT2", "STAT3", "STAT4", 
           "STAT5A", "STAT5B", "IRF1", "IRF7", "IRF9", "IFNAR1", "IFNAR2", "IFNGR1", 
           "IFNGR2", "IL6", "IL10", "IL12A", "IL12B", "IL18", "TNF", "TNFRSF1A", 
           "TNFRSF1B", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4", "CCL5", 
           "CCL8", "CCL19", "CCL21", "CCL22", "CCL23", "CCL24", "CCL26", 
           "CCL27", "CCL28", "CXCL9", "CXCL12", "CXCL13", "CXCL14", "CXCL15", 
           "CXCL16", "CXCL17", "CXCL18", "CXCL19", "CXCL20", "CXCL21", "CXCL22", 
           "CXCL23", "CXCL24", "CXCL25", "CXCL26", "CXCL27", "CXCL28", "CXCL29"),
  
  M3.4 = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "OAS1", "OAS2", 
           "OAS3", "MX1", "MX2", "RSAD2", "STAT1", "STAT2", "STAT3", "STAT4", 
           "STAT5A", "STAT5B", "IRF1", "IRF7", "IRF9", "IFNAR1", "IFNAR2", "IFNGR1", 
           "IFNGR2", "IL6", "IL10", "IL12A", "IL12B", "IL18", "TNF", "TNFRSF1A", 
           "TNFRSF1B", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4", "CCL5", 
           "CCL8", "CCL19", "CCL21", "CCL22", "CCL23", "CCL24", "CCL26", 
           "CCL27", "CCL28", "CXCL9", "CXCL12", "CXCL13", "CXCL14", "CXCL15", 
           "CXCL16", "CXCL17", "CXCL18", "CXCL19", "CXCL20", "CXCL21", "CXCL22", 
           "CXCL23", "CXCL24", "CXCL25", "CXCL26", "CXCL27", "CXCL28", "CXCL29"),
  
  M5.12 = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "OAS1", "OAS2", 
            "OAS3", "MX1", "MX2", "RSAD2", "STAT1", "STAT2", "STAT3", "STAT4", 
            "STAT5A", "STAT5B", "IRF1", "IRF7", "IRF9", "IFNAR1", "IFNAR2", "IFNGR1", 
            "IFNGR2", "IL6", "IL10", "IL12A", "IL12B", "IL18", "TNF", "TNFRSF1A", 
            "TNFRSF1B", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4", "CCL5", 
            "CCL8", "CCL19", "CCL21", "CCL22", "CCL23", "CCL24", "CCL26", 
            "CCL27", "CCL28", "CXCL9", "CXCL12", "CXCL13", "CXCL14", "CXCL15", 
            "CXCL16", "CXCL17", "CXCL18", "CXCL19", "CXCL20", "CXCL21", "CXCL22", 
            "CXCL23", "CXCL24", "CXCL25", "CXCL26", "CXCL27", "CXCL28", "CXCL29")
)


names(molecular_signatures) = paste0(names(molecular_signatures), ".", lapply(molecular_signatures, length))


cells_rankings_vst <- AUCell_buildRankings(m.only@assays$RNA$counts)

cells_AUC.IFN <- AUCell_calcAUC(molecular_signatures, cells_rankings_vst, 
                                aucMaxRank=nrow(cells_rankings_vst)*0.25, verbose = T)


m.only[["AUC_IFN"]] <- CreateAssayObject(cells_AUC.IFN@assays@data$AUC)

saveRDS(m.only, "m.only.with.AUC.IFN.rds")

m.only@assays$AUC_IFN




m.only$cell.type.jx.1.m.binned.factor

VlnPlot(m.only, group.by = "cell.type.jx.1.m.factor", features = c(names(molecular_signatures)), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")


VlnPlot(subset(b.only.jx.anno, cell.type.jx.0.5 %in% c("IgA+/IgG+ PC", "PB (Cycling)")), group.by = "cell.type.jx.0.5", features = c(names(molecular_signatures)), 
        stack = T, flip = T, fill.by = "ident", cols = custom_palette, split.by = "timepoint")

pdf("IFN.gene.set.vln.over.time.pdf")
VlnPlot(m.only, group.by = "cell.type.jx.1.m.binned.factor", features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")

VlnPlot(subset(b.only.jx.anno, !cell.type.jx.0.5 %in% c("IgA+/IgG+ PC", "PB (Cycling)") & timepoint %in% c("Base", "M3")), group.by = "cell.type.jx.0.5", features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")

dev.off()


pdf("IFN.gene.set.vln.over.time.pdf")
VlnPlot(subset(b.only.jx.anno, cell.type.jx.0.5 %in% c("IgA+/IgG+ PC", "PB (Cycling)")), group.by = "cell.type.jx.0.5", features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")

VlnPlot(subset(b.only.jx.anno, !cell.type.jx.0.5 %in% c("IgA+/IgG+ PC", "PB (Cycling)") & timepoint %in% c("Base", "M3")), group.by = "cell.type.jx.0.5", features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")

VlnPlot(subset(b.only.jx.anno, timepoint %in% c("Base", "M3")), group.by = "cell.type.jx.0.5", features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")
dev.off()


pdf("IFN.gene.set.vln.over.dose.pdf")
VlnPlot(subset(b.only.jx.anno, cell.type.jx.0.5 %in% c("IgA+/IgG+ PC", "PB (Cycling)")), group.by = "dose", features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint") + ggtitle("reduction in IFN in plasmablasts by dose")

dev.off()




#Analyze Neutrophils only
DimPlot(m.only, group.by = "RNA_snn_res.1", label = T)
m.only$cell.type.jx.1.m[m.only$RNA_snn_res.1 == "19"] = "Platelet/Neut Agg"
m.only$cell.type.jx.1.m.factor = factor(m.only$cell.type.jx.1.m, levels = c("Monocyte", "SIGLEC10+/TCF7L2+ Monocyte", "ZEB2-hi Monocyte", "Neut", "MMP8+ CEACAM8+ Neut", "DC2", "PDC" ,"Mast", "Megakaryocyte", "Granulocyte", "Platelet", "Platelet/Neut Agg",  "RBC"))
DimPlot(m.only, group.by = "cell.type.jx.1.m.factor", label = T)

neut.only = subset(m.only, cell.type.jx.1.m.factor %in% c("Neut", "MMP8+ CEACAM8+ Neut"))

neut.only = neut.only %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:20)
neut.only = neut.only %>% FindNeighbors() %>% FindClusters(resolution = c(0.5, 1, 2))

DimPlot(neut.only, group.by = "RNA_snn_res.1", label = T) + NoLegend()+
DimPlot(neut.only, group.by = "cell.type.jx.1.m.factor", label = T) + NoLegend() +
  DimPlot(neut.only, group.by = "timepoint", label = T) + NoLegend() +
  DimPlot(neut.only, group.by = "patient", label = T) + NoLegend()

marker_genes <- list(
  Nh0 = c("MMP9", "ITGAM", "FCN1", "CAMP", "CYBB", "CST3", "VIM", "TXN", "S100A6", "S100A8", "S100A9", "S100A11", "S100A12", "MME"),
  Nh1 = c("AIF1", "CXCR2", "TXNIP"),
  Nh2 = c("MALAT1", "NEAT1", "CSF3R"),
  Nh3 = c("HERC5", "IFI16", "IFIT1", "IFIT2", "IFITM2", "IFITM3", "ISG15")
)

marker_genes = marker_genes %>% unlist()

FeaturePlot(neut.only, marker_genes, ) + NoLegend()
DimPlot(neut.only, group.by = "RNA_snn_res.0.5", label = T) + NoLegend() +
VlnPlot(neut.only, marker_genes, group.by = "RNA_snn_res.0.5", flip = T, stack = T) + NoLegend()

c9.markers = FindMarkers(neut.only, ident.1 = "9", ident.2 = "8", group.by = "RNA_snn_res.0.5")


c9.markers <- c9.markers %>% as_tibble(rownames = "gene")
c9.markers <- c9.markers %>% filter(p_val_adj < 0.05)

pdf("Progenitor.vs.immature.pdf")
ggplot(c9.markers, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = c9.markers %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = c9.markers %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-100), max.overlaps = 50) +
  ggtitle("Progenitor vs Immature") +
  scale_x_continuous(breaks = seq(-15, 15, by = 1))  # Adjust the range as needed
dev.off()

neut.only@meta.data <- 
  neut.only@meta.data %>% mutate(cell.type.jx.0.5.neut = 
                                     case_when(
                                       RNA_snn_res.0.5 == 0 ~ "IFI-hi",
                                       RNA_snn_res.0.5 == 1 ~ "IFI-mid",
                                       RNA_snn_res.0.5 == 2 ~ "Canonical",
                                       RNA_snn_res.0.5 == 3 ~ "Maturing",
                                       RNA_snn_res.0.5 == 4 ~ "IFI-hi",
                                       RNA_snn_res.0.5 == 5 ~ "Maturing",
                                       RNA_snn_res.0.5 == 6 ~ "IFI-mid",
                                       RNA_snn_res.0.5 == 7 ~ "IFI-mid",
                                       RNA_snn_res.0.5 == 8 ~ "Immature",
                                       RNA_snn_res.0.5 == 9 ~ "DEFA3/4/1B/1+ FCER1G-low Neut Prog"
                                     ))


neut.only@meta.data <- 
  neut.only@meta.data %>% mutate(cell.type.jx.0.5.neut.binned = 
                                   case_when(
                                     RNA_snn_res.0.5 == 0 ~ "IFI+",
                                     RNA_snn_res.0.5 == 1 ~ "IFI+",
                                     RNA_snn_res.0.5 == 2 ~ "Canonical",
                                     RNA_snn_res.0.5 == 3 ~ "Canonical",
                                     RNA_snn_res.0.5 == 4 ~ "IFI+",
                                     RNA_snn_res.0.5 == 5 ~ "Canonical",
                                     RNA_snn_res.0.5 == 6 ~ "IFI+",
                                     RNA_snn_res.0.5 == 7 ~ "IFI+",
                                     RNA_snn_res.0.5 == 8 ~ "Immature",
                                     RNA_snn_res.0.5 == 9 ~ "DEFA3/4/1B/1+ FCER1G-low Neut Prog"
                                   ))


neut.only$cell.type.jx.0.5.neut.factor = factor(neut.only$cell.type.jx.0.5.neut, levels = c("DEFA3/4/1B/1+ FCER1G-low Neut Prog", 
                                                                                                          "Immature", "Maturing",
                                                                                                          "Canonical",
                                                                                                          "IFI-mid", 
                                                                                                          "IFI-hi"))

neut.only$cell.type.jx.0.5.neut.binned.factor = factor(neut.only$cell.type.jx.0.5.neut.binned, levels = c("DEFA3/4/1B/1+ FCER1G-low Neut Prog", 
  "Immature", 
  "Canonical",
  "IFI+"))

neut.colors = RColorBrewer::brewer.pal(n = 6, name = "Set3")
neut.colors2 = RColorBrewer::brewer.pal(n = 4, name = "Paired")

pdf("neut_UMAP.expression.pdf")
DimPlot(neut.only, group.by = "cell.type.jx.0.5.neut.factor", label = T, raster = T, cols = neut.colors, raster.dpi = c(600,600)) + NoLegend()
DimPlot(neut.only, group.by = "cell.type.jx.0.5.neut.binned", label = T, raster = T, cols = neut.colors2, raster.dpi = c(600,600)) + NoLegend()
  DimPlot(neut.only, group.by = "patient", label = T, raster = T, raster.dpi = c(600,600)) + NoLegend()
DimPlot(neut.only, group.by = "timepoint", label = T, raster = T, raster.dpi = c(600,600)) + NoLegend()
VlnPlot(neut.only, c("DEFA3", "DEFA4", "FCER1G", "S100A12", marker_genes, "MT.genes", "nCount_RNA"), cols = neut.colors ,group.by = "cell.type.jx.0.5.neut.factor", flip = T, stack = T, fill.by = "ident") + NoLegend()
VlnPlot(neut.only,  c("DEFA3", "DEFA4", "FCER1G", "S100A12", marker_genes, "MT.genes", "nCount_RNA"), cols = neut.colors2,  group.by = "cell.type.jx.0.5.neut.binned.factor", flip = T, stack = T, fill.by = "ident") + NoLegend()
dev.off()

pdf('neut.by.pt.time.pdf')

DimPlot(neut.only, group.by = "patient", label = T, raster = T, raster.dpi = c(600,600))
DimPlot(neut.only, group.by = "timepoint", label = T, raster = T, raster.dpi = c(600,600)) 
dev.off()
# Generate UMAP plot over time for neutrophils
pdf("UMAP.over.time.neut.pdf", height = 8, width = 24)
DimPlot(neut.only, group.by = "cell.type.jx.0.5.neut", label = TRUE, cols = myeloid.colors, split.by = "timepoint")
dev.off()

# Summarize statistics for the binned cell types for neutrophils
neut.summ.stats.res.1 <- neut.only@meta.data %>%
  group_by(timepoint, cell.type.jx.0.5.neut.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(neut.only@meta.data %>%
              group_by(timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

neut.summ.stats.res.2 <- neut.only@meta.data %>%
  group_by(timepoint, cell.type.jx.0.5.neut.binned.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(neut.only@meta.data %>%
              group_by(timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

# Create bar plot for changes over time for neutrophils
pdf("neut.bar.changes.over.time.neut.pdf")
ggplot(neut.summ.stats.res.1, aes(x = timepoint, y = prop, fill = cell.type.jx.0.5.neut.factor)) + 
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = neut.colors)+
  geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) + coord_flip() + theme(legend.position = "none")


ggplot(neut.summ.stats.res.2, aes(x = timepoint, y = prop, fill = cell.type.jx.0.5.neut.binned.factor)) + 
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = neut.colors2)+
  geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) + coord_flip() + theme(legend.position = "none")


dev.off()

Idents(neut.only) = neut.only$cell.type.jx.0.5.neut.binned.factor 

# For Canonical
m3.vs.base.Canonical <- FindMarkers(neut.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Canonical", 
                                    max.cells.per.ident = 5000)
m3.vs.base.Canonical <- m3.vs.base.Canonical %>% as_tibble(rownames = "gene")
m3.vs.base.Canonical <- m3.vs.base.Canonical %>% filter(p_val_adj < 0.05)

pdf("Canonical.m3.vs.base.pdf")
ggplot(m3.vs.base.Canonical, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Canonical %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Canonical %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50), max.overlaps = 50) +
  ggtitle("Canonical Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

VlnPlot(subset(neut.only, cell.type.jx.0.5.neut.binned.factor %in% c("Canonical", "IFI+")),  
        c("IL1R2", "TLR1", "TLR2", "TLR4", "CXCR4", "CEBPD", "ZFP36","CCND3", "FKBP5", "TPST1", "SLC8A1", "SMAP2", "TAGLN2"), 
        cols = c("#1f77b4", "#ff7f0e", "#2ca02c"),
        group.by = "cell.type.jx.0.5.neut.binned.factor", flip = T, stack = T, fill.by = "ident", split.by = "timepoint")

VlnPlot(subset(neut.only, cell.type.jx.0.5.neut.binned.factor %in% c("Canonical")),  
        c("IL1R2", "TLR1", "TLR2", "TLR4", "CXCR4", "CEBPD", "ZFP36","CCND3", "FKBP5", "TPST1", "SLC8A1", "SMAP2", "TAGLN2", "CXCL8"), 
        group.by = "cell.type.jx.0.5.neut.binned.factor", flip = T, stack = T, cols = c("#1f77b4", "#ff7f0e", "#2ca02c"), split.by = "timepoint")
dev.off()


# For IFI+
m3.vs.base.IFI <- FindMarkers(neut.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "IFI+", max.cells.per.ident = 10000)
m3.vs.base.IFI <- m3.vs.base.IFI %>% as_tibble(rownames = "gene")
m3.vs.base.IFI <- m3.vs.base.IFI %>% filter(p_val_adj < 0.05)


pdf("IFI_b_vs_m3.pdf")
ggplot(m3.vs.base.IFI, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.IFI %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.IFI %>% filter(abs(avg_log2FC) > 0.4 & p_val_adj < 1e-60 & !gene %in% isg.consensus), max.overlaps = 50) +
  ggtitle("IFI+ Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  geom_point(data = m3.vs.base.IFI %>% filter(gene %in% isg.consensus), color = "red") + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.IFI %>% filter(gene %in% isg.consensus, p_val_adj < 1e-50), color= "red")# Adjust the range as needed


VlnPlot(subset(neut.only, cell.type.jx.0.5.neut.binned.factor %in% c("IFI+")),  
        c("IL1R2", "CEBPD", "FCGR3A", "TLR2", "TLR4", "CXCR4", "FP671120.4", "CCND3", "SMAP2", "MX1", "RSAD2", "ISG15"), 
        cols = c("#1f77b4", "#ff7f0e", "#2ca02c"),
        group.by = "cell.type.jx.0.5.neut.binned.factor", flip = T, stack = T, fill.by = "ident", split.by = "timepoint")

dev.off()

# For DEFA3/4/1B/1+ FCER1G-low Neut Prog
m3.vs.base.DEFA <- FindMarkers(neut.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "DEFA3/4/1B/1+ FCER1G-low Neut Prog")
m3.vs.base.DEFA <- m3.vs.base.DEFA %>% as_tibble(rownames = "gene")
m3.vs.base.DEFA <- m3.vs.base.DEFA %>% filter(p_val_adj < 0.05)

ggplot(m3.vs.base.DEFA, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.DEFA %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.DEFA %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-10 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  ggtitle("DEFA3/4/1B/1+ FCER1G-low Neut Prog Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed




# For Immature
m3.vs.base.Immature <- FindMarkers(neut.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Immature")
m3.vs.base.Immature <- m3.vs.base.Immature %>% as_tibble(rownames = "gene")
m3.vs.base.Immature <- m3.vs.base.Immature %>% filter(p_val_adj < 0.05)

ggplot(m3.vs.base.Immature, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.Immature %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.Immature %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-10 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  ggtitle("Immature Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed



# For MMP8+ CEACAM8+ Neut
m3.vs.D28.MMP8_CEACAM8_Neut <- FindMarkers(m.only, ident.1 = "Base", ident.2 = "D28", group.by = "timepoint", subset.ident = "MMP8+ CEACAM8+ Neut")
m3.vs.D28.MMP8_CEACAM8_Neut <- m3.vs.D28.MMP8_CEACAM8_Neut %>% as_tibble(rownames = "gene")
m3.vs.D28.MMP8_CEACAM8_Neut <- m3.vs.D28.MMP8_CEACAM8_Neut %>% filter(p_val_adj < 0.05)
ggplot(m3.vs.D28.MMP8_CEACAM8_Neut, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.D28.MMP8_CEACAM8_Neut %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.D28.MMP8_CEACAM8_Neut %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 1e-5 | p_val_adj < 1e-20, !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.D28.MMP8_CEACAM8_Neut %>% filter(gene %in% isg.consensus), color = "red") + 
  ggtitle("MMP8+ CEACAM8+ Neut Base vs D28") +
  geom_text_repel(aes(label = gene), data = m3.vs.D28.MMP8_CEACAM8_Neut %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

VlnPlot(neut.only, features = c("FP671120.4", "FP236383.3", "MKI67", "MT.genes"), group.by = "cell.type.jx.0.5.neut.binned.factor", stack = T, flip = T, split.by = "timepoint")


cycle3 = readLines('/Users/jasonxu/Library/CloudStorage/Box-Box/ETP_SingleCell_Project/GCB_Courses/536_kim/traf6/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt')
s.genes = cycle3[1:43]
g2m.genes = cycle3[44:97]

neut.only = neut.only %>% CellCycleScoring(s.features = s.genes,
                                           g2m.features = g2m.genes, set.ident = F)
VlnPlot(neut.only, features = c("S.Score", "G2M.Score"), group.by = "cell.type.jx.0.5.neut.binned.factor", stack = T, flip = T, split.by = "timepoint")


cycling_cells_proportion <- neut.only@meta.data %>%
  # Filter for cycling cells (assuming cycling cells are labeled as "G1", "S", or "G2M" in the Phase column)
  filter(Phase %in% c("G1", "S", "G2M")) %>%
  group_by(patient, timepoint, cell.type.jx.0.5.neut.binned.factor, .drop = F) %>%
  summarise(
    total_cells = n(),
    cycling_cells = sum(Phase %in% c("S")),
    proportion_cycling = cycling_cells / total_cells
  ) %>%
  ungroup()



saveRDS(neut.only, "neut.only.withIFN.rds")



#monocyte only
DimPlot(m.only, group.by = "cell.type.jx.1.m.factor", label = T)


mono.only = subset(m.only, cell.type.jx.1.m.factor %in% c("Monocyte", "SIGLEC10+/TCF7L2+ Monocyte", "ZEB2-hi Monocyte", "Granulocyte"))


mono.only = mono.only %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25) %>% RunUMAP(dims = 1:10)
mono.only = mono.only %>% FindNeighbors() %>% FindClusters(resolution = c(0.5, 1, 2))
mono.only = mono.only %>% FindClusters(resolution = c(0.1, 0.25))
mono.only = mono.only %>% FindClusters(resolution = c(0.15, 0.2))
mono.only = mono.only %>% FindClusters(resolution = c(0.3, 0.4))

DimPlot(mono.only, group.by = "RNA_snn_res.0.5", label = T) + NoLegend()+
  DimPlot(mono.only, group.by = "cell.type.jx.1.m.factor", label = T) + NoLegend() +
  DimPlot(mono.only, group.by = "timepoint", label = T) + NoLegend() +
  DimPlot(mono.only, group.by = "patient", label = T) + NoLegend()

DimPlot(mono.only, group.by = "RNA_snn_res.0.25", label = T)
DimPlot(mono.only, group.by = "RNA_snn_res.0.5", label = T)

DimPlot(mono.only, group.by = "RNA_snn_res.1", label = T)


Idents(mono.only) = mono.only$RNA_snn_res.0.25
mono.markers = FindAllMarkers(mono.only)

mono.markers %>%
  group_by(cluster) %>%
  slice_head(n = 30) %>%
  ungroup() %>% print(n=300)




FeaturePlot(mono.only, "IFI44L")

monocyte_marker_genes_hgnc <- c(
  "CD14", "CD68", "CSF1R", "IL1B", "TNF", "CCL2", "CX3CR1", "S100A8", "S100A9", "FCGR1A",  # Classical Monocytes
  "FCGR3A", "HLA-DRA", "CD86", "IL6", "CCL3", "CCL4", "CXCL10",                          # Intermediate Monocytes
  "CCR2", "IL10", "IL12B", "FCGR1"                                                      # Non-Classical Monocytes
) %>% unique()
FeaturePlot(mono.only, monocyte_marker_genes)

# Classical Monocytes (CD14++CD16−)
classical_monocytes_genes <- c("CD14", "CD68", "CSF1R", "IL1B", "TNF", "CCL2", "CX3CR1", "S100A8", "S100A9", "FCGR1A")

# Intermediate Monocytes (CD14++CD16+)
intermediate_monocytes_genes <- c("CD14", "FCGR3A", "HLA-DRB1", "CD86", "IL6", "TNF", "CCL3", "CCL4", "CXCL10", "FCGR3A")

# Non-Classical Monocytes (CD14+CD16+)
non_classical_monocytes_genes <- c("CD14", "CD16", "FCGR3A", "CCR2", "CX3CR1", "S100A8", "S100A9", "IL10", "IL12B", "TNF", "CD64")

FeaturePlot(mono.only, non_classical_monocytes_genes)
FeaturePlot(mono.only, classical_monocytes_genes)

FeaturePlot(mono.only, c("CXCL5", "CXCL3", "CXCL1", "CCL2", "PPBP", "TGFBI", "FN1", "RNASE1", "S100A8"))
FeaturePlot(mono.only, c("MX1", "CD14", "FCGR3A", "C1QA", "IL1B", "CXCL8", "ISG15", "IFITM3"))

FeaturePlot(mono.only, c("HLA.genes", "MT.genes", "nCount_RNA"))



DimPlot(mono.only, group.by = "RNA_snn_res.0.25", label = T)

mono.only@meta.data <- 
  mono.only@meta.data %>% mutate(cell.type.jx.0.25.mono = 
                                   case_when(
                                     RNA_snn_res.0.25 == 0 ~ "IFI+",
                                     RNA_snn_res.0.25 == 1 ~ "HLA/MARCO-hi",
                                     RNA_snn_res.0.25 == 2 ~ "Classical",
                                     RNA_snn_res.0.25 == 3 ~ "CD16+ Non-Classical",
                                     RNA_snn_res.0.25 == 4 ~ "CCL8/CCL2/CXCL10+",
                                     RNA_snn_res.0.25 == 5 ~ "CCL3L1/CCL3+/CD83+",
                                     RNA_snn_res.0.25 == 6 ~ "Granulocyte",
                                     RNA_snn_res.0.25 == 7 ~ "MT-high/ZEB2-hi/Actin-hi",
                                   ))


mono.only@meta.data <- 
  mono.only@meta.data %>% mutate(cell.type.jx.0.25.mono.binned = 
                                   case_when(
                                     RNA_snn_res.0.25 == 0 ~ "IFI+",
                                     RNA_snn_res.0.25 == 1 ~ "HLA/MARCO-hi",
                                     RNA_snn_res.0.25 == 2 ~ "Classical",
                                     RNA_snn_res.0.25 == 3 ~ "CD16+ Non-Classical",
                                     RNA_snn_res.0.25 == 4 ~ "IFI+",
                                     RNA_snn_res.0.25 == 5 ~ "IFI+",
                                     RNA_snn_res.0.25 == 6 ~ "Granulocyte",
                                     RNA_snn_res.0.25 == 7 ~ "MT-high/ZEB2-hi/Actin-hi",
                                   ))



mono.only@meta.data <- 
  mono.only@meta.data %>% mutate(cell.type.jx.0.25.mono.binned2 = 
                                   case_when(
                                     RNA_snn_res.0.25 == 0 ~ "IFI+",
                                     RNA_snn_res.0.25 == 1 ~ "Classical",
                                     RNA_snn_res.0.25 == 2 ~ "Classical",
                                     RNA_snn_res.0.25 == 3 ~ "CD16+ Non-Classical",
                                     RNA_snn_res.0.25 == 4 ~ "IFI+",
                                     RNA_snn_res.0.25 == 5 ~ "IFI+",
                                     RNA_snn_res.0.25 == 6 ~ "Granulocyte",
                                     RNA_snn_res.0.25 == 7 ~ "Classical",
                                   ))



#mono.only = mono.only %>% RunUMAP(dims = 1:10)

DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono", label = T) +
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned", label = T) +
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned2", label = T)

pdf("Mono.umap.pdf", width = 5, height = 6)
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono", label = T, raster = T, raster.dpi = c(500,500)) + theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned", label = T, raster = T , raster.dpi = c(500,500))+ theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned2", label = T, raster = T, raster.dpi = c(500,500))+ theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "patient", label = T, raster = T, raster.dpi = c(500,500))+ theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "timepoint", label = T, raster = T, raster.dpi = c(500,500))+ theme(legend.position = "bottom")
dev.off()

pdf("Mono.Expression.pdf", width = 4, height = 6)

VlnPlot(mono.only, 
        features = c("FCGR3A", "C1QA", "FCGR3B",  "FCGR1A", "CD14", "LYZ", "S100A8", "S100A9",   "NAMPT", 
                     "SELL",  "MME","HLA.genes", "MARCO", "LDHB",  "CX3CR1",
                     "MX1", "IFI44L", "IFI6", "SIGLEC1", "IFI44", "ISG15", "IFITM3", "LY6E", "IL17RA",  
                      "ZEB2", "MT.genes","SLC8A1", "ARHGAP26","AOAH", "DOCK2",  
                      "nCount_RNA"),
        flip = T, group.by = "cell.type.jx.0.25.mono.binned", stack = T, fill.by = "ident") + NoLegend()

VlnPlot(mono.only, 
        features = c("FCGR3A",  "CX3CR1", "C1QA", "FCGR3B",  "SELL", "MME", "CD14", "LYZ", "S100A8", "S100A9",   "NAMPT", 
                    "MX1", "IFI44L", "IFI6", "SIGLEC1", "IFI44", "ISG15", "IFITM3", "LY6E", "IL17RA"),
        flip = T, group.by = "cell.type.jx.0.25.mono.binned2", stack = T, fill.by = "ident")  + NoLegend()

dev.off()

DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono", label = T, raster = T, raster.dpi = c(500,500)) + theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned", label = T, raster = T , raster.dpi = c(500,500))+ theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned2", label = T, raster = T, raster.dpi = c(500,500))+ theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "patient", label = T, raster = T, raster.dpi = c(500,500))+ theme(legend.position = "bottom")
DimPlot(mono.only, group.by = "timepoint", label = T, raster = T, raster.dpi = c(500,500))+ theme(legend.position = "bottom")
dev.off()



FeaturePlot(mono.only, c("MX1", "CCL8", "CCL2", "CXCL10", "MARCO", "CCL3L1", "CCL3", "CD83"))


# Summarize statistics for the binned cell types for monocytes
mono.summ.stats.res.1 <- mono.only@meta.data %>%
  group_by(timepoint, cell.type.jx.0.25.mono, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(mono.only@meta.data %>%
              group_by(timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

DimPlot(mono.only, group.by = "cell.type.jx.0.25.mono", label = T)

mono.summ.stats.res.2 <- mono.only@meta.data %>%
  group_by(timepoint, cell.type.jx.0.25.mono.binned2, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(mono.only@meta.data %>%
              group_by(timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)



  ggplot(mono.summ.stats.res.2, aes(x = timepoint, y = prop, fill = cell.type.jx.0.25.mono.binned2)) + 
  geom_col() + 
  theme_bw() + 
  #scale_fill_manual(values = mono.colors) +  # Ensure you have defined mono.colors
  geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) + 
  coord_flip()

  
  
ggplot(mono.summ.stats.res.1, aes(x = timepoint, y = prop, fill = cell.type.jx.0.25.mono)) + 
  geom_col() + 
  theme_bw() + 
  #scale_fill_manual(values = mono.colors) +  # Ensure you have defined mono.colors
  geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5)) + 
  coord_flip()




# Summarize statistics for monocytes
mono.summ.stats <- mono.only@meta.data %>%
  group_by(patient, timepoint, cell.type.jx.0.25.mono, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(mono.only@meta.data %>%
              group_by(patient, timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

# Create plots for monocyte cells
pdf("m.summ.stats.mono.2.pdf", height = 5, width = 15)
ggplot(mono.summ.stats, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) + 
  geom_line(aes(group = patient)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~cell.type.jx.0.25.mono, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3"))) +
  theme_bw() + 
  ggtitle("with unpaired p-values")

dev.off()



coord_flip()



mono.only@meta.data$cell.type.jx.0.25.mono.binned2.factor = factor(mono.only@meta.data$cell.type.jx.0.25.mono.binned2, levels = c("IFI+", "Classical", "CD16+ Non-Classical", "Granulocyte"))
# Summarize statistics for monocytes
mono.summ.stats2 <- mono.only@meta.data %>%
  group_by(patient, timepoint, cell.type.jx.0.25.mono.binned2.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(mono.only@meta.data %>%
              group_by(patient, timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

# Create plots for monocyte cells
pdf("m.summ.stats.mono.2.pdf", height = 5, width = 8)

ggplot(mono.summ.stats2, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) + 
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~cell.type.jx.0.25.mono.binned2.factor, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3")), paired = T) +
  theme_bw() + 
  ggtitle("with unpaired p-values")

ggplot(mono.summ.stats2, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) + 
  geom_line(aes(group = patient)) +
  geom_point() +
  facet_wrap(~cell.type.jx.0.25.mono.binned2.factor, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3")), paired = T) +
  theme_bw() + 
  ggtitle("with unpaired p-values")

dev.off()



mono.only@meta.data$cell.type.jx.0.25.mono.binned.factor = factor(mono.only@meta.data$cell.type.jx.0.25.mono.binned, levels = c("IFI+", "HLA/MARCO-hi", "MT-high/ZEB2-hi/Actin-hi", "Classical", "CD16+ Non-Classical", "Granulocyte"))
# Summarize statistics for monocytes
mono.summ.stats3 <- mono.only@meta.data %>%
  group_by(patient, timepoint, cell.type.jx.0.25.mono.binned.factor, .drop = FALSE) %>%
  summarize(n = length(patient), .drop = FALSE) %>%
  left_join(mono.only@meta.data %>%
              group_by(patient, timepoint, .drop = FALSE) %>%
              summarize(n.total = length(patient))) %>%
  mutate(prop = n / n.total)

# Create plots for monocyte cells
pdf("m.summ.stats.mono.3.pdf", height = 5, width = 8)

ggplot(mono.summ.stats3, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) + 
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~cell.type.jx.0.25.mono.binned.factor, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3")), paired = T) +
  theme_bw() + 
  ggtitle("with unpaired p-values")

ggplot(mono.summ.stats3, aes(x = timepoint, y = prop, fill = timepoint)) +
  geom_boxplot(outlier.size = -1) + 
  geom_line(aes(group = patient)) +
  geom_point() +
  facet_wrap(~cell.type.jx.0.25.mono.binned.factor, nrow = 1, scales = "free") +
  stat_compare_means(comparisons = list(c("Base", "D28"), c("D28", "M3"), c("Base", "M3")), paired = T) +
  theme_bw() + 
  ggtitle("with unpaired p-values")

dev.off()





# For Immature
Idents(mono.only) = mono.only$cell.type.jx.0.25.mono.binned2.factor
m3.vs.base.IFI.mono <- FindMarkers(mono.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "IFI+")
m3.vs.base.IFI.mono <- m3.vs.base.IFI.mono %>% as_tibble(rownames = "gene")
m3.vs.base.IFI.mono <- m3.vs.base.IFI.mono %>% filter(p_val_adj < 0.05)

pdf("IFI.mono.M3.vs.Base.pdf", height = 10, width = 10)
ggplot(m3.vs.base.IFI.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.IFI.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.IFI.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-25 & !gene %in% isg.consensus | abs(avg_log2FC) > 0.25 & p_val_adj < 1e-100 & !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.IFI.mono %>% filter(gene %in% isg.consensus, abs(avg_log2FC) > 0.25), color = "red") + 
  ggtitle("IFI mono Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.IFI.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()

VlnPlot(mono.only, 
        features = c("ISG15", "IFI27", "SIGLEC1", "FLT3", "ZBTB16", "IL1R2", "CXCL8", "DDIT4", "ADAMTS2", "EGR1", "S100A8", "S100A9", "FOX", "NAMPT", "ACSL1", "LGALS1", "HLA.genes", "MT.genes", "RB.genes"),
        flip = T, group.by = "cell.type.jx.0.25.mono.binned2", split.by = "timepoint", stack = T, fill.by = "ident")  + NoLegend()


# For Immature
m3.vs.base.classical.mono <- FindMarkers(mono.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "Classical")
m3.vs.base.classical.mono <- m3.vs.base.classical.mono %>% as_tibble(rownames = "gene")
m3.vs.base.classical.mono <- m3.vs.base.classical.mono %>% filter(p_val_adj < 0.05)


pdf("DEG.classical.mono.M3.vs.Base.pdf", height = 9, width = 9)

ggplot(m3.vs.base.classical.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.classical.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.classical.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-25 & !gene %in% isg.consensus | abs(avg_log2FC) > 0.25 & p_val_adj < 1e-100 & !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.classical.mono %>% filter(gene %in% isg.consensus, abs(avg_log2FC) > 0.25), color = "red") + 
  ggtitle("Classical mono Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.classical.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed

ggplot(m3.vs.base.classical.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.classical.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.classical.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-50 & !gene %in% isg.consensus | abs(avg_log2FC) > 0.5 & p_val_adj < 1e-100 & !gene %in% isg.consensus), max.overlaps = 50) +
  ggtitle("Classical mono Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()



# For Immature
m3.vs.base.cd16.nonclassical.mono <- FindMarkers(mono.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "CD16+ Non-Classical")
m3.vs.base.cd16.nonclassical.mono <- m3.vs.base.cd16.nonclassical.mono %>% as_tibble(rownames = "gene")
m3.vs.base.cd16.nonclassical.mono <- m3.vs.base.cd16.nonclassical.mono %>% filter(p_val_adj < 0.05)


pdf("DEG.nonclassical.mono.M3.vs.Base.pdf", height = 9, width = 9)

ggplot(m3.vs.base.cd16.nonclassical.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.cd16.nonclassical.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.cd16.nonclassical.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-25 & !gene %in% isg.consensus | abs(avg_log2FC) > 0.25 & p_val_adj < 1e-100 & !gene %in% isg.consensus), max.overlaps = 50) +
  geom_point(data = m3.vs.base.cd16.nonclassical.mono %>% filter(gene %in% isg.consensus, abs(avg_log2FC) > 0.25), color = "red") + 
  ggtitle("NonClassical mono Base vs M3") +
  geom_text_repel(aes(label = gene), data = m3.vs.base.cd16.nonclassical.mono %>% filter(gene %in% isg.consensus), color= "red") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed


ggplot(m3.vs.base.cd16.nonclassical.mono, aes(x = avg_log2FC, y = -log(p_val_adj, base = 10))) + 
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = m3.vs.base.cd16.nonclassical.mono %>% filter(p_val_adj < 0.05)) + 
  theme_bw() + 
  geom_text_repel(aes(label = gene), data = m3.vs.base.cd16.nonclassical.mono %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 1e-25 & !gene %in% isg.consensus | abs(avg_log2FC) > 0.25 & p_val_adj < 1e-50 & !gene %in% isg.consensus), max.overlaps = 50) +
  ggtitle("NonClassical mono Base vs M3") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1))  # Adjust the range as needed
dev.off()


VlnPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned2.factor", 
        features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")

VlnPlot(mono.only, group.by = "cell.type.jx.0.25.mono.binned.factor", 
        features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")


ifn.plot = cbind(mono.only@meta.data, t(mono.only@assays$AUC_IFN$counts[c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"),]))

ifn.plot = ifn.plot %>% pivot_longer(cols = banchereau.ifn.38:Landolt.5, names_to = "score.name", values_to = "score")


pdf("IFN.gene.sets.mono.only.pdf")
ifn.plot.sum = ifn.plot %>% group_by(cell.type.jx.0.25.mono.binned2.factor, timepoint, score.name) %>% summarize(mean.ifn = mean(score))
ggplot(ifn.plot.sum  %>% filter(timepoint != "D28"), aes(x = timepoint, y = mean.ifn, color = cell.type.jx.0.25.mono.binned2.factor, group = cell.type.jx.0.25.mono.binned2.factor)) + 
  geom_line() + geom_point() + facet_wrap(~score.name) + theme_bw()+ stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)


ifn.plot.sum2 = ifn.plot %>% group_by(cell.type.jx.0.25.mono.binned.factor, timepoint, score.name) %>% summarize(mean.ifn = mean(score))
ggplot(ifn.plot.sum2 %>% filter(timepoint != "D28"), aes(x = timepoint, y = mean.ifn, color = cell.type.jx.0.25.mono.binned.factor, group = cell.type.jx.0.25.mono.binned.factor)) + 
  geom_line() + geom_point() + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)
dev.off()


#plot for whole seurat
ifn.plot.myeloid = cbind(m.only@meta.data, t(m.only@assays$AUC_IFN$counts[c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"),]))

ifn.plot.myeloid = ifn.plot.myeloid %>% pivot_longer(cols = banchereau.ifn.38:Landolt.5, names_to = "score.name", values_to = "score")

pdf("IFN.gene.sets.myeloid.all.pdf")
ifn.plot.myeloid.sum = ifn.plot.myeloid %>% group_by(cell.type.jx.1.m.binned.factor, timepoint, score.name) %>% summarize(mean.ifn = mean(score), n.cells = length(score))
ggplot(ifn.plot.myeloid.sum %>% filter(timepoint != "D28"), aes(x = timepoint, y = mean.ifn, color = cell.type.jx.1.m.binned.factor, group = cell.type.jx.1.m.binned.factor)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)

ggplot(ifn.plot.myeloid.sum, aes(x = timepoint, y = mean.ifn, color = cell.type.jx.1.m.binned.factor, group = cell.type.jx.1.m.binned.factor)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)

ggplot(ifn.plot.myeloid.sum %>% filter(cell.type.jx.1.m.binned.factor %in% c("Monocyte", "Neut")),  
       aes(x = timepoint, y = mean.ifn, color = cell.type.jx.1.m.binned.factor, group = cell.type.jx.1.m.binned.factor)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)


dev.off()


pdf("neut.ifn.lines.pdf")
ifn.plot.neut = cbind(neut.only@meta.data, t(neut.only@assays$AUC_IFN$counts[c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"),]))

ifn.plot.neut = ifn.plot.neut %>% pivot_longer(cols = banchereau.ifn.38:Landolt.5, names_to = "score.name", values_to = "score")
ifn.plot.neut.sum = ifn.plot.neut %>% group_by(cell.type.jx.0.5.neut.binned.factor, timepoint, score.name) %>% summarize(mean.ifn = mean(score), n.cells = length(score))
ggplot(ifn.plot.neut.sum %>% filter(timepoint != "D28"), aes(x = timepoint, y = mean.ifn, color = cell.type.jx.0.5.neut.binned.factor, group = cell.type.jx.0.5.neut.binned.factor)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)


ifn.plot.neut.sum2 = ifn.plot.neut %>% group_by(cell.type.jx.0.5.neut.binned.factor, timepoint, score.name, patient, .drop = F) %>% summarize(mean.ifn = mean(score), n.cells = length(score))
ggplot(ifn.plot.neut.sum2 %>% filter(timepoint != "D28", cell.type.jx.0.5.neut.binned.factor == "IFI+"), 
       aes(x = timepoint, y = mean.ifn, color = patient, group = patient)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)

ggplot(ifn.plot.neut.sum2 %>% filter(timepoint != "D28", cell.type.jx.0.5.neut.binned.factor == "IFI+", patient != "S05"), 
       aes(x = timepoint, y = mean.ifn, color = patient, group = patient)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T)

dev.off()


pdf("IFN.gene.sets.neut.only.pdf")
ggplot(ifn.plot.neut.sum, aes(x = timepoint, y = mean.ifn, color = cell.type.jx.0.5.neut.binned.factor, group = cell.type.jx.0.5.neut.binned.factor)) + 
  geom_line() + geom_point(aes(size = n.cells)) + facet_wrap(~score.name) + theme_bw() + stat_compare_means(comparisons = list(c("Base", "M3")), paired = T, method = 't.test') + ggtitle('t.test')

VlnPlot(neut.only, group.by = "cell.type.jx.0.5.neut.binned.factor", 
        features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint")

VlnPlot(neut.only, group.by = "cell.type.jx.0.5.neut.binned.factor", 
        features = c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5"), 
        stack = T, flip = T, fill.by = "ident", split.by = "timepoint", idents = "IFI+")

gt::gt(m3.vs.base.neut.ifi.pathway)
gt::gt(m3.vs.base.neut.ifi.pathway %>% filter(pathway %in% c("banchereau.ifn.38", "M1.2.67", "Feng.5", "Yao.22", "Higgs.5", "Landolt.5")))


dev.off()

m3.vs.base.neut.ifi.pathway <- FindMarkers(neut.only, ident.1 = "Base", ident.2 = "M3", group.by = "timepoint", subset.ident = "IFI+", assay = "AUC_IFN")
m3.vs.base.neut.ifi.pathway <- m3.vs.base.neut.ifi.pathway %>% as_tibble(rownames = "pathway")

saveRDS(mono.only, "mono.only.final.anno.with.ifn.rds")
saveRDS(mono.only@meta.data, "mono.only.metadata.rds")
saveRDS(neut.only@meta.data, "neut.only.metadata.rds")

saveRDS(m.only@meta.data, "m.only.metadata.rds")
saveRDS(m.only, "m.only.167158.cells.rds")


save.image("Oct_19_myeloid.image.with.mono.all.rdata")

