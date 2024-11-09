 ## Author: Subinoy Biswas
library(Seurat)         ## Version 3
library(Matrix)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(sctransform)
library(Cairo)


set.seed(123)                                                               # Setting seed value for reproducibility
BC_tum_1246 <- readRDS("Aug_7_integration/BC_12_4_BC_6_aug_190807.rds")
BC_tum_1246_new <-UpdateSeuratObject(BC_tum_1246)

BC_tum_1246_new$tumor_2 <- rownames(BC_tum_1246_new@meta.data)


rm(BC_tum_1246)
BC_tum_1246 <- BC_tum_1246_new


# store mitochondrial percentage in object meta data
BC_tum_1246 <- PercentageFeatureSet(BC_tum_1246, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BC_tum_1246 <- SCTransform(BC_tum_1246, vars.to.regress = "percent.mt", verbose = TRUE)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BC_tum_1246), 10)

# These are now standard steps in the Seurat workflow for visualization and clustering
BC_tum_1246 <- RunPCA(BC_tum_1246, verbose = FALSE)


BC_tum_1246 <- RunUMAP(BC_tum_1246, dims = 1:30, verbose = FALSE)



tum_meta <- as.data.frame(BC_tum_1246@meta.data)
head(tum_meta)
tail(tum_meta)
tum_meta$tumor <- rownames(tum_meta)

tum_meta$tumor <- gsub(x=tum_meta$tumor, pattern= "(BC\\d).*", replacement = "\\1", perl=TRUE)

BC_tum_1246@meta.data <- tum_meta


head(BC_tum_1246@meta.data)
head(BC_tum_1246@meta.data$samples)
tail(rownames(BC_tum_1246@meta.data))
tail(BC_tum_1246@meta.data)

DimPlot(BC_tum_1246, label = TRUE)    # NoLegend

BC_tum_1246 <- FindNeighbors(BC_tum_1246, dims = 1:10, verbose = FALSE)

BC_tum_1246 <- FindClusters(BC_tum_1246, resolution = 0.5)


BC_tum_1246 <- RunTSNE(object = BC_tum_1246)

DimPlot(object = BC_tum_1246, reduction = "tsne",label = TRUE)

DimPlot(object = BC_tum_1246, reduction = "umap",label = TRUE)

BC_tum.markers_wt <- FindAllMarkers(BC_tum_1246)

DimPlot(object = BC_tum_1246, reduction = "tsne",label = TRUE)

write.csv(BC_tum.markers_wt, file="tum_6_Aug_9/BC_tum.markers_clust_17_aug_9.csv")

# Top 10 markers
top10_markers <- BC_tum.markers_wt %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(top10_markers , file="tum_6_Aug_9/top10_markers_clust_17.csv")


top20_markers <- BC_tum.markers_wt %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(top20_markers, file="tum_6_Aug_9/top20_markers_clust_17.csv")

saveRDS(BC_tum_1246, file="tum_6_Aug_9/BC_tum_1246_17.rds")
saveRDS(BC_tum.markers_wt, file="tum_6_Aug_9/BC_tum.markers.rds")
saveRDS(top10_markers, file="tum_6_Aug_9/top10_markers.rds")
## =============================================>>

DoHeatmap(BC_tum_1246, features = top10_markers$gene)


jpeg("tum_6_Aug_9/Cluster_b_cells.jpg", height = 8, width =12, units = 'in', res= 300)
FeaturePlot(BC_tum_1246, features = c("IGKC", "IGLC2", "IGLC3", "IGHG1",
                                       "CD79A", "IGHG3", "MS4A1", "IGHA1",
                                       "IGHG2", "IGHM" ), pt.size = 0.2, reduction = "tsne", ncol = 3)

dev.off()

# 2
jpeg("tum_6_Aug_9/Cluster_cd4.tconv.jpg", height = 8, width =12, units = 'in', res= 300 )
FeaturePlot(BC_tum_1246, features = c("IL7R", "LTB", "LDHB", "NOSIP",
                                       "MAL", "CCR7", "LEF1", "RPS6",
                                       "TCF7", "RPS12" ), pt.size = 0.2, reduction = "tsne", ncol = 3)

dev.off()

# 3
jpeg("tum_6_Aug_9/Cluster_cd4.treg.jpg", height = 8, width =12, units = 'in', res= 300)
FeaturePlot(BC_tum_1246, features = c("TNFRSF4", "TNFRSF18", "IL2RA", "CTLA4",
                                       "AC133644.2", "BATF", "IL32", "ICOS",
                                       "ISG15", "IFI6" ), pt.size = 0.2, reduction = "tsne", ncol = 3)

dev.off()
# 4
jpeg("tum_6_Aug_9/Cluster_cd8.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("CCL5", "GZMA", "GZMK", "CCL4",
                                       "CXCL13", "CD8A", "GZMB", "CCL4L2",
                                       "GZMH", "AC092580.4" ), pt.size = 0.2,reduction = "tsne", ncol = 3)

dev.off()

# 5
jpeg("tum_6_Aug_9/Cluster_cd14.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("S100A9", "S100A8", "LYZ", "S100A12",
                                       "RP11-1143G9.4", "CXCL8", "VCAN", "SPP1",
                                       "IL1B", "CD14" ), pt.size = 0.2,reduction = "tsne", ncol = 3)

dev.off()

# 6
jpeg("tum_6_Aug_9/Cluster_cd16.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("LST1", "AIF1", "FCGR3A", "IFITM3",
                                       "CFD", "SERPINA1", "CST3", "FCER1G",
                                       "MS4A7", "LINC01272" ), pt.size = 0.2,reduction = "tsne", ncol = 3)

dev.off()


# 7
jpeg("tum_6_Aug_9/markers_dc.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("CST3","HLA-DPB1","HLA-DPA1","HLA-DRA","HLA-DQA1",
                                      "HLA-DRB1","HLA-DQB1","HLA-DRB5","TXN","CD74"),
                                        pt.size = 0.2,reduction = "tsne", ncol = 3)

dev.off()


# 8
jpeg("tum_6_Aug_9/Cluster_mast.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("TPSAB1",
                                       "TPSB2",
                                       "CPA3",
                                       "CLU",
                                       "HPGDS",
                                       "LTC4S",
                                       "MS4A2",
                                       "CD9",
                                       "HPGD",
                                       "AREG"
                                    ), pt.size = 0.2, reduction = "tsne", ncol = 3)

dev.off()


# 9
jpeg("tum_6_Aug_9/Cluster_nk.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("FGFBP2",
                                       "GNLY",
                                       "NKG7",
                                       "SPON2",
                                       "KLRF1",
                                       "PRF1",
                                       "PTGDS",
                                       "KLRB1",
                                       "CTSW",
                                       "CLIC3"
                                    ), pt.size = 0.2, reduction = "tsne", ncol = 3)

dev.off()


# 10
jpeg("tum_6_Aug_9/markers_pdc.cells.jpg", height = 8, width =12, units = 'in', res= 300)

FeaturePlot(BC_tum_1246, features = c("PTGDS",
                                       "JCHAIN",
                                       "IRF7",
                                       "LILRA4",
                                       "PLD4",
                                       "GZMB",
                                       "IRF8",
                                       "ITM2C",
                                       "PPP1R14B",
                                       "SERPINF1"
                                    ), pt.size = 0.2, reduction = "tsne", ncol = 3)

dev.off()
