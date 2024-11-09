library(Seurat)
library(Matrix)
library(dplyr)
library(gplots)
library(RColorBrewer)

BC_tum_12469 <- readRDS("BC_124_6_BC_9.rds")

BC_tum_12469 <- NormalizeData(object = BC_tum_12469, 
                 normalization.method = "LogNormalize", 
                 scale.factor = 1e4)
# 
# BC_tum_12469 <- ScaleData(object = BC_tum_12469,
#                        vars.to.regress = c("nUMI", "percent.mito"))

## **** Since not all the samples have percent.mito information 
## I had to run ScaleData as below  
BC_tum_12469 <- ScaleData(object = BC_tum_12469,
                          vars.to.regress = c("nUMI", "nGene"))


BC_tum_12469 <- FindVariableGenes(object = BC_tum_12469,
                               mean.function = ExpMean,
                               dispersion.function = LogVMR,
                               x.low.cutoff = 0.0125,
                               x.high.cutoff = 3,
                               y.cutoff = 0.8)

# BC_tum_12469 <- ScaleData(object = BC_tum_12469,
#                        vars.to.regress = c("nUMI", "percent.mito"))


BC_tum_12469 <- RunPCA(object = BC_tum_12469,
                    pc.genes = BC_tum_12469@var.genes,
                    do.print = TRUE,
                    pcs.print = 1:5,
                    genes.print = 5)
 
PCAPlot(object = BC_tum_12469, dim.1 = 1, dim.2 = 2)



BC_tum_12469 <- RunTSNE(object = BC_tum_12469,
                     dims.use = 1:10,
                     do.fast = TRUE)


TSNEPlot(object = BC_tum_12469, do.label = TRUE)



BC_tum_12469 <- FindClusters(object = BC_tum_12469, reduction.type = "pca", dims.use = 1:10, 
                          resolution = c(0.3, 0.5, 0.7, 1), print.output = 0, save.SNN = TRUE)

BC_tum_12469 <- SetAllIdent(BC_tum_12469, id="res.0.5")
levels(BC_tum_12469@ident)

pdf("tsne_lot.pdf")

TSNEPlot(object = BC_tum_12469, do.label = TRUE)

dev.off()

# First lets stash our identities for later
BC_tum_12469 <- StashIdent(object = BC_tum_12469, save.name = "ClusterNames_0.5")

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
BC_tum_12469.markers <- FindAllMarkers(object = BC_tum_12469, only.pos = TRUE, min.pct = 0.25, 
                                    thresh.use = 0.25)

write.csv(BC_tum_12469.markers, file="BC_tum_12469.markers.csv")


BC_tum_12469.markers_20 <-BC_tum_12469.markers  %>% group_by(cluster) %>% top_n(20, avg_logFC)

write.csv(BC_tum_12469.markers_20, file="BC_tum_12469_cluster_genes_20.csv")

# Visualize a heatmap of markers
# DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting
# the top 20 markers (or all markers if less than 20) for each cluster.

BC_tum_12469.markers_top10 <- BC_tum_12469.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

pdf("Heatmap_BC_tum_12469.pdf")

DoHeatmap(object = BC_tum_12469, genes.use = BC_tum_12469.markers_top10$gene, slim.col.label = TRUE, remove.key = TRUE)


dev.off()

saveRDS(BC_tum_12469, file = "BC_tum_12469.rds")