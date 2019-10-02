## 
library(scater)
library(Seurat)
library(Matrix)
library(dplyr)
library(gplots)
library(RColorBrewer)
## -------------------------------------
## Function to read and create Seurat object
sample_process <- function(sample_name, sample_class, ID){
 # sample_class <- read.csv(sample_name, sep = ",", row.names = 1)
  sample_class <- read.delim(sample_name, sep=",", header=T)
  sample_class <- sample_class %>% distinct(X,.keep_all = T)

  rownames(sample_class) <- sample_class$X
  sample_class <- apply(sample_class ,2,function(x) replace(x,is.na(x),0))
  #head(sample_class[1:5, 1:5])

  sample_class <- Matrix(t(sample_class))
  #sample_class_obj <- paste0("class_name", "obj")
  sample_class <- CreateSeuratObject(raw.data=sample_class,
                                     min.cells=)
  #data.to.add_obj <- paste0(sample_class, "obj")
  data.to.add<- rep(ID,length=length(sample_class@cell.names))
  names(data.to.add) <- sample_class@cell.names

  sample_class <- AddMetaData(sample_class, metadata=data.to.add, col.name="samples")

}



## ***************************************************
## BC_1

BC_1_tum_1_data <- sample_process("GSM3148591_BC01_TUMOR1_counts.csv",BC_1_tum_1, "BC_1_tum_1" )
BC_1_tum_1_data
# 13386 genes across 1587 samples

BC_1_tum_2_data <- sample_process("GSM3148592_BC01_TUMOR2_counts.csv",BC_1_tum_2, "BC_1_tum_2" )
BC_1_tum_2_data
## 12852 genes across 1521 samples

BC_1_tum_3_data <- sample_process("GSM3148593_BC01_TUMOR3_counts.csv",BC_1_tum_3, "BC_1_tum_3" )
BC_1_tum_3_data
# 13809 genes across 1542 samples

BC_1_tum_4_data <- sample_process("GSM3148594_BC01_TUMOR4_counts.csv",BC_1_tum_4, "BC_1_tum_4" )
BC_1_tum_4_data
# 15126 genes across 2244 samples

BC_1_tum_1_2 <- MergeSeurat(object1 = BC_1_tum_1_data,
                            object2 = BC_1_tum_2_data,
                            add.cell.id1 = "BC_1_tum_1",
                            add.cell.id2 = "BC_1_tum_2",
                            project = "BC_1_tum_1_2")

BC_1_tum_1_2

BC_1_tum_12_3<- MergeSeurat(object1 = BC_1_tum_1_2,
                            object2 = BC_1_tum_3_data,
                            #add.cell.id1 = "BC_1_tum_1",
                            add.cell.id2 = "BC_1_tum_3",
                            project = "BC_1_tum_12_3")

head(BC_1_tum_12_3@meta.data)
tail(BC_1_tum_12_3@meta.data)


BC_1_tum_123_4<- MergeSeurat(object1 = BC_1_tum_12_3,
                             object2 = BC_1_tum_4_data,
                             min.cells = 3,
                             min.genes = 200,
                             #add.cell.id1 = "BC_1_tum_1",
                             add.cell.id2 = "BC_1_tum_4",
                             project = "BC_1_tum_123_4")

BC_1_tum_123_4

# 18036 genes across 6894 samples
# After min cell n genes filtering
# 12660 genes across 1868 samples ******

tail(BC_1_tum_123_4@meta.data)
saveRDS(BC_1_tum_123_4, file="BC_1_tum_123_4.rds")

# summary of total expression per single cell
summary(colSums(BC_1_tum_1_data@raw.data))

# check how many genes have at least one transcript in each cell
at_least_one <- apply(pbmc.data, 2, function(x) sum(x>0))

hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(pbmc.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

hist(colSums(BC_1_tum_1_data@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")
## *************************
## BC_2

BC_2_tum_1_data <- sample_process("GSM3148604_BC02_TUMOR1_counts.csv",BC_2_tum_1, "BC_2_tum_1" )
BC_2_tum_1_data
#16263 genes across 1359 samples

BC_2_tum_2_data <- sample_process("GSM3148605_BC02_TUMOR2_counts.csv",BC_2_tum_2, "BC_2_tum_2" )
BC_2_tum_2_data
# 16105 genes across 1406 samples

BC_2_tum_3_data <- sample_process("GSM3148606_BC02_TUMOR3_counts.csv",BC_2_tum_3, "BC_2_tum_3" )
BC_2_tum_3_data
# 17164 genes across 1950 samples

BC_2_tum_4_data <- sample_process("GSM3148607_BC02_TUMOR4_counts.csv",BC_2_tum_4, "BC_2_tum_4" )
BC_2_tum_4_data
# 12275 genes across 708 samples

BC_2_tum_1_2 <- MergeSeurat(object1 = BC_2_tum_1_data,
                            object2 = BC_2_tum_2_data,
                            add.cell.id1 = "BC_2_tum_1",
                            add.cell.id2 = "BC_2_tum_2",
                            project = "BC_2_tum_1_2")

BC_2_tum_1_2
# 18219 genes across 2765 samples
head(BC_2_tum_1_2@meta.data)
tail(BC_2_tum_1_2@meta.data)


BC_2_tum_12_3 <- MergeSeurat(object1 = BC_2_tum_1_2,
                            object2 = BC_2_tum_3_data,
                            #add.cell.id1 = "BC_2_tum_1",
                            add.cell.id2 = "BC_2_tum_3",
                            project = "BC_2_tum_12_3")
BC_2_tum_12_3
# 19947 genes across 4715 samples
tail(BC_2_tum_12_3@meta.data)


BC_2_tum_123_4 <- MergeSeurat(object1 = BC_2_tum_12_3,
                             object2 = BC_2_tum_4_data,
                             min.cells = 3,
                             min.genes = 200,
                             #add.cell.id1 = "BC_2_tum_1",
                             add.cell.id2 = "BC_2_tum_4",
                             project = "BC_2_tum_123_4")

BC_2_tum_123_4
# 20208 genes across 5423 samples
# After min cell n genes filtering
# 15241 genes across 2439 samples  *****

tail(BC_2_tum_123_4@meta.data)
saveRDS(BC_2_tum_123_4, file="BC_2_tum_123_4.rds")

## **************************
# BC_4

BC_4_tum_1_data <- sample_process("GSM3148621_BC04_TUMOR1_counts.csv",BC_4_tum_1, "BC_4_tum_1" )
BC_4_tum_1_data
# 21067 genes across 1643 samples

BC_4_tum_2_data <- sample_process("GSM3148622_BC04_TUMOR2_counts.csv",BC_4_tum_2, "BC_4_tum_2" )
BC_4_tum_2_data
# 24580 genes across 2481 samples

BC_4_tum_3_data <- sample_process("GSM3148623_BC04_TUMOR3_counts.csv",BC_4_tum_3, "BC_4_tum_3" )
BC_4_tum_3_data
# 23281 genes across 2429 samples

BC_4_tum_5_data <- sample_process("GSM3148625_BC04_TUMOR5_counts.csv",BC_4_tum_5, "BC_4_tum_5" )
BC_4_tum_5_data
# 22796 genes across 2426 samples

BC_4_tum_6_data <- sample_process("GSM3148626_BC04_TUMOR6_counts.csv",BC_4_tum_6, "BC_4_tum_6" )
BC_4_tum_6_data
# 22826 genes across 2052 samples

BC_4_tum_1_2 <- MergeSeurat(object1 = BC_4_tum_1_data,
                            object2 = BC_4_tum_2_data,
                            add.cell.id1 = "BC_4_tum_1",
                            add.cell.id2 = "BC_4_tum_2",
                            project = "BC_4_tum_1_2")

BC_4_tum_1_2
# 25939 genes across 4124 samples

head(BC_4_tum_1_2@meta.data)
tail(BC_4_tum_1_2@meta.data)



BC_4_tum_12_3 <- MergeSeurat(object1 = BC_4_tum_1_2,
                            object2 = BC_4_tum_3_data,
                            #add.cell.id1 = "BC_4_tum_1",
                            add.cell.id2 = "BC_4_tum_3",
                            project = "BC_4_tum_12_3")
BC_4_tum_12_3
# 27583 genes across 6553 samples

BC_4_tum_123_5 <- MergeSeurat(object1 = BC_4_tum_12_3,
                             object2 = BC_4_tum_5_data,
                             #add.cell.id1 = "BC_4_tum_1",
                             add.cell.id2 = "BC_4_tum_5",
                             project = "BC_4_tum_123_5")

BC_4_tum_123_5
# 28670 genes across 8979 samples
tail(BC_4_tum_123_5@meta.data)


BC_4_tum_1235_6 <- MergeSeurat(object1 = BC_4_tum_123_5,
                              object2 = BC_4_tum_6_data,
                              min.cells = 3,
                              min.genes = 200,
                              #add.cell.id1 = "BC_4_tum_1",
                              add.cell.id2 = "BC_4_tum_6",
                              project = "BC_4_tum_1235_6")

BC_4_tum_1235_6
# 29553 genes across 11031 samples
# After min cell n genes filtering
# 23730 genes across 8322 samples

tail(BC_4_tum_1235_6@meta.data)

saveRDS(BC_4_tum_1235_6, file="BC_4_tum_1235_6.rds")

##
## -----------------------------------------
## BC_6

BC_6_tum_1_data <- sample_process("GSM3148631_BC06_TUMOR1_counts.csv",BC_6_tum_1, "BC_6_tum_1")
head(BC_6_tum_1_data@meta.data)
#12502 genes across 479 samples


BC_6_tum_2_data <- sample_process("GSM3148632_BC06_TUMOR2_counts.csv", BC_6_tum_2, "BC_6_tum_2")
head(BC_6_tum_2_data@meta.data)
#12737 genes across 815 samples

BC_6_tum_3_data <- sample_process("GSM3148633_BC06_TUMOR3_counts.csv", BC_6_tum_3, "BC_6_tum_3")
head(BC_6_tum_3_data@meta.data)
BC_6_tum_3_data
# 20749 genes across 3441 samples


BC_6_tum_1_2 <- MergeSeurat(object1 = BC_6_tum_1_data,
                            object2 = BC_6_tum_2_data,
                            add.cell.id1 = "BC_6_tum_1",
                            add.cell.id2 = "BC_6_tum_2",
                            project = "BC_6_tum_1_2")

BC_6_tum_1_2
#14473 genes across 1294 samples


BC_6_tum_12_3 <- MergeSeurat(object1 = BC_6_tum_1_2,
                             object2 = BC_6_tum_3_data,
                             min.cells = 3,
                             min.genes = 200,
                             #add.cell.id1 = "BC_6_tum_1",
                             add.cell.id2 = "BC_6_tum_3",
                             project = "BC_6_tum_12_3")

head(BC_6_tum_12_3@meta.data)
tail(BC_6_tum_12_3@meta.data)
BC_6_tum_12_3
#21318 genes across 4735 samples
# After min cell n genes filtering
#16417 genes across 2876 samples *******
saveRDS(BC_6_tum_12_3, file="BC_6_tum_12_3.rds")

## *******************************************

BC_9_tum_1 <- Read10X(data.dir = "BC09_tum_1/filtered_gene_bc_matrices/hg19")

# summary of total expression per single cell
summary(colSums(BC_9_tum_1))

# check how many genes have at least one transcript in each cell
at_least_one <- apply(BC_9_tum_1, 2, function(x) sum(x>0))

hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(BC_9_tum_1),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")


BC_9_tum_1_data <- CreateSeuratObject(raw.data = BC_9_tum_1,
                          min.cells = 3,
                          min.genes = 200,
                          project = "BC_9_tum_1",
                          names.delim = "-",
                          names.field = 2)

BC_9_tum_1_data
# 15883 genes across 7096 samples\
BC_9_tum_1_samples <- rep("BC_9_tum_1",length(BC_9_tum_1_data@cell.names))
BC_9_tum_1_samples

BC_9_tum_1_data@meta.data$samples <- BC_9_tum_1_samples

mito.genes_BC_9_tum_1 <- grep(pattern = "^MT", x = rownames(x = BC_9_tum_1_data@data), value = TRUE)
head(mito.genes_BC_9_tum_1)

#length(BC_9_tum_1_data@meta.data$nGene)

percent.mito_BC_9_tum_1 <- Matrix::colSums(BC_9_tum_1_data@raw.data[mito.genes_BC_9_tum_1, ])/Matrix::colSums(BC_9_tum_1_data@raw.data)

BC_9_tum_1_data <- AddMetaData(object = BC_9_tum_1_data, metadata = percent.mito_BC_9_tum_1, col.name = "percent.mito")

VlnPlot(object = BC_9_tum_1_data,
        features.plot = c("nGene", "nUMI", "percent.mito"), cols.use = palette(brewer.pal(n = 8, name = "Set2")), nCol = 3)

ggplot(BC_9_tum_1_data@meta.data, aes(nUMI, percent.mito)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.09, linetype = "dashed", colour = "red")


head(BC_9_tum_1_data@meta.data)
length(BC_9_tum_1_data@meta.data)

# BC_9_tum_1_data_100 <- FilterCells(object = BC_9_tum_1_data,
#                                subset.names = c("nGene", "percent.mito"),
#                                low.thresholds = c(100, -Inf),
#                                high.thresholds = c(2500, 0.05))
#
# length(BC_9_tum_1_data@meta.data$nGene)
# length(BC_1_tum_1_data@meta.data$nGene)
#
# head(BC_1_tum_1_data@meta.data)
# head(BC_9_tum_1_data@meta.data)
# manual check; I already know all cells have >200 genes
table(BC_9_tum_1_data@meta.data$percent.mito < 0.05 & BC_9_tum_1_data@meta.data$nGene<2500)


BC_9_tum_2 <- Read10X(data.dir = "BC09_tum_2/filtered_gene_bc_matrices/hg19")

BC_9_tum_2_data <- CreateSeuratObject(raw.data = BC_9_tum_2,
                                      min.cells = 3,
                                      min.genes = 200,
                                      project = "BC_9_tum_2",
                                      names.delim = "-",
                                      names.field = 2)

BC_9_tum_2_data
# 15863 genes across 6889 samples
mito.genes <- grep(pattern = "^MT", x = rownames(x = mydata@data), value = TRUE)

BC_9_tum_2_samples <- rep("BC_9_tum_2",length(BC_9_tum_2_data@cell.names))
BC_9_tum_2_samples

BC_9_tum_2_data@meta.data$samples <- BC_9_tum_2_samples

mito.genes_BC_9_tum_2 <- grep(pattern = "^MT", x = rownames(x = BC_9_tum_2_data@data), value = TRUE)
head(mito.genes_BC_9_tum_2)

#length(BC_9_tum_1_data@meta.data$nGene)

percent.mito_BC_9_tum_2 <- Matrix::colSums(BC_9_tum_2_data@raw.data[mito.genes_BC_9_tum_2, ])/Matrix::colSums(BC_9_tum_2_data@raw.data)

BC_9_tum_2_data <- AddMetaData(object = BC_9_tum_2_data, metadata = percent.mito_BC_9_tum_2, col.name = "percent.mito")

head(BC_9_tum_2_data@meta.data)

BC_9_1_2 <- MergeSeurat(object1 = BC_9_tum_1_data,
                         object2 = BC_9_tum_2_data,
                         min.cells = 3,
                         min.genes = 200,
                         add.cell.id1 = "BC_9_tum_1",
                         add.cell.id2 = "BC_9_tum_2",
                         project = "BC_9_tum_1_2")

BC_9_1_2
# 16506 genes across 13985 samples

head(BC_9_1_2@meta.data)
tail(BC_9_1_2@meta.data)

saveRDS(BC_9_1_2, file="BC_9_1_2.rds")
tail(BC_9_1_2@meta.data)

# ____----------------
BC_1_tum_123_4
BC_2_tum_123_4
BC_4_tum_1235_6
BC_6_tum_12_3
BC_9_1_2

BC_1_BC_2 <- MergeSeurat(object1 = BC_1_tum_123_4,
                        object2 = BC_2_tum_123_4,
                        min.cells = 3,
                        min.genes = 200,
                        #add.cell.id1 = "BC_9_tum_1",
                        #add.cell.id2 = "BC_9_tum_2",
                        project = "BC_1_BC_2")

head(BC_1_BC_2@meta.data)
tail(BC_1_BC_2@meta.data)
#15730 genes across 4285 samples

BC_1_2_BC_4 <- MergeSeurat(object1 = BC_1_BC_2,
                           object2 = BC_4_tum_1235_6,
                           min.cells = 3,
                           min.genes = 200,
                           #add.cell.id1 = "BC_9_tum_1",
                           #add.cell.id2 = "BC_9_tum_2",
                           project = "BC_1_2_BC_4")

BC_1_2_BC_4
#23879 genes across 12606 samples

BC_12_4_BC_6<- MergeSeurat(object1 = BC_1_2_BC_4,
                           object2 = BC_6_tum_12_3,
                           min.cells = 3,
                           min.genes = 200,
                           #add.cell.id1 = "BC_9_tum_1",
                           #add.cell.id2 = "BC_9_tum_2",
                           project = "BC_12_4_BC_6")

BC_12_4_BC_6
# 24321 genes across 15482 samples

BC_124_6_BC_9<- MergeSeurat(object1 = BC_12_4_BC_6,
                            object2 = BC_9_1_2,
                            min.cells = 3,
                            min.genes = 200,
                            #add.cell.id1 = "BC_9_tum_1",
                            #add.cell.id2 = "BC_9_tum_2",
                            project = "BC_124_6_BC_9")
BC_124_6_BC_9
#27911 genes across 29467 samples

saveRDS(BC_124_6_BC_9, file="BC_124_6_BC_9.rds")
head(BC_124_6_BC_9@meta.data)
tail(BC_124_6_BC_9@meta.data)
samplename_all <- BC_124_6_BC_9@meta.data$samples
table(samplename_all)

# BC_1_tum_1 BC_1_tum_2 BC_1_tum_3 BC_1_tum_4 BC_2_tum_1 BC_2_tum_2 BC_2_tum_3 BC_2_tum_4 BC_4_tum_1 BC_4_tum_2
# 326        172        278       1072        663        616        912        246       1165       2237
# BC_4_tum_3 BC_4_tum_5 BC_4_tum_6 BC_6_tum_1 BC_6_tum_2 BC_6_tum_3 BC_9_tum_1 BC_9_tum_2
# 1496       1830       1593        135        194       2547       7096       6889
