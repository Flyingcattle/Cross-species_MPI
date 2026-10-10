

library(Seurat)
library(ggplot2)
library(dplyr)
library(MetaNeighbor)
library(scRNAseq)
library(scater)
library(patchwork)

#### load datasets
pig<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\MetaNeighbor\\All_cells\\Splited_seu_anno\\pig_seurat.rds")
mouse<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\MetaNeighbor\\All_cells\\Splited_seu_anno\\mouse_seurat.rds")
monkey<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\MetaNeighbor\\All_cells\\Splited_seu_anno\\monkey_seurat.rds")
human<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\MetaNeighbor\\All_cells\\Splited_seu_anno\\human_seurat.rds")
pig  #22141 features across 18956 samples
mouse #21608 features across 29984 samples
monkey #18414 features across 9700 samples
human #28820 features across 49043 samples

## downsample 
pig.subsampled <- pig[, sample(colnames(pig), size =1000, replace=F)]
mouse.subsampled <- mouse[, sample(colnames(mouse), size =1000, replace=F)]
monkey.subsampled <- monkey[, sample(colnames(monkey), size =1000, replace=F)]
human.subsampled <- human[, sample(colnames(human), size =1000, replace=F)]

rm(pig, mouse, monkey, human)
# convert seurat2sce
#human.sce <- as.SingleCellExperiment(human)
#monkey.sce <- as.SingleCellExperiment(monkey)
#mouse.sce <- as.SingleCellExperiment(mouse)
#pig.sce <- as.SingleCellExperiment(pig,Assays )

############## pig  #############
mycount<-as.data.frame(pig.subsampled@assays$RNA[,])
dim(mycount) #29325  9700
mygroup <- as.data.frame(pig.subsampled@meta.data)
dim(mygroup)
# create a SingleCellExperiment object
pig.sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(mycount)),
  colData = mygroup
)
############## pig  #############

############## mouse  #############
mycount<-as.data.frame(mouse.subsampled@assays$RNA[,])
dim(mycount) #29325  9700
mygroup <- as.data.frame(mouse.subsampled@meta.data)
dim(mygroup)
# create a SingleCellExperiment object
mouse.sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(mycount)),
  colData = mygroup
)
############## mouse  #############

############## monkey  #############
mycount<-as.data.frame(monkey.subsampled@assays$RNA[,])
dim(mycount) #29325  9700
mygroup <- as.data.frame(monkey.subsampled@meta.data)
dim(mygroup)
# create a SingleCellExperiment object
monkey.sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(mycount)),
  colData = mygroup
)
############## monkey  #############


############## human  #############
mycount<-as.data.frame(human.subsampled@assays$RNA[,])
dim(mycount) #29325  9700
mygroup <- as.data.frame(human.subsampled@meta.data)
dim(mygroup)
# create a SingleCellExperiment object
human.sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(mycount)),
  colData = mygroup
)
############## human  #############
human.sce
monkey.sce
mouse.sce
pig.sce

my_data<-list(human.sce,monkey.sce,mouse.sce,pig.sce)
my_data

lapply(my_data, function(x) head(rownames(x), 3))
lapply(my_data, function(x) colnames(colData(x)))
lapply(my_data, function(x) names(assays(x)))

library(MetaNeighbor)
fused_data = mergeSCE(my_data)
dim(fused_data)  #29325  4000
head(colData(fused_data))

saveRDS(fused_data,file="E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\MetaNeighbor\\merged_pancreas.rds")

pancreas_data<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\MetaNeighbor\\merged_pancreas.rds")

table(pancreas_data$finalcelltype, pancreas_data$species)


global_hvgs = variableGenes(dat = pancreas_data,
                            exp_labels = pancreas_data$species)
length(global_hvgs)  #2598

aurocs = MetaNeighborUS(var_genes = global_hvgs[1:500],
                        dat = pancreas_data,
                        study_id = pancreas_data$species,
                        cell_type = pancreas_data$finalcelltype,
                        fast_version = TRUE)
##下面这里有个报错，需要，原因是plotHeatmap函数重复了
detach("package:scater", unload = TRUE)
plotHeatmap(aurocs, cex = 0.5)  

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"magma"))(100))
plotHeatmap(aurocs, cex = 0.5,col = cols)  

cols <- rev(gray.colors(100, start = 0, end = 1)) 



topHits(aurocs, dat = pancreas_data, study_id = pancreas_data$species,
        cell_type = pancreas_data$finalcelltype, threshold = 0.9)

level1_split = splitClusters(aurocs, k = 2)
level1_split
first_split = level1_split[[2]]


full_labels = makeClusterName(pancreas_data$species,
                              pancreas_data$finalcelltype)
subdata = pancreas_data[, full_labels %in% first_split]
dim(subdata)

var_genes = variableGenes(dat = subdata, exp_labels = subdata
                          $species)
aurocs = MetaNeighborUS(var_genes = var_genes,
                        dat = subdata, fast_version = TRUE,
                        study_id = subdata$species,
                        cell_type = subdata$finalcelltype)
plotHeatmap(aurocs, cex = 0.7)


best_hits = MetaNeighborUS(var_genes = global_hvgs,
                           dat = pancreas_data,
                           study_id = pancreas_data$species,
                           cell_type = pancreas_data$finalcelltype,
                           fast_version = TRUE,
                           one_vs_best = TRUE, symmetric_output = FALSE)
plotHeatmap(best_hits, cex = 0.5)
