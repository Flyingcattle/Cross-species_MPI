# 加载Seurat包
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readxl)

###########################################################################################################   
#################################### load PLI score database ##############################

setwd("E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure3/Human_seprmatogenesis_test/pLI_score")
#BiocManager::install("AnnotationDbi")
# Loading pLI percentile
Human_pLI <- readRDS(file="Human.pLI_genesymbol.rds")
head(Human_pLI)
dim(Human_pLI)    # 14798     2

# 1. 删除第二列（保留第1和第3列）
Human_pLI <- Human_pLI[, -2, drop = FALSE]  # drop=FALSE确保结果仍是数据框
# 2. 将第一列列名从GeneSymbol改为Human_ID
colnames(Human_pLI)[1] <- "Human_ID"

#################################### load PLI score database ##############################
###########################################################################################################   


###########################################################################################################   
#################################### load Seurat object ##############################

combined<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\Figure3\\myTAI\\Seurat_sex_combined.rds")
cell_total<-combined
# 定义细胞类型顺序
celltype_order <- c("Mito", "Pre_L", "Lepto", 
                    "Zygot", "Pachy", "Diplo",
                    "Dicty","Post_meiotic")
table(cell_total@meta.data$sample_id)
table(cell_total@meta.data$sex)

DimPlot(combined, reduction = "umap", label=T,group.by = "celltype_myTAI")

#################################### load Seurat object ##############################
###########################################################################################################   

seu<-cell_total
Human_counts <- seu@assays$RNA@counts
Human_counts[1:5, 1:5]

Human_meta_data <- seu@meta.data

meta <- data.frame(
  CELL = rownames(Human_meta_data), 
  CELL_TYPE = Human_meta_data$celltype_myTAI,  
  check.names = FALSE,  
  stringsAsFactors = FALSE 
)
head(meta)
Human_meta_data<-meta

# Gathering and joining
Human_counts %>% 
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(Human_ID = rowname) %>% 
  inner_join( Human_pLI ) %>% 
  gather(CELL,count,-pLI_perc,-Human_ID) %>% 
  inner_join(Human_meta_data) %>%
  filter(count != 0) %>%
  group_by(CELL,CELL_TYPE) %>% 
  summarise( MEDIAN = median(pLI_perc),
             MEAN = mean(pLI_perc) )  -> Human_gathered

# Ordering cell types
Human_gathered$CELL_TYPE <- 
  factor(Human_gathered$CELL_TYPE, 
         levels = c("Mito", "Pre_L", "Lepto", 
                    "Zygot", "Pachy", "Diplo",
                    "Dicty","Post_meiotic") )

# 绘图（使用因子水平作为标签，确保顺序一致）
Human_gathered %>%
  ggplot(aes(CELL_TYPE, MEAN, fill = CELL_TYPE)) +
  geom_boxplot(notch = T) +
  theme_cowplot()  +
  scale_x_discrete(labels = levels(Human_gathered$CELL_TYPE)) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) +
  ylab("pLI percentile") + scale_y_continuous(limits = c(32, 45))

