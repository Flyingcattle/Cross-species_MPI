library(tidyverse)
library(cowplot)
library(tidyverse)
library(cowplot)
library(Seurat)

## load databse
Human_PSG <- readRDS(file="E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure3/Pleiotropy/Human.tau_genesymbol.rds")
head(Human_PSG)


###########################################################################################################   
#################################### load Seurat object ##############################
cell_total<-readRDS(file = "E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure2/germ_sub_updated.rds")
my_color3=c("#1F8A42", "#832688", "#F47D2C", "#FEE800", 
            "#6E4C9F", "#0C737B", "#282E6C", "#91D5E5",
            "#89C75F", "#E4989C", "#7B68EE", "#DB2228",
            "#D7A764", "#87A0D1", "#C06DAD", "#3DBBA7","#2CA029","#FD7F0D")
DimPlot(cell_total, cols = my_color3,
        reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(cell_total, cols = my_color3,
        reduction = "umap", label = T, group.by = "celltype")
table(cell_total@meta.data$celltype)
table(cell_total@meta.data$species)
# 细胞类型顺序
celltype_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")
# 去掉 NC_xxx
cell_total$celltype <- gsub("^NC_(6|7|8|14)$", "NC", cell_total$celltype)

#cell_total <- subset(cell_total, subset = species == "human" & celltype %in% celltype_order)
cell_total <- subset(x = cell_total,subset = celltype != "NC")  # 筛选 celltype 列中不等于 NC 的细胞

table(cell_total$celltype)

###########################################################################################################   
#################################### load Seurat object ##############################
###############################################
### 2. 提取 Seurat 数据
###############################################
seu <- cell_total

# Counts matrix
Human_counts <- seu@assays$RNA@counts

# meta 信息：CELL + CELL_TYPE
Human_meta_data <- data.frame(
  CELL = rownames(seu@meta.data),
  CELL_TYPE = seu@meta.data$celltype,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

###############################################
### 3. 数据整合（关键步骤）
###############################################
Human_counts %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENE") %>%                  # 基因作为第一列
  inner_join(Human_PSG, by = "GENE") %>%                # 合并基因 tau 数据
  gather(CELL, count, -GENE, -TissueTau) %>%            # 转长格式
  inner_join(Human_meta_data, by = "CELL") %>%          # 匹配细胞元信息
  filter(count > 0) %>%                                 # 只保留有表达的细胞
  filter(!is.na(TissueTau)) %>% 
  group_by(CELL, CELL_TYPE) %>%
  summarise(
    MEAN = mean(TissueTau),
    MEDIAN = median(TissueTau)
  ) %>% 
  ungroup() -> Human_gathered

###############################################
### 4. 排序（按你示例保持一致）
###############################################
Human_gathered$CELL_TYPE <- factor(
  Human_gathered$CELL_TYPE,
  levels = c("Mito", "Pre_L", "Lepto",
             "Zygot", "Pachy", "Diplo", 
             "Dicty")
)


###############################################
### 5. 绘图（与示例一致）
###############################################
ggplot(Human_gathered) +
  geom_boxplot(aes(CELL_TYPE, MEDIAN, fill = CELL_TYPE), notch = TRUE) +
  scale_fill_manual(values=c(
    "Mito" = "#99AEDA",
    "Pre_L" = "#5978BA",
    "Lepto" = "#B6D786",
    "Zygot"   = "#649744",
    "Pachy"       = "#2E572C",
    "Diplo"        = "#E6B648",
    "Dicty"         = "#DD7F38")) +
  geom_vline(xintercept = 2.5, color = "red") +
  theme_cowplot() +
  theme(
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("Tissue Tau")
