# 加载Seurat包
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(dplyr)
library(tibble)

cell_total<-readRDS(file = "E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure2/germ_sub_updated.rds")

#########  publication ready plot 
my_color3=c("#1F8A42", "#832688", "#F47D2C", "#FEE800", 
            "#6E4C9F", "#0C737B", "#282E6C", "#91D5E5",
            "#89C75F", "#E4989C", "#7B68EE", "#DB2228",
            "#D7A764", "#87A0D1", "#C06DAD", "#3DBBA7","#2CA029","#FD7F0D")

DimPlot(cell_total, cols = my_color3,
        reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(cell_total, cols = my_color3,
        reduction = "umap", label = T, group.by = "celltype")

# Step 1: 创建组合标签
cell_total$group <- paste(cell_total$celltype, cell_total$species, sep = "_")

# Step 2: 获取每个 group 的平均表达
#avg_exp <- AverageExpression(cell_total, group.by = "group", slot = "data")$RNA  # 获取 RNA slot 的平均表达

avg_exp <- AverageExpression(cell_total, group.by = "group", assays = "integrated", slot = "scale.data")$integrated
head(avg_exp)
# Step 3: 转置表达矩阵（行为样本，列为基因）
exp_mat <- t(as.matrix(avg_exp))

# Step 4: 去除方差为0的基因
var_genes <- apply(exp_mat, 2, var)
exp_mat_filtered <- exp_mat[, var_genes > 0]

# Step 5: PCA 分析
pca_res <- prcomp(exp_mat_filtered,center = TRUE, scale. = TRUE)
PCA<-pca_res
# Getting PCA variance
summary(PCA)$importance %>%
  as.data.frame() -> Variance_PCA


# Getting 10 first components
Variance_PCA[2,] %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>% 
  dplyr::rename(components = rowname) %>% 
  dplyr::rename(Proportion_of_Variance = "Proportion of Variance") %>% 
  filter(components == "PC1" |
           components == "PC2" |
           components == "PC3" |
           components == "PC4" |
           components == "PC5" |
           components == "PC6" |
           components == "PC7" |
           components == "PC8" |
           components == "PC9" |
           components == "PC10" ) -> Variance_PCA

# Ordering PCs
Variance_PCA$components <- factor( Variance_PCA$components , levels=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10") )

# Plotting variance
Variance_PCA %>% 
  ggplot() +
  geom_bar( aes(components,Proportion_of_Variance*100),stat = "identity" )+
  geom_text(aes(x=components,y=Proportion_of_Variance*100,label=sprintf("%0.1f", round(Proportion_of_Variance*100, digits = 1))),vjust=-0.5) +
  theme_cowplot()


# Getting species and cell types
PCA$x %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  mutate(Species = case_when( grepl("human",rowname) ~ "human",
                              grepl("monkey",rowname) ~ "monkey",
                              grepl("mouse",rowname) ~ "mouse",
                              grepl("pig",rowname) ~ "pig")) %>% 
  mutate(Cell_type = case_when( grepl("Dicty",rowname) ~ "Dicty",
                                grepl("Diplo",rowname) ~ "Diplo",
                                grepl("Lepto",rowname) ~ "Lepto",
                                grepl("NC_6",rowname) ~ "NC_6",
                                grepl("NC_7",rowname) ~ "NC_7",
                                grepl("NC_8",rowname) ~ "NC_8",
                                grepl("NC_14",rowname) ~ "NC_14",
                                grepl("Pachy",rowname) ~ "Pachy",
                                grepl("Mito",rowname) ~ "Mito",
                                grepl("Zygot",rowname) ~ "Zygot",
                                grepl("Pre_L",rowname) ~ "Pre_L") ) -> PCA_PC
# Ordering cell types
PCA_PC$Cell_type <- factor(PCA_PC$Cell_type,
                           levels=c("Mito",
                                    "Pre_L",
                                    "Lepto",
                                    "Zygot",
                                    "Pachy",
                                    "Diplo",
                                    "Dicty",
                                    "NC_6",
                                    "NC_7",
                                    "NC_8",
                                    "NC_14"))

# Ordering species
PCA_PC$Species <- factor(PCA_PC$Species,
                         levels=c("human",
                                  "monkey",
                                  "mouse",
                                  "pig")) 


# Plotting PCA
PCA_PC %>% 
  ggplot() +
  geom_point( aes(PC1,PC2,color=Cell_type,shape= Species),size=2,stroke = 1) +
  scale_shape_manual(values=1:nlevels(PCA_PC$Species)) +
  scale_color_manual(values= c("Mito"="#6E4C9F",
                               "Pre_L"="#E4989C",
                               "Lepto"="#C06DAD",
                               "Zygot"="#FEE800",
                               "Pachy"="#7B68EE",
                               "Diplo"="#F47D2C",
                               "Dicty"="#89C75F",
                               "NC_6"="#050505",
                               "NC_7"="#050505",
                               "NC_8"="#050505",
                               "NC_14"="#050505") ) +
  theme_cowplot()


PCA_PC %>% 
  ggplot() +
  geom_point( aes(PC1,PC2,color=Cell_type,shape= Species),size=2,stroke = 1) +
  scale_shape_manual(values=1:nlevels(PCA_PC$Species)) +
  scale_color_manual(values= c("Mito"="#6E4C9F",
                               "Pre_L"="#E4989C",
                               "Lepto"="#C06DAD",
                               "Zygot"="#FEE800",
                               "Pachy"="#7B68EE",
                               "Diplo"="#F47D2C",
                               "Dicty"="#89C75F",
                               "NC_6"="red",
                               "NC_7"="blue",
                               "NC_8"="black",
                               "NC_14"="green") ) +
  theme_cowplot()


