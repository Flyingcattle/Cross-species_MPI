# 加载Seurat包
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readxl)

setwd("/home/lab401_B/scRNA_downstream/scRNA_species/postive_selection")
###########################################################################################################   
#################################### load PLI score database ##############################
Human_PSG <- readRDS(file="/home/lab401_B/scRNA_downstream/scRNA_species/postive_selection/Human.PSG_genesymbol.rds")
head(Human_PSG)
dim(Human_PSG)    # 14798     2
#################################### load PLI score database ##############################
###########################################################################################################   

########################### load male female combined Seurat data ##############################
###########################################################################################################   
combined <- readRDS("/home/lab401_B/scRNA_downstream/scRNA_species/pLI/Seurat_sex_combined.rds")
seu <- combined
table(seu$species)

seu <- subset(combined, subset = sex == "female")
# 再检查过滤是否成功
table(seu$sex)

# extract count & meta
Human_counts <- seu@assays$RNA@counts
Human_meta_raw <- seu@meta.data
table(Human_meta_raw$species)
# format meta — 加入 species 信息（四物种比较关键）
Human_meta_data <- data.frame(
  CELL      = rownames(Human_meta_raw),
  CELL_TYPE = Human_meta_raw$celltype_myTAI,
  species   = Human_meta_raw$species,       # ★★★ 四物种比较关键字段
  check.names = FALSE,
  stringsAsFactors = FALSE
)

head(Human_meta_data)
table(Human_meta_data$CELL_TYPE)
########################### load male female combined Seurat data ##############################
###########################################################################################################   


# STEP 1 — Gather counts into long format
Human_counts %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(GENE = rowname) %>%
  gather(CELL, COUNTS, -GENE) %>%
  filter(COUNTS != 0) -> Human_Gathered

head(Human_Gathered)
head(Human_Gathered_PSG)

# STEP 2 — Join PSG
Human_Gathered %>%
  inner_join(
    Human_PSG %>% 
      rename(GENE = Tested_For_Postive_selection) %>% 
      mutate(Tested = 1),
    by = "GENE"
  ) -> Human_Gathered_PSG

# STEP 3 — Join meta & compute fraction
Human_Gathered_PSG %>%
  inner_join(Human_meta_data, by = "CELL") %>%
  group_by(CELL, CELL_TYPE, species) %>%     # ★★★ 用 species 来分组
  summarise(Frac_PSG = sum(PSG) / sum(Tested)) ->
  Human_Gathered_PSG_Summary

# STEP 4 — Set cell type order
Human_Gathered_PSG_Summary$CELL_TYPE <- factor(
  Human_Gathered_PSG_Summary$CELL_TYPE,
  levels = c("Mito","Pre_L","Lepto","Zygot","Pachy","Diplo","Dicty","Post_meiotic")
)

###########################################################################################################   
########################################### 物种颜色配置 #################################################

species_colors <- c(
  "human"  = "#BD4E93",
  "pig"    = "#DF7628",
  "monkey" = "#6BAD2B",
  "mouse"  = "#5879BD"
)

###########################################################################################################   
############################################ 绘图：四物种比较 #############################################

ggplot(Human_Gathered_PSG_Summary,
       aes(x = CELL_TYPE, y = Frac_PSG * 100, fill = species)) +
  
  geom_boxplot(
    alpha = 0.7,
    size = 0.3,
    position = position_dodge(width = 0.85),
    outlier.shape = 16,
    outlier.size = 0.5,
    outlier.colour = "black"
  ) +
  
  scale_fill_manual(values = species_colors) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title  = element_text(color = "black"),
    legend.text  = element_text(color = "black"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  
  labs(
    x = "",
    y = "Positive selection (%)",
    title = "Comparing PSG Fraction Across Four Species"
  )


