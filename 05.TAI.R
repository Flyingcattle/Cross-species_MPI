# 加载Seurat包
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(myTAI)
library(readxl)

setwd("E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure3/Human_seprmatogenesis_test/pLI_score")

############### load sperm data 
seu<-readRDS(file = "./human_sperm_celltype3.rds")
DimPlot(seu, reduction = "umap", label = T, group.by = "celltype3")
table(seu@meta.data$celltype3)
head(seu@meta.data)

# 添加myTAI注释信息
seu@meta.data$celltype_myTAI<-seu@meta.data$celltype3
seu@meta.data$gender<-c("male")
seu@meta.data$group_myTAI<-seu@meta.data$sample_id

###############  load oocyte data 
cell_total<-readRDS(file = "E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure2/germ_sub_updated.rds")
DimPlot(cell_total, reduction = "umap", label = T, group.by = "celltype")
# 定义细胞类型顺序
celltype_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")
# 简化细胞类型标签
cell_total$celltype <- gsub("^NC_(6|7|8|14)$", "NC", cell_total$celltype)
# 过滤指定细胞类型
cell_total <- subset(cell_total, subset = celltype %in% celltype_order)
DimPlot(cell_total, reduction = "umap", label = T, group.by = "celltype")
head(cell_total@meta.data)

# 添加myTAI注释信息
cell_total@meta.data$celltype_myTAI<-cell_total@meta.data$celltype
cell_total@meta.data$gender<-c("female")
cell_total@meta.data$group_myTAI<-cell_total@meta.data$PCW_time

seu_male<-seu
seu_female<-cell_total
head(seu_male)
head(seu_female)

table(seu_male@meta.data$group_myTAI)
table(seu_female@meta.data$group_myTAI)

######################### combined dataset ##############
seu_male$sex <- "male"
seu_female$sex <- "female"

seu_list <- list(male = seu_male, female = seu_female)

# 1. Normalize + HVG
seu_list <- lapply(seu_list, function(x){
  NormalizeData(x)
  FindVariableFeatures(x)
})

# 2. 找 anchors
anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:30)

# 3. 整合
combined <- IntegrateData(anchorset = anchors, dims = 1:30)


# 6️⃣ 设置默认 assay
DefaultAssay(combined) <- "integrated"

# 7️⃣ 后续分析（可选）
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)

# 8️⃣ 查看整合结果
DimPlot(combined, reduction = "umap", group.by = "sex")
DimPlot(combined, reduction = "umap", group.by = "group_myTAI")
DimPlot(combined, reduction = "umap", group.by = "celltype_myTAI")

DimPlot(combined, reduction = "umap", label=T,group.by = "celltype_myTAI")

saveRDS(combined,file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\Figure3\\myTAI\\Seurat_sex_combined.rds")
######################### combined dataset ##############

combined<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\Figure3\\myTAI\\Seurat_sex_combined.rds")
cell_total<-combined
# 定义细胞类型顺序
celltype_order <- c("Mito", "Pre_L", "Lepto", 
                    "Zygot", "Pachy", "Diplo",
                    "Dicty","Post_meiotic")

table(cell_total@meta.data$sample_id)

# -----------------------------
# 1) AverageExpression by sex + celltype
# -----------------------------
avg_exp_sex_celltype <- AverageExpression(
  cell_total,
  group.by = c("sex", "celltype_myTAI"),
  assays = "RNA",
  slot = "data"
)$RNA
dim(avg_exp_sex_celltype)   # 36988    14
# -----------------------------
# 2) Split male / female
# -----------------------------
# male
avg_exp_male <- avg_exp_sex_celltype[, grepl("^male_", colnames(avg_exp_sex_celltype))]
colnames(avg_exp_male) <- sub("^male_", "", colnames(avg_exp_male))
avg_exp_male <- data.frame(GeneID = rownames(avg_exp_male), avg_exp_male, stringsAsFactors = FALSE)

# female
avg_exp_female <- avg_exp_sex_celltype[, grepl("^female_", colnames(avg_exp_sex_celltype))]
colnames(avg_exp_female) <- sub("^female_", "", colnames(avg_exp_female))
avg_exp_female <- data.frame(GeneID = rownames(avg_exp_female), avg_exp_female, stringsAsFactors = FALSE)

# -----------------------------
# 3) Merge phylorank
# -----------------------------
setwd("E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure3/myTAI/1.phylomap")
phylorank_df <- read_excel("hg19_ver73_age_symbol_filtered.xlsx")
# male
avg_exp_male <- avg_exp_male %>%
  left_join(phylorank_df %>% select(GeneID, phylorank), by = "GeneID") %>%
  select(phylorank, GeneID, everything()) %>%
  filter(!is.na(phylorank))

# female
avg_exp_female <- avg_exp_female %>%
  left_join(phylorank_df %>% select(GeneID, phylorank), by = "GeneID") %>%
  select(phylorank, GeneID, everything()) %>%
  filter(!is.na(phylorank))

# -----------------------------
# 4) 按各自存在的细胞类型顺序排列
# -----------------------------
# male 的顺序（没有 Dicty）
male_order <- c("Mito", "Pre_L", "Lepto", 
                "Zygot", "Pachy", "Diplo", "Post_meiotic")
male_cols_exist <- intersect(male_order, colnames(avg_exp_male))
avg_exp_male <- avg_exp_male[, c("phylorank", "GeneID", male_cols_exist)]

# female 的顺序（没有 Post_meiotic）
female_order <- c("Mito", "Pre_L", "Lepto", 
                  "Zygot", "Pachy", "Diplo", "Dicty")
female_cols_exist <- intersect(female_order, colnames(avg_exp_female))
avg_exp_female <- avg_exp_female[, c("phylorank", "GeneID", female_cols_exist)]

# -----------------------------
# 5) 构建 PhyloExpressionSet
# -----------------------------
male_phyex <- myTAI::BulkPhyloExpressionSet_from_df(avg_exp_male, name = "male_TAI")
female_phyex <- myTAI::BulkPhyloExpressionSet_from_df(avg_exp_female, name = "female_TAI")

# -----------------------------
# 6) 绘图
# -----------------------------
myTAI::plot_signature(male_phyex)
myTAI::plot_signature(female_phyex)


################################  bootstrap TAI for male & female ##################################

# -----------------------------
# 1) 提取表达矩阵和 phylorank
# -----------------------------
expr_male <- male_phyex@expression
phy_male  <- male_phyex@strata

expr_female <- female_phyex@expression
phy_female  <- female_phyex@strata

# -----------------------------
# 2) 定义 bootstrap TAI 函数
# -----------------------------
bootstrap_TAI <- function(expr_mat, phylorank, n_boot = 1000) {
  
  phy_num <- as.numeric(as.character(phylorank))
  gene_n <- nrow(expr_mat)
  celltypes <- colnames(expr_mat)
  
  boot_res <- matrix(NA, nrow = n_boot, ncol = length(celltypes))
  colnames(boot_res) <- celltypes
  
  for (b in 1:n_boot) {
    # 基因 bootstrap sampling
    idx <- sample(1:gene_n, gene_n, replace = TRUE)
    
    expr_b <- expr_mat[idx, , drop = FALSE]
    phy_b  <- phy_num[idx]
    
    numerator <- colSums(expr_b * phy_b, na.rm = TRUE)
    denominator <- colSums(expr_b, na.rm = TRUE)
    denominator[denominator == 0] <- NA  # 避免除 0
    
    boot_res[b, ] <- numerator / denominator
  }
  
  boot_res <- as.data.frame(boot_res)
  boot_res[sapply(boot_res, is.nan)] <- NA
  return(boot_res)
}

# -----------------------------
# 3) 设置随机种子
# -----------------------------
set.seed(100)

# male bootstrap
tai_boot_male <- bootstrap_TAI(expr_male, phy_male, n_boot = 1000)
tai_boot_male_long <- tai_boot_male %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "celltype",
                      values_to = "TAI")
# factor 顺序
male_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Post_meiotic")
tai_boot_male_long$celltype <- factor(tai_boot_male_long$celltype,
                                      levels = male_order)

# female bootstrap
tai_boot_female <- bootstrap_TAI(expr_female, phy_female, n_boot = 1000)
tai_boot_female_long <- tai_boot_female %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "celltype",
                      values_to = "TAI")
# factor 顺序
female_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")
tai_boot_female_long$celltype <- factor(tai_boot_female_long$celltype,
                                        levels = female_order)

# -----------------------------
# 4) 绘图
# -----------------------------
library(ggplot2)

# male TAI
ggplot(tai_boot_male_long, aes(x = celltype, y = TAI)) +
  geom_boxplot(outlier.shape = NA, fill = "steelblue", alpha = 0.7) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = "Cell type", y = "Transcriptome Age Index (TAI)",
       title = "Male: Bootstrap-based TAI (1,000 replicates)")

# female TAI
ggplot(tai_boot_female_long, aes(x = celltype, y = TAI)) +
  geom_boxplot(outlier.shape = NA, fill = "tomato", alpha = 0.7) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = "Cell type", y = "Transcriptome Age Index (TAI)",
       title = "Female: Bootstrap-based TAI (1,000 replicates)")


############################### split by sex #############################
################################  bootstrap TAI for male & female (facet plot) ##################################

# -----------------------------
# 1) 定义男女各自的细胞类型顺序
# -----------------------------
male_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Post_meiotic")
female_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")

# -----------------------------
# 2) 准备表达矩阵和 phylorank
# -----------------------------
# male
avg_exp_male <- avg_exp_male[, c("phylorank", male_order)]
expr_mat_male <- avg_exp_male %>%
  dplyr::select(-phylorank) %>%
  dplyr::mutate(across(everything(), as.numeric))

# female
avg_exp_female <- avg_exp_female[, c("phylorank", female_order)]
expr_mat_female <- avg_exp_female %>%
  dplyr::select(-phylorank) %>%
  dplyr::mutate(across(everything(), as.numeric))

# -----------------------------
# 3) 定义 bootstrap TAI 函数
# -----------------------------
bootstrap_TAI <- function(expr_mat, phylorank, n_boot = 1000) {
  phy_num <- as.numeric(as.character(phylorank))
  gene_n <- nrow(expr_mat)
  celltypes <- colnames(expr_mat)
  
  boot_res <- matrix(NA, nrow = n_boot, ncol = length(celltypes))
  colnames(boot_res) <- celltypes
  
  for (b in 1:n_boot) {
    idx <- sample(1:gene_n, gene_n, replace = TRUE)
    expr_b <- expr_mat[idx, , drop = FALSE]
    phy_b  <- phy_num[idx]
    
    numerator <- colSums(expr_b * phy_b, na.rm = TRUE)
    denominator <- colSums(expr_b, na.rm = TRUE)
    denominator[denominator == 0] <- NA
    
    boot_res[b, ] <- numerator / denominator
  }
  
  boot_res <- as.data.frame(boot_res)
  boot_res[sapply(boot_res, is.nan)] <- NA
  return(boot_res)
}

# -----------------------------
# 4) 计算 bootstrap TAI
# -----------------------------
set.seed(100)
tai_boot_male <- bootstrap_TAI(expr_mat_male, avg_exp_male$phylorank, n_boot = 1000)
tai_boot_female <- bootstrap_TAI(expr_mat_female, avg_exp_female$phylorank, n_boot = 1000)

# -----------------------------
# 5) 转长格式并添加 sex 标签
# -----------------------------
tai_boot_male_long <- tai_boot_male %>%
  tidyr::pivot_longer(cols = everything(), names_to = "celltype", values_to = "TAI") %>%
  mutate(sex = "male")

tai_boot_female_long <- tai_boot_female %>%
  tidyr::pivot_longer(cols = everything(), names_to = "celltype", values_to = "TAI") %>%
  mutate(sex = "female")

tai_boot_long <- bind_rows(tai_boot_male_long, tai_boot_female_long)

# -----------------------------
# 6) 设置 factor 顺序
# -----------------------------
tai_boot_long$celltype <- factor(tai_boot_long$celltype,
                                 levels = unique(c(male_order, female_order)))

tai_boot_long$sex <- factor(tai_boot_long$sex, levels = c("male", "female"))

# -----------------------------
# 7) 绘图：Box plot，按 sex 分面
# -----------------------------
# -----------------------------
# 7) 绘图：Box plot，sex 并列显示每个细胞类型
# -----------------------------
ggplot(tai_boot_long, aes(x = celltype, y = TAI, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell type",
    y = "Transcriptome Age Index (TAI)",
    title = ""
  )
ggplot(tai_boot_long, aes(x = celltype, y = TAI, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("male" = "#E46729", "female" = "#83589C")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.position = "top",           # 图例放上方
    legend.direction = "horizontal",   # 图例水平排列
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "",
    y = "TAI index",
    title = ""
  )

#######################  折线图  ##########################
# 计算每个 celltype × sex 的平均值和 95% CI
tai_ci <- tai_boot_long %>%
  group_by(celltype, sex) %>%
  summarise(
    mean_TAI = mean(TAI, na.rm = TRUE),
    lower = quantile(TAI, probs = 0.025, na.rm = TRUE),
    upper = quantile(TAI, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# 绘制折线图 + 阴影
ggplot(tai_ci, aes(x = celltype, y = mean_TAI, group = sex, color = sex, fill = sex)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("male" = "#E46729", "female" = "#83589C")) +
  scale_fill_manual(values = c("male" = "#E46729", "female" = "#83589C")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "",
    y = "TAI index",
    title = ""
  )


ggplot(tai_ci, aes(x = celltype, y = mean_TAI, group = sex)) +
  # 灰色置信区间
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  
  # 黑色折线
  geom_line(color = "black", size = 1) +
  
  # 不同性别用不同点形状区分
  geom_point(aes(shape = sex), size = 3, color = "black") +
  
  scale_shape_manual(values = c("male" = 17, "female" = 16)) +  # 16圆点，17三角
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "Cell type",
    y = "Transcriptome Age Index (TAI)",
    title = ""
  )

#########################  reanalysis boxplot  #####################################

# 自定义顺序
new_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo","Dicty", "Post_meiotic")

tai_boot_long$celltype <- factor(tai_boot_long$celltype, levels = new_order)

table(tai_boot_long$celltype)

ggplot(tai_boot_long, aes(x = celltype, y = TAI, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("male" = "#E46729", "female" = "#83589C")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "",
    y = "TAI index",
    title = ""
  )




ggplot(tai_boot_long, aes(x = celltype, y = TAI, fill = sex)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 0.5,
    outlier.colour = "black",
    alpha = 0.7,
    size = 0.1,
    position = position_dodge(width = 0.85)
  ) +
  
  scale_fill_manual(values = c(
    "male"   = "#E46729",
    "female" = "#83589C"
  )) +
  
  theme_bw(base_size = 14) +
  theme(
    # —— 去掉刻度文字 —— #
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    
    # —— 保留刻度线 —— #
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    
    # —— 去掉坐标轴标题 —— #
    axis.title = element_blank(),
    
    # —— 去掉 legend —— #
    legend.position = "none",
    
    # —— 背景与边框统一 —— #
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(title = "")


################################# 折线图 ###############################
library(ggplot2)
ggplot(tai_ci, aes(x = celltype, y = mean_TAI, group = sex, color = sex, fill = sex)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("male" = "#E46729", "female" = "#83589C")) +
  scale_fill_manual(values = c("male" = "#E46729", "female" = "#83589C")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_blank(),       # 去掉 X 轴刻度文字
    axis.text.y = element_blank(),       # 去掉 Y 轴刻度文字
    axis.title.x = element_blank(),      # 去掉 X 轴标题
    axis.title.y = element_blank(),      # 去掉 Y 轴标题
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none",            # 去掉 legend
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(title = "")
