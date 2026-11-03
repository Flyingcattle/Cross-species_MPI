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

seu@meta.data$celltype_myTAI<-seu@meta.data$celltype3
seu@meta.data$gender<-c("male")
seu@meta.data$group_myTAI<-seu@meta.data$sample_id

###############  load oocyte data 
cell_total<-readRDS(file = "E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure2/germ_sub_updated.rds")
DimPlot(cell_total, reduction = "umap", label = T, group.by = "celltype")
celltype_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")
cell_total$celltype <- gsub("^NC_(6|7|8|14)$", "NC", cell_total$celltype)
cell_total <- subset(cell_total, subset = celltype %in% celltype_order)
DimPlot(cell_total, reduction = "umap", label = T, group.by = "celltype")
head(cell_total@meta.data)

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

seu_list <- lapply(seu_list, function(x){
  NormalizeData(x)
  FindVariableFeatures(x)
})

anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:30)

combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)

DimPlot(combined, reduction = "umap", group.by = "sex")
DimPlot(combined, reduction = "umap", group.by = "group_myTAI")
DimPlot(combined, reduction = "umap", group.by = "celltype_myTAI")

DimPlot(combined, reduction = "umap", label=T,group.by = "celltype_myTAI")

saveRDS(combined,file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\Figure3\\myTAI\\Seurat_sex_combined.rds")
######################### combined dataset ##############

combined<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\1.All_species_Seurat\\Figure3\\myTAI\\Seurat_sex_combined.rds")
cell_total<-combined
celltype_order <- c("Mito", "Pre_L", "Lepto", 
                    "Zygot", "Pachy", "Diplo",
                    "Dicty","Post_meiotic")

table(cell_total@meta.data$sample_id)

avg_exp_sex_celltype <- AverageExpression(
  cell_total,
  group.by = c("sex", "celltype_myTAI"),
  assays = "RNA",
  slot = "data"
)$RNA
dim(avg_exp_sex_celltype)   # 36988    14

avg_exp_male <- avg_exp_sex_celltype[, grepl("^male_", colnames(avg_exp_sex_celltype))]
colnames(avg_exp_male) <- sub("^male_", "", colnames(avg_exp_male))
avg_exp_male <- data.frame(GeneID = rownames(avg_exp_male), avg_exp_male, stringsAsFactors = FALSE)

avg_exp_female <- avg_exp_sex_celltype[, grepl("^female_", colnames(avg_exp_sex_celltype))]
colnames(avg_exp_female) <- sub("^female_", "", colnames(avg_exp_female))
avg_exp_female <- data.frame(GeneID = rownames(avg_exp_female), avg_exp_female, stringsAsFactors = FALSE)

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

male_order <- c("Mito", "Pre_L", "Lepto", 
                "Zygot", "Pachy", "Diplo", "Post_meiotic")
male_cols_exist <- intersect(male_order, colnames(avg_exp_male))
avg_exp_male <- avg_exp_male[, c("phylorank", "GeneID", male_cols_exist)]

# female 的顺序（没有 Post_meiotic）
female_order <- c("Mito", "Pre_L", "Lepto", 
                  "Zygot", "Pachy", "Diplo", "Dicty")
female_cols_exist <- intersect(female_order, colnames(avg_exp_female))
avg_exp_female <- avg_exp_female[, c("phylorank", "GeneID", female_cols_exist)]

male_phyex <- myTAI::BulkPhyloExpressionSet_from_df(avg_exp_male, name = "male_TAI")
female_phyex <- myTAI::BulkPhyloExpressionSet_from_df(avg_exp_female, name = "female_TAI")

myTAI::plot_signature(male_phyex)
myTAI::plot_signature(female_phyex)


################################  bootstrap TAI for male & female ##################################

expr_male <- male_phyex@expression
phy_male  <- male_phyex@strata

expr_female <- female_phyex@expression
phy_female  <- female_phyex@strata

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
    denominator[denominator == 0] <- NA  # 避免除 0
    boot_res[b, ] <- numerator / denominator
  }
  boot_res <- as.data.frame(boot_res)
  boot_res[sapply(boot_res, is.nan)] <- NA
  return(boot_res)
}

set.seed(100)

tai_boot_male <- bootstrap_TAI(expr_male, phy_male, n_boot = 1000)
tai_boot_male_long <- tai_boot_male %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "celltype",
                      values_to = "TAI")
male_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Post_meiotic")
tai_boot_male_long$celltype <- factor(tai_boot_male_long$celltype,
                                      levels = male_order)

tai_boot_female <- bootstrap_TAI(expr_female, phy_female, n_boot = 1000)
tai_boot_female_long <- tai_boot_female %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "celltype",
                      values_to = "TAI")
female_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")
tai_boot_female_long$celltype <- factor(tai_boot_female_long$celltype,
                                        levels = female_order)

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

