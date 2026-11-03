# ================== 0. 加载包 ==================
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

library(GO.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

cell_human <- subset(cell_total, species == "pig")
cell_human$celltype <- factor(
  cell_human$celltype,
  levels = c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty"))
modules_nucleolus <- list(
  rDNA_transcription = c("UBTF", "RRN3", "POLR1A", "POLR1D", "TAF1D"),
  pre_rRNA_processing = c("FBL", "NOP56", "NOP58", "NOP10", "DKC1"),
  ribosome_assembly = c("RRS1", "RPF2", "NMD3", "SBDS", "BOP1"),
  nucleolar_organization = c("TCOF1", "NCL", "NPM1", "NOP2", "WDR12"))

modules_nucleolus <- list(
  rDNA_transcription = c("UBTF", "RRN3", "POLR1A", "POLR1B", "POLR1D", "TAF1D"),
  pre_rRNA_processing = c("FBL", "NOP56", "NOP58", "NOP10", "DKC1", "PES1", "BOP1", "WDR12"),
  ribosome_assembly = c("RRS1", "RPF2", "NMD3", "SBDS", "LSG1", "NVL"),
  nucleolar_organization = c("TCOF1", "NCL", "NPM1", "NOP2", "SURF6"))

cell_human <- AddModuleScore(
  object = cell_human,
  features = modules_nucleolus,
  name = "Nucleolus")

df_module <- FetchData(
  cell_human,
  vars = c(
    "celltype",
    "Nucleolus1",
    "Nucleolus2",
    "Nucleolus3",
    "Nucleolus4"))

df_avg <- df_module %>%
  group_by(celltype) %>%
  summarise(
    rDNA_tx      = mean(Nucleolus1),
    pre_rRNA     = mean(Nucleolus2),
    ribo_asm     = mean(Nucleolus3),
    nucleolar_org= mean(Nucleolus4)
  ) %>%
  as.data.frame()

rownames(df_avg) <- df_avg$celltype
df_avg$celltype <- NULL
mat <- t(as.matrix(df_avg))
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[
  c("rDNA_tx", "pre_rRNA", "ribo_asm", "nucleolar_org"),
  c("Mito","Pre_L","Lepto","Zygot","Pachy","Diplo","Dicty")]

annotation_col <- data.frame(
  Stage = factor(colnames(mat_scaled),  
                 levels = c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")))
rownames(annotation_col) <- colnames(mat_scaled) 

ann_colors <- list(
  Stage = c(
    "Mito"   = "#6E4C9F",
    "Pre_L"  = "#E4989C",
    "Lepto"  = "#C06DAD",
    "Zygot"  = "#FEE800",
    "Pachy"  = "#7B68EE",
    "Diplo"  = "#F47D2C",
    "Dicty"  = "#89C75F"))

stage_col_counts <- table(annotation_col$Stage)
gaps_col <- if(ncol(mat_scaled) > 1) {
  cumsum(stage_col_counts)[-length(stage_col_counts)]
} else {
  NULL
}

pheatmap(
  mat = mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  annotation_fontsize = 10,
  gaps_col = gaps_col,
  color = colorRampPalette(
    c("white", "grey90", "grey60", "grey30", "black")
  )(100),
  border_color = NA,
  fontsize_row = 12,   
  fontsize_col = 12,   
  fontsize = 10,       
  main = "",
  legend = TRUE,
  legend_breaks = c(-1, 0, 1, 2),  
  legend_labels = c("-1", "0", "1", "2"),
  treeheight_row = 0,  
  treeheight_col = 0)
