# ================== 0. 加载包 ==================
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

library(GO.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ================== 1. 仅保留 monkey（注意：你代码中变量名是cell_human，但实际筛选的是monkey） ==================
cell_human <- subset(cell_total, species == "pig")

cell_human$celltype <- factor(
  cell_human$celltype,
  levels = c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")
)

# ================== 2. 定义核仁相关基因模块 ==================
modules_nucleolus <- list(
  rDNA_transcription = c(
    "UBTF", "RRN3", "POLR1A", "POLR1D", "TAF1D"
  ),
  
  pre_rRNA_processing = c(
    "FBL", "NOP56", "NOP58", "NOP10", "DKC1"
  ),
  
  ribosome_assembly = c(
    "RRS1", "RPF2", "NMD3", "SBDS", "BOP1"
  ),
  
  nucleolar_organization = c(
    "TCOF1", "NCL", "NPM1", "NOP2", "WDR12"
  )
)

modules_nucleolus <- list(
  # 1. rDNA转录：包含RNA聚合酶I核心亚基及起始因子
  rDNA_transcription = c(
    "UBTF", "RRN3", "POLR1A", "POLR1B", "POLR1D", "TAF1D"
  ),
  
  # 2. pre-rRNA加工：包含snoRNP复合物及PeB复合物（负责剪切与修饰）
  pre_rRNA_processing = c(
    "FBL", "NOP56", "NOP58", "NOP10", "DKC1", 
    "PES1", "BOP1", "WDR12"  # PeB Complex，核心加工单元
  ),
  
  # 3. 核糖体组装：侧重于亚基的成熟、组装与核输出衔接
  ribosome_assembly = c(
    "RRS1", "RPF2", "NMD3", "SBDS", "LSG1", "NVL"
  ),
  
  # 4. 核仁组织与结构：维持核仁相分离和颗粒区（GC）稳定的骨架蛋白
  nucleolar_organization = c(
    "TCOF1", "NCL", "NPM1", "NOP2", "SURF6"
  )
)

# ================== 3. 计算模块评分 ==================
cell_human <- AddModuleScore(
  object = cell_human,
  features = modules_nucleolus,
  name = "Nucleolus"
)

df_module <- FetchData(
  cell_human,
  vars = c(
    "celltype",
    "Nucleolus1",
    "Nucleolus2",
    "Nucleolus3",
    "Nucleolus4"
  )
)

# ================== 4. 按细胞类型计算平均评分 ==================
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

# ================== 5. 整理矩阵并标准化 ==================
mat <- t(as.matrix(df_avg))

mat_scaled <- t(scale(t(mat)))

# 按指定顺序排列行和列
mat_scaled <- mat_scaled[
  c("rDNA_tx", "pre_rRNA", "ribo_asm", "nucleolar_org"),
  c("Mito","Pre_L","Lepto","Zygot","Pachy","Diplo","Dicty")
]

# ================== 6. 定义发育阶段注释及配色（核心新增部分） ==================
# 1) 构建列注释数据框（行名与热图列名完全一致）
annotation_col <- data.frame(
  Stage = factor(colnames(mat_scaled),  # 列名直接是发育阶段
                 levels = c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty"))
)
rownames(annotation_col) <- colnames(mat_scaled)  # 强制对齐行名和热图列名

# 2) 定义Stage配色（沿用你指定的颜色）
ann_colors <- list(
  Stage = c(
    "Mito"   = "#6E4C9F",
    "Pre_L"  = "#E4989C",
    "Lepto"  = "#C06DAD",
    "Zygot"  = "#FEE800",
    "Pachy"  = "#7B68EE",
    "Diplo"  = "#F47D2C",
    "Dicty"  = "#89C75F"
  )
)

# ================== 7. 绘制带注释的热图（修改后的核心代码） ==================
pheatmap(
  mat = mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  # 新增注释相关参数
  annotation_col = annotation_col,    # 添加列注释（Stage）
  annotation_colors = ann_colors,     # 应用Stage配色
  annotation_names_col = TRUE,        # 显示注释标题（Stage）
  annotation_legend = TRUE,           # 显示Stage配色图例
  # 原有参数保留
  color = colorRampPalette(
    c("white", "grey90", "grey60", "grey30", "black")
  )(100),
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "",
  # 可选优化：调整注释字体大小
  annotation_fontsize = 10,
  # 可选：隐藏图例（如果不需要）
  # legend = FALSE
)


# ================== 7. 计算列间隙位置（核心美化） ==================
# 计算每个Stage对应的列数（你的数据中每个Stage对应1列）
stage_col_counts <- table(annotation_col$Stage)
# 计算间隙位置（在每个Stage列之后添加间隙，最后一列除外）
gaps_col <- if(ncol(mat_scaled) > 1) {
  cumsum(stage_col_counts)[-length(stage_col_counts)]
} else {
  NULL
}

# ================== 8. 最终美化版热图 ==================
pheatmap(
  mat = mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  # 注释相关
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  annotation_fontsize = 10,
  # 列间隙（核心美化）
  gaps_col = gaps_col,
  # 颜色优化
  color = colorRampPalette(
    c("white", "grey90", "grey60", "grey30", "black")
  )(100),
  # 边框与间距
  border_color = NA,
 # cellwidth = 40,  # 自定义列宽，让列更宽更清晰
#  cellheight = 30, # 自定义行高
  # 字体优化
  fontsize_row = 12,   # 行名（模块名）字体放大
  fontsize_col = 12,   # 列名（Stage）字体放大
  fontsize = 10,       # 整体字体大小
  # 标题与图例
  main = "",
  legend = TRUE,
  legend_breaks = c(-1, 0, 1, 2),  # 自定义图例刻度
  legend_labels = c("-1", "0", "1", "2"),
  # 其他美化
  treeheight_row = 0,  # 关闭行聚类树（本就关闭，显式指定更稳妥）
  treeheight_col = 0   # 关闭列聚类树
)
