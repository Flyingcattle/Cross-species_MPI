# 加载必要的包
library(reshape2)
library(dplyr)
library(tidyr)

setwd("E:/Lab of Germ Cell Bio/10xGenomics_species/1.All_species_Seurat/Figure3/conserved_gene_characterazation")

# 1. 定义文件路径和物种名称
species_files <- list(
  "human" = "all_genes_avg_matrix_human.csv",
  "monkey" = "all_genes_avg_matrix_monkey.csv",
  "pig" = "all_genes_avg_matrix_pig.csv",
  "mouse" = "all_genes_avg_matrix_mouse.csv"
)

# 2. 定义发育阶段顺序
stage_order <- c("Mito", "Pre_L", "Lepto", "Zygot", "Pachy", "Diplo", "Dicty")

# 3. 数据标准化函数（z-score标准化：(x-mean)/sd）
normalize_data <- function(matrix_data) {
  # 对每行（每个基因）进行标准化
  t(apply(matrix_data, 1, function(row) {
    row_mean <- mean(row, na.rm = TRUE)
    row_sd <- sd(row, na.rm = TRUE)
    # 避免标准差为0导致的除以零问题
    if (row_sd == 0) {
      return(rep(0, length(row)))
    } else {
      return((row - row_mean) / row_sd)
    }
  }))
}

# 4. 读取、过滤并标准化四个物种的数据
species_data <- list()
for (sp in names(species_files)) {
  # 读取数据
  data <- read.csv(species_files[[sp]], row.names = 1, check.names = FALSE)
  
  # 按统一时间序列排序
  data_ordered <- data[, stage_order, drop = FALSE]
  
  # 过滤1：移除所有阶段表达值均为0的基因
  non_zero_rows <- rowSums(data_ordered) > 0
  data_filtered <- data_ordered[non_zero_rows, , drop = FALSE]
  cat("在", sp, "数据中过滤掉", sum(!non_zero_rows), "个全零基因，保留", nrow(data_filtered), "个基因\n")
  
  # 标准化：对每个基因进行z-score标准化
  data_normalized <- normalize_data(as.matrix(data_filtered))
  # 将标准化后的数据转回数据框
  data_normalized_df <- as.data.frame(data_normalized, check.names = FALSE)
  colnames(data_normalized_df) <- stage_order
  
  # 提取基因名（去掉物种前缀）
  gene_names <- sub(paste0("^", sp, "_"), "", rownames(data_filtered))
  
  # 存储处理后的数据（包含基因名）
  species_data[[sp]] <- data.frame(
    Gene = gene_names,
    data_normalized_df,
    row.names = rownames(data_filtered),
    check.names = FALSE
  )
}

# 5. 获取所有物种共有的非零基因
all_genes <- lapply(species_data, function(x) unique(x$Gene))
common_genes <- Reduce(intersect, all_genes)
cat("\n四个物种共有的非全零基因数量：", length(common_genes), "\n")
 ## 四个物种共有的非全零基因数量： 7952

# 6. 计算保守性系数（使用标准化后的数据）
conservation_results <- data.frame(
  Gene = character(),
  Conservation_Score = numeric(),
  stringsAsFactors = FALSE
)

for (gene in common_genes) {
  expr_list <- list()
  valid <- TRUE
  
  for (sp in names(species_data)) {
    sp_gene_row <- which(species_data[[sp]]$Gene == gene)
    if (length(sp_gene_row) == 0) {
      valid <- FALSE
      break
    }
    
    # 提取标准化后的表达值
    expr_values <- as.numeric(species_data[[sp]][sp_gene_row, stage_order])
    expr_list[[sp]] <- expr_values
  }
  
  if (!valid) next
  
  # 构建表达矩阵（阶段×物种）
  expr_matrix <- do.call(cbind, expr_list)
  rownames(expr_matrix) <- stage_order
  
  # 计算皮尔逊相关系数
  cor_matrix <- cor(expr_matrix, method = "pearson")
  
  # 计算平均两两相关性
  cor_values <- cor_matrix[lower.tri(cor_matrix)]
  conservation_score <- mean(cor_values, na.rm = TRUE)
  
  conservation_results <- rbind(conservation_results,
                                data.frame(Gene = gene,
                                           Conservation_Score = conservation_score))
}

# 7. 排序并输出结果
conservation_results <- conservation_results %>%
  arrange(desc(Conservation_Score))

# 过滤出保守性系数大于0.5的基因
filtered_conservation <- conservation_results %>%
  filter(Conservation_Score > 0.5)
dim(filtered_conservation)  ## 1104

# 按保守性系数从高到低排序（保持原有排序逻辑）
filtered_conservation <- filtered_conservation %>%
  arrange(desc(Conservation_Score))


# 仅保存保守性系数>0.5的基因数据
write.csv(filtered_conservation, "gene_conservation_scores_above0.5_4species.csv", row.names = FALSE)
###########################################################################




