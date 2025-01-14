## ====================================================================
## DAY 3. R语言数据清洗-TCGA&GEO
## 
## 描述：
##   从GEO数据库下载阿尔茨海默病相关数据集，并进行预处理和标准化
##
## 数据集：
##   - GSE48350: AD和对照组海马体样本
##   - GSE5281: AD和对照组海马体样本
##
## 输出：
##   - raw_geo_data.Rda: 原始GEO数据
##   - normalized_expression_data.Rda: 标准化后的表达矩阵
##   - batch_effect_plots.pdf: 批次效应校正前后的可视化图
##
## 注意事项：
##   - 需要稳定的网络连接
##   - 数据下载可能需要较长时间
##   - 建议预留足够的内存空间
## ====================================================================

#####################################################
## Step1：环境初始化
## 
## 描述：初始化R环境，设置参数，加载必要的包
#####################################################

# 清空工作环境
rm(list = ls())

# 基础设置
options(
    scipen = 20,      # 避免科学计数法
    timeout = 1000    # 下载超时设置(秒)
)
setwd("./D3-GEO数据下载预处理")   # 设置工作目录

# 加载必要的R包
library(GEOquery)    # GEO数据下载和处理
library(tidyverse)   # 数据处理和管道操作
library(limma)       # 表达数据标准化
library(sva)         # 批次效应校正
library(ggplot2)     # 数据可视化
library(patchwork)   # 图形组合工具

#####################################################
## Step2：数据获取
## 
## 描述：从GEO数据库下载数据集并保存
#####################################################

# 下载GSE48350数据集
gse_1 <- getGEO(
    "GSE48350",
    destdir = "./input",
    AnnotGPL = FALSE,
    getGPL = FALSE
)

# 下载GSE5281数据集
gse_2 <- getGEO(
    "GSE5281",
    destdir = "./input",
    AnnotGPL = FALSE,
    getGPL = FALSE
)

# 保存原始数据
save(gse_1, gse_2, file = 'raw_geo_data.Rda')

#####################################################
## Step3：数据预处理
## 
## 描述：提取海马体样本数据并进行初步处理
#####################################################

# 处理GSE48350数据集
data1 <- gse_1[[1]] %>%
    pData() %>%
    filter(grepl("hippocampus", title, ignore.case = TRUE)) %>%
    rownames() %>%
    {exprs(gse_1[[1]])[, colnames(exprs(gse_1[[1]])) %in% .]}

# 处理GSE5281数据集
data2 <- gse_2[[1]] %>%
    pData() %>%
    filter(grepl("HIP", title)) %>%
    rownames() %>%
    {exprs(gse_2[[1]])[, colnames(exprs(gse_2[[1]])) %in% .]}

# 处理平台注释信息
platform1 <- Table(getGEO(unique(pData(gse_1[[1]])$"platform_id"), destdir = "./input")) %>%
    select("ID", "Gene Symbol") %>%
    filter(!grepl("///", .$"Gene Symbol")) %>%
    rename(GENE_SYMBOL = "Gene Symbol") %>%
    filter(GENE_SYMBOL != "")

platform2 <- Table(getGEO(unique(pData(gse_2[[1]])$"platform_id"), destdir = "./input")) %>%
    select("ID", "Gene Symbol") %>%
    filter(!grepl("///", .$"Gene Symbol")) %>%
    rename(GENE_SYMBOL = "Gene Symbol") %>%
    filter(GENE_SYMBOL != "")

# 合并探针数据为基因表达值
expr1 <- data1 %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(platform1, by = "ID") %>%
    filter(!is.na(GENE_SYMBOL)) %>%
    group_by(GENE_SYMBOL) %>%
    summarise(across(where(is.numeric), mean)) %>%
    column_to_rownames("GENE_SYMBOL")

expr2 <- data2 %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(platform2, by = "ID") %>%
    filter(!is.na(GENE_SYMBOL)) %>%
    group_by(GENE_SYMBOL) %>%
    summarise(across(where(is.numeric), mean)) %>%
    column_to_rownames("GENE_SYMBOL")

# 找到共同基因并过滤数据
common_genes <- intersect(rownames(expr1), rownames(expr2))
expr1 <- expr1[common_genes, ]
expr2 <- expr2[common_genes, ]

#####################################################
## Step4：数据标准化和批次效应校正
## 
## 描述：标准化数据并进行批次效应校正
#####################################################
#' 数据标准化和对数转换函数
#' @param exp_data 表达矩阵
#' @return 标准化后的矩阵
normalize_expression <- function(exp_data) {
  # 负值处理
  if (range(exp_data)[1] < 0) {
    exp_data <- exp_data + abs(min(exp_data))
  }
  
  # 计算分位数并判断是否需要对数转换
  qx <- quantile(exp_data, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE) %>% 
    as.numeric()
  
  need_log <- (qx[5] > 100) || 
    (qx[6] - qx[1] > 50) || 
    (qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  
  # 对数转换
  if (need_log) {
    exp_data <- log2(exp_data + 1)
  }
  
  # 标准化
  normalizeBetweenArrays(exp_data)
}

# 对两个数据集分别进行标准化
data1_norm <- normalize_expression(expr1)
data2_norm <- normalize_expression(expr2)

# 批次效应校正
batch <- factor(c(rep("2", ncol(data2_norm)), rep("1", ncol(data1_norm))))
adjusted_matrix <- ComBat(cbind(data2_norm, data1_norm), batch = batch)
adjusted_matrix <- normalizeBetweenArrays(adjusted_matrix)

# 保存标准化数据
save(adjusted_matrix, file = "normalized_expression_data.Rda")

# 创建批次效应校正前后的可视化图
# 1. PCA图 - 校正前
pca_before <- prcomp(t(cbind(data2_norm, data1_norm)), scale. = TRUE)
pca_before_df <- data.frame(
    PC1 = pca_before$x[,1],
    PC2 = pca_before$x[,2],
    Batch = factor(c(rep("GSE5281", ncol(data2_norm)), 
                    rep("GSE48350", ncol(data1_norm))))
)

p1 <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3) +
    stat_ellipse() +
    theme_bw() +
    labs(title = "Before Batch Correction") +
    theme(plot.title = element_text(hjust = 0.5))

# PCA图 - 校正后
pca_after <- prcomp(t(adjusted_matrix), scale. = TRUE)
pca_after_df <- data.frame(
    PC1 = pca_after$x[,1],
    PC2 = pca_after$x[,2],
    Batch = factor(c(rep("GSE5281", ncol(data2_norm)), 
                    rep("GSE48350", ncol(data1_norm))))
)

p2 <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3) +
    stat_ellipse() +
    theme_bw() +
    labs(title = "After Batch Correction") +
    theme(plot.title = element_text(hjust = 0.5))

# 2. 箱线图
# 准备数据
boxplot_data_before <- data.frame(
    Expression = c(as.vector(data2_norm), as.vector(data1_norm)),
    Batch = factor(rep(c("GSE5281", "GSE48350"), 
                      c(length(data2_norm), length(data1_norm))))
)

boxplot_data_after <- data.frame(
    Expression = as.vector(adjusted_matrix),
    Batch = factor(rep(c("GSE5281", "GSE48350"), 
                      c(ncol(data2_norm) * nrow(data2_norm), 
                        ncol(data1_norm) * nrow(data1_norm))))
)

p3 <- ggplot(boxplot_data_before, aes(x = Batch, y = Expression, fill = Batch)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(title = "Expression Distribution Before Correction",
         y = "Expression Value") +
    theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(boxplot_data_after, aes(x = Batch, y = Expression, fill = Batch)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(title = "Expression Distribution After Correction",
         y = "Expression Value") +
    theme(plot.title = element_text(hjust = 0.5))

# 3. 密度图
p5 <- ggplot(boxplot_data_before, aes(x = Expression, fill = Batch)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    labs(title = "Density Plot Before Correction",
         x = "Expression Value") +
    theme(plot.title = element_text(hjust = 0.5))

p6 <- ggplot(boxplot_data_after, aes(x = Expression, fill = Batch)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    labs(title = "Density Plot After Correction",
         x = "Expression Value") +
    theme(plot.title = element_text(hjust = 0.5))

# 组合所有图形并保存
batch_effect_plots <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
    plot_layout(heights = c(1, 1, 1))

# 保存图形
ggsave("plots/pca_plot.pdf", p1, width = 6, height = 5, dpi = 300)
ggsave("plots/pca_after_plot.pdf", p2, width = 6, height = 5, dpi = 300)
ggsave("plots/boxplot_plot.pdf", p3, width = 6, height = 5, dpi = 300)
ggsave("plots/boxplot_after_plot.pdf", p4, width = 6, height = 5, dpi = 300)
ggsave("plots/density_plot.pdf", p5, width = 6, height = 5, dpi = 300)
ggsave("plots/density_after_plot.pdf", p6, width = 6, height = 5, dpi = 300)
ggsave("plots/batch_effect_plots.pdf", batch_effect_plots,width = 12, height = 15, dpi = 300)
