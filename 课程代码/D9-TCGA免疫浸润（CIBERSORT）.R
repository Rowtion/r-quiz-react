## ====================================================================
## DAY 9. TCGA免疫浸润分析（CIBERSORT）
## 
## 功能：使用CIBERSORT算法对TCGA数据进行免疫细胞浸润分析和可视化
## 输入：
##   - TCGA-BRCA_tpm.Rda: TCGA表达矩阵数据
##   - TCGA_BRCA_cli_clean.Rda: 清洗后的临床信息，包含OLR1高低表达组分组信息
## 输出：
##   - immune_infiltration_results.Rda: 免疫浸润分析结果
## ====================================================================

#####################################################
## 环境准备
## - 清理工作环境
## - 加载必要的R包
## - 导入数据
#####################################################

# 清空环境
rm(list = ls())

#设置工作路径
setwd("./D9-TCGA免疫浸润（CIBERSORT）")

# 加载必要的包
library(CIBERSORT)    # 免疫浸润分析
library(tidyverse)    # 数据处理和可视化
library(ggpubr)      # 统计检验和可视化
library(reshape2)     # 数据重构
library(org.Hs.eg.db) # 基因ID转换
library(aplot)        # 图表插入

#####################################################
## Step1：数据预处理
## - 导入原始数据
## - 筛选肿瘤样本
## - 统一样本标识
## - 基因ID转换
#####################################################

## ------------ 1. 数据导入 ------------
load("TCGA-BRCA_tpm.Rda")  # 表达矩阵数据
load("TCGA_BRCA_cli_clean.Rda")  # 清洗后的临床信息，包含OLR1高低分组

## ------------ 2. 样本筛选-----------
# 提取OLR1高低分组信息
cluster <- tumor_data_final %>%
  dplyr::select(sample, group) %>%
  arrange(group)

#筛选疾病组样本
tpm_data <- tpm_data[, cluster$sample] %>% dplyr::select(tumor_data_final$sample)

## ------------ 3. 基因ID转换 ------------
# 获取Ensembl ID
ensembl_ids <- rownames(tpm_data)

# 转换为Gene Symbol
gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
)

# 移除未转换成功的基因
valid_genes <- !is.na(gene_symbols)
tpm_data <- tpm_data[valid_genes,]
gene_symbols <- gene_symbols[valid_genes]

## ------------ 4. 处理重复基因 ------------
# 将表达矩阵转换为数据框并添加基因信息
tpm_data_df <- as.data.frame(tpm_data)
tpm_data_df$symbol <- gene_symbols
tpm_data_df$mean_expr <- rowMeans(tpm_data)

# 对重复的基因名保留表达量最高的
tpm_data_df <- tpm_data_df %>%
    group_by(symbol) %>%
    slice_max(order_by = mean_expr, n = 1) %>%
    ungroup()

# 重构表达矩阵
tpm_data <- as.matrix(tpm_data_df[, !colnames(tpm_data_df) %in% c("symbol", "mean_expr")])
rownames(tpm_data) <- tpm_data_df$symbol

#####################################################
## Step2：CIBERSORT免疫浸润分析
## - 运行CIBERSORT算法
## - 提取免疫细胞浸润结果
#####################################################

## ------------ 1. 运行CIBERSORT ------------
dat_res_cibersort <- cibersort(
    LM22,           # 免疫细胞特征基因集
    tpm_data,       # 表达矩阵
    perm = 100,     # 置换检验次数
    QN = TRUE       # 分位数标准化
) %>% as.data.frame()

# 提取22种免疫细胞的浸润丰度结果
dat_res_cibersort <- dat_res_cibersort[, 1:22]

# 去除样本中完全不存在的免疫细胞类型
dat_res_cibersort <- dat_res_cibersort[, colSums(dat_res_cibersort) > 0]

#####################################################
## Step3：可视化分析
## - 准备数据
## - 绘制堆叠柱状图
#####################################################
## ------------ 1. 准备数据 ------------
# 构建堆叠柱状图数据
dat_cibersort_stack <- dat_res_cibersort
dat_cibersort_stack <- dat_cibersort_stack[rownames(dat_cibersort_stack) %in% cluster$sample, ]
dat_cibersort_stack$sample <- rownames(dat_cibersort_stack)
dat_cibersort_stack <- reshape2::melt(dat_cibersort_stack)
colnames(dat_cibersort_stack) <- c('sample', 'cell', 'score')

# 添加分组信息并按组排序
dat_cibersort_stack <- dat_cibersort_stack %>%
  left_join(cluster, by = "sample") %>%
  mutate(sample = factor(sample, levels = cluster$sample))

## ------------ 2. 绘制堆叠柱状图 ------------
# 设置免疫细胞类型的颜色
immune_colors <- c(
    '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2',
    '#DC0000', '#7E6148', '#B09C85', '#0072B5', '#BC3C29',
    '#E18727', '#20854E', '#7876B1', '#6F99AD', '#FFDC91',
    '#EE4C97', '#3B4992', '#EE0000', '#008B45', '#631879',
    '#008280', '#BB0021'
)

# 创建主图
p <- (
    ggplot(dat_cibersort_stack, aes(sample, score, fill = cell)) +
    geom_col(position = 'fill') +
    theme_classic() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    labs(
        x = element_blank(),
        y = 'Relative Percentage',
        fill = 'Immune Cell'
    ) +
    scale_fill_manual(values = immune_colors)
) %>%
    # 添加分组注释条
    aplot::insert_bottom(
        ggplot(cluster, aes(sample, 'Group', fill = group)) +
        geom_tile() +
        theme_classic() +
        theme(
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank()
        ) +
        labs(x = element_blank(),
             y = element_blank(),
             fill = 'Group') +
        scale_fill_manual(values = c('#4DBBD5', '#E64B35')),
        height = 0.03
    )

## ------------ 3. 保存图片 ------------
ggsave(
    filename = 'output/CIBERSORT_stackplot.pdf',
    plot = p,
    width = 10,
    height = 6
)

#保存数据
save(dat_res_cibersort, file = "immune_infiltration_results.Rda")
