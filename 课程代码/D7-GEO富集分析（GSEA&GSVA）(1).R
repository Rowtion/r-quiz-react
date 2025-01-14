## ====================================================================
## DAY 7. R语言GESA&GSVA功能富集分析
## 
## 功能：进行基因集富集分析(GSEA)和基因集变异分析(GSVA)
## 输入：
##   - DEGs_data.Rda: 差异分析结果
##   - normalized_expression_data.Rda: 标准化表达矩阵
##   - c2.all.v2024.1.Hs.symbols.gmt: MSigDB基因集数据
## 输出：
##   - gsva_result.Rda: GSVA分析结果
##   - gsea_result.Rda： GSEA分析结果
## ====================================================================

#####################################################
## Step1：环境准备
## - 加载必要的R包
#####################################################
# 设置工作路径
setwd("./D7-GEO富集分析（GSEA&GSVA）")

# 加载必要的包
library(clusterProfiler)  # GSEA分析
library(GSVA)            # 基因集变异分析
library(GSEABase)        # 基因集处理
library(GseaVis)         # GSEA结果可视化
library(pheatmap)        # 热图绘制
library(ggplot2)         # 绘图基础包
library(patchwork)       # 图形组合
library(tidyr)          # 数据整理
library(dplyr)          # 数据处理
library(tidyverse)      # 数据处理

#####################################################
## Step2：数据准备
## - 读取差异分析和表达矩阵数据
## - 准备GSEA分析所需的排序基因列表
## - 读取MSigDB基因集数据
#####################################################
# 读取数据
load("DEGs_data.Rda")  # 差异分析结果
load("normalized_expression_data.Rda")  # 标准化的表达矩阵

# 准备GSEA所需的基因列表
gene_list <- DEGs$logFC
names(gene_list) <- rownames(DEGs)
gene_list <- sort(gene_list, decreasing = TRUE)

# 读取基因集
gmt <- read.gmt("c2.all.v2024.1.Hs.symbols.gmt")

#####################################################
## Step3：GSEA分析
## - 运行GSEA分析：
##   * 最小基因集大小: 10
##   * 最大基因集大小: 500
##   * p值阈值: 0.05
## - 绘制GSEA富集图：
##   * 展示第一个显著通路
##   * 添加统计信息
#####################################################
# 运行GSEA分析
gsea_result <- GSEA(
  geneList = gene_list,
  TERM2GENE = gmt,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# 绘制GSEA富集图
pdf("gsea_enrichment_plots.pdf", width = 12, height = 8)

# 获取第一个通路的信息
pathway_stats <- gsea_result@result[1, ]
pathway_name <- pathway_stats$Description

# 使用gseaNb绘图
gseaNb(object = gsea_result, 
       geneSetID = pathway_name,       
       addPval = T,
       pCol = "steelblue",
       pvalX = 0.65, # 位置
       pvalY = 0.7,
       pHjust = 0, # 对齐方式
       nesDigit = 4, # 小数点位数
       pDigit = 4)

dev.off()

# 保存GSEA分析结果
save(gsea_result, file = "gsea_result.Rda")

#####################################################
## Step4：GSVA分析
## - 准备基因集和表达矩阵
## - 设置GSVA参数
## - 运行GSVA分析
## - 选择前20个基因集进行可视化
#####################################################
# GSVA分析
# 准备基因集和表达矩阵
geneSets <- split(gmt$gene, gmt$term)
expr_matrix <- as.matrix(adjusted_matrix)

#准备GSVA分析对象
gsvaParamObj <- gsvaParam(adjusted_matrix, geneSets, kcdf = "Gaussian")

# 直接运行GSVA分析
gsva_result <- gsva(gsvaParamObj)

# 选择前30个基因集进行热图展示
top_pathways <- head(rownames(gsva_result), 20)
gsva_matrix <- gsva_result[top_pathways, ]

# 1. 绘制热图
p1 <- pheatmap(gsva_matrix,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "GSVA Score Heatmap",
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 8,
         scale = "row",
         show_colnames = TRUE,
         show_rownames = TRUE,
         border_color = NA,
         angle_col = 45)

# 保存热图
ggsave("gsva_heatmap.pdf", p1, width = 8, height = 5)

# 保存GSVA分析结果
save(gsva_result, file = "gsva_result.Rda")
