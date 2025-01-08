## ====================================================================
## DAY 8. R语言相关性分析
## 
## 功能：分析和可视化基因间的相关性关系
## 输入：
##   - normalized_expression_data.Rda: 标准化的表达数据
##   - Common_genes.csv: 交集基因列表
## 输出：
##   - 可视化结果
## ====================================================================

#####################################################
## Step1：环境准备
## - 加载必要的R包：
#####################################################
#设置工作路径
setwd("./D8-GEO相关性分析")

# 加载必要的R包
library(circlize)        # 和弦图绘制
library(ComplexHeatmap)  # 复杂热图绘制
library(pheatmap)        # 基础热图绘制
library(corrplot)        # 相关性图绘制

#####################################################
## Step2：数据导入与预处理
## - 读取标准化的表达数据
## - 读取交集基因列表
## - 提取目标基因的表达数据
## - 计算基因间相关性矩阵
#####################################################
# 读取标准化的表达数据
load("normalized_expression_data.Rda")

# 读取交集基因列表
common_genes <- read.csv("Common_genes.csv", stringsAsFactors = FALSE)

# 提取交集基因的表达数据
gene_expr <- adjusted_matrix[common_genes$Gene, ]

# 计算相关性矩阵
cor_matrix <- cor(t(gene_expr))

#####################################################
## Step3：可视化
## Step3.1：和弦图
## - 设置相关性颜色范围（-1到1）
## - 绘制基因间相关性和弦图：
##   * 使用彩虹色区分基因
##   * 相关性强度用颜色深浅表示
##   * 添加基因名称标签
##   * 添加相关性图例
#####################################################
pdf("gene_correlation_chord.pdf", width = 8, height = 8)
# 设置相关性的颜色范围（-1到1）
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

# 绘制和弦图
chordDiagram(
    cor_matrix,
    grid.col = rainbow(nrow(cor_matrix)),  # 设置每个基因的颜色
    col = col_fun,                         # 设置相关性的颜色
    symmetric = TRUE,                      # 对称矩阵
    directional = FALSE,                   # 无方向性
    annotationTrack = "grid",             # 添加网格线
    preAllocateTracks = list(
        track.height = mm_h(4),
        track.margin = c(mm_h(2), 0)
    )
)

# 添加基因名称标签
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# 在右下角添加简单的图例
lgd = Legend(col_fun = col_fun, 
            title = "Correlation",
            at = c(-1, 0, 1),
            direction = "horizontal",
            legend_width = unit(4, "cm"))  # 增加图例的长度
draw(lgd, x = unit(0.95, "npc"), y = unit(0.05, "npc"), just = c("right", "bottom"))
dev.off()

#####################################################
## Step3.2: 共表达热图
## - 对表达数据进行标准化
## - 绘制基因表达热图：
##   * 使用层次聚类
##   * 蓝白红配色方案
##   * 隐藏样本名称
#####################################################
pdf("gene_expression_heatmap.pdf", width = 10, height = 8)
# 对表达数据进行标准化
scaled_expr <- t(scale(t(gene_expr)))
# 绘制热图
pheatmap(scaled_expr,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Gene Expression Heatmap",
         fontsize_row = 10,
         border_color = NA)
dev.off()

#####################################################
## Step3.3：相关性热图
## - 使用corrplot绘制相关性热图：
##   * 上三角矩阵展示
##   * 层次聚类排序
##   * 显示相关系数
##   * 45度角旋转标签
## - 打印相关性矩阵（保留两位小数）
#####################################################
pdf("correlation_heatmap.pdf", width = 8, height = 8)
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200))
dev.off()

# 打印相关性矩阵（保留两位小数）
print(round(cor_matrix, 2))
