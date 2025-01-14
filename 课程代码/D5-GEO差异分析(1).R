## ====================================================================
## DAY 5. TCGA/GEO差异分析
## 
## 功能：对标准化后的GEO数据进行差异基因分析和可视化
## 输入：
##   - normalized_expression_data.Rda: 标准化后的表达矩阵
##   - raw_geo_data.Rda: 原始GEO数据
## 输出：
##   - DEGs.csv: 差异基因列表
##   - group_GEO.Rda: 样本分组情况
##   - DEGs_data.Rda: 差异分析结果
## ====================================================================

#####################################################
## 环境准备
## - 清理工作环境
## - 加载必要的R包
## - 导入预处理数据
## - 设置通用绘图主题
#####################################################
rm(list = ls())
setwd("./D5-GEO差异分析")

# 加载必要的R包
library(GEOquery)    # GEO数据处理
library(limma)       # 差异分析
library(tidyverse)   # 数据处理和基础可视化
library(ggrepel)     # 标签标注
library(ggforce)     # 椭圆添加
library(ComplexHeatmap) # 热图绘制
library(circlize)    # ComplexHeatmap依赖包
library(ggpubr)      # ggplot2拓展
library(patchwork)   # 拼图

# 加载预处理后的数据
load("./normalized_expression_data.Rda")
load("./raw_geo_data.Rda")

# 通用可视化主题设置
plot_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.border = element_rect(linewidth = 0.2),
    panel.grid = element_blank()
  )


#####################################################
## Step1：样本分组信息处理
## - 处理GSE122063样本分组（海马体样本）
## - 处理GSE5281样本分组（海马体样本）
## - 合并两个数据集的分组信息
#####################################################
# GSE122063样本分组
sample_group1 <- pData(gse_1[[1]]) %>%
  filter(grepl("hippocampus", title, ignore.case = TRUE)) %>% 
  transmute(
    ID = rownames(.),
    group = ifelse(grepl("indiv", title, ignore.case = TRUE), 'Normal', 'AD')
  )

# GSE5281样本分组
sample_group2 <- pData(gse_2[[1]]) %>% 
  filter(grepl("HIP", title)) %>% 
  transmute(
    ID = rownames(.),
    group = if_else(str_detect(title, "affected"), "AD", "Normal")
  )

# 合并分组信息
combined_groups <- rbind(sample_group1, sample_group2) %>% 
  arrange(group)

#####################################################
## Step2：表达矩阵预处理
## - 过滤低表达基因
## - 标准：在75%以上样本中表达的基因被保留
#####################################################
# 过滤低表达基因（在75%样本中表达）
filtered_exp <- adjusted_matrix[, combined_groups$ID] %>% 
  .[apply(., 1, function(x) sum(x > 0) > ncol(.) * 0.75), ]

#####################################################
## Step3：降维分析
## - 数据准备和转置
## - 主成分分析(PCA)
## - 可视化PCA结果
## - 保存PCA图
#####################################################
# 准备数据
expr_matrix <- t(filtered_exp)
rownames(expr_matrix) <- colnames(filtered_exp)

# PCA分析
pca_result <- prcomp(expr_matrix, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[,1:2])
pca_df$Group <- combined_groups$group
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)

p4 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_mark_ellipse(aes(fill = Group), alpha = 0.1, linetype = 2, show.legend = FALSE) +
  scale_color_manual(values = c("AD" = "#d95f02", "Normal" = "#1b9e77")) +
  scale_fill_manual(values = c("AD" = "#d95f02", "Normal" = "#1b9e77")) +
  plot_theme +
  theme(legend.position = "right") +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    title = "PCA Plot",
    tag = "A"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = 0.3))

# 保存单图
ggsave("plot/pca_plot.pdf", p4, width = 6, height = 5, dpi = 300)

#####################################################
## Step4：差异表达分析
## - 构建实验设计矩阵
## - 使用limma进行差异分析
## - 设置阈值筛选差异基因：
##   * |logFC| > 0.5
##   * adj.P.Val < 0.05
## - 标记上调/下调基因
#####################################################
# 定义一个函数来进行差异分析，明确指定实验组和对照组
run_differential_analysis <- function(exp_matrix, group_info, 
                                    experiment_group = "AD", 
                                    control_group = "Normal",
                                    lfc_threshold = 0.5,
                                    p_threshold = 0.05) {
  # 确保分组信息正确
  if (!all(c(experiment_group, control_group) %in% unique(group_info$group))) {
    stop("指定的实验组或对照组名称不存在于数据中！")
  }
  
  # 创建因子，明确指定水平顺序，实验组在前，对照组在后
  group_factor <- factor(group_info$group, 
                        levels = c(experiment_group, control_group))
  
  # 构建设计矩阵
  design_matrix <- model.matrix(~ 0 + group_factor)
  colnames(design_matrix) <- levels(group_factor)
  rownames(design_matrix) <- colnames(exp_matrix)
  
  # 构建对比矩阵：实验组 vs 对照组
  contrast_matrix <- makeContrasts(
    contrasts = paste0(experiment_group, "-", control_group),
    levels = design_matrix
  )
  
  # 差异分析
  fit <- lmFit(exp_matrix, design_matrix)
  fit2 <- contrasts.fit(fit, contrast_matrix) %>% eBayes()
  
  # 提取差异基因
  DEGs <- topTable(fit2, coef = 1, n = Inf) %>% 
    na.omit() %>%
    mutate(
      change = case_when(
        abs(logFC) < lfc_threshold | adj.P.Val >= p_threshold ~ "Not",
        logFC >= lfc_threshold ~ "Up",
        logFC <= -lfc_threshold ~ "Down"
      ),
      logP = -log10(adj.P.Val)
    )
  
  # 返回结果
  return(list(
    DEGs = DEGs,
    design = design_matrix,
    contrast = contrast_matrix,
    fit = fit2
  ))
}

# 使用新函数进行差异分析
diff_results <- run_differential_analysis(
  exp_matrix = filtered_exp,
  group_info = combined_groups,
  experiment_group = "AD",    # 指定实验组
  control_group = "Normal",   # 指定对照组
  lfc_threshold = 0.5,        # 设置logFC阈值
  p_threshold = 0.05          # 设置P值阈值
)

# 获取差异分析结果
DEGs <- diff_results$DEGs

# 获取显著差异基因
sig_DEGs <- DEGs %>% filter(change != 'Not')

#####################################################
## Step5：标注差异基因
## - 筛选显著差异基因
## - 提取上下调TOP5基因
#####################################################
# 获取上下调前五的基因
top_up <- DEGs %>% 
  filter(change == "Up") %>% 
  arrange(desc(logFC)) %>% 
  head(5)

top_down <- DEGs %>% 
  filter(change == "Down") %>% 
  arrange(logFC) %>% 
  head(5)

top_DEGs <- rbind(top_up, top_down)
top_DEGs$label <- rownames(top_DEGs)

#####################################################
## Step6：差异分析结果可视化
## 6.1 统计差异基因数量
## 6.2 创建可视化图形：
##    - 火山图：展示差异基因分布
##    - 箱线图：展示重要基因表达差异
##    - 热图：展示差异基因表达模式
#####################################################
## 6.1 统计差异基因数量
deg_stats <- data.frame(
  change = c("Up", "Down"),
  n = c(
    sum(DEGs$change == "Up"),
    sum(DEGs$change == "Down")
  )
)

## 6.2 创建各个图形对象
# 准备上下调各TOP5差异基因的表达数据
top_genes_data <- filtered_exp[rownames(top_DEGs), ] %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Group = combined_groups$group) %>% 
  pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")

# 火山图
p1 <- ggplot(DEGs, aes(x = logFC, y = logP)) +
  geom_point(aes(color = change), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Down" = "#1f77b4", "Not" = "grey90", "Up" = "#e31a1c")) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(data = top_DEGs, size = 2, shape = 21, fill = NA, color = "black") +
  geom_text_repel(data = top_DEGs,
                  aes(label = label),
                  size = 3,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  force = 10,
                  max.overlaps = Inf) +
  plot_theme +
  theme(legend.position = "right") +
  labs(
    y = "-Log10 (adj.P.Val)",
    x = "Log2 (Fold Change)",
    title = "Volcano Plot",
    tag = "A"
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  scale_x_continuous(expand = expansion(mult = 0.2))

# Top5上下调差异基因表达箱线图
p2 <- ggplot(top_genes_data, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 0.5, width = 0.7) +
  scale_fill_manual(values = c("AD" = "#d95f02", "Normal" = "#1b9e77")) +
  stat_compare_means(aes(group = Group), 
                    method = "t.test",
                    label = "p.signif",
                    label.y = max(top_genes_data$Expression) + 0.5) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "right"
  ) +
  labs(
    x = "Gene",
    y = "Expression Level",
    title = "Expression of Top 10 DEGs",
    tag = "B"
  )

# 差异基因数量统计
p3 <- ggplot(deg_stats, aes(x = change, y = n, fill = change)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8,
           color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Down" = "#1f77b4", "Up" = "#e31a1c")) +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  plot_theme +
  theme(legend.position = "none") +
  labs(
    x = "Regulation",
    y = "Number of DEGs",
    title = "DEG Statistics",
    tag = "C"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

## 6.3 热图可视化
# 准备热图数据
# 1. 选择显著差异基因并标准化
sig_matrix <- filtered_exp[rownames(sig_DEGs), ]
sig_matrix_scaled <- t(scale(t(sig_matrix)))

# 2. 计算基因间相关性用于聚类
gene_cors <- cor(t(sig_matrix_scaled), method = "pearson")
gene_dist <- as.dist(1 - gene_cors)
gene_clust <- hclust(gene_dist, method = "complete")

# 3. 按照分组排序样本
sample_order <- order(combined_groups$group)
sig_matrix_scaled <- sig_matrix_scaled[, sample_order]

# 4. 准备注释数据
ha <- HeatmapAnnotation(
  Group = combined_groups$group[sample_order],
  col = list(Group = c("AD" = "#d95f02", "Normal" = "#1b9e77")),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Group = list(
      title = "Group",
      at = c("AD", "Normal"),
      labels = c("AD", "Normal")
    )
  ),
  height = unit(0.3, "cm")
)

# 5. 设置颜色
col_fun <- colorRamp2(
  breaks = c(min(sig_matrix_scaled), 0, max(sig_matrix_scaled)),
  colors = c("#1f77b4", "white", "#e31a1c")
)

# 6. 绘制热图
set.seed(2025)  # 保证结果可重复
heatmap_plot <- Heatmap(
  sig_matrix_scaled,
  name = "Expression",  # 颜色图例标题
  
  # 聚类设置
  cluster_rows = gene_clust,
  cluster_columns = FALSE,
  
  # 注释设置
  top_annotation = ha,
  
  # 颜色设置
  col = col_fun,
  
  # 显示设置
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  # 图例设置
  heatmap_legend_param = list(
    title = "Z-score",
    legend_height = unit(4, "cm"),
    grid_width = unit(0.5, "cm")
  ),
  
  # 其他视觉参数
  border = FALSE,
  
  # 尺寸设置
  width = unit(8, "cm"),
  height = unit(12, "cm")
)

# 7. 保存热图
pdf("plot/heatmap.pdf", width = 8, height = 12)
draw(heatmap_plot, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

## 6.4 保存单图
ggsave("plot/volcano_plot.pdf", p1, width = 6, height = 5, dpi = 300)
ggsave("plot/expression_boxplot.pdf", p2, width = 7, height = 5, dpi = 300)
ggsave("plot/deg_stats.pdf", p3, width = 5, height = 5, dpi = 300)

## 6.5 拼图
deg_plot <- p1 + p2 + p3 +
  plot_layout(guides = "keep")

ggsave("plot/deg_analysis.pdf", deg_plot, width = 17, height = 5, dpi = 300)

# 保存差异基因结果
write.csv(sig_DEGs, file = 'DEGs.csv')
save(combined_groups, file = "group_GEO.Rda")
save(DEGs,file = 'DEGs_data.Rda')
