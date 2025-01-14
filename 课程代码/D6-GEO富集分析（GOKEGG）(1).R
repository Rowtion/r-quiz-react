## ====================================================================
## DAY 6. R语言GOKEGG功能富集分析
## 
## 功能：对差异表达基因进行GO和KEGG通路富集分析
## 输入：
##   - DEGs.csv: 差异表达基因列表
##   - PANoptosis_genes.csv: PANoptosis相关基因列表
## 输出：
##   - Common_genes.csv: 交集基因列表
## ====================================================================

#####################################################
## Step1：环境初始化
## - 清理工作环境
## - 设置工作目录
## - 加载R包
#####################################################
# 清空工作环境
rm(list = ls())

# 设置工作目录
setwd("./D6-GEO富集分析（GOKEGG）")

# 加载必要的R包
library(clusterProfiler) # GO和KEGG富集分析
library(org.Hs.eg.db)    # 基因ID转换
library(VennDiagram)     # 韦恩图绘制
library(grid)            # 韦恩图布局
library(gridExtra)       # 手动给图片添加文字注释

#####################################################
## Step2：差异基因与表型基因取交集
## - 读取差异基因和PANoptosis基因数据
## - 计算基因集交集
## - 生成基因统计信息
## - 创建韦恩图可视化：
##   * 设置图形参数
##   * 绘制基本韦恩图
##   * 添加数值标注
## - 转换基因ID（Symbol到Entrez）
#####################################################
# 读取差异表达基因和PANoptosis基因数据
data <- read.csv("DEGs.csv", row.names = 1)
PANoptosis_genes <- read.csv("PANoptosis_genes.csv", stringsAsFactors = FALSE)

# 提取基因列表并找到交集
deg_genes <- rownames(data)
PANoptosis_genes <- PANoptosis_genes$gene
common_genes <- intersect(deg_genes, PANoptosis_genes)

# 打印基因数量
cat("\n基因数量统计：\n")
cat("DEGs数量:", length(deg_genes), "\n")
cat("PANoptosis基因数量:", length(PANoptosis_genes), "\n")
cat("交集基因数量:", length(common_genes), "\n")

# 保存交集基因到CSV文件
write.csv(data.frame(Gene = common_genes), "Common_genes.csv", row.names = FALSE)

# 创建韦恩图
# 设置颜色
venn_colors <- c("#87CEEB", "#FFB6C1")  # 浅蓝色和浅粉色

# 创建韦恩图
venn_plot <- draw.pairwise.venn(
  area1 = 100,
  area2 = 100,
  cross.area = 30,
  category = c("DEGs", "PANoptosis genes"),
  lwd = 0,  # 移除边框
  lty = "solid",
  col = "white",  # 边框颜色设为白色
  fill = adjustcolor(venn_colors, alpha = 0.5),
  cex = 0,
  cat.cex = 1.2,
  cat.pos = c(0, 0),
  cat.dist = 0.05,
  cat.col = "black",
  scaled = FALSE,  # 设置为FALSE以强制使用标准韦恩图
  euler.d = FALSE,  # 设置为FALSE以强制使用标准韦恩图
  sep.dist = 0.04,
  rotation.degree = 0,
  ind = FALSE,
  fontfamily = "sans"
)

# 创建画布
pdf("plot/Venn_diagram.pdf", width = 6, height = 6)

# 绘制韦恩图
grid.arrange(gTree(children = venn_plot))

# 添加实际数值
grid.text(length(deg_genes), x=0.2, y=0.5, gp=gpar(fontsize=20))
grid.text(length(PANoptosis_genes), x=0.8, y=0.5, gp=gpar(fontsize=20))
grid.text(length(common_genes), x=0.5, y=0.5, gp=gpar(fontsize=20))

# 关闭设备
dev.off()

# 将基因符号转换为ENTREZ ID用于富集分析
gene_ids <- bitr(common_genes, 
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

#####################################################
## Step3：GO富集分析
## - 生物学过程(BP)富集分析
##   * p值阈值: 0.05
##   * q值阈值: 0.25
## - 分子功能(MF)富集分析
## - 细胞组分(CC)富集分析
#####################################################
# 生物学过程(BP)富集分析
go_bp <- enrichGO(gene = gene_ids$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID",
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.25)

# 分子功能(MF)富集分析
go_mf <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "MF",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.25)

# 细胞组分(CC)富集分析
go_cc <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "CC",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.25)

#####################################################
## Step4：KEGG通路富集分析
## - 使用KEGG数据库进行通路富集
## - 设置参数：
##   * 物种: 人类(hsa)
##   * p值阈值: 0.05
##   * q值阈值: 0.25
#####################################################
# 进行KEGG通路富集分析
kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.25)

#####################################################
## Step5：富集分析结果可视化
## - 数据预处理函数：
##   * 转换基因比例
##   * 截断过长描述
##   * 选择前10个显著通路
## - 处理各类富集结果
## - 创建组合气泡图：
##   * 展示通路富集程度
##   * 显示基因数量
##   * 按类别分面展示
#####################################################
# 处理各个富集分析结果数据
process_for_plot <- function(enrichment_result, category) {
  if(is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    return(NULL)
  }
  
  result_df <- as.data.frame(enrichment_result)
  
  # 确保所需的列都存在
  required_cols <- c("Description", "GeneRatio", "p.adjust", "Count", "geneID")
  if(!all(required_cols %in% colnames(result_df))) {
    return(NULL)
  }
  
  # 转换GeneRatio为数值
  result_df$GeneRatio <- sapply(strsplit(result_df$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
  # 只保留需要的列
  result_df <- result_df[, required_cols]
  
  # 按照p.adjust值排序并只取前10个
  result_df <- result_df[order(result_df$p.adjust), ]
  if(nrow(result_df) > 10) {
    result_df <- result_df[1:10, ]
  }
  
  # 截断过长的Description
  result_df$Description <- sapply(result_df$Description, function(x) {
    if(nchar(x) > 50) {
      paste0(substr(x, 1, 47), "...")
    } else {
      x
    }
  })
  
  # 添加类别标签
  result_df$Category <- category
  return(result_df)
}

# 处理每个类别的结果
bp_plot <- process_for_plot(go_bp, "BP")
mf_plot <- process_for_plot(go_mf, "MF")
cc_plot <- process_for_plot(go_cc, "CC")
kegg_plot <- process_for_plot(kegg, "KEGG")

# 合并有效的结果
plot_results_list <- list(bp_plot, mf_plot, cc_plot, kegg_plot)
valid_plot_results <- plot_results_list[!sapply(plot_results_list, is.null)]

if(length(valid_plot_results) > 0) {
  # 合并所有有效结果
  plot_data <- do.call(rbind, valid_plot_results)
  
  # 创建分面气泡图
  bubble_plot <- ggplot(plot_data, 
                       aes(x = GeneRatio, 
                           y = reorder(Description, GeneRatio), 
                           size = Count, 
                           color = -log10(p.adjust))) +
    geom_point(alpha = 0.7) + 
    facet_grid(Category ~ ., scales = "free_y", space = "free") +
    scale_size_continuous(name = "Count", 
                         range = c(2, 8),
                         breaks = c(50, 75, 100, 125, 150)) +
    scale_color_gradient(name = "P.adj", 
                        low = "#4393C3", 
                        high = "#D6604D",
                        breaks = c(5, 10, 15)) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 10, face = "bold", color = "white"),
      strip.background = element_rect(fill = "gray40", color = "gray40"),
      panel.spacing = unit(0.2, "cm"),
      panel.grid.major = element_line(color = "gray95"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.title = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      panel.border = element_rect(color = "gray80")
    ) +
    labs(x = "GeneRatio", 
         y = NULL) +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x),
                      breaks = seq(0, max(plot_data$GeneRatio), length.out = 5))
  
  # 保存气泡图
  ggsave("GOKEGG_bubble_facet.pdf", 
         plot = bubble_plot, 
         width = 8,  
         height = max(6, 1.5 + nrow(plot_data) * 0.25),  
         dpi = 300)
}

#####################################################
## Step6：结果统计和保存
## - 打印各类基因集和富集结果的统计信息
## - 保存富集分析结果为CSV文件
#####################################################
# 打印富集结果的统计信息
cat("GO BP富集通路数量:", nrow(as.data.frame(go_bp)), "\n")
cat("GO MF富集通路数量:", nrow(as.data.frame(go_mf)), "\n")
cat("GO CC富集通路数量:", nrow(as.data.frame(go_cc)), "\n")
cat("KEGG富集通路数量:", nrow(as.data.frame(kegg)), "\n")

# 保存富集分析结果为CSV文件
write.csv(as.data.frame(go_bp), "GO_BP_results.csv")
write.csv(as.data.frame(go_mf), "GO_MF_results.csv")
write.csv(as.data.frame(go_cc), "GO_CC_results.csv")
write.csv(as.data.frame(kegg), "KEGG_results.csv")