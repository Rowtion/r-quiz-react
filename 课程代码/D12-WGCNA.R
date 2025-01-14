## ====================================================================
## WGCNA基因共表达网络分析
##
## 描述：
##   使用WGCNA包构建基因共表达网络，识别共表达模块，并进行模块-性状关联分析
##
## 输入：
##   - normalized_expression_data.Rda: 标准化后的表达矩阵
##   - group_GEO.Rda: 样本分组信息
##
## 输出：
##   - output/network_construction.RData: 网络构建结果
##   - output/merged_modules.RData: 合并后的模块信息
##   - output/module_trait_correlation.csv: 模块-性状相关性
##   - output/module_trait_pvalue.csv: 模块-性状p值
##   - output/module_genes.csv: 关键模块基因列表
##   - 多个可视化PDF文件
##
## 注意事项：
##   - 需要较大的内存空间
##   - 建议在多核CPU上运行
## ====================================================================

# 清空工作环境
rm(list = ls())

#####################################################
## Step1：环境初始化
##
## 描述：初始化R环境，加载必要的包，设置工作目录
#####################################################

# 加载必要的R包
if (!require("WGCNA")) install.packages("WGCNA")
if (!require("tidyverse")) install.packages("tidyverse")
library(WGCNA)
library(tidyverse)

# 设置工作目录和输出路径
setwd("./D12")
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# 加载输入数据
load("input/normalized_expression_data.Rda")
load("input/group_GEO.Rda")

#####################################################
## Step2：基因筛选
##
## 描述：选择高变异基因用于后续分析
## 方法：选择表达变异系数在前20%的基因
#####################################################

# 基因筛选 - 选择高变异基因
print("正在进行基因筛选...")
vars_all_gene <- apply(adjusted_matrix, 1, var)
gene0 <- names(vars_all_gene[vars_all_gene >= quantile(vars_all_gene, probs = 0.8)])
adjusted_matrix <- adjusted_matrix[gene0, ] %>% t() %>% as.data.frame()

#####################################################
## Step3：数据预处理
##
## 描述：对表达数据进行质量控制和标准化处理
## 步骤：
##   1. 检查基因和样本质量
##   2. 样本聚类分析
##   3. 数据分布检查
#####################################################

# 1. 数据预处理
print("正在进行数据预处理...")
enableWGCNAThreads()  # 启用多线程计算

# 1.1 检查数据质量
gsg <- goodSamplesGenes(adjusted_matrix, verbose = 3)
if (!gsg$allOK) {
  # 记录并保存被移除的基因和样本
  removed_genes <- names(adjusted_matrix)[!gsg$goodGenes]
  removed_samples <- rownames(adjusted_matrix)[!gsg$goodSamples]
  
  if (length(removed_genes) > 0) {
    print(paste("移除基因:", length(removed_genes), "个"))
    write.csv(data.frame(Removed_Genes=removed_genes),
              file.path(output_dir, "removed_genes.csv"))
  }
  
  if (length(removed_samples) > 0) {
    print(paste("移除样本:", length(removed_samples), "个"))
    write.csv(data.frame(Removed_Samples=removed_samples),
              file.path(output_dir, "removed_samples.csv"))
  }
  
  adjusted_matrix <- adjusted_matrix[gsg$goodSamples, gsg$goodGenes]
}

# 1.2 样本聚类检查
pdf(file.path(output_dir, "sample_clustering.pdf"), width=10, height=6)
sampleTree <- hclust(dist(adjusted_matrix), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")
dev.off()

# 1.3 数据分布检查
pdf(file.path(output_dir, "data_distribution.pdf"), width=10, height=6)
par(mfrow=c(1,2))
hist(adjusted_matrix, breaks=100, main="Expression Distribution", xlab="Expression Value")
boxplot(adjusted_matrix, main="Expression Boxplot")
dev.off()

#####################################################
## Step4：网络构建
##
## 描述：构建基因共表达网络并识别模块
## 步骤：
##   1. 选择软阈值
##   2. 构建拓扑重叠矩阵(TOM)
##   3. 识别共表达模块
##   4. 合并相似模块
#####################################################

# 2. 网络构建
print("正在构建基因共表达网络...")

# 2.1 选择软阈值
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(adjusted_matrix, powerVector = powers, RsquaredCut = 0.9, verbose = 5)

# 2.2 绘制scale-free拓扑拟合指数图
pdf(file.path(output_dir, "scale_free_topology.pdf"), width=10, height=5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.90, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity",
     type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

# 2.3 选择软阈值
if (is.na(sft$powerEstimate)) {
  print("无法自动选择软阈值，请查看output/scale_free_topology.pdf并手动选择合适软阈值")
  print("建议选择使R^2接近0.9的最小软阈值")
  power <- readline(prompt = "请输入手动选择的软阈值: ")
  power <- as.numeric(power)
} else {
  power <- sft$powerEstimate
  print(paste("选择的软阈值功率为:", power))
}

# 2.4 构建网络
print("正在构建基因共表达网络...")
net <- blockwiseModules(adjusted_matrix, 
                       power = power,
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       pamRespectsDendrogram = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = file.path(output_dir, "TOM"),
                       verbose = 3)

# 2.5 模块合并
print("正在合并相似模块...")
me_diss_threshold <- 0.25
merged_modules <- mergeCloseModules(adjusted_matrix, 
                                   net$colors,
                                   cutHeight = me_diss_threshold,
                                   verbose = 3)

moduleColors <- merged_modules$colors
MEs <- merged_modules$newMEs

# 保存合并后的模块信息
save(moduleColors, MEs, file = file.path(output_dir, "merged_modules.RData"))

#####################################################
## Step5：结果保存与可视化
##
## 描述：保存分析结果并生成可视化图表
## 步骤：
##   1. 保存模块信息
##   2. 模块-性状关联分析
##   3. 生成可视化图表
##   4. 保存会话信息
#####################################################

# 3. 结果保存与可视化
print("正在保存结果和生成可视化...")

# 3.1 保存模块信息
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

save(moduleLabels, moduleColors, MEs, geneTree, 
     file = file.path(output_dir, "network_construction.RData"))

# 构建分组信息
combined_groups$group <- as.factor(combined_groups$group)
dat_trait <- model.matrix( ~ 0 + combined_groups$group) %>% as.data.frame()
colnames(dat_trait) <- levels(combined_groups$group)
rownames(dat_trait) <- combined_groups$ID

# 3.2 模块与性状关联分析
mat_module_trait_cor <- WGCNA::cor(MEs, dat_trait, method = 'spearman', use = 'p')
moduleTraitPvalue <- corPvalueStudent(mat_module_trait_cor, nrow(adjusted_matrix))

# 保存关联分析结果
write.csv(mat_module_trait_cor, file.path(output_dir, "module_trait_correlation.csv"))
write.csv(moduleTraitPvalue, file.path(output_dir, "module_trait_pvalue.csv"))

# 3.3 可视化
# 模块树状图
pdf(file.path(output_dir, "dendrogram.pdf"), width=12, height=9)
plotDendroAndColors(geneTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 相关性和p值标注
mat_text <- paste0(
  'Correlation = ',
  round(mat_module_trait_cor, 2),
  '\npvalue = ',
  signif(moduleTraitPvalue, 3)
)
dim(mat_text) <- dim(mat_module_trait_cor)

# 模块-性状关系热图
pdf(file.path(output_dir, "module_trait_heatmap.pdf"), width=8, height=8)
par(mar = c(5, 10, 3, 5))
labeledHeatmap(
  Matrix = mat_module_trait_cor,
  xLabels = colnames(dat_trait),
  xLabelsAngle = 0,
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  colors = colorRampPalette(colors = c('#4DBBD5', 'white', '#E64B35'))(50),
  textMatrix = mat_text,
  setStdMargins = F,
  zlim = c(-1, 1),
  main = 'Module-Trait Relationship'
)

dev.off()

# 关闭多线程计算
disableWGCNAThreads()

# 3.4 保存会话信息
save.image(file.path(output_dir, "WGCNA_session.RData"))

#####################################################
## Step6：模块基因分析
##
## 描述：提取关键模块中的基因并进行后续分析
## 步骤：
##   1. 选择感兴趣的关键模块
##   2. 提取模块中的基因
##   3. 保存模块基因列表
#####################################################

# 4. 模块基因分析
select_modules <- c('ME3', 'ME9', 'ME12')  # 可根据实际情况调整
module_genes <- gene0[moduleLabels %in% select_modules]
write.csv(data.frame(Gene=module_genes), 
          file.path(output_dir, "module_genes.csv"), 
          row.names = FALSE)
