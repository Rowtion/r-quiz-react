export const quizData = {
    'week1': {
        title: "R代码讲席营第一周周末测验",
        description: "本测验为包含D1-D6课程内容的综合测试。",
        questions: [
            {
                category: 'R基础与包管理',
                difficulty: 'easy',
                question: '在R中，以下哪个包集合最适合用于数据处理和可视化？',
                code: `# 示例代码
library(tidyverse)
library(ggplot2)
library(dplyr)`,
                options: [
                    'base R + stats',
                    'tidyverse + ggplot2',
                    'data.table + lattice',
                    'reshape2 + plotly'
                ],
                correctAnswer: 1,
                explanation: 'tidyverse是一个包集合，它包含了数据处理(dplyr)、可视化(ggplot2)等多个强大的包，是现代R编程的标准工具集。'
            },
            {
                category: 'TCGA数据分析',
                difficulty: 'medium',
                question: '在处理TCGA数据时，将counts数据转换为TPM的正确步骤是什么？',
                code: `# 转换步骤示例
counts2fpkm <- function(counts, gene_length) {
    N <- sum(counts)
    return(counts * 1e9 / gene_length / N)
}

fpkm2tpm <- function(fpkm) {
    return(fpkm * 1e6 / sum(fpkm))
}`,
                options: [
                    '直接将counts除以总reads数',
                    '先转换为RPKM，再标准化',
                    '先转换为FPKM，再转换为TPM',
                    '直接将counts乘以一个常数'
                ],
                correctAnswer: 2,
                explanation: '正确的转换流程是：先考虑基因长度和测序深度转换为FPKM，然后将FPKM标准化为TPM。这样可以消除基因长度和测序深度的影响。'
            },
            {
                category: '数据可视化',
                difficulty: 'medium',
                question: '使用ggplot2创建复杂可视化图形时，下列哪个语句是正确的？',
                code: `# 示例代码
ggplot(climate_data, 
       aes(x = temperature, 
           y = rainfall, 
           size = year, 
           color = city)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis(discrete = TRUE)`,
                options: [
                    'geom_point()必须在aes()之前',
                    'color必须在geom_point()中设置',
                    '可以用+号添加多个图层',
                    'scale_color_viridis()不能用于离散变量'
                ],
                correctAnswer: 2,
                explanation: 'ggplot2基于图形语法，使用+号逐层添加图形元素是其核心特征。每个图层都可以独立设置属性，并且按照添加顺序绘制。'
            },
            {
                category: 'GEO数据分析',
                difficulty: 'hard',
                question: '在处理来自不同GEO数据集的数据时，如何正确处理批次效应？',
                code: `# 批次效应处理示例
library(sva)
batch <- c(rep(1, ncol(data1)), rep(2, ncol(data2)))
adjusted_matrix <- ComBat(
    dat = cbind(data1, data2),
    batch = batch
)`,
                options: [
                    '直接合并数据集即可',
                    '使用ComBat进行批次效应校正',
                    '分别对每个数据集进行标准化',
                    '删除批次效应明显的样本'
                ],
                correctAnswer: 1,
                explanation: 'ComBat是一种广泛使用的批次效应校正方法，它可以有效去除非生物学变异，同时保留生物学信号。在合并多个数据集时，批次效应校正是必要的步骤。'
            },
            {
                category: '数据预处理',
                difficulty: 'hard',
                question: '在进行差异表达分析之前，需要进行哪些必要的数据预处理步骤？',
                code: `# 数据预处理示例
# 1. 数据标准化
normalized_data <- normalizeBetweenArrays(raw_data)

# 2. 质量控制
boxplot(normalized_data)
hist(normalized_data)

# 3. 批次效应检查
pca_plot <- plotPCA(normalized_data)`,
                options: [
                    '仅需要标准化数据',
                    '标准化和质量控制',
                    '标准化、质量控制和批次效应检查',
                    '直接进行差异分析'
                ],
                correctAnswer: 2,
                explanation: '完整的预处理流程应包括：数据标准化、质量控制和批次效应检查。这些步骤能确保后续分析的可靠性。'
            },
            {
                category: '差异表达分析',
                difficulty: 'medium',
                question: '在使用limma包进行差异表达分析时，以下哪个步骤是正确的？',
                code: `# 差异分析示例
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(
    AD-Normal,
    levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEGs <- topTable(fit2, coef = 1, n = Inf)`,
                options: [
                    '直接使用t检验比较两组',
                    '不需要设计矩阵，直接用eBayes',
                    '先构建设计矩阵，再用lmFit和eBayes',
                    '只需要用wilcox.test就够了'
                ],
                correctAnswer: 2,
                explanation: 'limma包的差异分析流程需要：1)构建设计矩阵；2)使用lmFit拟合线性模型；3)使用contrasts.fit设置对比；4)使用eBayes进行贝叶斯估计。这样可以得到更稳健的结果。'
            },
            {
                category: '差异表达分析',
                difficulty: 'hard',
                question: '在筛选差异表达基因时，通常使用什么标准？',
                code: `# 差异基因筛选示例
DEGs <- topTable(fit2, coef = 1, n = Inf) %>% 
    filter(abs(logFC) > 1 & adj.P.Val < 0.05)

# 添加上下调标记
DEGs <- DEGs %>% 
    mutate(change = case_when(
        logFC > 1 & adj.P.Val < 0.05 ~ "Up",
        logFC < -1 & adj.P.Val < 0.05 ~ "Down",
        TRUE ~ "NS"
    ))`,
                options: [
                    '只看p值小于0.05',
                    '只看fold change大于2',
                    'p值小于0.05且|logFC|大于1',
                    '随便选择一些基因就行'
                ],
                correctAnswer: 2,
                explanation: '差异基因的筛选通常需要同时考虑统计显著性(p值或校正后的p值)和生物学意义(fold change)。通常的标准是adj.P.Val < 0.05且|logFC| > 1，这样可以确保筛选出的基因既具有统计显著性，又有足够大的表达差异。'
            },
            {
                category: '功能富集分析',
                difficulty: 'medium',
                question: '使用clusterProfiler进行GO富集分析时，enrichGO函数的关键参数是什么？',
                code: `# GO富集分析示例
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
)`,
                options: [
                    '只需要提供基因列表',
                    '需要基因列表和物种数据库',
                    '基因列表、物种数据库和ID类型都需要',
                    '不需要设置任何参数'
                ],
                correctAnswer: 2,
                explanation: 'enrichGO函数需要设置多个关键参数：1)gene：基因列表；2)OrgDb：物种注释数据库；3)keyType：基因ID类型；4)ont：GO分类(BP/MF/CC)；5)pAdjustMethod：p值校正方法。这些参数都很重要，确保富集分析的准确性。'
            },
            {
                category: '功能富集分析',
                difficulty: 'hard',
                question: '在进行KEGG富集分析之前，为什么需要进行基因ID转换？',
                code: `# 基因ID转换示例
gene_list <- bitr(
    gene = gene_list,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

# KEGG富集分析
ekegg <- enrichKEGG(
    gene = gene_list$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05
)`,
                options: [
                    '不需要转换，直接用基因名就可以',
                    'KEGG数据库只认ENTREZ ID',
                    '转换是可选的',
                    '任何ID类型都可以'
                ],
                correctAnswer: 1,
                explanation: 'KEGG数据库使用ENTREZ ID作为基因标识符。因此，如果我们的基因列表使用其他类型的ID（如SYMBOL、ENSEMBL等），就需要先转换为ENTREZ ID才能进行KEGG富集分析。'
            },
            {
                category: '数据可视化',
                difficulty: 'medium',
                question: '在使用pheatmap绘制热图时，如何正确处理数据标准化？',
                code: `# 热图数据处理和绘制
# 数据标准化
scaled_data <- t(scale(t(expr_matrix)))

# 绘制热图
pheatmap(
    scaled_data,
    scale = "none",
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE
)`,
                options: [
                    '不需要标准化',
                    '只需要log2转换',
                    '需要进行Z-score标准化',
                    '直接使用原始值'
                ],
                correctAnswer: 2,
                explanation: '在绘制热图之前，通常需要对数据进行Z-score标准化，使不同基因之间的表达值具有可比性。这可以通过scale函数实现，它会将数据转换为均值为0、标准差为1的标准正态分布。'
            }
        ]
    },
    'd7': {
        title: "D7作业 - GSEA & GSVA分析",
        description: "本测验涵盖基因集富集分析(GSEA)和基因集变异分析(GSVA)的相关内容。",
        questions: [
            {
                category: 'GSEA分析',
                difficulty: 'medium',
                question: "在进行GSEA分析时，我们需要准备排序后的基因列表(gene_list)。以下哪个语句正确描述了gene_list的准备过程？",
                code: `# 准备gene_list
res <- read.csv("deg_results.csv")
gene_list <- res$logFC
names(gene_list) <- res$SYMBOL
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA分析
gsea_result <- GSEA(
    gene_list,
    TERM2GENE = pathway2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05
)`,
                options: [
                    "gene_list是根据基因名字的字母顺序排序的",
                    "gene_list是根据logFC值从大到小排序的",
                    "gene_list是根据p值从小到大排序的",
                    "gene_list是随机排序的"
                ],
                correctAnswer: 1,
                explanation: "在代码中可以看到：gene_list <- sort(gene_list, decreasing = TRUE)，这表明gene_list是根据logFC值从大到小排序的。这个排序对于GSEA分析很重要，因为它决定了基因的排名。"
            },
            {
                category: 'GSEA分析',
                difficulty: 'medium',
                question: "在脚本中设置GSEA分析参数时，以下哪个参数设置是错误的？",
                code: `# GSEA参数设置
gsea_result <- GSEA(
    gene_list,
    TERM2GENE = pathway2gene,
    minGSSize = 10,    # 最小基因集大小
    maxGSSize = 500,   # 最大基因集大小
    pvalueCutoff = 0.05,
    by = "fgsea"       # 使用fgsea算法
)`,
                options: [
                    "最小基因集大小(minGSSize)设置为10",
                    "最大基因集大小(maxGSSize)设置为500",
                    "p值阈值(pvalueCutoff)设置为0.05",
                    "以上设置都是正确的"
                ],
                correctAnswer: 3,
                explanation: "脚本中的GSEA参数设置都是合理的：minGSSize = 10（避免过小的基因集），maxGSSize = 500（避免过大的基因集），pvalueCutoff = 0.05（标准的显著性水平）。"
            },
            {
                category: 'GSVA分析',
                difficulty: 'hard',
                question: "关于GSVA分析，以下哪个说法是正确的？",
                code: `# GSVA分析示例
library(GSVA)
gsva_matrix <- gsva(
    expr = expr_matrix,
    gset.idx.list = pathway_list,
    method = "gsva",
    kcdf = "Gaussian",
    verbose = TRUE
)`,
                options: [
                    "GSVA分析不需要分组信息就能进行",
                    "GSVA分析只能用于两组样本的比较",
                    "GSVA分析必须使用原始计数矩阵",
                    "GSVA分析只能用于人类基因组数据"
                ],
                correctAnswer: 0,
                explanation: "GSVA（Gene Set Variation Analysis）是一种无监督的方法，它不需要预先的分组信息就能计算基因集富集分数。这也是它与GSEA的一个重要区别。"
            },
            {
                category: '结果可视化',
                difficulty: 'medium',
                question: "在脚本中绘制GSVA热图时使用了哪些参数来优化可视化效果？",
                code: `# 热图绘制
pheatmap(
    gsva_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8
)`,
                options: [
                    "只使用了默认参数",
                    "使用了自定义的颜色方案和聚类方法",
                    "只改变了图形大小",
                    "没有进行热图可视化"
                ],
                correctAnswer: 1,
                explanation: "脚本中使用了多个参数来优化热图：自定义颜色方案(colorRampPalette)、聚类距离(euclidean)、聚类方法(complete)、字体大小设置等。"
            },
            {
                category: 'GSEA可视化',
                difficulty: 'medium',
                question: "在进行GSEA结果可视化时，使用了gseaNb函数。以下哪个参数设置是正确的？",
                code: `# GSEA结果可视化
gseaNb(
    object = gsea_result,
    geneSetID = "HALLMARK_APOPTOSIS",
    addPval = TRUE,
    pvalX = 0.65,
    pvalY = 0.7
)`,
                options: [
                    "设置了pCol参数为'red'",
                    "设置了pvalX和pvalY参数来调整p值的位置",
                    "没有添加任何统计信息",
                    "只显示了基因名称"
                ],
                correctAnswer: 1,
                explanation: "在代码中，gseaNb函数使用了pvalX = 0.65和pvalY = 0.7来调整p值的显示位置，这样可以更好地展示统计信息。"
            }
        ]
    },
    'gene_correlation': {
        title: "D8作业-基因相关性分析",
        description: "本测验主要考察基因相关性分析的方法和可视化技巧。",
        questions: [
            {
                category: '基因相关性分析',
                difficulty: 'medium',
                question: '在进行基因相关性分析时，以下哪种可视化方法最适合展示大量基因之间的相互关系并同时显示相关性强度？',
                code: `# 相关性分析代码示例
# 计算相关性矩阵
cor_matrix <- cor(t(gene_expr))

# 方法1：和弦图
chordDiagram(cor_matrix,
    grid.col = rainbow(nrow(cor_matrix)),
    symmetric = TRUE)

# 方法2：热图
pheatmap(cor_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(100))

# 方法3：散点图矩阵
pairs(t(gene_expr))`,
                options: [
                    '散点图矩阵 - 因为它可以直观显示每对基因的具体关系',
                    '和弦图 - 因为它能够优雅地展示基因间的连接关系',
                    '热图 - 因为它既能展示层次聚类结构又能显示相关性强度',
                    '线图 - 因为它能显示基因表达量的变化趋势'
                ],
                correctAnswer: 2,
                explanation: '热图是展示大规模基因相关性最理想的方式，因为它不仅能通过颜色深浅直观地展示相关性强度，还能通过层次聚类展示基因间的关系模式。和弦图虽然美观，但在基因数量较多时可能会显得混乱；散点图矩阵在基因数量大时会占用过多空间；而线图则不适合展示相关性关系。'
            },
            {
                category: '数据预处理',
                difficulty: 'medium',
                question: '在计算基因相关性之前，为什么需要对表达数据进行标准化？',
                code: `# 数据标准化示例
# 方法1：Z-score标准化
scaled_expr <- t(scale(t(gene_expr)))

# 方法2：Min-Max标准化
min_max_scale <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}
normalized_expr <- t(apply(gene_expr, 1, min_max_scale))`,
                options: [
                    '标准化会改变基因间的相关性',
                    '使不同基因的表达值在相同尺度上便于比较',
                    '标准化是可选的，不会影响结果',
                    '只需要对部分基因进行标准化'
                ],
                correctAnswer: 1,
                explanation: '标准化是必要的预处理步骤，因为不同基因的表达水平可能差异很大。通过标准化（如Z-score标准化），可以消除量级差异的影响，使基因表达值具有可比性。这样计算出的相关性才能真实反映基因表达模式的相似度。'
            },
            {
                category: '相关性分析',
                difficulty: 'hard',
                question: '在R中计算基因相关性时，以下哪个方法最适合处理基因表达数据？',
                code: `# 不同相关性计算方法
# Pearson相关系数
pearson_cor <- cor(gene1_expr, gene2_expr, 
                  method = "pearson")

# Spearman相关系数
spearman_cor <- cor(gene1_expr, gene2_expr, 
                   method = "spearman")

# Kendall相关系数
kendall_cor <- cor(gene1_expr, gene2_expr, 
                   method = "kendall")`,
                options: [
                    'Pearson相关系数 - 因为它最常用',
                    'Spearman相关系数 - 因为它对异常值不敏感',
                    'Kendall相关系数 - 因为它计算最准确',
                    '三种方法都需要尝试比较'
                ],
                correctAnswer: 1,
                explanation: 'Spearman相关系数是处理基因表达数据的优选方法，因为：1)它不要求数据呈正态分布；2)对异常值不敏感；3)能够捕捉非线性关系。基因表达数据通常包含噪声和异常值，使用Spearman相关系数可以得到更稳健的结果。'
            },
            {
                category: '相关性分析',
                difficulty: 'medium',
                question: '在进行基因共表达网络分析时，为什么要计算基因表达的相关系数矩阵？',
                code: `# 计算相关系数矩阵
cor_matrix <- cor(t(expr_matrix), 
                 method = "spearman",
                 use = "pairwise.complete.obs")

# 设置相关性阈值
threshold <- 0.7
cor_matrix[abs(cor_matrix) < threshold] <- 0`,
                options: [
                    '相关系数矩阵可以反映基因间的表达模式相似性',
                    '只是为了数据可视化',
                    '相关系数计算是可选的步骤',
                    '相关系数只用于检测异常值'
                ],
                correctAnswer: 0,
                explanation: '相关系数矩阵是构建基因共表达网络的基础，它反映了基因间表达模式的相似性。高相关性表明基因可能参与相似的生物学过程或受相同的调控机制控制。通过设置相关性阈值，可以筛选出显著的基因对关系。'
            },
            {
                category: '相关性分析',
                difficulty: 'hard',
                question: '在分析基因表达数据的相关性时，Pearson相关系数和Spearman相关系数的选择依据是什么？',
                code: `# 比较Pearson和Spearman相关系数
pearson_cor <- cor(gene1_expr, gene2_expr, 
                  method = "pearson")
spearman_cor <- cor(gene1_expr, gene2_expr, 
                   method = "spearman")

# 可视化两个基因的表达关系
plot(gene1_expr, gene2_expr)
abline(lm(gene2_expr ~ gene1_expr))`,
                options: [
                    'Pearson更好，因为计算速度快',
                    'Spearman更好，因为考虑了秩次关系',
                    '根据数据分布特征选择：正态分布用Pearson，非正态或有异常值用Spearman',
                    '两种方法没有区别，随便选择'
                ],
                correctAnswer: 2,
                explanation: '选择相关系数类型需要考虑数据特征：1) Pearson相关系数适用于呈正态分布且具有线性关系的数据；2) Spearman相关系数基于秩次，对异常值不敏感，适用于非正态分布或存在异常值的数据。在基因表达数据分析中，如果数据经过良好的预处理且呈正态分布，可以使用Pearson；如果数据分布不确定或存在异常值，建议使用Spearman。'
            }
        ]
    },
    'D9': {
        title: "D9作业-免疫浸润分析",
        description: "本测验主要考察免疫浸润分析相关知识点。",
        questions: [
            {
                category: '免疫浸润分析',
                difficulty: 'medium',
                question: '在使用CIBERSORT进行免疫浸润分析时，为什么需要对基因表达数据进行标准化处理？',
                code: `# 数据标准化示例
# 1. TPM标准化
tpm_matrix <- counts2tpm(counts_matrix, gene_length)

# 2. log2转换
log2_tpm <- log2(tpm_matrix + 1)

# 3. 进行CIBERSORT分析
results <- CIBERSORT(sig_matrix, mixture_data)`,
                options: [
                    '标准化不是必需的，直接用原始数据即可',
                    '为了消除测序深度和基因长度的影响，使不同样本可比',
                    '标准化只是为了数据可视化',
                    '标准化会降低分析的准确性'
                ],
                correctAnswer: 1,
                explanation: '在进行CIBERSORT分析之前，需要对数据进行适当的标准化处理：1)将counts数据转换为TPM可以消除测序深度和基因长度的影响；2)log2转换可以使数据分布更接近正态分布；3)标准化后的数据更适合进行跨样本比较，提高免疫浸润估计的准确性。'
            },
            {
                category: '免疫浸润分析',
                difficulty: 'hard',
                question: '在解释CIBERSORT结果时，以下哪个说法是正确的？',
                code: `# CIBERSORT结果分析
cibersort_results <- read.table("cibersort_results.txt")

# 检查P值
significant_samples <- cibersort_results$P.value < 0.05

# 计算免疫细胞比例
cell_proportions <- cibersort_results[, 1:22]
boxplot(cell_proportions, las = 2)`,
                options: [
                    'P值大于0.05的样本结果更可靠',
                    'P值小于0.05表示免疫细胞组成估计结果可靠',
                    '不需要考虑P值，所有结果都同样可靠',
                    'P值只与样本数量有关'
                ],
                correctAnswer: 1,
                explanation: 'CIBERSORT的P值反映了免疫细胞组成估计的可靠性：1)P值<0.05表示估计结果统计学显著，更可靠；2)应该优先考虑P值显著的样本进行后续分析；3)P值的计算基于置换检验，反映了估计结果的置信度。'
            },
            {
                category: '免疫浸润分析',
                difficulty: 'medium',
                question: '在比较不同组织或疾病状态下的免疫浸润差异时，应该使用什么统计方法？',
                code: `# 免疫浸润差异分析
# 方法1：t检验
t_test_results <- apply(cell_proportions, 2, function(x) {
    t.test(x ~ group)$p.value
})

# 方法2：Wilcoxon秩和检验
wilcox_results <- apply(cell_proportions, 2, function(x) {
    wilcox.test(x ~ group)$p.value
})

# 多重检验校正
adjusted_pvals <- p.adjust(wilcox_results, method = "BH")`,
                options: [
                    '只用t检验就够了',
                    '只需要看细胞比例的平均值',
                    '应该结合非参数检验和多重检验校正',
                    '不需要统计检验'
                ],
                correctAnswer: 2,
                explanation: '比较免疫浸润差异时的统计考虑：1)由于免疫细胞比例数据通常不符合正态分布，建议使用非参数检验如Wilcoxon秩和检验；2)因为同时比较多个免疫细胞类型，需要进行多重检验校正以控制假阳性率；3)可以结合箱线图等可视化方法展示差异。'
            },
            {
                category: '免疫浸润分析',
                difficulty: 'hard',
                question: '如何评估和解释不同免疫细胞类型之间的相关性？',
                code: `# 免疫细胞相关性分析
# 计算相关性矩阵
cor_matrix <- cor(cell_proportions, 
                 method = "spearman")

# 可视化相关性
library(corrplot)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45)`,
                options: [
                    '正相关意味着细胞类型完全相同',
                    '负相关表示细胞类型之间没有关系',
                    '相关性反映了细胞类型在微环境中的共存或拮抗关系',
                    '不同免疫细胞之间不存在相关性'
                ],
                correctAnswer: 2,
                explanation: '免疫细胞相关性分析的解释：1)正相关可能反映细胞类型在功能上的协同作用或共同调控机制；2)负相关可能表示细胞类型之间的拮抗关系或相互抑制；3)相关性分析有助于理解肿瘤微环境中免疫细胞的相互作用网络；4)应结合生物学知识解释相关性结果。'
            },
            {
                category: '免疫浸润分析',
                difficulty: 'hard',
                question: '在进行免疫浸润分析时，如何评估结果的生物学意义？',
                code: `# 免疫浸润结果与临床相关性分析
# 生存分析
library(survival)
library(survminer)

# 根据免疫细胞比例分组
median_value <- median(cell_proportions$CD8T)
groups <- ifelse(cell_proportions$CD8T > median_value, 
                "High", "Low")

# KM曲线
fit <- survfit(Surv(time, status) ~ groups)
ggsurvplot(fit, 
           pval = TRUE,
           risk.table = TRUE)`,
                options: [
                    '只需要看免疫细胞比例的高低',
                    '直接用P值判断重要性',
                    '需要结合临床特征、生存预后等多个方面综合分析',
                    '不同样本之间不能比较'
                ],
                correctAnswer: 2,
                explanation: '评估免疫浸润结果的生物学意义需要多个层面：1)结合临床特征分析免疫细胞组成与疾病进展的关系；2)通过生存分析评估免疫细胞比例与预后的关联；3)整合其他分子特征数据，如基因突变、表达谱等；4)参考已有文献和数据库中的相关研究结果；5)考虑样本类型和疾病背景的特异性。'
            }
        ]
    },
    'D10': {
        title: "R代码讲席营D10课程测验",
        description: "本测验涵盖阿尔茨海默病(AD)诊断模型构建与评估的相关内容。",
        questions: [
            {
                category: '诊断模型构建',
                difficulty: 'medium',
                question: '在构建AD诊断模型时，ROC分析中的AUC值代表什么？',
                code: `# ROC分析示例
roc_res <- roc(group ~ expression,
               data = dat_expr,
               auc = TRUE,    
               ci = TRUE)     

auc_text <- sprintf("AUC = %.3f (%.3f-%.3f)",
                   roc_res$auc,
                   roc_res$ci[1],
                   roc_res$ci[3])`,
                options: [
                    'ROC曲线下的面积，表示模型的整体诊断效能',
                    '真阳性率和假阳性率的比值',
                    '模型的准确率',
                    '敏感性和特异性的平均值'
                ],
                correctAnswer: 0,
                explanation: 'AUC（Area Under Curve）是ROC曲线下的面积，取值范围在0-1之间。AUC越接近1，表示模型的诊断效能越好。AUC=0.5表示模型的诊断效能等同于随机猜测。'
            },
            {
                category: '特征选择',
                difficulty: 'hard',
                question: '在使用LASSO进行特征选择时，lambda参数的作用是什么？',
                code: `# LASSO回归示例
cv_fit <- cv.glmnet(x = as.matrix(expr_selected),
                    y = as.factor(group_data),
                    family = "binomial",
                    alpha = 1,
                    nfolds = 5)`,
                options: [
                    'lambda值越大，模型越简单，选择的特征越少',
                    'lambda值越大，模型越复杂，选择的特征越多',
                    'lambda值不影响特征选择',
                    'lambda值只影响模型的收敛速度'
                ],
                correctAnswer: 0,
                explanation: 'LASSO回归中的lambda是正则化参数，控制L1惩罚项的强度。lambda值越大，惩罚越强，更多特征的系数会被压缩为0，从而实现特征选择。这有助于避免过拟合并构建更简约的模型。'
            },
            {
                category: '模型评估',
                difficulty: 'medium',
                question: '在构建诊断列线图(Nomogram)时，为什么需要进行校正曲线(Calibration Curve)分析？',
                code: `# 校正曲线分析示例
cal <- calibrate(log_fit_lrm, method = "boot", B = 100)
plot(cal, xlab = "Predicted Probability", 
     ylab = "Observed Probability")`,
                options: [
                    '评估模型的区分度',
                    '评估预测概率与实际观察概率的一致性',
                    '计算模型的敏感性和特异性',
                    '确定最佳截断值'
                ],
                correctAnswer: 1,
                explanation: '校正曲线用于评估模型预测的概率与实际观察到的概率之间的一致性。理想情况下，预测概率应该与观察概率相近，即校正曲线应该接近45度对角线。这是评估诊断模型可靠性的重要工具。'
            },
            {
                category: '统计分析',
                difficulty: 'medium',
                question: '在进行Logistic回归分析时，OR值(优势比)大于1代表什么？',
                code: `# Logistic回归分析示例
fit0 <- glm(group ~ gene_expr, 
            data = expr_data, 
            family = binomial)
OR <- exp(coef(fit0)[2])  # 计算OR值`,
                options: [
                    '该变量与疾病风险呈负相关',
                    '该变量与疾病风险无关',
                    '该变量与疾病风险呈正相关',
                    '无法判断相关性'
                ],
                correctAnswer: 2,
                explanation: 'OR值(优势比)大于1表示该变量与疾病风险呈正相关，即该变量每增加一个单位，患病的几率会增加。具体增加的倍数等于OR值。这是评估基因表达与疾病关联性的重要指标。'
            }
        ]
    }
}
