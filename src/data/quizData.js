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
    }
}
