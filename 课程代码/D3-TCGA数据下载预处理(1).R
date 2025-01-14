## ====================================================================
## DAY 3. R语言数据清洗-TCGA&GEO
## 
## 功能：下载和预处理TCGA数据，包括表达谱和临床数据
## 输出：
##   - [cancer_type]_counts.Rda: 原始计数矩阵
##   - [cancer_type]_cli.Rda: 临床数据
##   - [cancer_type]_tpm.Rda: TPM标准化表达矩阵
## ====================================================================

#####################################################
## Step1：环境初始化
## - 清理工作环境
## - 设置下载超时参数
#####################################################
rm(list = ls())
options(timeout = 100000000)
setwd("./D3-TCGA数据下载预处理")

# 加载必要的R包
library(TCGAbiolinks)   # TCGA数据下载和处理
library(tidyverse)      # 数据处理和管道操作
library(GenomicFeatures)# 基因组特征注释

#####################################################
## Step2：TCGA数据获取
## - 设置癌症类型
## - 下载RNA-seq表达数据：
##   * 数据类型：基因表达定量
##   * 工作流：STAR - Counts
## - 下载临床数据：
##   * 格式：BCR XML
##   * 类型：患者信息
#####################################################
cancer_type <- 'TCGA-BRCA' # 设置癌症类型为乳腺浸润癌

query_exp <- GDCquery(
    project = cancer_type,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

GDCdownload(query_exp, 
    directory = "GDCdata")

counts_data <- GDCprepare(query_exp) %>% 
    SummarizedExperiment::assay()

query_cli <- GDCquery(
    project = cancer_type,
    data.category = "Clinical",
    data.format = "BCR XML"
)

GDCdownload(query_cli)
cli_data <- GDCprepare_clinic(query_cli, clinical.info = 'patient')

#####################################################
## Step3：基因组注释处理
## - 下载最新版本的Gencode注释文件
## - 创建转录组数据库
## - 计算基因长度：
##   * 合并重叠外显子
##   * 计算总长度
#####################################################
gtf_url <- 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v47.annotation.gtf.gz'
gtf_file <- 'gencode.v47.annotation.gtf.gz'
download.file(gtf_url, gtf_file)

txdb <- makeTxDbFromGFF(gtf_file, format = 'gtf')
gene_length <- exonsBy(txdb, by = 'gene') %>% 
  lapply(function(x) sum(width(reduce(x)))) %>% 
  unlist() %>% 
  data.frame(gene_length = .)

#####################################################
## Step4：表达数据标准化
## - 基因ID处理：
##   * 移除版本号
##   * 筛选共同基因
## - 表达量转换：
##   * counts -> FPKM：考虑基因长度和测序深度
##   * FPKM -> TPM：标准化到每百万转录本
#####################################################
process_ids <- function(x) {
    x <- sub("\\.[0-9]+.*$", "", x)
    return(x)
}
rownames(counts_data) <- process_ids(rownames(counts_data))
rownames(gene_length) <- process_ids(rownames(gene_length))
gene_length <- as.matrix(gene_length)

common_genes <- intersect(rownames(counts_data), rownames(gene_length))
counts_data <- counts_data[common_genes, ]
gene_length <- gene_length[common_genes, ]

counts2fpkm <- function(counts, gen_len) {
    N <- sum(counts)
    counts * 1e9 / gen_len / N
}

fpkm2tpm <- function(fpkm) {
    fpkm / sum(fpkm) * 1e6
}

fpkm_data <- apply(counts_data, 2, counts2fpkm, 
                  gen_len = gene_length) %>% as.data.frame()
tpm_data <- apply(fpkm_data, 2, fpkm2tpm) %>% as.data.frame()

#####################################################
## Step5：结果保存
## - 保存原始计数矩阵
## - 保存临床数据
## - 保存TPM标准化表达矩阵
#####################################################
save(counts_data, file = paste0(cancer_type, '_counts.Rda'))
save(cli_data, file = paste0(cancer_type, '_cli.Rda'))
save(tpm_data, file = paste0(cancer_type, '_tpm.Rda'))
