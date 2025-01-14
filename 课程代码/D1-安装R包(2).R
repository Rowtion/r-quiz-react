## ====================================================================
## DAY 1. R语言基础入门
## 
## 描述：
##   安装本课程中所需的所有R包，包括CRAN包和Bioconductor包
##
## 注意事项：
##   - 如遇到安装错误，请检查R版本和镜像设置
## ====================================================================

#####################################################
## Step1：CRAN包安装
## 
## 描述：安装来自CRAN的R包
#####################################################

# 定义CRAN包列表
cran_packages <- c(
    "tidyverse",     # 数据处理和管道操作
    "pheatmap",      # 热图绘制
    "circlize",      # 和弦图绘制
    "VennDiagram",   # 韦恩图绘制
    "corrplot",      # 相关性图绘制
    "ggforce",       # 椭圆添加
    "ggrepel",       # 标签标注
    "ggpubr",        # 可视化拓展包
    "gridExtra",      # 添加图片信息
    "patchwork"      # 拼图
)
# 安装CRAN包
for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        cat(sprintf("已安装包：%s\n", pkg))
    } else {
        cat(sprintf("包已存在：%s\n", pkg))
    }
}

#####################################################
## Step2：Bioconductor包安装
## 
## 描述：安装来自Bioconductor的R包
#####################################################

# 安装BiocManager
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    cat("已安装BiocManager\n")
}

# 定义Bioconductor包列表
bioc_packages <- c(
    "GEOquery",        # GEO数据处理
    "limma",           # 差异分析
    "TCGAbiolinks",    # TCGA数据下载和处理
    "GenomicFeatures", # 基因组特征注释
    "clusterProfiler", # 富集分析
    "org.Hs.eg.db",    # 人类基因注释数据库
    "GSVA",            # 基因集变异分析
    "GSEABase",        # 基因集处理
    "ComplexHeatmap",  # 复杂热图绘制
    "GseaVis"          # GSEA结果可视化
)

# 安装Bioconductor包
for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg)
        cat(sprintf("已安装包：%s\n", pkg))
    } else {
        cat(sprintf("包已存在：%s\n", pkg))
    }
}

# 打印安装完成信息
cat("\n所有包安装完成！\n")
