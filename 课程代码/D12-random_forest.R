## ====================================================================
## 随机森林基因筛选
##
## 描述：
##   使用随机森林算法筛选common_gene中的重要基因
##
## 输入：
##   - Common_genes.csv: 候选基因列表
##   - combined_groups.Rda: 样本分组信息
##   - normalized_expression_data.Rda: 标准化表达矩阵
##
## 输出：
##   - important_genes.csv: 重要基因列表
##   - lollipop_plot.pdf: 重要基因棒棒糖图
##
## 注意事项：
##   - 需要足够的内存空间
##   - 建议在服务器上运行
## ====================================================================

#####################################################
## Step1：环境初始化
##
## 描述：初始化R环境，设置参数，加载必要的包
#####################################################

# 清空工作环境
rm(list = ls())

# 基础设置
options(
    scipen = 20,      # 避免科学计数法
    stringsAsFactors = FALSE
)

# 加载必要的R包
library(randomForest)  # 随机森林算法
library(caret)         # 模型评估
library(tidyverse)     # 数据处理
library(pROC)          # ROC曲线分析
library(ggpubr)        # 绘图

#####################################################
## Step2：数据加载
##
## 描述：加载输入数据并进行预处理
#####################################################

# 加载基因列表
common_genes <- read.csv("input/Common_genes.csv", header = TRUE)

# 加载分组信息
load("input/group_GEO.Rda")

# 加载标准化表达数据
load("input/normalized_expression_data.Rda")

# 数据预处理
expression_data <- as.data.frame(t(adjusted_matrix))
expression_data <- expression_data[, colnames(expression_data) %in% common_genes$Gene]

# 合并表达数据和分组信息
if(is.data.frame(combined_groups)) {
  combined_groups <- combined_groups[,2]
}

data <- cbind(group = as.factor(combined_groups), expression_data)
colnames(data)[1] <- "group"  # 确保第一列名为group

#####################################################
## Step3：随机森林分析
##
## 描述：构建随机森林模型并筛选重要基因
#####################################################

# 设置随机种子保证结果可重复
set.seed(123)

# 构建随机森林模型
rf_model <- randomForest(
    group ~ .,
    data = data,
    importance = TRUE,
    ntree = 1000,  # 增加树的数量
    mtry = sqrt(ncol(data) - 1),  # 使用平方根特征
    nodesize = 5,  # 调整节点大小
    sampsize = c(min(table(data$group)), min(table(data$group)))  # 平衡样本
)

# 提取基因重要性
importance <- importance(rf_model, type = 1)
importance <- data.frame(
    Gene = rownames(importance),
    Importance = importance[, 1]
) %>%
    arrange(desc(Importance))

# 选择top 10重要基因
top_genes <- importance[1:10, ]

#####################################################
## Step4：结果可视化
##
## 描述：生成棒棒糖图展示重要基因
#####################################################

# 生成棒棒糖图
pdf("output/lollipop_plot.pdf", width = 10, height = 6)
ggdotchart(
    top_genes,
    x = "Gene",
    y = "Importance",
    sorting = "descending",
    add = "segments",
    rotate = TRUE,
    dot.size = 6,
    label = round(top_genes$Importance, 2),
    font.label = list(color = "white", size = 9, vjust = 0.5),
    ggtheme = theme_pubr()
) +
labs(
    title = "Top 10 Important Genes",
    x = "Gene",
    y = "Importance Score"
)
dev.off()

#####################################################
## Step5：结果保存
##
## 描述：保存重要基因列表
#####################################################

# 保存重要基因列表
write.csv(
    top_genes,
    file = "output/important_genes.csv",
    row.names = FALSE
)

