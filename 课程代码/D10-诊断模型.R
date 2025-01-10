# ==============================================================================
# 阿尔茨海默病(AD)诊断模型构建与评估
# 
# 功能：构建和评估AD诊断模型，包括：
#   1. ROC分析
#   2. Logistic回归分析
#   3. 诊断列线图
#   4. LASSO特征选择
# 
# 输入文件：
#   - normalized_expression_data.Rda: 标准化表达矩阵
#   - group_GEO.Rda: 样本分组信息
#   - Common_genes.csv: 交集基因列表
# 
# 输出文件：
#   - 1-ROC.pdf: ROC曲线图
#   - 2-Logistic_forestplot.pdf: Logistic回归森林图
#   - 3-Nomogram.pdf: 诊断列线图
#   - 4-CalibrationCurve.pdf: 校正曲线
#   - 5-LASSO_lambda.pdf: LASSO系数变化图
#   - 6-LASSO_coef.pdf: LASSO交叉验证结果图
#   - selected_genes.csv: LASSO筛选的基因列表
# ==============================================================================


# ==============================================================================
# Step 1: 环境准备
# - 清空环境变量
# - 加载必要的R包
# - 设置统一的主题和配色
# - 导入分析所需数据
# ==============================================================================

#------------------------------------
# 1.1 清空环境
#------------------------------------
# 清除所有变量，确保环境干净
rm(list = ls())

#------------------------------------
# 1.2 加载必要的R包
#------------------------------------
# 数据处理和可视化相关包
library(tidyverse)  # 数据处理和管道操作
library(ggplot2)    # 绘图系统

# 模型构建和评估相关包
library(pROC)       # ROC曲线分析和AUC计算
library(rms)        # 回归建模和诊断图
library(Hmisc)      # rms包的依赖包，提供数据处理函数
library(survival)   # rms包的依赖包，生存分析相关
library(glmnet)     # LASSO回归实现
library(forestplot) # 森林图可视化
library(ggDCA)      # 决策曲线分析

#------------------------------------
# 1.3 设置配色
#------------------------------------
# 定义统一的配色方案，确保图形风格一致性
color_palette <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488")

#------------------------------------
# 1.4 导入数据
#------------------------------------
# 导入标准化后的基因表达数据
load('./input/normalized_expression_data.Rda')
# 导入样本分组信息
load('./input/group_GEO.Rda')
# 导入经过筛选的共同基因列表
genes <- read.csv("./input/Common_genes.csv")


# ==============================================================================
# Step 2: ROC分析
# - 以PSMD3为例进行ROC分析
# - 计算AUC值和置信区间
# - 绘制ROC曲线并保存结果
# ==============================================================================

#------------------------------------
# 2.1 创建ROC分析函数
#------------------------------------
# 函数说明：
# - 输入参数：
#   gene_expr: 单个基因的表达值向量
#   group_data: 样本分组信息
#   gene_name: 基因名称
# - 返回值：
#   包含ROC对象和ggplot2图形对象的列表
perform_roc_analysis <- function(gene_expr, group_data, gene_name) {
  # 整理数据为数据框格式
  dat_expr <- data.frame(
    expression = gene_expr,
    group = group_data
  )
  
  # 计算ROC曲线及其置信区间
  roc_res <- roc(group ~ expression,
                 data = dat_expr,
                 auc = TRUE,    # 计算AUC值
                 ci = TRUE)     # 计算95%置信区间
  
  # 格式化AUC结果文本
  auc_text <- sprintf("AUC = %.3f (%.3f-%.3f)",
                     roc_res$auc,
                     roc_res$ci[1],
                     roc_res$ci[3])
  
  # 创建ROC曲线图
  p <- ggroc(roc_res,
             legacy.axes = TRUE,     # 使用传统坐标轴
             col = color_palette[1], # 使用预设配色
             lwd = 1.2) +
    # 添加ROC曲线下方的阴影区域
    geom_polygon(data = data.frame(
      x = c(1, 1-roc_res$specificities, 0),
      y = c(0, roc_res$sensitivities, 0)
    ), aes(x = x, y = y), fill = color_palette[1], alpha = 0.2) +
    labs(title = paste("ROC Curve for", gene_name),
         x = "1-Specificity (FPR)",
         y = "Sensitivity (TPR)") +
    # 添加AUC值标注
    annotate("text",
             x = 0.75,
             y = 0.25,
             label = auc_text,
             size = 4,
             color = color_palette[1]) +
    coord_fixed() +  # 保持坐标轴比例1:1
    # 添加对角参考线
    geom_abline(
      slope = 1,
      intercept = 0,
      lty = "dashed",
      color = "gray50"
    ) +
    theme_bw()

  return(list(roc = roc_res, plot = p))
}

#------------------------------------
# 2.2 执行ROC分析
#------------------------------------
# 提取PSMD3基因表达数据并转换格式
dat_expr_PSMD3 <- data.frame(
  PSMD3 = as.numeric(adjusted_matrix["PSMD3", ]),  # 提取PSMD3基因表达值
  group = combined_groups$group                     # 对应的分组信息
)

# 对PSMD3基因进行ROC分析
roc_results_PSMD3 <- perform_roc_analysis(
  dat_expr_PSMD3$PSMD3, 
  dat_expr_PSMD3$group,
  "PSMD3"
)

# 保存ROC曲线图为高质量PDF文件
ggsave(
  file = "./output/1-ROC.pdf",
  roc_results_PSMD3$plot,
  height = 6,
  width = 6,
  dpi = 300  # 设置较高分辨率
)


# ==============================================================================
# Step 3: Logistic回归分析
# - 包括单因素和多因素Logistic回归
# - 用于评估基因表达与AD的关联性
# - 计算OR值及其置信区间
# ==============================================================================

# 准备分析数据
expr_genes <- adjusted_matrix[genes$Gene, ] %>% 
  na.omit() %>%                # 移除缺失值
  t() %>%                      # 转置矩阵以适应建模需求
  as.data.frame()
expr_genes$group <- factor(combined_groups$group, 
                          levels = c('Normal', 'AD'))  # 设置对照组为参考水平

# 单因素Logistic回归分析
log_uni <- data.frame()

# 对每个基因进行单独的Logistic回归分析
for (gene in genes$Gene) {
  # 构建并拟合单因素模型
  fit0 <- glm(as.formula(paste0('group ~ ', gene)), 
              data = expr_genes, 
              family = binomial)  # 使用logit连接函数
  
  # 提取模型结果
  fit_sum <- summary(fit0)
  conf_int <- exp(confint(fit0))  # 计算OR值的95%置信区间
  
  # 整理单个基因的分析结果
  result <- data.frame(
    gene = gene,
    OR = exp(fit_sum$coef[2, 'Estimate']),  # 转换为优势比
    OR1 = conf_int[2, 1],                   # 置信区间下限
    OR2 = conf_int[2, 2],                   # 置信区间上限
    pvalue = fit_sum$coef[2, 'Pr(>|z|)']    # p值
  )
  
  log_uni <- rbind(log_uni, result)
}

# 筛选显著性基因(p < 0.5)用于多因素分析
genes_selected <- log_uni %>% 
  dplyr::filter(pvalue < 0.5) %>% 
  .$gene

# 构建多因素Logistic回归模型
formula_mul <- as.formula(paste0('group ~ ', paste(genes_selected, collapse = ' + ')))
log_mul_fit <- glm(formula_mul, 
                   data = expr_genes, 
                   family = 'binomial')  # 多因素logistic回归

# 提取多因素回归结果
mul_sum <- summary(log_mul_fit)
df2 <- confint(log_mul_fit)  # 计算置信区间

# 准备森林图数据
estimate <- mul_sum$coefficients[-1, "Estimate"]      # 回归系数
ci_lower <- df2[-1, "2.5 %"]                         # 置信区间下限
ci_upper <- df2[-1, "97.5 %"]                        # 置信区间上限
p_values <- mul_sum$coefficients[-1, "Pr(>|z|)"]     # p值

# 创建OR和CI的文本表示
or_text <- sprintf("%.2f (%.2f-%.2f)", 
                  exp(estimate), 
                  exp(ci_lower), 
                  exp(ci_upper))

# 创建P值的文本表示
p_text <- sprintf("%.3f", p_values)

# 创建森林图数据矩阵
tabletext <- cbind(
  c("Gene", names(estimate)),
  c("OR (95% CI)", or_text),
  c("P-value", p_text)
)

# 绘制森林图
pdf(file = './output/2-Logistic_uni_forestplot.pdf', width = 10, height = max(6, length(estimate) * 0.4),onefile = F)
forestplot(
  tabletext,
  mean = c(NA, exp(estimate)),
  lower = c(NA, exp(ci_lower)),
  upper = c(NA, exp(ci_upper)),
  is.summary = c(TRUE, rep(FALSE, length(estimate))),
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(6, "mm"),  
  colgap = unit(4, "mm"),  
  col = fpColors(
    box = color_palette[1],
    line = color_palette[1],
    summary = color_palette[1]
  ),
  xlab = "Odds Ratio (95% CI)",
  title = NULL,
  xticks = seq(0, max(exp(ci_upper)) + 1, by = 1),
  clip = c(0, max(exp(ci_upper)) + 1),
  graphwidth = unit(4, "inches"),
  txt_gp = fpTxtGp(
    xlab = gpar(cex = 1),
    ticks = gpar(cex = 0.9)
  ),
  hrzl_lines = list(
    "2" = gpar(lwd = 1, col = "#444444"),
    "1" = gpar(lwd = 1, col = "#444444")
  ),
  fn.ci_norm = fpDrawNormalCI,
  ci.vertices.height = 0.1
)
dev.off()

# ==============================================================================
#
# Step 4: 构建诊断列线图（Nomogram）
# - 创建列线图用于评估个体患AD的风险
# - 列线图绘制和校准曲线分析
# 
# ==============================================================================

# 准备数据
model_data <- expr_genes
model_data$group <- as.numeric(model_data$group == "AD")

# 构建列线图
dd <- datadist(model_data)
options(datadist = "dd")

# 构建logistic回归模型
formula_str <- as.formula(paste0('group ~ ', paste(genes_selected, collapse = ' + ')))
log_fit_lrm <- lrm(formula_str,
                   data = model_data,
                   x = TRUE,
                   y = TRUE)

# 创建列线图
nom <- nomogram(log_fit_lrm,
                fun = plogis,
                funlabel = "Risk of AD",
                lp = FALSE)

# 保存列线图
pdf(width = 12, height = 8, file = "./output/3-Nomogram.pdf", onefile = FALSE)  
plot(nom, xfrac = 0.3, cex.axis = 0.8, cex.var = 0.8, lmgp = 0)  
dev.off()

# 绘制校准曲线
model_data$prob <- predict(log_mul_fit, type = 'response')

pdf(file = "./output/4-CalibrationCurve.pdf", width = 9, height = 8.5)
par(mar = c(5, 5, 4, 2))
val.prob(model_data$prob,
         model_data$group,
         cex = 1,
         xlab = "Predicted Probability",
         ylab = "Observed Probability")
title("Calibration Curve")
dev.off()

# 清理环境
options(datadist = NULL)
rm(dd)

# ==============================================================================
#
# Step 5: LASSO特征选择
# - 目的：通过LASSO回归筛选重要的预测变量（基因）
# - 使用交叉验证确定最优的惩罚参数lambda
# - 输出筛选后的基因列表
# 
# ==============================================================================

#------------------------------------
# 5.1 创建LASSO分析函数
#------------------------------------
# 函数说明：
# - 输入参数：
#   expr_data: 基因表达数据
#   group_data: 样本分组信息
#   seed: 随机种子
# - 返回值：
#   筛选到的基因列表
perform_lasso_analysis <- function(expr_data, group_data, seed = 2025) {
  set.seed(seed)
  
  # LASSO回归
  fit_lasso <- glmnet(as.matrix(expr_data),
                      group_data,
                      family = 'binomial',
                      alpha = 1)
  
  # 交叉验证
  fit_cv <- cv.glmnet(as.matrix(expr_data),
                      group_data,
                      nfolds = 10,
                      family = 'binomial')
  
  # 绘制lambda图
  pdf(file = './output/LASSO_lambda.pdf', width = 8, height = 6)
  plot(fit_lasso, xvar = 'lambda', label = TRUE)
  dev.off()
  
  # 绘制系数图
  pdf(file = './output/LASSO_coef.pdf', width = 8, height = 6)
  plot(fit_cv)
  dev.off()
  
  # 提取筛选到的基因
  mycoef <- coef(fit_cv, s = 'lambda.min')
  lasso_genes <- mycoef@Dimnames[[1]][which(mycoef != 0)]
  
  return(lasso_genes)
}

#------------------------------------
# 5.2 执行LASSO分析并保存结果
#------------------------------------
lasso_selected_genes <- perform_lasso_analysis(expr_genes, combined_groups$group)
write.csv(lasso_selected_genes, 
          "./output/selected_genes.csv", 
          row.names = FALSE)
