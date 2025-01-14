## ====================================================================
## 1.R - 乳腺癌生存分析
##
## 描述：
##   本脚本对TCGA-BRCA数据进行全面的生存分析，包括：
##   1. 差异表达分析：识别高低表达组间的差异基因
##   2. 时间依赖ROC分析：评估基因表达对生存时间的预测能力
##   3. KM生存分析：比较高低表达组的生存差异
##   4. Cox回归分析：识别影响生存的独立危险因素
##
## 数据集：
##   - TCGA-BRCA_counts.Rda: 原始counts数据，包含所有样本的基因表达量
##   - TCGA-BRCA_tpm.Rda: TPM标准化数据，用于生存分析
##   - TCGA_BRCA_cli_clean.Rda: 临床数据，包含生存时间、生存状态等信息
##
## 输出：
##   - 1-ROC.pdf: 时间依赖ROC曲线，展示1年、3年、5年生存预测能力
##   - 2-KM.pdf: KM生存曲线，比较高低表达组的生存差异
##   - 3-Cox_uni_forestplot.pdf: 单因素Cox回归森林图，展示各因素的风险比
##   - Cox_uni_res.csv: 单因素Cox回归结果，包含HR值和p值
##   - Cox_multi_res.csv: 多因素Cox回归结果，识别独立危险因素
##   - dds.Rda: DESeq2分析结果对象
##   - sig_genes.Rda: 差异表达基因列表
##
## 注意事项：
##   - 确保输入数据格式正确：
##     * counts数据：行为基因，列为样本
##     * 临床数据：包含OS.time（生存时间）和OS（生存状态）列
## ====================================================================

#####################################################
## Step1：环境初始化
## 
## 描述：
##   1. 清空当前工作环境
##   2. 加载分析所需的R包
##   3. 设置工作目录
##   4. 加载输入数据
##
## 关键参数：
##   - rm(list = ls()): 清空当前环境变量
##   - library(): 加载所需R包
##   - load(): 加载输入数据文件
##
## 注意事项：
##   - 确保所有依赖包已正确安装
##   - 检查数据文件路径是否正确
##   - 建议在分析开始前保存当前工作环境
#####################################################
rm(list = ls())
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(timeROC)
library(ggplot2)

load('input/TCGA-BRCA_counts.Rda')
load('input/TCGA-BRCA_tpm.Rda')
load('input/TCGA_BRCA_cli_clean.Rda')

#####################################################
## Step2：差异表达分析
## 
## 描述：
##   使用DESeq2包进行差异表达分析，主要步骤包括：
##   1. 构建DESeqDataSet对象
##   2. 进行标准化和差异分析
##   3. 提取显著差异基因
##   4. 基因ID转换（ENSEMBL to SYMBOL）
##   5. 保存分析结果
##
## 输入：
##   - TCGA-BRCA_counts.Rda: 原始counts数据
##   - TCGA_BRCA_cli_clean.Rda: 临床数据（用于分组）
##
## 输出：
##   - dds.Rda: DESeq2分析结果对象
##   - sig_genes.Rda: 显著差异基因列表（SYMBOL格式）
##
## 关键参数：
##   - design = ~ conditions: 实验设计公式
##   - contrast = c('conditions', c('High', 'Low')): 对比组设置
##   - logFC_cutoff = 1: 差异倍数阈值
##   - padj < 0.05: 校正p值阈值
##
## 注意事项：
##   - 确保counts数据为非负整数
##   - 检查分组信息是否正确
##   - 基因ID转换时注意处理重复ID
##   - 建议保存中间结果以便复查
#####################################################

# 2.1 数据加载与预处理
expr <- counts_data[,tumor_data_final$sample]

conditions <- data.frame(conditions = factor(tumor_data_final$group))

dds <- DESeqDataSetFromMatrix(countData = expr, colData = conditions, design = ~ conditions)
dds <- DESeq(dds)

contrast <- c('conditions', c('High', 'Low'))

res_diff <- results(dds, contrast)
summary(res_diff)

DEGs <- as.data.frame(res_diff) %>% na.omit()

logFC_cutoff <- 1

DEGs$change <- ifelse(dplyr::between(DEGs$log2FoldChange, -logFC_cutoff, logFC_cutoff) | DEGs$padj >= 0.05,
                      'Not', ifelse(DEGs$log2FoldChange >= logFC_cutoff, 'Up', 'Down'))

#提取上下调Top10差异基因
DEGs_sig <- filter(DEGs, change != 'Not')

DEGs_sig <- order(DEGs_sig$log2FoldChange, decreasing = T) %>% DEGs_sig[., ] %>%
  .[c(1:10, (nrow(DEGs_sig) -9):(nrow(DEGs_sig))), ]

#转换基因id
genelist <- bitr(rownames(DEGs_sig), fromType = 'ENSEMBL', toType = 'SYMBOL', 
                 OrgDb = 'org.Hs.eg.db')
genelist <- genelist[!duplicated(genelist$ENSEMBL), ]

sig_genes <- genelist$SYMBOL
save(sig_genes, file = 'output/sig_genes.Rda')

#####################################################
## Step3：时间依赖ROC分析
## 
## 描述：
##   使用timeROC包进行时间依赖ROC分析，评估基因表达对生存时间的预测能力，
##   主要步骤包括：
##   1. 数据准备：整合表达数据和临床数据
##   2. 计算1年、3年、5年ROC曲线
##   3. 绘制ROC曲线图
##   4. 计算AUC值及其置信区间
##
## 输入：
##   - TCGA-BRCA_tpm.Rda: TPM标准化数据
##   - TCGA_BRCA_cli_clean.Rda: 临床数据
##
## 输出：
##   - 1-ROC.pdf: 时间依赖ROC曲线图，包含：
##     * 1年、3年、5年ROC曲线
##     * 各时间点AUC值及95%置信区间
##
## 关键参数：
##   - times = c(1,3,5) * 365.25: 评估时间点（转换为天数）
##   - cause = 1: 主要事件类型
##   - weighting = 'marginal': 权重计算方法
##
## 注意事项：
##   - 确保生存时间单位为天
##   - 检查事件类型编码（1=事件发生，0=删失）
##   - 建议保存原始ROC计算结果以便后续分析
#####################################################

# 3.1 数据准备
rm(list = ls())
library(tidyverse)
library(survival)
library(timeROC)
library(ggplot2)
load('input/TCGA-BRCA_tpm.Rda')
load('input/TCGA_BRCA_cli_clean.Rda')

# 整理分组数据
dat_expr_OLR1 <- tumor_data_final[,c(5,4,2)]

# ROC分析
time_roc <- timeROC(T = dat_expr_OLR1$OS.time, delta = dat_expr_OLR1$OS, marker = dat_expr_OLR1$OLR1,
                    cause = 1, weighting = 'marginal', times = c(1,3,5) * 365.25,
                    ROC = T,
                    iid = T)

time_roc_df <- data.frame(
  TP_1 = time_roc$TP[, 1],
  FP_1 = time_roc$FP[, 1],
  TP_3 = time_roc$TP[, 2],
  FP_3 = time_roc$FP[, 2],
  TP_5 = time_roc$TP[, 3],
  FP_5 = time_roc$FP[, 3]
)

pdf("output/1-ROC.pdf", width = 8, height = 8)
ggplot(time_roc_df) +
  geom_line(aes(FP_1, TP_1), lwd = 1, color = '#4DBBD5') +
  geom_line(aes(FP_3, TP_3), lwd = 1, color = '#E64B35') +
  geom_line(aes(FP_5, TP_5), lwd = 1, color = '#00A087') +
  geom_abline(slope = 1, intercept = 0, lty = 'dashed') +
  labs(subtitle = 'OLR1', 
       x = '1 - Specificities (FPR)',
       y = 'Sensitivities (TPR)') +
  coord_fixed() +
  theme_bw() +
  annotate(
    'text',
    label = paste0('AUC at 1 year = ', round(time_roc$AUC[[1]], 3), ' ( ',
                   round(confint(time_roc, level = 0.95)$CI_AUC[1, 1] / 100, 3), ' - ',
                   round(confint(time_roc, level = 0.95)$CI_AUC[1, 2] / 100, 3), ' )'),
    x = 0.75, y = 0.1, color = '#4DBBD5'
  ) +
  annotate(
    'text',
    label = paste0('AUC at 3 year = ', round(time_roc$AUC[[2]], 3), ' ( ',
                   round(confint(time_roc, level = 0.95)$CI_AUC[2, 1] / 100, 3), ' - ',
                   round(confint(time_roc, level = 0.95)$CI_AUC[2, 2] / 100, 3), ' )'),
    x = 0.75, y = 0.05, color = '#E64B35'
  ) +
  annotate(
    'text',
    label = paste0('AUC at 5 year = ', round(time_roc$AUC[[3]], 3), ' ( ',
                   round(confint(time_roc, level = 0.95)$CI_AUC[3, 1] / 100, 3), ' - ',
                   round(confint(time_roc, level = 0.95)$CI_AUC[3, 2] / 100, 3), ' )'),
    x = 0.75, y = 0, color = '#00A087'
  )
dev.off()

#####################################################
## Step4：KM生存分析
## 
## 描述：
##   使用survival和survminer包进行Kaplan-Meier生存分析，主要步骤包括：
##   1. 数据准备：整合表达数据和临床数据
##   2. 根据基因表达水平分组（高/低表达）
##   3. 拟合生存曲线
##   4. 绘制KM曲线图
##   5. 计算log-rank检验p值
##
## 输入：
##   - TCGA-BRCA_tpm.Rda: TPM标准化数据
##   - TCGA_BRCA_cli_clean.Rda: 临床数据
##
## 输出：
##   - 2-KM.pdf: KM生存曲线图，包含：
##     * 高低表达组生存曲线
##     * 风险表
##     * log-rank检验p值
##
## 关键参数：
##   - palette = c('#4DBBD5', '#E64B35'): 曲线颜色
##   - pval = T: 显示log-rank检验p值
##   - conf.int = T: 显示置信区间
##   - risk.table = T: 显示风险表
##
## 注意事项：
##   - 确保分组标准一致
##   - 检查生存时间单位
##   - 建议保存生存分析结果对象以便后续分析
#####################################################

# 4.1 数据准备
rm(list = ls())
library(survival)
library(survminer)
load('input/TCGA-BRCA_tpm.Rda')
load('input/TCGA_BRCA_cli_clean.Rda')


#整理数据
dat_expr_OLR1 <- tumor_data_final[,c(5,4,3)]

fit_surv_OLR1 <- surv_fit(Surv(OS.time, OS) ~ group, data = dat_expr_OLR1)

pdf("output/2-KM.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(
  fit_surv_OLR1,
  data = dat_expr_OLR1,
  pval = T,
  linetype = 'solid',
  palette = c('#4DBBD5', '#E64B35'),
  title = 'OLR1',
  legend.title = 'OLR1',
  legend = c(0.7, 0.9),
  legend.labs = c('High', 'Low'),
  conf.int = T,
  conf.int.style = 'ribbon',
  conf.int.alpha = 0.1,
  risk.table = T
)
dev.off()

#####################################################
## Step5：Cox回归分析
## 
## 描述：
##   使用survival包进行Cox比例风险回归分析，包括：
##   1. 单因素Cox回归：
##     * 对每个变量单独进行回归分析
##     * 计算风险比（HR）和p值
##   2. 多因素Cox回归：
##     * 对单因素分析中显著的变量进行多因素分析
##     * 识别独立危险因素
##   3. 结果可视化：
##     * 绘制森林图展示单因素分析结果
##
## 输入：
##   - TCGA-BRCA_tpm.Rda: TPM标准化数据
##   - TCGA_BRCA_cli_clean.Rda: 临床数据
##   - sig_genes.Rda: 差异表达基因
##
## 输出：
##   - 3-Cox_uni_forestplot.pdf: 单因素Cox回归森林图
##   - Cox_uni_res.csv: 单因素Cox回归结果，包含：
##     * 变量名称
##     * HR值及95%置信区间
##     * p值
##   - Cox_multi_res.csv: 多因素Cox回归结果，包含：
##     * 独立危险因素
##     * 调整后的HR值及95%置信区间
##     * p值
##
## 关键参数：
##   - pvalue_cutoff = 0.1: 单因素到多因素筛选阈值
##   - clip = c(0, 5): 森林图HR值显示范围
##   - boxsize = 0.15: 森林图点大小
##
## 注意事项：
##   - 检查变量类型（连续/分类）
##   - 处理分类变量时注意参考水平设置
##   - 检查比例风险假设
##   - 建议保存完整的Cox模型对象以便后续分析
#####################################################

# 5.1 单因素Cox回归
rm(list = ls())
library(tidyverse)
library(survival)
library(rms)
library(forestplot)
library(ggDCA)
library(clusterProfiler)
library(org.Hs.eg.db)

load('input/TCGA-BRCA_tpm.Rda')
load('output/sig_genes.Rda')
load('input/TCGA_BRCA_cli_clean.Rda')

#转换基因id
dat_expr <- tpm_data
dat_expr$ENSEMBL<-rownames(dat_expr)
geneid<-bitr(dat_expr$ENSEMBL,
             fromType = "ENSEMBL",
             toType = "SYMBOL",
             OrgDb = "org.Hs.eg.db")
dat_expr<-left_join(geneid,dat_expr,by = "ENSEMBL") %>% .[,-1]

#对Symbol去重
dat_expr<-dat_expr[!duplicated(dat_expr$SYMBOL),]
rownames(dat_expr)<-dat_expr$SYMBOL
dat_expr<-dat_expr[,-1]

#整理临床数据
rownames(tumor_data_final) <- tumor_data_final$sample
dat_cli <- tumor_data_final[, -c(1:3,11)]
dat_expr <- dat_expr[, rownames(dat_cli)]

# 整理数据
sig_genes <- c(sig_genes, 'OLR1')
expr_genes <- dat_expr[sig_genes, ] %>% na.omit() %>% t() %>% as.data.frame()
expr_genes_cli <- cbind(dat_cli, expr_genes)

expr_genes_cli$pathologic_stage <- ifelse(
  expr_genes_cli$pathologic_stage %in% c('Stage_I', 'Stage_II'), 'Stage I/II',
  ifelse(expr_genes_cli$pathologic_stage %in% c('Stage_III', 'Stage_IV'), 'Stage III/IV', NA))

expr_genes_cli$T_stage <- ifelse(
  expr_genes_cli$T_stage %in% c('T1', 'T2'), 'T1/T2',
  ifelse(expr_genes_cli$T_stage %in% c('T3', 'T4'), 'T3/T4', NA))

expr_genes_cli$N_stage <- ifelse(
  expr_genes_cli$N_stage %in% c('N0', 'N1'), 'N0/N1',
  ifelse(expr_genes_cli$N_stage %in% c('N2', 'N3'), 'N2/N3', NA))

expr_genes_cli$pathologic_stage <- factor(expr_genes_cli$pathologic_stage, levels = c('Stage I/II', 'Stage III/IV'))
expr_genes_cli$T_stage <- factor(expr_genes_cli$T_stage, levels = c('T1/T2', 'T3/T4'))
expr_genes_cli$N_stage <- factor(expr_genes_cli$N_stage, levels = c('N0/N1', 'N2/N3'))
expr_genes_cli$M_stage <- factor(expr_genes_cli$M_stage, levels = c('M0', 'M1'))
expr_genes_cli <- expr_genes_cli[,-c(8,14)]


# 单因素Cox回归
# 初始化单因素Cox回归列表
cox_uni <- data.frame()
# 批量对每一个基因或临床变量做单因素Cox回归
for(i in 3:ncol(expr_genes_cli)) {
  cox <- coxph(as.formula(paste0('Surv(OS.time, OS) ~ `', colnames(expr_genes_cli)[i], '`')), data = expr_genes_cli)
  cox_sum <- summary(cox)
  # 整理结果
  dat_res0 <- cbind(
    characteristics = colnames(expr_genes_cli)[i],
    HR = cox_sum$conf.int[, 'exp(coef)'],
    HR1 = cox_sum$conf.int[, 'lower .95'],
    HR2 = cox_sum$conf.int[, 'upper .95'],
    pvalue = cox_sum$coefficients[, 'Pr(>|z|)']
  )
  cox_uni <- rbind(cox_uni, dat_res0)
  # 移除中间变量
  rm(i, cox, cox_sum, dat_res0)
}


for (i in 2:ncol(cox_uni)) {
  cox_uni[,i] <- as.numeric(cox_uni[,i])
}

# 输出单因素Cox结果
cox_uni_output <- data.frame(
  characteristic = cox_uni$characteristics,
  HR = paste0(
    round(cox_uni$HR, 2),
    ' (',
    round(cox_uni$HR1, 2),
    '-',
    round(cox_uni$HR2, 2),
    ')'
  ),
  pvalue = signif(cox_uni$pvalue, 3)
)

write.csv(cox_uni_output, file = 'output/Cox_uni_res.csv', row.names = F)

# 5.2 多因素Cox回归

# 单因素进行多因素的p值阈值可以放宽到0.1
pvalue_cutoff <- 0.1
characteristic_CCRDEGs_cox <- cox_uni %>%
  dplyr::filter(pvalue < pvalue_cutoff) %>%
  .$characteristic

cox_mul_fit <- coxph(as.formula(paste0('Surv(OS.time,OS) ~ ', paste(characteristic_CCRDEGs_cox, collapse = ' + '))),
                     data = expr_genes_cli)

cox_mul_sum <- summary(cox_mul_fit)

cox_mul <- cbind(
  characteristic = characteristic_CCRDEGs_cox,
  HR = cox_mul_sum$conf.int[, 'exp(coef)'],
  HR1 = cox_mul_sum$conf.int[, 'lower .95'],
  HR2 = cox_mul_sum$conf.int[, 'upper .95'],
  pvalue = cox_mul_sum$coefficients[, 'Pr(>|z|)']
) %>% as.data.frame()

for (i in 2:ncol(cox_mul)) {
  cox_mul[,i] <- as.numeric(cox_mul[,i])
}

# 输出多因素Cox结果
cox_mul_output <- data.frame(
  characteristic = cox_mul$characteristic,
  HR = paste0(
    round(cox_mul$HR, 2),
    ' (',
    round(cox_mul$HR1, 2),
    '-',
    round(cox_mul$HR2, 2),
    ')'
  ),
  pvalue = signif(cox_mul$pvalue, 3)
)
#cox_mul_output$characteristic[1:3] <- c('T Stage', 'N Stage', 'M Stage')
write.csv(cox_mul_output, file = 'output/Cox_multi_res.csv', row.names = F)

# 5.3 结果可视化 - 森林图

# 森林图数据
dat_forest_cox_mul <- rbind(colnames(cox_mul_output), cox_mul_output)

# 森林图
pdf(file = 'output/3-Cox_uni_forestplot.pdf', width = 8, height = 6,
    # 不加这一句的话pdf会出两页
    onefile = F
)
forestplot(
  labeltext = as.matrix(dat_forest_cox_mul),
  mean = c(NA, cox_mul$HR),
  lower = c(NA, cox_mul$HR1),
  upper = c(NA, cox_mul$HR2),
  ci.vertices = T,
  # 森林图在第几列``
  graph.pos = 3,
  # 点型和颜色
  fn.ci_norm = fpDrawCircleCI,
  col = fpColors(
    box = '#4DBBD5',
    lines = 'black',
    zero = 'black'
  ),
  clip = c(0, 5),
  boxsize = 0.15,
  zero = 1,
  xlab = 'Hazard Ratio',
  # 三线表的三条线
  hrzl_lines = rlang::set_names(list(
    gpar(col = 'black'), gpar(col = 'black'), gpar(col = 'black')
    # 线分别在第1行，第2行和最后一行
  ), c(1, 2, nrow(cox_mul) + 2))
)
dev.off()
