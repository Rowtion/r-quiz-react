## ====================================================================
## DAY 4. R语言临床数据整合
## 
## 功能：整理TCGA临床数据，并按照OLR1单基因分高低表达组
## 
## 输入:
##   - TCGA-BRCA_tpm.Rda: TCGA-BRCA TPM表达矩阵
##   - TCGA-BRCA_cli.Rda: TCGA-BRCA 临床数据
## 输出:
##   - TCGA_BRCA_cli_clean.Rda: TCGA-BRCA 整理好的临床数据
##   - TCGA_cluster.Rda: OLR1高低表达组分组信息
## ====================================================================

#--------------------1. 数据准备--------------------
# 加载必要的包
library(tidyverse)
library(ggpubr)
library(stringr)
library(UCSCXenaTools)
library(ggplot2)
library(tidyr)

# 加载数据
setwd("./D4-TCGA临床信息整理")
load("TCGA-BRCA_tpm.Rda")
load("TCGA-BRCA_cli.Rda")

#--------------------2. 样本筛选--------------------
# 处理样本分组信息（正常vs肿瘤）
group_TCGA <- data.frame(
    sample_id = colnames(tpm_data),
    group = ifelse(substr(colnames(tpm_data), 14, 14) == "1", "Normal", "Tumor"),
    row.names = colnames(tpm_data)
)

#--------------------3. OLR1表达量提取--------------------
# 提取OLR1表达数据（包含所有样本）
OLR1_exp <- data.frame(
      sample = colnames(tpm_data),
      OLR1 = as.numeric(tpm_data["ENSG00000173391",]),
      group = group_TCGA$group
  )

#--------------------4. 生存数据获取--------------------
# 从XENA下载TCGA-BRCA生存数据
xe <- XenaGenerate(subset = XenaCohorts == "GDC TCGA Breast Cancer (BRCA)") %>% 
    XenaFilter(filterDatasets = "survival") %>% 
    XenaQuery()
XenaDownload(xe, destdir = './')

# 处理生存数据
survival_clean <- read.table("TCGA-BRCA.survival.tsv.gz", header = TRUE, sep = "\t") %>% 
    dplyr::select(sample = sample, OS.time, OS)

#合并OLR1表达数据和生存数据
OLR1_exp_survival <- OLR1_exp %>% 
    mutate(sample_short = substr(sample, 1, 16)) %>%
    left_join(survival_clean, by = c("sample_short" = "sample")) %>% 
    filter(!is.na(OS.time), !is.na(OS)) %>%
    dplyr::select(-sample_short)

#--------------------5. 临床信息整合--------------------
# 处理临床信息
clinical_clean <- cli_data %>%
    # 提取基本临床信息
    dplyr::select(
        sample = bcr_patient_barcode,
        age = age_at_initial_pathologic_diagnosis,
        pathologic_stage = stage_event_pathologic_stage,
        tnm = stage_event_tnm_categories
    ) %>% 
    mutate(
        sample = substr(sample, 1, 12),
        # 处理病理分期，去掉最后的字母
        pathologic_stage = str_extract(pathologic_stage, "Stage [I|V|X]+") %>% 
            str_replace("Stage ", "Stage_"),
        T_stage = str_extract(tnm, "T\\d[a-z]?") %>% str_replace("[a-z]$", ""),
        N_stage = str_extract(tnm, "N\\d[a-z]?") %>% str_replace("[a-z]$", ""),
        M_stage = str_extract(tnm, "M\\d[a-z]?") %>% str_replace("[a-z]$", ""),
        # 添加年龄分组
        age_group = case_when(
            age <= 60 ~ "≤60",
            age > 60 ~ ">60"
        )
    ) %>% 
  filter(
    !is.na(age), 
    !is.na(pathologic_stage), 
    !is.na(T_stage), !is.na(N_stage), !is.na(M_stage)
  ) %>% 
    dplyr::select(-tnm)

# 分别处理正常组和肿瘤组
# 正常组数据
normal_data <- OLR1_exp_survival %>% 
    filter(group == "Normal") %>% 
    dplyr::select(sample, OLR1, group)

# 肿瘤组数据（合并临床信息）
tumor_data <- OLR1_exp_survival %>% 
    filter(group == "Tumor") %>% 
    mutate(sample_short = substr(sample, 1, 12)) %>%
    left_join(clinical_clean, by = c("sample_short" = "sample")) %>% 
    filter(
        !is.na(age), 
        !is.na(pathologic_stage), 
        !is.na(T_stage), !is.na(N_stage), !is.na(M_stage)
    ) %>%
    dplyr::select(-sample_short)

#--------------------6. 定义绘制箱型图的函数--------------------
# 自定义主题
my_theme <- function() {
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "right"
    )
}

# 绘制箱型图的函数
plot_boxplot <- function(data, group_var, title = NULL) {
    ggplot(data, aes(x = !!sym(group_var), y = OLR1, fill = !!sym(group_var))) +
        geom_boxplot(
            alpha = 0.7,
            outlier.shape = 16,
            outlier.size = 0.5,
            width = 0.7
        ) +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
        labs(
            title = title,
            x = group_var,
            y = "OLR1 Expression (log2 TPM)"
        ) +
        my_theme() +
        scale_fill_brewer(palette = "Set2")
}

#--------------------7. 绘制各组箱型图--------------------
# 检查数据
print("正常样本数量：")
print(sum(OLR1_exp$group == "Normal"))
print("肿瘤样本数量：")
print(sum(OLR1_exp$group == "Tumor"))

# 准备绘图数据
# 获取正常样本数据
normal_data_for_plot <- OLR1_exp %>% 
    filter(group == "Normal") %>% 
    mutate(
        sample = sample,
        OLR1 = OLR1,
        age_group = "Normal",
        pathologic_stage = "Normal",
        T_stage = "Normal",
        N_stage = "Normal",
        M_stage = "Normal"
    )

# 准备肿瘤样本数据
tumor_data_for_plot <- tumor_data %>% 
    dplyr::select(sample, OLR1, age_group, pathologic_stage, T_stage, N_stage, M_stage)

# 合并正常和肿瘤样本数据
plot_data <- bind_rows(
    normal_data_for_plot,
    tumor_data_for_plot
)

# 设置因子水平，确保Normal总是第一个
plot_data <- plot_data %>% 
    mutate(
        age_group = factor(age_group, levels = c("Normal", "≤60", ">60")),
        pathologic_stage = factor(pathologic_stage, levels = c("Normal", sort(unique(tumor_data_for_plot$pathologic_stage)))),
        T_stage = factor(T_stage, levels = c("Normal", sort(unique(tumor_data_for_plot$T_stage)))),
        N_stage = factor(N_stage, levels = c("Normal", sort(unique(tumor_data_for_plot$N_stage)))),
        M_stage = factor(M_stage, levels = c("Normal", sort(unique(tumor_data_for_plot$M_stage))))
    )

# 1. 正常vs肿瘤
normal_tumor_plot <- ggplot(OLR1_exp, aes(x = group, y = OLR1, fill = group)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.shape = 16,
        outlier.size = 0.5,
        width = 0.7
    ) +
    labs(
        title = "Normal vs Tumor",
        x = "Group",
        y = "OLR1 Expression (log2 TPM)"
    ) +
    stat_compare_means(method = "wilcox.test",
                      comparisons = list(c("Normal", "Tumor")),
                      label = "p.signif",
                      label.y = max(OLR1_exp$OLR1) * 1.2) +
    my_theme() +
    scale_fill_brewer(palette = "Set2")

# 2. 年龄分组
age_plot <- ggplot(plot_data, aes(x = age_group, y = OLR1, fill = age_group)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.shape = 16,
        outlier.size = 0.5,
        width = 0.7
    ) +
    labs(
        title = "Age Groups",
        x = "Age Group",
        y = "OLR1 Expression (log2 TPM)"
    ) +
    stat_compare_means(method = "wilcox.test",
                      comparisons = list(
                          c("Normal", "≤60"),
                          c("Normal", ">60"),
                          c("≤60", ">60")
                      ),
                      label = "p.signif",
                      label.y = c(max(plot_data$OLR1) * 1.2,
                                max(plot_data$OLR1) * 1.3,
                                max(plot_data$OLR1) * 1.1)) +
    my_theme() +
    scale_fill_brewer(palette = "Set2")

# 3. 病理分期
stage_plot <- ggplot(plot_data, aes(x = pathologic_stage, y = OLR1, fill = pathologic_stage)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.shape = 16,
        outlier.size = 0.5,
        width = 0.7
    ) +
    labs(
        title = "Pathologic Stage",
        x = "Stage",
        y = "OLR1 Expression (log2 TPM)"
    ) +
    stat_compare_means(method = "wilcox.test",
                      comparisons = list(
                          c("Normal", "Stage I"),
                          c("Normal", "Stage II"),
                          c("Normal", "Stage III"),
                          c("Normal", "Stage IV"),
                          c("Stage I", "Stage II"),
                          c("Stage I", "Stage III"),
                          c("Stage I", "Stage IV"),
                          c("Stage II", "Stage III"),
                          c("Stage II", "Stage IV"),
                          c("Stage III", "Stage IV")
                      ),
                      label = "p.signif",
                      label.y = seq(max(plot_data$OLR1) * 1.1,
                                  max(plot_data$OLR1) * 2.0,
                                  length.out = 10)) +
    my_theme() +
    scale_fill_brewer(palette = "Set2")

# 4. T分期
t_plot <- ggplot(plot_data, aes(x = T_stage, y = OLR1, fill = T_stage)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.shape = 16,
        outlier.size = 0.5,
        width = 0.7
    ) +
    labs(
        title = "T Stage",
        x = "T Stage",
        y = "OLR1 Expression (log2 TPM)"
    ) +
    stat_compare_means(method = "wilcox.test",
                      comparisons = list(
                          c("Normal", "T1"),
                          c("Normal", "T2"),
                          c("Normal", "T3"),
                          c("Normal", "T4"),
                          c("T1", "T2"),
                          c("T1", "T3"),
                          c("T1", "T4"),
                          c("T2", "T3"),
                          c("T2", "T4"),
                          c("T3", "T4")
                      ),
                      label = "p.signif",
                      label.y = seq(max(plot_data$OLR1) * 1.1,
                                  max(plot_data$OLR1) * 2.0,
                                  length.out = 10)) +
    my_theme() +
    scale_fill_brewer(palette = "Set2")

# 5. N分期
n_plot <- ggplot(plot_data, aes(x = N_stage, y = OLR1, fill = N_stage)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.shape = 16,
        outlier.size = 0.5,
        width = 0.7
    ) +
    labs(
        title = "N Stage",
        x = "N Stage",
        y = "OLR1 Expression (log2 TPM)"
    ) +
    stat_compare_means(method = "wilcox.test",
                      comparisons = list(
                          c("Normal", "N0"),
                          c("Normal", "N1"),
                          c("Normal", "N2"),
                          c("Normal", "N3"),
                          c("N0", "N1"),
                          c("N0", "N2"),
                          c("N0", "N3"),
                          c("N1", "N2"),
                          c("N1", "N3"),
                          c("N2", "N3")
                      ),
                      label = "p.signif",
                      label.y = seq(max(plot_data$OLR1) * 1.1,
                                  max(plot_data$OLR1) * 2.0,
                                  length.out = 10)) +
    my_theme() +
    scale_fill_brewer(palette = "Set2")

# 6. M分期
m_plot <- ggplot(plot_data, aes(x = M_stage, y = OLR1, fill = M_stage)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.shape = 16,
        outlier.size = 0.5,
        width = 0.7
    ) +
    labs(
        title = "M Stage",
        x = "M Stage",
        y = "OLR1 Expression (log2 TPM)"
    ) +
    stat_compare_means(method = "wilcox.test",
                      comparisons = list(
                          c("Normal", "M0"),
                          c("Normal", "M1"),
                          c("M0", "M1")
                      ),
                      label = "p.signif",
                      label.y = c(max(plot_data$OLR1) * 1.2,
                                max(plot_data$OLR1) * 1.4,
                                max(plot_data$OLR1) * 1.6)) +
    my_theme() +
    scale_fill_brewer(palette = "Set2")

# 保存单个图片
ggsave("plots/normal_tumor_plot.pdf", normal_tumor_plot, width = 6, height = 5)
ggsave("plots/age_plot.pdf", age_plot, width = 6, height = 5)
ggsave("plots/stage_plot.pdf", stage_plot, width = 7, height = 5)
ggsave("plots/t_plot.pdf", t_plot, width = 6, height = 5)
ggsave("plots/n_plot.pdf", n_plot, width = 6, height = 5)
ggsave("plots/m_plot.pdf", m_plot, width = 6, height = 5)

# 使用patchwork组合图片
library(patchwork)

# 创建2x3布局的组合图
combined_plot <- (normal_tumor_plot + age_plot + stage_plot) /
                 (t_plot + n_plot + m_plot) +
                 plot_annotation(
                     title = "OLR1 Expression Across Different Clinical Features",
                     theme = theme(
                         plot.title = element_text(size = 16, hjust = 0.5)
                     )
                 )

# 保存组合图
ggsave("plots/combined_plot.pdf", combined_plot, width = 18, height = 12)

# 保存tumor_data
tumor_data_final <- tumor_data %>% 
    mutate(group = ifelse(OLR1 > median(OLR1), "High", "Low")) %>% 
    mutate(group = factor(group, levels = c("Low", "High"))) %>%
    arrange(group)

# 保存为Rda文件
save(tumor_data_final, file = "TCGA_BRCA_cli_clean.Rda")
