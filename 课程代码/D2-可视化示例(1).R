# ============================================================================
# DAY 2. R语言数据可视化
# ============================================================================
# 加载必要的包并准备数据
# ============================ 1. 环境准备 ============================

# 加载必要的包
library(ggplot2)
library(dplyr)
library(viridis)
library(hrbrthemes)


#设置工作路径
setwd("./D2-可视化示例")

# 创建模拟数据：四个城市近年来的温度变化
set.seed(123)
cities <- c("北京", "上海", "广州", "成都")
years <- 2015:2023

climate_data <- expand.grid(
  city = cities,
  year = years
) %>%
  as.data.frame() %>%
  mutate(
    temperature = case_when(
      city == "北京" ~ rnorm(n(), mean = 12, sd = 2),
      city == "上海" ~ rnorm(n(), mean = 16, sd = 1.5),
      city == "广州" ~ rnorm(n(), mean = 22, sd = 1),
      city == "成都" ~ rnorm(n(), mean = 16, sd = 1.8)
    ),
    rainfall = case_when(
      city == "北京" ~ rnorm(n(), mean = 600, sd = 50),
      city == "上海" ~ rnorm(n(), mean = 1200, sd = 100),
      city == "广州" ~ rnorm(n(), mean = 1800, sd = 150),
      city == "成都" ~ rnorm(n(), mean = 900, sd = 80)
    )
  )

# ============================ 2. 创建ggplot2图片 ============================

# 创建折线图
p1 <- ggplot(climate_data, aes(x = year, y = temperature, color = city, group = city)) +
  # 添加平滑线
  geom_smooth(aes(fill = city), alpha = 0.2, size = 1.2) +
  # 添加数据点
  geom_point(size = 3, alpha = 0.6) +
  # 使用viridis配色方案
  scale_color_viridis(discrete = TRUE, option = "D", 
                      labels = c("Beijing", "Shanghai", "Guangzhou", "Chengdu")) +
  scale_fill_viridis(discrete = TRUE, option = "D",
                     labels = c("Beijing", "Shanghai", "Guangzhou", "Chengdu")) +
  # 使用现代化主题
  theme_ipsum_rc() +
  # 自定义主题元素
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    legend.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) +
  # 添加标签
  labs(
    title = "Temperature Trends in Major Chinese Cities",
    subtitle = "Temperature Changes from 2015 to 2023 (Simulated Data)",
    x = "Year",
    y = "Temperature (°C)",
    color = "City",
    fill = "City"
  )

# 创建气泡图
p2 <- ggplot(climate_data, 
       aes(x = temperature, y = rainfall, size = year, color = city)) +
  # 添加气泡
  geom_point(alpha = 0.7) +
  # 使用viridis配色
  scale_color_viridis(discrete = TRUE, option = "C",
                      labels = c("Beijing", "Shanghai", "Guangzhou", "Chengdu")) +
  scale_size_continuous(range = c(3, 12)) +
  # 使用现代化主题
  theme_ipsum_rc() +
  # 自定义主题
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) +
  # 添加标签
  labs(
    title = "Temperature-Rainfall Relationship",
    subtitle = "Bubble size represents year",
    x = "Temperature (°C)",
    y = "Annual Rainfall (mm)",
    color = "City",
    size = "Year"
  )

# 创建箱线图
p3 <- ggplot(climate_data, aes(x = city, y = temperature, fill = city)) +
  # 添加箱线图
  geom_boxplot(alpha = 0.7, width = 0.5) +
  # 添加抖动的数据点
  geom_jitter(color = "grey40", size = 1, alpha = 0.4, width = 0.2) +
  # 使用viridis配色
  scale_fill_viridis(discrete = TRUE, option = "E",
                     labels = c("Beijing", "Shanghai", "Guangzhou", "Chengdu")) +
  # 修改x轴标签
  scale_x_discrete(labels = c("Beijing", "Shanghai", "Guangzhou", "Chengdu")) +
  # 使用现代化主题
  theme_ipsum_rc() +
  # 自定义主题
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) +
  # 添加标签
  labs(
    title = "Temperature Distribution by City",
    subtitle = "Boxplot showing temperature distribution characteristics",
    y = "Temperature (°C)"
  )

# ============================ 3. 图片预览和保存 ============================

# 预览图表
print(p1)
print(p2)
print(p3)

# 保存图表
ggsave("图1.png", p1, width = 12, height = 8, dpi = 300)
ggsave("图2.png", p2, width = 12, height = 8, dpi = 300)
ggsave("图3.png", p3, width = 10, height = 7, dpi = 300)

