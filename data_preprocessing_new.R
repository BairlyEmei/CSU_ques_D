library(readxl)
library(dplyr)
library(writexl)
library(entropy)

# 读取数据
data <- read_excel("D题附件：法医物证多人身份鉴定问题数据集/附件2：不同混合比例的STR图谱数据_processed.xlsx")

# 构建增强的特征列名
enhanced_char <- c("gene_sum","height_sum","height_mean","height_std","height_max",
                   "height_cv", "skewness", "kurtosis", "effective_peaks",
                   "triplet_ratio", "phr", "low_phr", "high_freq_cluster",
                   "height_main_ratio", "peak_entropy")  # 新增特征

marker <- c('D8S1179', 'D21S11', 'D7S820', 'CSF1PO', 'D3S1358', 'TH01', 'D13S317',
            'D16S539', 'D2S1338', 'D19S433', 'vWA', 'TPOX', 'D18S51', 'AMEL', 'D5S818', 'FGA')

marker_enhanced <- paste(rep(marker, each=length(enhanced_char)), enhanced_char, sep="_")
global_features <- c("prop_max", "prop_min", "prop_sd", "prop_entropy",
                     "avg_gene_sum", "sd_effective_peaks", "avg_phr", "var_phr",
                     "total_low_phr", "avg_main_ratio")  # 全局特征

colnames <- c("sample_file","gene_from","people","proportion", "gene_total","height_total",
              "cv_std", global_features, marker_enhanced)

# 构建数据框
df <- data.frame(matrix(ncol=length(colnames),nrow=50))
colnames(df) <- colnames

# 填充基础信息
unique_data <- data %>%
  group_by(`Sample File`) %>%
  summarise(
    gene_from = first(`基因来源（贡献者）`),
    people = first(`人数`),
    proportion = first(`比例`)
  )

df$sample_file <- unique_data$`Sample File`
df$gene_from <- unique_data$gene_from
df$people <- unique_data$people
df$proportion <- unique_data$proportion

# 混合比例解析函数
parse_proportion <- function(prop_str) {
  if (is.na(prop_str)) return(list(prop_max=NA, prop_min=NA, prop_sd=NA, prop_entropy=NA))
  ratios <- as.numeric(strsplit(prop_str, ":")[[1]])
  total <- sum(ratios)
  proportions <- ratios / total
  return(list(
    prop_max = max(proportions),
    prop_min = min(proportions),
    prop_sd = ifelse(length(proportions)>=2, sd(proportions), 0),
    prop_entropy = entropy(proportions)
  ))
}

# 计算近似熵的函数
calculate_approx_entropy <- function(max_height, total_height) {
  p <- max_height / (total_height + 1e-6)
  q <- 1 - p
  entropy <- - (p * log2(p + 1e-6) + q * log2(q + 1e-6))
  return(ifelse(is.finite(entropy), entropy, NA))
}

# 特征计算主循环
height_cols <- paste("Height", 1:100, sep = " ")

for (i in 1:nrow(df)) {
  current_sample <- df$sample_file[i]
  sample_data <- data %>% filter(`Sample File` == current_sample)

  # 解析混合比例特征
  prop_features <- parse_proportion(df$proportion[i])
  df[i, c("prop_max", "prop_min", "prop_sd", "prop_entropy")] <- prop_features

  # 初始化全局特征
  gene_sums <- numeric(length(marker))
  effective_peaks <- numeric(length(marker))
  phr_values <- numeric(length(marker))
  main_ratios <- numeric(length(marker))

  for (m_idx in seq_along(marker)) {
    m <- marker[m_idx]
    marker_data <- sample_data %>% filter(Marker == m)

    if (nrow(marker_data) == 0) next

    heights <- as.numeric(marker_data[1, height_cols])
    valid_heights <- heights[!is.na(heights) & heights > 0]

    # 动态噪声过滤：结合最大峰高和总峰高
    if (length(valid_heights) > 0) {
      total_height <- sum(valid_heights)
      threshold <- max(0.1*max(valid_heights), 0.01*total_height)
      valid_heights <- valid_heights[valid_heights > threshold]
    }

    # 计算基本特征
    gene_sum <- length(valid_heights)
    height_sum <- sum(valid_heights)
    height_mean <- ifelse(gene_sum > 0, mean(valid_heights), NA)
    height_std <- ifelse(gene_sum > 1, sd(valid_heights), NA)
    height_max <- ifelse(gene_sum > 0, max(valid_heights), NA)
    height_cv <- ifelse(gene_sum > 0 & height_mean > 0, height_std/height_mean, NA)

    # 新增特征计算
    sorted_heights <- sort(valid_heights, decreasing = TRUE)
    phr <- ifelse(length(sorted_heights) >= 2, sorted_heights[2]/sorted_heights[1], NA)
    height_main_ratio <- ifelse(height_sum > 0, height_max / height_sum, NA)
    peak_entropy <- ifelse(gene_sum > 0, entropy::entropy(valid_heights), NA)

    # 保存到全局特征数组
    gene_sums[m_idx] <- gene_sum
    effective_peaks[m_idx] <- sum(valid_heights > 0.05 * height_sum)
    phr_values[m_idx] <- phr
    main_ratios[m_idx] <- height_main_ratio

    approx_entropy <- calculate_approx_entropy(height_max, height_sum)
    df[i, paste(m, "approx_entropy", sep="_")] <- round(approx_entropy, 4)

    # 填充标记级特征
    df[i, paste(m, "gene_sum", sep="_")] <- gene_sum
    df[i, paste(m, "height_sum", sep="_")] <- height_sum
    df[i, paste(m, "height_mean", sep="_")] <- round(height_mean, 2)
    df[i, paste(m, "height_std", sep="_")] <- round(height_std, 2)
    df[i, paste(m, "height_max", sep="_")] <- height_max
    df[i, paste(m, "height_cv", sep="_")] <- round(height_cv, 4)
    df[i, paste(m, "phr", sep="_")] <- round(phr, 4)
    df[i, paste(m, "height_main_ratio", sep="_")] <- round(height_main_ratio, 4)
    df[i, paste(m, "peak_entropy", sep="_")] <- round(peak_entropy, 4)
  }

  # 计算全局特征
  df$avg_gene_sum[i] <- mean(gene_sums, na.rm=TRUE)
  df$sd_effective_peaks[i] <- sd(effective_peaks, na.rm=TRUE)
  df$avg_phr[i] <- mean(phr_values, na.rm=TRUE)
  df$var_phr[i] <- var(phr_values, na.rm=TRUE)
  df$avg_main_ratio[i] <- mean(main_ratios, na.rm=TRUE)
  entropy_cols <- grep("_approx_entropy$", names(df))
  df$global_entropy <- rowMeans(df[, entropy_cols], na.rm=TRUE)
  # 添加交互特征
  df$entropy_ratio <- df$global_entropy * df$prop_entropy
  df$height_balance <- df$height_total / (df$gene_total + 1e-6)
}

# 保存结果
write_xlsx(df, "data/enhanced_processed_Q2_data.xlsx")