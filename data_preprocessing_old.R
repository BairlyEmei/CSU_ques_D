library(readxl)
library(dplyr)
library(writexl)

#读取数据
data<-read_excel("D题附件：法医物证多人身份鉴定问题数据集/附件2：不同混合比例的STR图谱数据_processed.xlsx")
# head(data)

#构建列名
marker<-c('D8S1179', 'D21S11', 'D7S820', 'CSF1PO', 'D3S1358', 'TH01', 'D13S317',
'D16S539', 'D2S1338', 'D19S433', 'vWA', 'TPOX', 'D18S51', 'AMEL', 'D5S818', 'FGA')
char<-c("gene_sum","height_sum","height_mean","height_std","height_max"
        ,"height_cv", "skewness", "kurtosis", "effective_peaks",
        "triplet_ratio", "phr", "low_phr", "high_freq_cluster",
        "peak_pair_ratio", "height_balance", "dominant_pair_ratio",
        "multi_peak_ratio", "height_dispersion", "peak_cluster_count"
)

marker_char<-paste(rep(marker,each=length(char)),char,sep="_")
colnames<-c("sample_file","gene_from","people","proportion", "gene_total","height_total",
            "cv_std","global_phr_mean","low_phr_count",
            marker_char)

#构建数据框
df<-data.frame(matrix(ncol=length(colnames),nrow=50))
colnames(df)<-colnames

#计算指标
df$sample_file<-unique(data[['Sample File']])
# df$people <- rep(c(2,3,4,5), times=c(15,16,11,9))
# df$proportion <- rep(c('1:1','1:1:1','1:1:1:1','1:1:1:1:1'), times=c(15,16,11,9))

# 按样本分组
unique_data <- data %>%
  group_by(`Sample File`) %>%
  summarise(
    gene_from = first(`基因来源（贡献者）`),
    people = first(`人数`),
    proportion = first(`比例`)
  )

# 赋值到df
df$sample_file <- unique_data$`Sample File`
df$gene_from <- unique_data$gene_from
df$people <- unique_data$people
df$proportion <- unique_data$proportion


# 获取所有Height列名
height_cols <- paste("Height", 1:100, sep = " ")

# 遍历每个样本
for (i in 1:nrow(df)) {
  current_sample <- df$sample_file[i]

  # 提取当前样本的所有数据
  sample_data <- data %>% filter(`Sample File` == current_sample)

  # 遍历每个marker
  for (m in marker) {
    # 获取当前marker的数据行
    marker_data <- sample_data %>% filter(Marker == m)

    # 检查是否存在该marker的数据
    if (nrow(marker_data) == 0) {
      warning(paste("Marker", m, "not found for sample", current_sample))
      next
    }

    # 提取所有Height值并转换为数值
    heights <- as.numeric(marker_data[1, height_cols])

    # 过滤有效Height值（>0且非NA）
    valid_heights <- heights[!is.na(heights) & heights > 0]

    #按阈值过滤噪声
    valid_heights <- valid_heights[valid_heights > 0.1*max(valid_heights)]


    # 计算统计量
    gene_sum <- length(valid_heights)
    height_sum <- sum(valid_heights)
    height_mean <- ifelse(gene_sum > 0, mean(valid_heights), NA)
    height_std <- ifelse(gene_sum > 1, sd(valid_heights), NA)
    height_max <- ifelse(gene_sum > 0, max(valid_heights), NA)
    height_cv <- ifelse(gene_sum > 0 & height_mean > 0, height_std/height_mean, NA)  # 新增变异系数
    # 计算峰高比(PHR) = 次高峰/最高峰
    sorted_heights <- sort(valid_heights, decreasing = TRUE)
    phr <- ifelse(length(sorted_heights) >= 2,
                 sorted_heights[2]/sorted_heights[1],
                 NA)
    triplet_ratio <- ifelse(length(sorted_heights) >= 3,
                           sorted_heights[3]/sorted_heights[1],
                           0)
    low_phr <- ifelse(!is.na(phr) & phr < 0.5, 1, 0)
    # 计算高频等位基因簇数量(峰高占比>15%)
    total_height <- sum(valid_heights)
    high_freq_cluster <- sum(valid_heights > 0.15 * total_height)  # 计算高频簇数量
    # 计算峰高分布偏度和峰度
    skewness <- ifelse(gene_sum > 2,
                      (sum((valid_heights - height_mean)^3)/gene_sum)/(height_std^3),
                      NA)
    kurtosis <- ifelse(gene_sum > 3,
                      (sum((valid_heights - height_mean)^4)/gene_sum)/(height_std^4) - 3,
                      NA)

    # 计算有效峰数(峰高超过总高度5%的峰)
    effective_peaks <- sum(valid_heights > 0.05 * total_height)

    # 在计算统计量的循环中添加新指标
    # 计算峰高分布均匀性 (1 - Gini系数)
    height_uniformity <- ifelse(gene_sum > 0,
                               1 - (sum(abs(valid_heights - mean(valid_heights)))/(2 * gene_sum * mean(valid_heights))),
                               NA)

    # 计算主峰占比
    dominant_peak_ratio <- ifelse(gene_sum > 0,
                                 max(valid_heights)/sum(valid_heights),
                                 NA)

    # 计算峰高分布熵值
    height_entropy <- ifelse(gene_sum > 0,
                            -sum((valid_heights/sum(valid_heights)) *
                                 log(valid_heights/sum(valid_heights))),
                            NA)
    # 在计算统计量的循环中添加新指标
    # 峰高对比例 (适用于1:1混合)
    peak_pair_ratio <- ifelse(length(sorted_heights) >= 2,
                             min(sorted_heights[1:2])/max(sorted_heights[1:2]),
                             NA)

    # 高度平衡指标 (1:1时接近1，非均衡时偏离)
    height_balance <- ifelse(length(sorted_heights) >= 2,
                           (sorted_heights[1] - sorted_heights[2])/sum(sorted_heights[1:2]),
                           NA)

    # 主峰对占比 (前两高峰占总和比例)
    dominant_pair_ratio <- ifelse(length(sorted_heights) >= 2,
                                 sum(sorted_heights[1:2])/sum(valid_heights),
                                 NA)
  # 多人峰比例 (前3高峰占总和比例)
    multi_peak_ratio <- ifelse(length(sorted_heights) >= 3,
                              sum(sorted_heights[1:3])/sum(valid_heights),
                              NA)

    # 高度离散度 (衡量峰高分布的离散程度)
    height_dispersion <- ifelse(gene_sum > 1,
                               sd(valid_heights)/mean(valid_heights),
                               NA)

    # 峰簇数量 (峰高超过平均高度一定比例的峰数)
    peak_cluster_count <- ifelse(gene_sum > 0,
                                sum(valid_heights > 0.3 * mean(valid_heights)),
                                NA)

    # 更新数据框
    df[i, paste(m, "gene_sum", sep = "_")] <- gene_sum
    df[i, paste(m, "height_sum", sep = "_")] <- height_sum
    df[i, paste(m, "height_mean", sep = "_")] <- round(height_mean, 2)
    df[i, paste(m, "height_std", sep = "_")] <- round(height_std, 2)
    df[i, paste(m, "height_max", sep = "_")] <- height_max
    df[i, paste(m, "height_cv", sep = "_")] <- round(height_cv, 4)
    df[i, paste(m, "phr", sep = "_")] <- round(phr, 4)
    df[i, paste(m, "low_phr", sep = "_")] <- low_phr
    df[i, paste(m, "triplet_ratio", sep = "_")] <- round(triplet_ratio, 4)
    df[i, paste(m, "high_freq_cluster", sep = "_")] <- high_freq_cluster
    df[i, paste(m, "skewness", sep = "_")] <- round(skewness, 4)
    df[i, paste(m, "kurtosis", sep = "_")] <- round(kurtosis, 4)
    df[i, paste(m, "effective_peaks", sep = "_")] <- effective_peaks
    df[i, paste(m, "height_uniformity", sep = "_")] <- round(height_uniformity, 4)
    df[i, paste(m, "dominant_peak_ratio", sep = "_")] <- round(dominant_peak_ratio, 4)
    df[i, paste(m, "height_entropy", sep = "_")] <- round(height_entropy, 4)
    df[i, paste(m, "peak_pair_ratio", sep = "_")] <- round(peak_pair_ratio, 4)
    df[i, paste(m, "height_balance", sep = "_")] <- round(height_balance, 4)
    df[i, paste(m, "dominant_pair_ratio", sep = "_")] <- round(dominant_pair_ratio, 4)
    df[i, paste(m, "multi_peak_ratio", sep = "_")] <- round(multi_peak_ratio, 4)
    df[i, paste(m, "height_dispersion", sep = "_")] <- round(height_dispersion, 4)
    df[i, paste(m, "peak_cluster_count", sep = "_")] <- peak_cluster_count

  }
}

# 计算总基因型和总高度
df$gene_total <- rowSums(df[, grep("_gene_sum$", names(df))], na.rm = TRUE)
df$height_total <- rowSums(df[, grep("_height_sum$", names(df))], na.rm = TRUE)

cv_cols <- grep("_height_cv$", names(df))
df$cv_std <- apply(df[, cv_cols], 1, sd, na.rm = TRUE)  # 新增变异系数标准差列

phr_cols <- grep("_phr$", names(df))
low_phr_cols <- grep("_low_phr$", names(df))

df$global_phr_mean <- rowMeans(df[, phr_cols], na.rm = TRUE)
df$low_phr_count <- rowSums(df[, low_phr_cols], na.rm = TRUE)

effective_peaks_cols <- grep("_effective_peaks$", names(df))
df$effective_peaks_mean <- rowMeans(df[, effective_peaks_cols], na.rm = TRUE)

skewness_cols <- grep("_skewness$", names(df))
df$skewness_mean <- rowMeans(df[, skewness_cols], na.rm = TRUE)

kurtosis_cols <- grep("_kurtosis$", names(df))
df$kurtosis_mean <- rowMeans(df[, kurtosis_cols], na.rm = TRUE)

# 计算峰数变异系数
gene_sum_cols <- grep("_gene_sum$", names(df))
df$gene_sum_cv <- apply(df[, gene_sum_cols], 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))

# 计算峰高均匀性均值
uniformity_cols <- grep("_height_uniformity$", names(df))
df$uniformity_mean <- rowMeans(df[, uniformity_cols], na.rm = TRUE)

# 计算主峰占比均值
dominant_ratio_cols <- grep("_dominant_peak_ratio$", names(df))
df$dominant_ratio_mean <- rowMeans(df[, dominant_ratio_cols], na.rm = TRUE)

# 计算熵值均值
entropy_cols <- grep("_height_entropy$", names(df))
df$entropy_mean <- rowMeans(df[, entropy_cols], na.rm = TRUE)

# 计算多人峰比例均值
multi_peak_cols <- grep("_multi_peak_ratio$", names(df))
df$multi_peak_mean <- rowMeans(df[, multi_peak_cols], na.rm = TRUE)

# 计算高度离散度均值
dispersion_cols <- grep("_height_dispersion$", names(df))
df$dispersion_mean <- rowMeans(df[, dispersion_cols], na.rm = TRUE)

# 计算峰簇数量均值
cluster_cols <- grep("_peak_cluster_count$", names(df))
df$cluster_mean <- rowMeans(df[, cluster_cols], na.rm = TRUE)

#保存结果
if(!dir.exists("data")) dir.create("data")
write_xlsx(df, path = "data/processed_Q2_data.xlsx")
