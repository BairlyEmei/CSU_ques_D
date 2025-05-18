library(readxl)
library(dplyr)
library(writexl)

#读取数据
data<-read_excel("D题附件：法医物证多人身份鉴定问题数据集/附件4：去噪后的STR图谱数据.xlsx")
# head(data)

#构建列名
marker<-c('D8S1179', 'D21S11', 'D7S820', 'CSF1PO', 'D3S1358', 'TH01', 'D13S317',
'D16S539', 'D2S1338', 'D19S433', 'vWA', 'TPOX', 'D18S51', 'AMEL', 'D5S818', 'FGA')
char<-c("gene_sum","height_sum","height_mean","height_std","height_max"
                ,"height_cv", "phr", "low_phr", "high_freq_cluster"
        )
marker_char<-paste(rep(marker,each=length(char)),char,sep="_")
colnames<-c("sample_file","people","proportion", "gene_total","height_total",
            "cv_std","global_phr_mean","low_phr_count",
            marker_char)

#构建数据框
df<-data.frame(matrix(ncol=length(colnames),nrow=51))
colnames(df)<-colnames

#计算指标
df$sample_file<-unique(data[['Sample File']])
df$people <- rep(c(2,3,4,5), times=c(15,16,11,9))
df$proportion <- rep(c('1:1','1:1:1','1:1:1:1','1:1:1:1:1'), times=c(15,16,11,9))

# 获取所有Height列名
height_cols <- paste("Height", 1:23, sep = " ")

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
    low_phr <- ifelse(!is.na(phr) & phr < 0.5, 1, 0)
    # 计算高频等位基因簇数量(峰高占比>15%)
    total_height <- sum(valid_heights)
    high_freq_cluster <- sum(valid_heights > 0.15 * total_height)  # 计算高频簇数量

    # 更新数据框
    df[i, paste(m, "gene_sum", sep = "_")] <- gene_sum
    df[i, paste(m, "height_sum", sep = "_")] <- height_sum
    df[i, paste(m, "height_mean", sep = "_")] <- round(height_mean, 2)
    df[i, paste(m, "height_std", sep = "_")] <- round(height_std, 2)
    df[i, paste(m, "height_max", sep = "_")] <- height_max
    df[i, paste(m, "height_cv", sep = "_")] <- round(height_cv, 4)
    df[i, paste(m, "phr", sep = "_")] <- round(phr, 4)
    df[i, paste(m, "low_phr", sep = "_")] <- low_phr
    df[i, paste(m, "high_freq_cluster", sep = "_")] <- high_freq_cluster

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

#保存结果
if(!dir.exists("data")) dir.create("data")
write_xlsx(df, path = "data/processed_Q4_data.xlsx")
