library(readxl)
data1<-read_excel("D题附件：法医物证多人身份鉴定问题数据集/附件1：不同人数的STR图谱数据.xlsx")
data2<-read_excel("D题附件：法医物证多人身份鉴定问题数据集/附件4：去噪后的STR图谱数据.xlsx")
height_raw <- unlist(data1[, paste("Height", 1:100, sep = " ")])
height_raw<-na.omit(height_raw)

height<-unlist(data2[, paste("Height", 1:23, sep = " ")])
height<-na.omit(height)

# 创建boxplot并保存
if(!dir.exists("images")) dir.create("images")
png("images/boxplot.png", width = 8, height = 6, units = "in", res = 300)  # 设置300dpi分辨率
box_stats<-boxplot(list(height_raw, height),
                  names = c("raw", "denoise"),
                  col = c("lightblue", "pink"),  # 注意是col不是color
                  main = "Boxplot of Height")
# 添加分位数标注
text(x = rep(1.1:2.1, each=5),
     y = box_stats$stats,  # 垂直偏移
     labels = round(box_stats$stats, 2),
     pos = c(1,3,3,3,3),  # 上下交替位置
     col = rep(c("blue", "red"), each=5),
     cex = 0.6)

dev.off()  # 关闭图形设备并保存

hist(height_raw,breaks = 100,main = "Histogram of Height_Raw",xlab = "Height",ylab = "Frequency")
hist(height,breaks = 100,main = "Histogram of Height",xlab = "Height",ylab = "Frequency")