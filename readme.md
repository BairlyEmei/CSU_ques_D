# 介绍
这是一个用于csu2025数模校赛的项目~~

队伍：渐变的金黄色宇宙队  
语言: python, R  
题目：[D题：法医物证多人身份鉴定问题.docx](论文材料/D题：法医物证多人身份鉴定问题.docx)
## 问题1
具体思路见[quesD_try.pdf](论文材料/quesD_try.pdf)
### 特征工程
1. 目标特征：
   - 混合人数
   - 混合比例
2. 每个基因座特征（*16）：
   ```R
   enhanced_char <- c("gene_sum","height_sum","height_mean","height_std","height_max",
                      "height_cv", "skewness", "kurtosis", "effective_peaks",
                      "triplet_ratio", "phr", "low_phr", "high_freq_cluster",
                      "height_main_ratio", "peak_entropy")
   ```
3. 全局特征
   ```R
   global_features <- c("sample_file","gene_from","people","proportion", "gene_total","height_total",
                        "cv_std","prop_max", "prop_min", "prop_sd", "prop_entropy",
                        "avg_gene_sum", "sd_effective_peaks", "avg_phr", "var_phr",
                        "total_low_phr", "avg_main_ratio")
   ```
4. 动态噪声过滤：结合最大峰高和总峰高
   ```R
    valid_heights <- heights[!is.na(heights)]
    if (length(valid_heights) > 0) {
      total_height <- sum(valid_heights)
      threshold <- max(0.1*max(valid_heights), 0.01*total_height)
      valid_heights <- valid_heights[valid_heights > threshold]
    }
    ```
5. 数据预处理代码:[data_preprocessing_new.R](data_preprocessing_new.R)
### 运行结果
1. 数据集：[enhanced_processed_Q1_data.xlsx](data/enhanced_processed_Q1_data.xlsx)，八二分
2. 模型：随机森林[random_forest.py](random_forest.py)
3. 交叉验证最佳分数: 0.9666666666666666
4. 最佳参数: {'max_depth': 5, 'max_features': 'sqrt', 'min_samples_leaf': 1, 'min_samples_split': 2, 'n_estimators': 200}
5. 模型评估指标：
   - 准确率: 1.0000
   - 精确率: 1.0000
   - 召回率: 1.0000
   - F1分数: 1.0000
   - AUC值: 1.0000
6. 模型性能可视化：
   ![特征重要性图](images/Q1/feature_importance.png)
   ![ROC曲线](images/Q1/roc_curve.png)
   ![混淆矩阵](images/Q1/confusion_matrix.png)

## 问题2
具体思路见[第二问建模2.pdf](论文材料/第二问建模2.pdf)
### 判断人数
数据预处理同上
1. 数据集：[enhanced_processed_Q2_data.xlsx](data/enhanced_processed_Q2_data.xlsx)
2. 模型：随机森林[random_forest_Q2.py](random_forest_Q2.py)
3. 交叉验证最佳分数: 0.5369047619047619
4. 最佳参数: {'max_depth': 5, 'max_features': 'sqrt', 'min_samples_leaf': 2, 'min_samples_split': 5, 'n_estimators': 300}
5. 模型评估指标:
   - 准确率: 0.8000
   - 精确率: 0.8542
   - 召回率: 0.8333
   - F1分数: 0.8310
   - AUC值: 0.9259
6. 模型性能可视化：
   ![特征重要性图](images/Q2/feature_importance.png)
   ![ROC曲线](images/Q2/roc_curve.png)
   ![混淆矩阵](images/Q2/confusion_matrix.png)
### 判断比例
（以2人的样本为例）
1. 设出两人的基因型为矩阵A1,A2=（a_ij)16*j,行为16个基因座，列为等位基因，a_ij={0,1,2}表示有该基因的数量。a_ij对应的具体基因已经约定。
2. 穷举a_ij,Ai=p1*A1+p2*A2，为理论混合数据（单位：基因个数）
3. 构建样本矩阵Bi=（b_ij)16*j，查找对应位置的基因对应的峰高，未查到的记为NA。为实际混合数据（单位：基因峰高）。并且按第一问的方法过滤峰高。
4. 输入矩阵Ci：对Bi按行做归一化。
5. 比较Ai与Ci（后续改为按行（基因座）进行比较求出比例，再进行加权平均）
6. 减小穷举量，让A根据过滤后的C来穷举，有值的位置才取0，1，2，没值的位置取0。
7. 并对A进行去重，防止相同的情况重复计算。
8. 具体实现步骤看我Q2代码，反正jupyter标注好了，懒得誊到这，一切以代码为准。
### 编程手表示想吐槽：
1. 这辈子不想碰深圳杯第二次
2. 经过我对数据的进一步研究和理解，对同一个基因，峰高可以反映其量的多少。但对于不同基因，峰高还受到基因长度的影响，不能完美反映比例。
尤其是AMEL基因座的性染色体基因，TMD你告诉我X比Y还要少，这对劲吗？？我本来以为是因为Y的长度显著小于X，所以Y的峰高高，结果一看它俩的size，MDY的size还比X大一点点，
而且差距也不大。我本来都假设了size差距不大的情况下峰高比例是有效的了，不是那你怎么给我保证其他常染色体的基因座上的峰高比例有效？
3. 难道我还需要再去研究研究不同基因的size对峰高的影响，每个拟合个函数出来，再去求解？
4. woc深圳杯怎么这么难啊
5. 不管了将就做吧跑个值出来就行哩管他对不对呢
6. 但是AMEL基因座你肯定是别想进入我的模型了👊
7. 算了你还是进来一起算吧因为第三问要算基因型
8. woc穷举法，给我的电脑感动坏了
### 运行结果
见result文件夹，详细基因型预测请翻阅每个样本的log文件  
请结合log文件和[问题2测试集](result/Q2_data_test.xlsx)进行对比，第三问请手动对照[附件3：各个贡献者对应的基因型数据.xlsx](D题附件：法医物证多人身份鉴定问题数据集/附件3：各个贡献者对应的基因型数据.xlsx)
1. [A05_RD14-0003-46_47_48-1;1;1-M3I22-0.189IP-Q0.9_001.5sec.fsa](result/A05_RD14-0003-46_47_48-1;1;1-M3I22-0.189IP-Q0.9_001.5sec.fsa.log)
2. [A06_RD14-0003-49_50_29-1;4;1-M3I15-0.186IP-Q0.9_001.5sec.fsa](result/A06_RD14-0003-49_50_29-1;4;1-M3I15-0.186IP-Q0.9_001.5sec.fsa.log)
3. [A05_RD14-0003-30_31_32_33_34-1;1;1;1;1-M3I22-0.315IP-Q1.3_001.5sec.fsa](result/A05_RD14-0003-30_31_32_33_34-1;1;1;1;1-M3I22-0.315IP-Q1.3_001.5sec.fsa.log)
4. [A07_RD14-0003-31_32_33_34_35-1;1;2;9;1-M3d-0.21IP-Q1.8_001.5sec.fsa](result/A07_RD14-0003-31_32_33_34_35-1;1;2;9;1-M3d-0.21IP-Q1.8_001.5sec.fsa.log)
5. [A09_RD14-0003-36_37_38_39_40-1;9;9;9;1-M2d-0.5IP-Q1.6_001.5sec.fsa](result/A09_RD14-0003-36_37_38_39_40-1;9;9;9;1-M2d-0.5IP-Q1.6_001.5sec.fsa.log)
6. [A02_RD14-0003-50_29_30_31-1;1;2;1-M3c-0.075IP-Q0.5_001.5sec.fsa](result/A02_RD14-0003-50_29_30_31-1;1;2;1-M3c-0.075IP-Q0.5_001.5sec.fsa.log)
7. [A04_RD14-0003-32_33_34_35-1;1;9;1-M2c-0.372IP-Q0.4_001.5sec.fsa](result/A04_RD14-0003-32_33_34_35-1;1;9;1-M2c-0.372IP-Q0.4_001.5sec.fsa.log)
8. [A10_RD14-0003-48_49_50_29-1;4;4;4-M2a-0.5IP-Q0.4_001.5sec.fsa](result/A10_RD14-0003-48_49_50_29-1;4;4;4-M2a-0.5IP-Q0.4_001.5sec.fsa.log)
9. [A09_RD14-0003-49_50_29-1;4;1-M3e-0.186IP-Q1.6_001.5sec.fsa](result/A09_RD14-0003-49_50_29-1;4;1-M3e-0.186IP-Q1.6_001.5sec.fsa.log)
10. [A04_RD14-0003-40_41-1;4-M3e-0.155IP-Q0.8_001.5sec.fsa](result/A04_RD14-0003-40_41-1;4-M3e-0.155IP-Q0.8_001.5sec.fsa.log)