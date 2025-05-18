import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
def plot_heatmap(matrix, len, wid, title="Heatmap", save_path="heatmap.png"):
    # 设置支持中文的字体
    # plt.rcParams['font.family'] = 'SimHei'

    # 解决负号显示问题
    plt.rcParams['axes.unicode_minus'] = False

    # 读取数据
    df = matrix

    # 绘制热力图
    plt.figure(figsize=(len, wid))
    sns.heatmap(
        df,
        annot=True,        # 显示数值
        fmt=".2f",         # 数值格式
        cmap="coolwarm",
        vmin=-1, vmax=1,   # 颜色范围固定为-1到1
        linewidths=0.5,
    )
    plt.title(title)
    plt.xticks(rotation=45, ha="right")  # 调整X轴标签角度
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()
    plt.close()

#读取数据
data=pd.read_excel("data/processed_Q4_data.xlsx")
data = data.drop(columns=['sample_file', 'proportion'])
data = data.iloc[:,0:15]

#绘制热图
correlation_matrix = data.corr()
plot_heatmap(correlation_matrix,10,8,title="Correlation Matrix Heatmap",save_path="./images/Q4/correlation_heatmap.png")
