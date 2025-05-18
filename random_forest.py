import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, roc_auc_score, roc_curve, confusion_matrix, ConfusionMatrixDisplay
)
import os
import sys

# 设置字体
if sys.platform == 'darwin':  # macOS
    font_path = '/System/Library/Fonts/STHeiti Light.ttc'
elif sys.platform == 'win32': # Windows
    plt.rcParams['font.sans-serif'] = ['SimHei']  # Windows系统自带黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
else: # Linux/其他系统
    font_path = '/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc'  # 文泉驿字体

# ======================
# 1. 数据读取与预处理
# ======================
def load_and_preprocess(data_path):
    """读取数据并进行预处理"""
    # 读取Excel文件
    df = pd.read_excel(data_path)

    #选取列
    df = df.drop(columns=['sample_file', 'proportion'])

    # 处理缺失值
    # 方法1：用列均值填充
    # df = df.fillna(df.mean())
    # 方法2：用0填充（适合你的高度数据）
    df = df.fillna(0)
    return df


# ======================
# 2. 随机森林模型构建
# ======================
from sklearn.model_selection import GridSearchCV

def train_random_forest(X_train, y_train):
    """训练随机森林分类器(带参数调优)"""
    # 定义基础模型
    rf = RandomForestClassifier(class_weight='balanced', random_state=42)

    # 定义参数网格
    param_grid = {
        'n_estimators': [100, 200, 300],
        'max_depth': [5, 8, 10, None],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4]
    }

    # 创建GridSearchCV对象
    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        cv=5,
        scoring='f1_macro',
        n_jobs=-1,
        verbose=1
    )

    # 执行网格搜索
    grid_search.fit(X_train, y_train)

    # 输出最佳参数
    print("最佳参数:", grid_search.best_params_)

    return grid_search.best_estimator_


# ======================
# 3. 可视化与评估
# ======================
def plot_feature_importance(model, feature_names, save_path):
    """绘制特征重要性图"""
    importance = model.feature_importances_
    indices = np.argsort(importance)[::-1]

    indices = indices[:20]
    importance = importance[indices]
    feature_names = [feature_names[i] for i in indices]

    plt.figure(figsize=(12, 8))
    plt.title("随机森林前20特征重要性排序", fontsize=14)
    plt.barh(range(len(indices)), importance, align='center')  # 修改这里
    plt.yticks(range(len(indices)), feature_names)  # 修改这里
    plt.gca().invert_yaxis()
    plt.xlabel("相对重要性", fontsize=12)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_roc_curve(y_true, y_proba, save_path):
    """绘制多类ROC曲线"""
    from sklearn.preprocessing import label_binarize
    from sklearn.metrics import roc_curve, auc

    # 二值化标签
    y_true_bin = label_binarize(y_true, classes=[2, 3, 4, 5])
    n_classes = y_true_bin.shape[1]

    # 计算每个类的ROC曲线和AUC
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_true_bin[:, i], y_proba[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # 绘制所有类的ROC曲线
    plt.figure()
    colors = ['darkorange', 'navy', 'green', 'red']
    for i, color in zip(range(n_classes), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=2,
                 label=f'类别 {i + 2} (AUC = {roc_auc[i]:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('假正率', fontsize=12)
    plt.ylabel('真正率', fontsize=12)
    plt.title('多类接收者操作特征曲线', fontsize=14)
    plt.legend(loc="lower right")
    plt.savefig(save_path, dpi=300)
    plt.show()
    plt.close()

# ======================
# 主流程
# ======================
if __name__ == "__main__":
    # 创建图片保存目录
    os.makedirs("images/Q1", exist_ok=True)
    os.makedirs("images/Q2", exist_ok=True)
    os.makedirs("images/Q3", exist_ok=True)
    os.makedirs("images/Q4", exist_ok=True)

    # 1. 数据加载与预处理
    df = load_and_preprocess(data_path = "data/processed_Q4_data.xlsx")

    X = df.drop('people', axis=1)
    y = df['people']

    # 2. 划分训练集/测试集
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.3,              # 测试集比例
        stratify=y,                 # 按类别分层抽样
        random_state=42)

    # 3. 模型训练
    rf_model = train_random_forest(X_train, y_train)

    # 4. 模型预测
    y_pred = rf_model.predict(X_test)
    y_proba = rf_model.predict_proba(X_test)  # 获取所有类别的概率

    # 5. 模型评估
    print("=" * 40)
    print("模型评估指标：")
    print(f"准确率: {accuracy_score(y_test, y_pred):.4f}")
    print(f"精确率: {precision_score(y_test, y_pred, average='macro'):.4f}")
    print(f"召回率: {recall_score(y_test, y_pred, average='macro'):.4f}")
    print(f"F1分数: {f1_score(y_test, y_pred, average='macro'):.4f}")
    print(f"AUC值: {roc_auc_score(y_test, y_proba, multi_class='ovo'):.4f}")

    # 6. 可视化输出
    # 特征重要性图
    plot_feature_importance(rf_model, X.columns.tolist(), "images/Q4/feature_importance.png")

    # ROC曲线
    plot_roc_curve(y_test, y_proba, "images/Q4/roc_curve.png")

    # 混淆矩阵
    cm = confusion_matrix(y_test, y_pred)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=rf_model.classes_)
    disp.plot(cmap=plt.cm.Blues)
    plt.title("混淆矩阵")
    plt.savefig("images/Q4/confusion_matrix.png", dpi=300)
    plt.show()
    plt.close()