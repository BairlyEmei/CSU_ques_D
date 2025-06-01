from random_forest_Q2 import *
from scipy.optimize import minimize
from itertools import product, combinations
from tqdm import tqdm
import numpy as np
from collections import defaultdict
import pandas as pd

# 生成样本矩阵C
def generate_C_matrix(sample_name, data_test, GENE_LOCI):
    """
    返回NumPy矩阵的最终版本

    参数:
    sample_name (str): 样本名称（需完全匹配）
    data_test (pd.DataFrame): 包含所有样本数据的DataFrame
    GENE_LOCI (dict): 各基因座的标准等位基因定义

    返回:
    np.ndarray: 16xN的归一化峰高比例矩阵，N为各基因座最大等位基因数
    """
    # 输入验证
    if not isinstance(data_test, pd.DataFrame):
        raise ValueError("输入数据必须为Pandas DataFrame")
    if sample_name not in data_test['Sample File'].values:
        raise ValueError(f"样本 {sample_name} 不存在")

    # 获取基因座顺序列表
    loci_order = list(GENE_LOCI.keys())

    # 预计算最大等位基因数量
    max_alleles = max(len(v) for v in GENE_LOCI.values())

    # 初始化矩阵容器
    matrix_data = np.zeros((len(GENE_LOCI), max_alleles), dtype=np.float32)

    # 处理每个基因座
    for locus_idx, (locus_name, alleles) in enumerate(GENE_LOCI.items()):
        # 转换标准等位基因格式
        if locus_name == 'AMEL':
            std_alleles = [str(a).upper() for a in alleles]  # 性别基因大写处理
        else:
            std_alleles = [f"{float(a):g}" for a in alleles]  # 数值型统一格式

        # 获取当前基因座数据
        sample_mask = (
            (data_test['Sample File'] == sample_name) &
            (data_test['Marker'].str.upper() == locus_name.upper())
        )
        locus_data = data_test[sample_mask].copy()

        if not locus_data.empty:
            allele_counts = defaultdict(float)

            # 遍历所有Allele列
            for col in [c for c in locus_data.columns if c.startswith('Allele')]:
                col_num = col.split()[-1]
                height_col = f'Height {col_num}'

                # 提取等位基因和峰高
                allele = str(locus_data[col].iloc[0]).strip()
                height = locus_data[height_col].iloc[0]

                # 过滤无效数据
                if pd.isna(allele) or pd.isna(height) or 'OL' in allele:
                    continue

                # 格式处理
                try:
                    processed_allele = (
                        allele.upper() if locus_name == 'AMEL'
                        else f"{float(allele):g}"
                    )
                except ValueError:
                    continue

                # 累加峰高值
                allele_counts[processed_allele] += float(height)

            # 应用动态阈值筛选
            valid_heights = [h for h in allele_counts.values() if h > 0]
            if valid_heights:
                max_height = max(valid_heights)
                total_height = sum(valid_heights)
                threshold = max(0.1 * max_height, 0.01 * total_height)
                filtered = {k:v for k,v in allele_counts.items() if v >= threshold}

                # 归一化处理
                valid_total = sum(filtered.get(a, 0) for a in std_alleles)
                if valid_total > 0:
                    for idx, allele in enumerate(std_alleles):
                        if idx >= max_alleles:  # 防止索引越界
                            break
                        matrix_data[locus_idx, idx] = filtered.get(allele, 0.0) / valid_total

    # 验证矩阵维度
    assert matrix_data.shape == (16, max_alleles), f"矩阵维度异常，应为(16, {max_alleles})"

    return matrix_data


def estimate_contributors_by_locus_simple_avg(k, C_matrix, valid_comb_dict, GENE_LOCI):
    """
    按基因座优化并使用简单平均的贡献者比例估计

    参数:
    k (int): 贡献者人数
    C_matrix (list): C矩阵，16个基因座向量
    valid_comb_dict (dict): 每个基因座所有合法组合的字典
    GENE_LOCI (dict): 基因座定义字典

    返回:
    final_p: 最终估计的贡献者比例
    locus_results: 每个基因座的优化结果
    """
    loci_order = list(GENE_LOCI.keys())
    n_loci = len(loci_order)

    # 存储每个基因座的结果
    locus_results = []

    print(f"开始按基因座优化，贡献者人数k={k}\n")
    for i, locus_name in enumerate(loci_order):
        C_vector = np.array(C_matrix[i])
        valid_combs = valid_comb_dict[locus_name]

        # 输出当前基因座信息
        print(f"处理基因座 '{locus_name}'")

        # 初始化最佳结果
        min_residual = float('inf')
        best_p = None

        # 如果该基因座没有数据或无效，使用均匀分布
        if not valid_combs:
            uniform_p = np.ones(k) / k
            locus_results.append(uniform_p)
            print(f"  无有效组合，使用均匀分布: {uniform_p.round(3)}\n")
            continue

        # 穷举所有可能的组合
        comb_product = list(product(valid_combs, repeat=k))

        print(f"  考虑{len(comb_product)}种组合...", end='')

        for comb in comb_product:
            # 构建设计矩阵 (n_alleles x k)
            A_matrix = np.array(comb).T

            # 目标函数：最小化残差平方和
            def objective(p):
                return np.linalg.norm(A_matrix @ p - C_vector)

            # 约束：比例和为1且非负
            constraints = (
                {'type': 'eq', 'fun': lambda p: np.sum(p) - 1},
                {'type': 'ineq', 'fun': lambda p: p}
            )

            # 初始猜测值（均匀分布）
            p0 = np.ones(k) / k

            # 求解优化问题
            res = minimize(objective, p0, method='SLSQP', constraints=constraints)

            if res.success and res.fun < min_residual:
                min_residual = res.fun
                best_p = res.x

        # 存储并输出该基因座的最佳结果
        if best_p is None:  # 如果没有找到有效解
            best_p = np.ones(k) / k
            min_residual = 0

        best_p = np.clip(best_p, 0, 1)  # 确保比例在[0,1]范围内
        best_p /= best_p.sum()  # 重新归一化

        locus_results.append(best_p)

        print(f"完成! 最小残差={min_residual:.4f}")
        print(f"  最优比例向量: {best_p.round(4)}")
        print()

    # 3. 使用简单平均替代加权平均
    print("\n所有基因座处理完成，开始简单平均")
    final_p = np.zeros(k)

    # 累加所有基因座的比例
    for i, locus_p in enumerate(locus_results):
        final_p += locus_p
        print(f"  基因座 '{loci_order[i]}' 贡献: {locus_p.round(4)}")

    # 计算简单平均
    final_p /= n_loci
    final_p /= final_p.sum()  # 最终归一化

    print(f"\n最终混合比例(简单平均): {final_p.round(4)}")
    print(f"比例总和: {final_p.sum():.4f}")

    return final_p, locus_results

def estimate_contributors_by_locus(k, C_matrix, locus_total_heights, valid_comb_dict, GENE_LOCI):
    """
    按基因座优化并加权平均的贡献者比例估计（带每个基因座结果输出）

    参数:
    k (int): 贡献者人数
    C_matrix (list): C矩阵，16个基因座向量
    locus_total_heights (list): 每个基因座的总峰高
    valid_comb_dict (dict): 每个基因座所有合法组合的字典
    GENE_LOCI (dict): 基因座定义字典

    返回:
    final_p: 最终估计的贡献者比例
    locus_results: 每个基因座的优化结果
    weights: 每个基因座的权重
    """
    from itertools import product
    from scipy.optimize import minimize
    import numpy as np

    loci_order = list(GENE_LOCI.keys())
    n_loci = len(loci_order)

    # 存储每个基因座的结果
    locus_results = []
    locus_weights = []

    # 1. 计算权重 (ω_i = √(Σc_ij))
    weights = [np.sqrt(th) if th > 0 else 0 for th in locus_total_heights]
    total_weight = sum(weights)
    normalized_weights = [w / total_weight if total_weight > 0 else 1/n_loci for w in weights]

    # 2. 处理每个基因座
    print(f"开始按基因座优化，贡献者人数k={k}\n")
    for i, locus_name in enumerate(loci_order):
        C_vector = np.array(C_matrix[i])
        valid_combs = valid_comb_dict[locus_name]

        # 输出当前基因座信息
        print(f"处理基因座 '{locus_name}'，总峰高={locus_total_heights[i]:.0f}")

        # 初始化最佳结果
        min_residual = float('inf')
        best_p = None
        best_combo = None

        # 如果该基因座没有数据，使用均匀分布
        if locus_total_heights[i] <= 0 or not valid_combs:
            uniform_p = np.ones(k) / k
            locus_results.append(uniform_p)
            locus_weights.append(normalized_weights[i])
            print(f"  无有效数据，使用均匀分布: {uniform_p.round(3)}\n")
            continue

        # 穷举所有可能的组合（k个贡献者的所有基因型组合）
        comb_product = list(product(valid_combs, repeat=k))
        total_combs = len(comb_product)

        print(f"  考虑{total_combs}种组合...", end='')

        for j, comb in enumerate(comb_product):
            # 构建设计矩阵 (n_alleles x k)
            A_matrix = np.array(comb).T

            # 目标函数：最小化残差平方和
            def objective(p):
                return np.linalg.norm(A_matrix @ p - C_vector)

            # 约束：比例和为1且非负
            constraints = (
                {'type': 'eq', 'fun': lambda p: np.sum(p) - 1},
                {'type': 'ineq', 'fun': lambda p: p}
            )

            # 初始猜测值（均匀分布）
            p0 = np.ones(k) / k

            # 求解优化问题
            res = minimize(objective, p0, method='SLSQP', constraints=constraints)

            if res.success and res.fun < min_residual:
                min_residual = res.fun
                best_p = res.x
                best_combo = comb

        # 存储并输出该基因座的最佳结果
        best_p = np.clip(best_p, 0, 1)  # 确保比例在[0,1]范围内
        best_p /= best_p.sum()  # 重新归一化

        locus_results.append(best_p)
        locus_weights.append(normalized_weights[i])

        print(f"完成! 最小残差={min_residual:.4f}")
        print(f"  最优比例向量: {best_p.round(4)}")
        print(f"  贡献者基因型组合:")

        # 解码最优组合
        alleles = GENE_LOCI[locus_name]
        for p_idx in range(k):
            combo_vec = best_combo[p_idx]
            allele_names = []
            for idx, count in enumerate(combo_vec):
                if count > 0:
                    allele_names.append(f"{alleles[idx]}×{int(count)}")
            print(f"    贡献者{p_idx+1}: {' + '.join(allele_names)}")
        print()

    # 3. 按权重进行加权平均 (p_j = Σ(ω_i′ * p_ji))
    print("\n所有基因座处理完成，开始加权平均")
    final_p = np.zeros(k)
    for i, locus_name in enumerate(loci_order):
        # 输出每个基因座的权重和贡献
        weight_contrib = normalized_weights[i] * locus_results[i]
        print(f"  基因座 '{locus_name}' (权重={normalized_weights[i]:.4f}) 贡献:")
        print(f"      比例向量: {locus_results[i].round(4)}")
        print(f"      加权贡献: {weight_contrib.round(4)}")

        # p_i 是该基因座下各贡献者的比例向量
        final_p += weight_contrib

    # 最终归一化处理
    final_p = np.clip(final_p, 0, 1)
    final_p /= final_p.sum()

    print(f"\n最终混合比例: {final_p.round(4)}")
    print(f"比例总和: {final_p.sum():.4f}")

    return final_p, locus_results, normalized_weights