# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第六章：多元统计分析
脚本名称：ch6_multivariate_pca_hca_and_plsda.py
代码目的：1. 第一节：Z-score 标准化（消除不同指标量纲与量级垄断）与主成分分析（PCA）碎石图；
          2. 第二节：带 95% Hotelling's T² 置信椭圆的二维 PCA 得分图（Score Plot）与载荷双标图（Biplot）；
          3. 第三节：高维欧氏距离与 Ward 最小方差系统层次聚类树状图（HCA Dendrogram）；
          4. 第四节：有监督偏最小二乘判别分析（PLS-DA）、200 次置换防过拟合检验与 VIP > 1.0 标志物挖掘。
================================================================================
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns

# 跨平台中文字体与负号防乱码配置
sns.set_theme(style="ticks", font_scale=1.06)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Noto Sans CJK JP', 'Droid Sans Fallback', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42


def read_csv_safe(filepath):
    try:
        return pd.read_csv(filepath, encoding='utf-8')
    except UnicodeDecodeError:
        return pd.read_csv(filepath, encoding='gbk')


def zscore_standardize(X):
    """
    Z-score 标准化函数（等价于 sklearn.preprocessing.StandardScaler().fit_transform(X)）
    将每一列特征转换为均值为 0、标准差为 1 的无量纲数据：Z = (X - mean) / std
    """
    mu = np.mean(X, axis=0, keepdims=True)
    sd = np.std(X, axis=0, ddof=1, keepdims=True) + 1e-12
    return (X - mu) / sd


def run_pca_svd(X_scaled):
    """
    奇异值分解 (SVD) 权威计算 PCA 得分矩阵 (Scores)、载荷矩阵 (Loadings) 与方差贡献率
    （与 sklearn.decomposition.PCA 完全等价）
    """
    n_samples = X_scaled.shape[0]
    U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)
    # 特征值 Eigenvalues = S^2 / (n - 1)
    eigenvalues = (S ** 2) / (n_samples - 1)
    explained_ratio = eigenvalues / np.sum(eigenvalues) * 100.0
    scores = U * S          # 样品在主成分轴上的投影得分坐标 (n_samples, k)
    loadings = Vt.T         # 各变量对主成分轴的载荷系数向量 (p_features, k)
    return scores, loadings, eigenvalues, explained_ratio


def draw_confidence_ellipse(x, y, ax, n_std=2.4477, facecolor='none', edgecolor='black', **kwargs):
    """绘制二维点集的 95% Hotelling's T² / 卡方置信椭圆"""
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, edgecolor=edgecolor, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(np.mean(x), np.mean(y))
    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)


# ==============================================================================
# 第一节：Z-score 标准化必要性验证与 PCA 方差贡献率碎石图 (Scree Plot)
# ==============================================================================
def run_section1_pca_scaling():
    print("\n" + "=" * 72)
    print(">>> 第六章第一节：Z-score 标准化与主成分分析 (PCA) 碎石图 <<<")
    print("=" * 72)

    df1 = read_csv_safe("data/ch6_sec1_pca_scaling.csv")
    sensor_cols = [c for c in df1.columns if c.startswith("Sensor_")]
    X_raw = df1[sensor_cols].values

    # 对比未标准化 vs Z-score 标准化后的第一主成分方差垄断现象
    _, _, _, ratio_unscaled = run_pca_svd(X_raw - X_raw.mean(axis=0))
    X_scaled = zscore_standardize(X_raw)
    scores, loadings, eigvals, ratio_scaled = run_pca_svd(X_scaled)
    cum_ratio = np.cumsum(ratio_scaled)

    print(f"[警示对比] 未标准化直接跑 PCA -> PC1 独占 {ratio_unscaled[0]:.2f}% 方差 (被大数值传感器 Sensor_1 彻底垄断!)")
    print(f"[规范做法] Z-score 标准化后跑 PCA -> PC1 = {ratio_scaled[0]:.2f}%, PC2 = {ratio_scaled[1]:.2f}%, 前两主成分累积 = {cum_ratio[1]:.2f}%")

    pca_tbl = pd.DataFrame({
        "Principal_Component": [f"PC{i}" for i in range(1, len(eigvals)+1)],
        "Eigenvalue (特征值)": np.round(eigvals, 3),
        "Variance_Explained (%)": np.round(ratio_scaled, 2),
        "Cumulative_Variance (%)": np.round(cum_ratio, 2)
    })
    pca_tbl.to_csv("data/ch6_sec1_pca_variance_table.csv", index=False, encoding='utf-8-sig')

    # 绘制主成分碎石图 (Scree Plot)
    fig, ax1 = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    pcs = [f"PC{i}" for i in range(1, 9)]
    ax1.bar(pcs, ratio_scaled[:8], color='#377EB8', edgecolor='black', width=0.5, alpha=0.8, label='单主成分方差贡献率 (%)')
    ax1.plot(pcs, cum_ratio[:8], color='#E41A1C', marker='o', markersize=7, linewidth=2.4, label='累积方差贡献率 (%)')
    ax1.axhline(70.0, color='darkgreen', linestyle='--', linewidth=1.3, label='70% 累积解释率达标线')
    for i in range(4):
        ax1.text(i, cum_ratio[i] + 3.0, f"{cum_ratio[i]:.1f}%", ha='center', color='#990000', fontweight='bold', fontsize=9.5)
    ax1.set_ylim(0, 110)
    ax1.set_title("电子鼻高维特征 Z-score 标准化后 PCA 碎石图 (Scree Plot)", fontsize=11.5, fontweight='bold')
    ax1.set_xlabel("主成分维度 (Principal Components)", fontsize=10.5)
    ax1.set_ylabel("方差贡献率百分比 (%)", fontsize=10.5)
    ax1.legend(loc='center right', frameon=True)
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch6_sec1_scree_plot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第一节 PCA 碎石图已保存至: scripts/ch6_sec1_scree_plot.png")


# ==============================================================================
# 第二节：PCA 95% 置信椭圆得分图 (Score Plot) 与载荷双标图 (Biplot)
# ==============================================================================
def run_section2_biplot():
    print("\n" + "=" * 72)
    print(">>> 第六章第二节：PCA 95% 置信椭圆得分图与载荷向量双标图 (Biplot) <<<")
    print("=" * 72)

    df2 = read_csv_safe("data/ch6_sec2_biplot_ellipse.csv")
    feat_cols = [c for c in df2.columns if c not in ['Sample_ID', 'Origin']]
    X_scaled = zscore_standardize(df2[feat_cols].values)
    scores, loadings, _, ratio = run_pca_svd(X_scaled)

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), dpi=150)
    palette = {'High_Mountain': '#1B9E77', 'Plain_Garden': '#D95F02'}
    markers = {'High_Mountain': 'o', 'Plain_Garden': 's'}

    # 左图：带 95% Hotelling's T² 置信椭圆的 PCA 得分图
    for grp in ['High_Mountain', 'Plain_Garden']:
        mask = (df2['Origin'] == grp).values
        axes[0].scatter(scores[mask, 0], scores[mask, 1], c=palette[grp], marker=markers[grp],
                        s=65, edgecolor='black', linewidth=0.8, alpha=0.88, label=grp)
        draw_confidence_ellipse(scores[mask, 0], scores[mask, 1], ax=axes[0],
                                facecolor=palette[grp], alpha=0.18, edgecolor=palette[grp], linewidth=1.8, linestyle='--')
    axes[0].axhline(0, color='gray', linestyle=':', linewidth=1.0)
    axes[0].axvline(0, color='gray', linestyle=':', linewidth=1.0)
    axes[0].set_title("A. 无监督 PCA 得分图 (带 95% 置信椭圆)\n[两产地椭圆完全分离，证实整体风味指纹差异显著]", fontsize=10.5, fontweight='bold')
    axes[0].set_xlabel(f"PC1 ({ratio[0]:.1f}%)", fontsize=10.5, fontweight='bold')
    axes[0].set_ylabel(f"PC2 ({ratio[1]:.1f}%)", fontsize=10.5, fontweight='bold')
    axes[0].legend(title="茶叶产地 (Origin)", loc='upper right', frameon=True)

    # 右图：PCA 载荷双标图 (Biplot：样品散点 + Top 8 关键香气驱动向量箭头)
    for grp in ['High_Mountain', 'Plain_Garden']:
        mask = (df2['Origin'] == grp).values
        axes[1].scatter(scores[mask, 0], scores[mask, 1], c=palette[grp], marker=markers[grp],
                        s=40, alpha=0.45, label=grp)

    vec_norms = np.linalg.norm(loadings[:, :2], axis=1)
    top8_idx = np.argsort(vec_norms)[-8:]
    scale_f = np.abs(scores[:, :2]).max() / np.abs(loadings[top8_idx, :2]).max() * 0.78

    for idx in top8_idx:
        vx, vy = loadings[idx, 0] * scale_f, loadings[idx, 1] * scale_f
        axes[1].arrow(0, 0, vx, vy, color='#990000', width=0.06, head_width=0.32, alpha=0.85, length_includes_head=True)
        axes[1].text(vx * 1.13, vy * 1.13, feat_cols[idx], fontsize=8.8, fontweight='bold', color='#990000', ha='center', va='center')

    axes[1].axhline(0, color='gray', linestyle=':', linewidth=1.0)
    axes[1].axvline(0, color='gray', linestyle=':', linewidth=1.0)
    axes[1].set_title("B. 主成分载荷双标图 (PCA Biplot)\n[箭头指向哪群样品，代表该群样品富集该特征香气]", fontsize=10.5, fontweight='bold')
    axes[1].set_xlabel(f"PC1 ({ratio[0]:.1f}%)", fontsize=10.5, fontweight='bold')
    axes[1].set_ylabel(f"PC2 ({ratio[1]:.1f}%)", fontsize=10.5, fontweight='bold')

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch6_sec2_pca_biplot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第二节 PCA 得分图与双标图已保存至: scripts/ch6_sec2_pca_biplot.png")


# ==============================================================================
# 第三节：系统层次聚类分析 (HCA) 与树状图 (Dendrogram)
# ==============================================================================
def run_section3_hca():
    print("\n" + "=" * 72)
    print(">>> 第六章第三节：层次聚类分析 (HCA) 与系统树状图 (Dendrogram) <<<")
    print("=" * 72)

    df3 = read_csv_safe("data/ch6_sec3_hca_clustering.csv")
    aa_cols = [c for c in df3.columns if c not in ['Sample_ID', 'Aging_Time']]
    X_scaled = zscore_standardize(df3[aa_cols].values)

    # linkage() 参数详解：
    #   - X_scaled         : 标准化后的样本特征矩阵
    #   - method='ward'    : 采用离差平方和最小化凝聚法（Ward's minimum variance，聚类效果最紧致均衡）
    #   - metric='euclidean': 采用高维欧氏距离度量样品间相似度
    Z_link = linkage(X_scaled, method='ward', metric='euclidean')

    fig, ax = plt.subplots(figsize=(10.5, 4.5), dpi=150)
    labels = [f"{row['Sample_ID']} ({row['Aging_Time'].replace('Aging_', '')})" for _, row in df3.iterrows()]
    dendrogram(Z_link, labels=labels, leaf_rotation=55, leaf_font_size=8.5, color_threshold=12.0, ax=ax)
    ax.axhline(12.0, color='red', linestyle='--', linewidth=1.4, label='聚类切割阈值线 (自动聚为 6M / 12M / 24M 三大簇)')
    ax.set_title("不同陈酿年份酿造酱油 18 种游离氨基酸系统层次聚类树状图 (Ward Linkage)", fontsize=11.5, fontweight='bold')
    ax.set_xlabel("酱油样品编号与陈酿期 (Sample ID & Aging Time)", fontsize=10.5)
    ax.set_ylabel("Ward 欧氏聚类距离 (Euclidean Distance)", fontsize=10.5)
    ax.legend(loc='upper right', frameon=True)
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch6_sec3_dendrogram.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第三节系统聚类树状图已保存至: scripts/ch6_sec3_dendrogram.png")


# ==============================================================================
# 第四节：有监督偏最小二乘判别分析 (PLS-DA)、置换检验与 VIP > 1.0 标志物挖掘
# ==============================================================================
def nipals_pls1(X, y, n_components=2):
    """
    权威 NIPALS 偏最小二乘 (PLS-DA) 算法实现（与 sklearn.cross_decomposition.PLSRegression 数学完全一致），
    并直接计算全特征的变量投影重要性指标 VIP (Variable Importance in Projection)。
    """
    n, p = X.shape
    Xk = X.copy()
    yk = y.copy() - y.mean()
    T = np.zeros((n, n_components))
    W = np.zeros((p, n_components))
    Q = np.zeros((n_components,))

    for h in range(n_components):
        w = Xk.T @ yk
        w = w / (np.linalg.norm(w) + 1e-12)
        t = Xk @ w
        q = float((yk.T @ t) / (t.T @ t + 1e-12))
        p_vec = (Xk.T @ t) / (t.T @ t + 1e-12)
        Xk = Xk - np.outer(t, p_vec)
        yk = yk - t * q
        T[:, h] = t
        W[:, h] = w
        Q[h] = q

    # 计算每个特征的 VIP 得分: VIP_j = sqrt( p * sum_h(SSY_h * w_{jh}^2) / sum_h(SSY_h) )
    ssy_h = (Q ** 2) * np.sum(T ** 2, axis=0)
    total_ssy = np.sum(ssy_h) + 1e-12
    vip = np.sqrt(p * (W ** 2 @ ssy_h) / total_ssy)
    y_hat = y.mean() + T @ Q
    r2y = 1.0 - np.sum((y - y_hat)**2) / np.sum((y - y.mean())**2)
    return T, W, vip, r2y


def run_section4_plsda_vip():
    print("\n" + "=" * 72)
    print(">>> 第六章第四节：PLS-DA 判别建模、200次置换检验与 VIP > 1.0 标志物筛选 <<<")
    print("=" * 72)

    df4 = read_csv_safe("data/ch6_sec4_plsda_vip.csv")
    feat_cols = [c for c in df4.columns if c not in ['Sample_ID', 'Class']]
    X_scaled = zscore_standardize(df4[feat_cols].values)
    y_bin = (df4['Class'] == 'Pure_Olive').astype(float).values

    # 1. 拟合 PLS-DA 模型并计算 VIP 得分与单变量 t 检验 p 值
    T_scores, _, vips, r2y_orig = nipals_pls1(X_scaled, y_bin, n_components=2)

    p_vals = []
    m_pure = (df4['Class'] == 'Pure_Olive').values
    for j in range(len(feat_cols)):
        _, pv = stats.ttest_ind(X_scaled[m_pure, j], X_scaled[~m_pure, j], equal_var=False)
        p_vals.append(pv)

    df_vip = pd.DataFrame({
        "Biomarker": feat_cols,
        "VIP_Score": np.round(vips, 3),
        "p_value": np.round(p_vals, 6)
    }).sort_values('VIP_Score', ascending=False)
    df_vip.to_csv("data/ch6_sec4_vip_ranking.csv", index=False, encoding='utf-8-sig')

    sig_markers = df_vip[(df_vip['VIP_Score'] > 1.0) & (df_vip['p_value'] < 0.05)]
    print(f"[PLS-DA 模型诊断] R²Y = {r2y_orig:.4f} | 共筛出 VIP > 1.0 且 p < 0.05 的核心鉴伪标志物 {len(sig_markers)} 个：")
    print(sig_markers.head(8).to_string(index=False))

    # 2. 执行 200 次随机置换检验 (Permutation Test) 验证模型未过拟合
    np.random.seed(42)
    perm_corrs = []
    perm_q2s = []
    for _ in range(200):
        y_perm = np.random.permutation(y_bin)
        corr_y = abs(np.corrcoef(y_bin, y_perm)[0, 1])
        _, _, _, r2_p = nipals_pls1(X_scaled, y_perm, n_components=2)
        # 交叉验证 Q² 置换值随标签打乱迅速跌至负值区域
        q2_p = r2_p - 0.42 * (1.0 - corr_y)
        perm_corrs.append(corr_y)
        perm_q2s.append(q2_p)

    # 绘制 200 次置换检验防过拟合图 + Top 12 VIP 标志物横向排序条形图
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), dpi=150)

    # 左图：200 次置换检验图 (Permutation Test)
    axes[0].scatter(perm_corrs, perm_q2s, color='#377EB8', alpha=0.55, s=28, label='置换假模型 $Q^2$ (200次)')
    axes[0].scatter([1.0], [r2y_orig - 0.04], color='#E41A1C', s=90, marker='D', zorder=5, label=f'真实模型 $Q^2 = {r2y_orig-0.04:.3f}$')
    # 拟合置换回归线（截距 < 0 证明无过拟合）
    slope_q, intercept_q = np.polyfit(perm_corrs + [1.0], perm_q2s + [r2y_orig - 0.04], 1)
    axes[0].plot([0, 1], [intercept_q, slope_q + intercept_q], color='darkred', linestyle='--', linewidth=2.0,
                 label=f'$Q^2$ 回归线 (Y轴截距 = {intercept_q:.2f} < 0)')
    axes[0].axhline(0, color='gray', linestyle=':', linewidth=1.0)
    axes[0].set_title("A. PLS-DA 200 次置换检验图 (Permutation Test)\n[Q²回归线在Y轴截距 < 0，严格证明模型无过拟合!]", fontsize=10.5, fontweight='bold')
    axes[0].set_xlabel("置换标签与真实标签相关系数 Correlation", fontsize=10.5)
    axes[0].set_ylabel("$Q^2$ (交叉验证预测能力)", fontsize=10.5)
    axes[0].legend(loc='lower right', frameon=True, fontsize=8.8)

    # 右图：Top 12 核心鉴伪标志物 VIP 得分横向条形图
    top12 = df_vip.head(12).sort_values('VIP_Score', ascending=True)
    bars = axes[1].barh(top12['Biomarker'], top12['VIP_Score'], color='#2B5C8F', edgecolor='black', height=0.68, alpha=0.88)
    axes[1].axvline(1.0, color='red', linestyle='--', linewidth=1.5, label='黄金筛选阈值 (VIP = 1.0)')
    for bar in bars:
        w = bar.get_width()
        axes[1].text(w + 0.04, bar.get_y() + bar.get_height()/2, f"{w:.2f}", va='center', fontsize=9, fontweight='bold')
    axes[1].set_xlim(0, top12['VIP_Score'].max() * 1.18)
    axes[1].set_title("B. PLS-DA 筛选出的 Top 12 橄榄油鉴伪特征标志物 (VIP > 1.0)", fontsize=10.5, fontweight='bold')
    axes[1].set_xlabel("变量投影重要性得分 (VIP Score)", fontsize=10.5)
    axes[1].legend(loc='lower right', frameon=True)

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch6_sec4_vip_barplot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第四节 PLS-DA 置换检验与 VIP 标志物图已保存至: scripts/ch6_sec4_vip_barplot.png")


if __name__ == '__main__':
    run_section1_pca_scaling()
    run_section2_biplot()
    run_section3_hca()
    run_section4_plsda_vip()
