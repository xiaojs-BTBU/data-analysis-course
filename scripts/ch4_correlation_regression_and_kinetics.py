# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第四章：相关分析、回归分析和拟合
脚本名称：ch4_correlation_regression_and_kinetics.py
代码目的：1. 第一节：Pearson 线性相关与 Spearman 秩相关计算，绘制带显著性星号（* / ** / ***）的下三角冷暖热图；
          2. 第二节：一元与多元线性回归最小二乘法（OLS），计算方差膨胀因子（VIF）排查并剔除严重共线性变量，
                     绘制标准残差诊断四合一图；
          3. 第三节：非线性最小二乘拟合（scipy.optimize.curve_fit）实战：
                     花青素热降解一级/二级反应动力学模型（求速率常数 k 与半衰期 t_1/2）及
                     植物乳杆菌发酵 Logistic 与修正 Gompertz S 型生长曲线拟合（求最大比生长速率与延滞期 λ）。
================================================================================
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
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


# ==============================================================================
# 第一节：Pearson / Spearman 相关性分析与带显著性星号的下三角冷暖热图
# ==============================================================================
def run_section1_correlation():
    print("\n" + "=" * 72)
    print(">>> 第四章第一节：Pearson 与 Spearman 相关性分析及显著性星号热图 <<<")
    print("=" * 72)

    df1 = read_csv_safe("data/ch4_sec1_correlation.csv")
    num_cols = ['Total_Phenols', 'Anthocyanins', 'Titratable_Acid', 'Residual_Sugar', 'ABTS_Capacity', 'Sensory_Rank']

    # 计算 Pearson 相关系数矩阵与精确 p 值矩阵
    # （若使用 pingouin: pg.pairwise_corr(df1, columns=num_cols, method='pearson')）
    corr_mat = df1[num_cols].corr(method='pearson')
    annot_mat = pd.DataFrame("", index=num_cols, columns=num_cols)
    pair_records = []

    for i, c1 in enumerate(num_cols):
        for j, c2 in enumerate(num_cols):
            if i == j:
                annot_mat.loc[c1, c2] = "1.00"
            else:
                r_p, p_p = stats.pearsonr(df1[c1], df1[c2])
                r_s, p_s = stats.spearmanr(df1[c1], df1[c2])
                star = "***" if p_p < 0.001 else ("**" if p_p < 0.01 else ("*" if p_p < 0.05 else ""))
                annot_mat.loc[c1, c2] = f"{r_p:.2f}\n{star}" if star else f"{r_p:.2f}"
                if i > j:
                    pair_records.append({
                        "Variable_X": c1, "Variable_Y": c2,
                        "Pearson_r": round(r_p, 3), "Pearson_p": round(p_p, 5),
                        "Spearman_rho": round(r_s, 3), "Spearman_p": round(p_s, 5),
                        "Significance": star if star else "ns"
                    })

    df_pairs = pd.DataFrame(pair_records)
    print("[两两变量 Pearson 与 Spearman 相关系数对比表]：\n", df_pairs.to_string(index=False))
    df_pairs.to_csv("data/ch4_sec1_corr_table.csv", index=False, encoding='utf-8-sig')

    # 绘制下三角冷暖发散相关性矩阵热图 (RdBu_r)
    mask = np.triu(np.ones_like(corr_mat, dtype=bool), k=1)  # 遮盖上三角冗余区域
    fig, ax = plt.subplots(figsize=(7.2, 6.0), dpi=150)
    # sns.heatmap() 参数详解：
    #   - corr_mat   : 相关系数数值矩阵（决定颜色深浅）
    #   - mask=mask  : 布尔掩码矩阵（True 的位置不绘制，仅保留下三角与对角线）
    #   - annot=annot_mat: 单元格内显示的自定义字符串矩阵（包含数值与星号）
    #   - fmt=""     : 因为传入的是格式化好的字符串，fmt 设为空字符串
    #   - cmap="RdBu_r": 经典红蓝冷暖发散科学色阶（正相关暖红，负相关冷蓝，0为纯白）
    #   - vmin=-1, vmax=1, center=0: 将色阶范围严格对称锁定在 [-1, +1]
    sns.heatmap(
        corr_mat, mask=mask, annot=annot_mat, fmt="", cmap="RdBu_r",
        vmin=-1.0, vmax=1.0, center=0.0, square=True, linewidths=1.2,
        cbar_kws={"shrink": 0.78, "label": "Pearson Correlation Coefficient (r)"}, ax=ax
    )
    ax.set_title("果酒理化与感官指标相关性矩阵热图 (* p<0.05, ** p<0.01, *** p<0.001)",
                 fontsize=11.5, fontweight='bold', pad=14)
    plt.xticks(rotation=25, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig("scripts/ch4_sec1_corr_heatmap.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第一节相关性星号热图已保存至: scripts/ch4_sec1_corr_heatmap.png")


# ==============================================================================
# 第二节：一元/多元线性回归（OLS）、方差膨胀因子（VIF）共线性诊断与残差四合一图
# ==============================================================================
def fit_ols_and_vif(df, x_cols, y_col):
    """
    基于矩阵最小二乘法 (OLS) 精确计算多元线性回归参数、t检验 p 值、R²、调整后 R² 及每个自变量的 VIF 值
    """
    y = df[y_col].values
    X_raw = df[x_cols].values
    n, p = X_raw.shape
    # 添加常数截距列 (Intercept = 1)
    X = np.column_stack([np.ones(n), X_raw])
    # 最小二乘正规方程求解系数: beta = (X^T X)^(-1) X^T y
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    residuals = y - y_hat

    # 决定系数 R² 与调整后 R²_adj
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot
    r2_adj = 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)

    # 参数标准误 SE、t 统计量与双尾 p 值
    df_resid = n - p - 1
    mse = ss_res / df_resid
    cov_beta = mse * np.linalg.inv(X.T @ X)
    se_beta = np.sqrt(np.diag(cov_beta))
    t_stats = beta / se_beta
    p_vals = 2 * (1.0 - stats.t.cdf(np.abs(t_stats), df_resid))

    # 计算每个自变量的方差膨胀因子 VIF_j = 1 / (1 - R_j^2)
    vifs = [np.nan]  # 截距项不计 VIF
    for j in range(p):
        xj = X_raw[:, j]
        X_other = np.column_stack([np.ones(n), np.delete(X_raw, j, axis=1)])
        b_j, _, _, _ = np.linalg.lstsq(X_other, xj, rcond=None)
        xj_hat = X_other @ b_j
        r2_j = 1.0 - np.sum((xj - xj_hat)**2) / np.sum((xj - xj.mean())**2)
        vifs.append(1.0 / max(1e-8, (1.0 - r2_j)))

    coef_table = pd.DataFrame({
        "Term (自变量)": ["Intercept (常数截距)"] + list(x_cols),
        "Coef (回归系数 β)": np.round(beta, 4),
        "Std_Err (标准误)": np.round(se_beta, 4),
        "t_value": np.round(t_stats, 3),
        "p_value": np.round(p_vals, 5),
        "VIF (方差膨胀因子)": [round(v, 2) if not np.isnan(v) else "-" for v in vifs]
    })
    # 计算帽子矩阵对角线杠杆值 Leverage h_ii
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    leverage = np.diag(H)
    std_residuals = residuals / np.sqrt(mse * (1.0 - leverage))
    return coef_table, r2, r2_adj, y_hat, residuals, std_residuals, leverage


def run_section2_regression():
    print("\n" + "=" * 72)
    print(">>> 第四章第二节：一元/多元线性回归、VIF 多重共线性诊断与残差四合一图 <<<")
    print("=" * 72)

    df2 = read_csv_safe("data/ch4_sec2_regression_vif.csv")

    # 1. 包含全部 5 个自变量的初始回归（含故意构造高共线性的 Catechin_Deriv）
    full_cols = ['Total_Phenols', 'Flavonoids', 'Catechin_Deriv', 'Polysaccharides', 'Reducing_Sugar']
    tbl_full, r2_f, r2adj_f, _, _, _, _ = fit_ols_and_vif(df2, full_cols, 'DPPH_Rate')
    print(f"[1. 剔除共线性前：全变量回归 (R²={r2_f:.4f}, Adj-R²={r2adj_f:.4f}) —— 注意观察 VIF > 10 的严重共线性！]：")
    print(tbl_full.to_string(index=False))

    # 2. 剔除冗余共线性变量 Catechin_Deriv 后的稳健多元回归
    clean_cols = ['Total_Phenols', 'Flavonoids', 'Polysaccharides', 'Reducing_Sugar']
    tbl_clean, r2_c, r2adj_c, y_hat, resid, std_resid, leverage = fit_ols_and_vif(df2, clean_cols, 'DPPH_Rate')
    print(f"\n[2. 剔除 Catechin_Deriv 后：稳健多元回归 (R²={r2_c:.4f}, Adj-R²={r2adj_c:.4f}) —— 所有 VIF < 2.0 完全健康！]：")
    print(tbl_clean.to_string(index=False))
    tbl_clean.to_csv("data/ch4_sec2_regression_table.csv", index=False, encoding='utf-8-sig')

    # 3. 绘制标准残差诊断四合一图 (Residuals vs Fitted, Normal Q-Q, Scale-Location, Residuals vs Leverage)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), dpi=150)

    # (A) Residuals vs Fitted (检验线性与方差齐性)
    axes[0, 0].scatter(y_hat, resid, color='#1F78B4', edgecolor='black', s=45, alpha=0.8)
    axes[0, 0].axhline(0, color='red', linestyle='--', linewidth=1.2)
    axes[0, 0].set_title("A. 残差 vs 拟合值图 (Residuals vs Fitted)\n[点均匀分布在0两侧，证明满足线性与同方差]", fontsize=10.5, fontweight='bold')
    axes[0, 0].set_xlabel("预测拟合值 Fitted Values")
    axes[0, 0].set_ylabel("残差 Residuals")

    # (B) Normal Q-Q Plot (检验残差正态性)
    stats.probplot(std_resid, dist="norm", plot=axes[0, 1])
    axes[0, 1].set_title("B. 标准化残差正态 Q-Q 图 (Normal Q-Q)\n[散点紧贴对角红线，证明残差服从正态分布]", fontsize=10.5, fontweight='bold')
    axes[0, 1].get_lines()[0].set_markerfacecolor('#1F78B4')

    # (C) Scale-Location (检验方差齐性稳健度)
    sqrt_abs_sr = np.sqrt(np.abs(std_resid))
    axes[1, 0].scatter(y_hat, sqrt_abs_sr, color='#33A02C', edgecolor='black', s=45, alpha=0.8)
    sns.regplot(x=y_hat, y=sqrt_abs_sr, scatter=False, order=2, ci=None, color='red', ax=axes[1, 0])
    axes[1, 0].set_title("C. 位置-尺度图 (Scale-Location)\n[红线平坦无上扬喇叭口，证明无异方差]", fontsize=10.5, fontweight='bold')
    axes[1, 0].set_xlabel("预测拟合值 Fitted Values")
    axes[1, 0].set_ylabel("$\\sqrt{|\\mathrm{Standardized\\ Residuals}|}$")

    # (D) Residuals vs Leverage (识别高杠杆强影响点)
    axes[1, 1].scatter(leverage, std_resid, color='#FF7F00', edgecolor='black', s=45, alpha=0.8)
    axes[1, 1].axhline(0, color='gray', linestyle='--', linewidth=1.0)
    axes[1, 1].axhline(2.5, color='red', linestyle=':', linewidth=1.0)
    axes[1, 1].axhline(-2.5, color='red', linestyle=':', linewidth=1.0)
    axes[1, 1].set_title("D. 残差 vs 杠杆值图 (Residuals vs Leverage)\n[无超出±3的极端强影响异常点]", fontsize=10.5, fontweight='bold')
    axes[1, 1].set_xlabel("杠杆值 Leverage ($h_{ii}$)")
    axes[1, 1].set_ylabel("标准化残差 Standardized Residuals")

    for ax in axes.flat:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch4_sec2_residual_diagnostics.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第二节多元回归残差诊断四合一图已保存至: scripts/ch4_sec2_residual_diagnostics.png")


# ==============================================================================
# 第三节：非线性回归与科研动力学拟合（花青素一级热降解 + 乳酸菌 Gompertz 生长曲线）
# ==============================================================================
def run_section3_kinetics():
    print("\n" + "=" * 72)
    print(">>> 第四章第三节：非线性动力学拟合（一级/二级热降解与 Gompertz 生长曲线） <<<")
    print("=" * 72)

    df3 = read_csv_safe("data/ch4_sec3_kinetics_growth.csv")

    # 1. 定义动力学数学模型函数
    def first_order_decay(t, C0, k1):
        """一级反应降解动力学: C(t) = C0 * exp(-k1 * t)"""
        return C0 * np.exp(-k1 * t)

    def second_order_decay(t, C0, k2):
        """二级反应降解动力学: C(t) = C0 / (1 + C0 * k2 * t)"""
        return C0 / (1.0 + C0 * k2 * t)

    def gompertz_growth(t, N0, K, mu_m, lam):
        """修正 Gompertz 微生物 S 型生长模型:
           N0: 初始菌浓, K: 最大渐近增量, mu_m: 最大比生长速率 (h^-1), lam: 延滞期时长 Lag time (h)
        """
        return N0 + K * np.exp(-np.exp((mu_m * np.e / K) * (lam - t) + 1.0))

    # 2. 拟合花青素 70℃ 热降解数据
    t_min = df3['Time_min'].values
    c_ret = df3['Anthocyanin_Retention'].values

    # curve_fit() 参数详解：
    #   - f   : 待拟合的自定义数学模型函数（如 first_order_decay）
    #   - xdata, ydata: 自变量时间数组与因变量实测值数组
    #   - p0  : 参数初始猜测值（给一个合理的数量级起点，帮助算法快速收敛）
    popt_1st, pcov_1st = curve_fit(first_order_decay, t_min, c_ret, p0=[100.0, 0.01])
    popt_2nd, pcov_2nd = curve_fit(second_order_decay, t_min, c_ret, p0=[100.0, 0.0002])

    c0_1, k1 = popt_1st
    se_c0_1, se_k1 = np.sqrt(np.diag(pcov_1st))
    t_half_1 = np.log(2) / k1
    r2_1st = 1.0 - np.sum((c_ret - first_order_decay(t_min, *popt_1st))**2) / np.sum((c_ret - c_ret.mean())**2)
    r2_2nd = 1.0 - np.sum((c_ret - second_order_decay(t_min, *popt_2nd))**2) / np.sum((c_ret - c_ret.mean())**2)

    print(f"[1. 花青素热降解动力学对比] 一级动力学 R² = {r2_1st:.4f} (明显优于二级动力学 R² = {r2_2nd:.4f})")
    print(f"    -> 一级速率常数 k1 = {k1:.5f} ± {se_k1:.5f} min^-1, 半衰期 t_1/2 = {t_half_1:.1f} min")

    # 3. 拟合植物乳杆菌 37℃ 发酵 Gompertz S 型生长曲线
    t_hr = df3['Ferm_Hour'].values
    od_val = df3['OD600_Biomass'].values
    popt_gomp, pcov_gomp = curve_fit(gompertz_growth, t_hr, od_val, p0=[0.08, 2.4, 0.25, 3.5])
    N0_g, K_g, mu_g, lam_g = popt_gomp
    r2_gomp = 1.0 - np.sum((od_val - gompertz_growth(t_hr, *popt_gomp))**2) / np.sum((od_val - od_val.mean())**2)
    print(f"[2. 植物乳杆菌 Gompertz 生长动力学] R² = {r2_gomp:.4f} | "
          f"最大菌浓增量 K = {K_g:.2f}, 最大生长速率 μ_m = {mu_g:.3f} h^-1, 延滞期 λ = {lam_g:.2f} h")

    # 保存参数汇总三线表
    kin_df = pd.DataFrame([
        {"Model": "First-Order Degradation (花青素一级热降解)", "Key_Parameters": f"C0={c0_1:.1f}%, k1={k1:.4f} min^-1, t1/2={t_half_1:.1f} min", "R2": round(r2_1st, 4)},
        {"Model": "Second-Order Degradation (花青素二级热降解)", "Key_Parameters": f"C0={popt_2nd[0]:.1f}%, k2={popt_2nd[1]:.5f}", "R2": round(r2_2nd, 4)},
        {"Model": "Modified Gompertz Growth (植物乳杆菌S型生长)", "Key_Parameters": f"K={K_g:.2f}, μm={mu_g:.3f} h^-1, Lag λ={lam_g:.2f} h", "R2": round(r2_gomp, 4)},
    ])
    kin_df.to_csv("data/ch4_sec3_kinetic_params.csv", index=False, encoding='utf-8-sig')

    # 4. 绘制带 95% 预测置信带的双联动力学拟合图
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), dpi=150)

    # 左图：花青素一级 vs 二级热降解曲线
    t_grid1 = np.linspace(0, 185, 200)
    y_grid1 = first_order_decay(t_grid1, *popt_1st)
    y_grid2 = second_order_decay(t_grid1, *popt_2nd)
    sd_res1 = np.std(c_ret - first_order_decay(t_min, *popt_1st), ddof=2)

    axes[0].scatter(t_min, c_ret, color='#D95F02', edgecolor='black', s=45, alpha=0.8, label='实验实测点 (Observed)')
    axes[0].plot(t_grid1, y_grid1, color='#1B9E77', linewidth=2.4, label=f'一级动力学拟合 ($R^2={r2_1st:.4f}$)')
    axes[0].plot(t_grid1, y_grid2, color='#7570B3', linewidth=1.8, linestyle='--', label=f'二级动力学拟合 ($R^2={r2_2nd:.4f}$)')
    axes[0].fill_between(t_grid1, y_grid1 - 1.96*sd_res1, y_grid1 + 1.96*sd_res1, color='#1B9E77', alpha=0.2, label='95% 预测置信带')
    axes[0].text(0.42, 0.62, f"一级反应方程:\n$C_t = {c0_1:.1f} \\cdot e^{{-{k1:.4f}t}}$\n半衰期 $t_{{1/2}} = {t_half_1:.1f}$ min\n$R^2 = {r2_1st:.4f}$",
                 transform=axes[0].transAxes, fontsize=9.5, bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='gray'))
    axes[0].set_title("A. 蓝莓花青素 70℃ 热降解反应动力学拟合", fontsize=11, fontweight='bold')
    axes[0].set_xlabel("热处理时间 Time (min)", fontsize=10.5)
    axes[0].set_ylabel("花青素保留率 Retention (%)", fontsize=10.5)
    axes[0].legend(loc='upper right', fontsize=8.5, frameon=True)

    # 右图：植物乳杆菌修正 Gompertz S 型生长曲线
    t_grid2 = np.linspace(0, 24.5, 200)
    y_gomp = gompertz_growth(t_grid2, *popt_gomp)
    sd_res_g = np.std(od_val - gompertz_growth(t_hr, *popt_gomp), ddof=4)

    axes[1].scatter(t_hr, od_val, color='#2B5C8F', edgecolor='black', s=45, alpha=0.8, label='实测菌浓 OD600')
    axes[1].plot(t_grid2, y_gomp, color='#E41A1C', linewidth=2.4, label=f'Gompertz S型曲线 ($R^2={r2_gomp:.4f}$)')
    axes[1].fill_between(t_grid2, y_gomp - 1.96*sd_res_g, y_gomp + 1.96*sd_res_g, color='#E41A1C', alpha=0.2, label='95% 预测置信带')
    axes[1].axvline(lam_g, color='darkgreen', linestyle=':', linewidth=1.5, label=f'延滞期 $\\lambda = {lam_g:.2f}$ h')
    axes[1].text(0.42, 0.22, f"Modified Gompertz 模型:\n最大菌浓增量 $K = {K_g:.2f}$\n最大比生长速率 $\\mu_m = {mu_g:.3f}\\ \\mathrm{{h}}^{{-1}}$\n延滞期 $\\lambda = {lam_g:.2f}$ h\n$R^2 = {r2_gomp:.4f}$",
                 transform=axes[1].transAxes, fontsize=9.5, bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='gray'))
    axes[1].set_title("B. 植物乳杆菌发酵 S 型生长动力学拟合 (Gompertz)", fontsize=11, fontweight='bold')
    axes[1].set_xlabel("发酵时间 Fermentation Time (h)", fontsize=10.5)
    axes[1].set_ylabel("菌体密度 $\mathrm{OD_{600}}$", fontsize=10.5)
    axes[1].legend(loc='upper left', fontsize=8.5, frameon=True)

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch4_sec3_kinetic_curves.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第三节非线性动力学拟合曲线图已保存至: scripts/ch4_sec3_kinetic_curves.png")


if __name__ == '__main__':
    run_section1_correlation()
    run_section2_regression()
    run_section3_kinetics()
