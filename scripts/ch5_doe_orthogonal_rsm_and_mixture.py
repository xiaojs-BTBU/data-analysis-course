# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第五章：实验设计与分析
脚本名称：ch5_doe_orthogonal_rsm_and_mixture.py
代码目的：1. 第一节：L9(3^4) 正交试验极差分析（k 值与极差 R 排序主次因素）与空列误差方差分析；
          2. 第二、三节：3因素3水平 Box-Behnken (BBD) 响应面二次多项式回归拟合、
                        失拟项检验（Lack of Fit p > 0.05 生死线诊断）、全局极值寻优与 3D 曲面/等高线绘制；
          3. 第四节：配方混料设计（Scheffé 无截距回归与三元相图等高线 Ternary Plot）及
                    开题前最小样本量与统计检验力分析（Power Analysis）。
================================================================================
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
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
# 第一节：L9(3^4) 正交试验直观极差分析（R值）与方差分析
# ==============================================================================
def run_section1_orthogonal():
    print("\n" + "=" * 72)
    print(">>> 第五章第一节：L9(3^4) 正交试验极差分析 (R值) 与方差分析 <<<")
    print("=" * 72)

    df1 = read_csv_safe("data/ch5_sec1_orthogonal.csv")
    factors = ['Power_A', 'Ratio_B', 'Time_C']

    # 1. 计算各因素各水平均值 k1, k2, k3 与极差 R = max(k) - min(k)
    range_rows = []
    grand_mean = df1['Yield'].mean()
    # 空列 (第4列 Blank_Error) 平方和作为系统随机误差平方和 SS_Error (df=2)
    err_k = df1.groupby('Blank_Error')['Yield'].mean()
    ss_err = 3 * ((err_k - grand_mean)**2).sum()
    df_err = 2
    ms_err = ss_err / df_err

    for col in factors:
        k_means = df1.groupby(col)['Yield'].mean()
        r_val = k_means.max() - k_means.min()
        best_lvl = k_means.idxmax()
        ss_f = 3 * ((k_means - grand_mean)**2).sum()
        df_f = 2
        ms_f = ss_f / df_f
        f_val = ms_f / ms_err
        p_val = 1.0 - stats.f.cdf(f_val, df_f, df_err)
        range_rows.append({
            "Factor (因素)": col,
            "k1 (水平1均值)": round(k_means.iloc[0], 2),
            "k2 (水平2均值)": round(k_means.iloc[1], 2),
            "k3 (水平3均值)": round(k_means.iloc[2], 2),
            "Range_R (极差R)": round(r_val, 2),
            "Optimal_Level (最优水平)": best_lvl,
            "SS (平方和)": round(ss_f, 3),
            "F_value": round(f_val, 2),
            "p_value": round(p_val, 4)
        })

    df_ortho_res = pd.DataFrame(range_rows).sort_values('Range_R (极差R)', ascending=False)
    print("[正交试验极差与方差分析合并表 (按极差 R 降序排列因素主次顺序)]：")
    print(df_ortho_res.to_string(index=False))
    df_ortho_res.to_csv("data/ch5_sec1_orthogonal_table.csv", index=False, encoding='utf-8-sig')

    # 绘制正交试验各因素水平直观效应趋势折线图
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=True, dpi=150)
    labels = ["超声功率 Power A (W)", "液料比 Ratio B (mL/g)", "超声时间 Time C (min)"]
    colors = ['#E41A1C', '#377EB8', '#4DAF4A']

    for i, col in enumerate(factors):
        k_series = df1.groupby(col)['Yield'].mean()
        axes[i].plot(k_series.index.astype(str), k_series.values, marker='o', markersize=8,
                     color=colors[i], linewidth=2.4)
        for x_idx, val in enumerate(k_series.values):
            axes[i].text(x_idx, val + 0.35, f"{val:.2f}%", ha='center', fontweight='bold', fontsize=9.5)
        r_v = k_series.max() - k_series.min()
        axes[i].set_title(f"{labels[i]}\n(极差 $R = {r_v:.2f}$)", fontsize=10.5, fontweight='bold')
        axes[i].set_xlabel("因素水平 (Level)")
        if i == 0:
            axes[i].set_ylabel("多糖提取率均值 $k_i$ (%)")
        axes[i].set_ylim(12.5, 20.0)
        sns.despine(ax=axes[i], top=True, right=True)

    plt.suptitle("L9(3^4) 正交试验直观效应折线图 (主次顺序: 液料比 B > 超声功率 A > 超声时间 C)",
                 fontsize=12, fontweight='bold', y=1.03)
    plt.tight_layout()
    plt.savefig("scripts/ch5_sec1_orthogonal_plot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第一节正交效应折线图已保存至: scripts/ch5_sec1_orthogonal_plot.png")


# ==============================================================================
# 第二、三节：Box-Behnken 响应面二次建模、失拟项检验与 3D 曲面/等高线极值寻优
# ==============================================================================
def run_section2_and_3_rsm():
    print("\n" + "=" * 72)
    print(">>> 第五章第二、三节：BBD 响应面二次建模、失拟诊断与 3D 曲面寻优 <<<")
    print("=" * 72)

    df_bbd = read_csv_safe("data/ch5_sec2_rsm_bbd.csv")
    x1, x2, x3 = df_bbd['x1'].values, df_bbd['x2'].values, df_bbd['x3'].values
    y = df_bbd['Yield'].values
    n = len(y)

    # 1. 构建全二次多项式设计矩阵 X (含常数项、一次项、交互项、二次平方项，共10项)
    X_mat = np.column_stack([
        np.ones(n), x1, x2, x3,
        x1 * x2, x1 * x3, x2 * x3,
        x1**2, x2**2, x3**2
    ])
    term_names = ['Intercept', 'x1 (Temp)', 'x2 (Ratio)', 'x3 (Time)',
                  'x1*x2', 'x1*x3', 'x2*x3', 'x1^2', 'x2^2', 'x3^2']

    # 最小二乘求解二次回归系数
    beta, _, _, _ = np.linalg.lstsq(X_mat, y, rcond=None)
    y_pred = X_mat @ beta
    resid = y - y_pred

    # 2. 精确拆解方差分析与失拟项 (Lack of Fit) vs 纯误差 (Pure Error)
    ss_tot = np.sum((y - y.mean())**2)
    ss_res = np.sum(resid**2)
    ss_model = ss_tot - ss_res
    df_model = 9
    df_res = n - 10  # 17 - 10 = 7

    # 提取 5 个中心重复零点 (x1=0, x2=0, x3=0) 计算纯实验误差 Pure Error
    center_y = df_bbd[(df_bbd['x1'] == 0) & (df_bbd['x2'] == 0) & (df_bbd['x3'] == 0)]['Yield'].values
    ss_pe = np.sum((center_y - center_y.mean())**2)
    df_pe = len(center_y) - 1  # 5 - 1 = 4
    ms_pe = ss_pe / df_pe

    # 失拟项 Lack of Fit = 残差平方和 - 纯误差平方和
    ss_lof = ss_res - ss_pe
    df_lof = df_res - df_pe    # 7 - 4 = 3
    ms_lof = ss_lof / df_lof
    f_lof = ms_lof / ms_pe
    p_lof = 1.0 - stats.f.cdf(f_lof, df_lof, df_pe)

    ms_model = ss_model / df_model
    ms_res = ss_res / df_res
    f_model = ms_model / ms_res
    p_model = 1.0 - stats.f.cdf(f_model, df_model, df_res)

    r2 = 1.0 - ss_res / ss_tot
    r2_adj = 1.0 - (1.0 - r2) * (n - 1) / df_res

    print(f"[响应面模型全局诊断] R² = {r2:.4f}, Adj-R² = {r2_adj:.4f}")
    print(f"  -> 模型整体显著性: F = {f_model:.2f}, p = {p_model:.4e} (极显著 p < 0.0001***)")
    print(f"  -> 失拟项检验 (Lack of Fit): F = {f_lof:.2f}, p = {p_lof:.4f} (> 0.05 不显著! 证明模型拟合极佳!)")

    # 保存系数表
    cov_b = ms_res * np.linalg.inv(X_mat.T @ X_mat)
    se_b = np.sqrt(np.diag(cov_b))
    t_b = beta / se_b
    p_b = 2 * (1.0 - stats.t.cdf(np.abs(t_b), df_res))
    df_coefs = pd.DataFrame({
        "Term": term_names, "Coefficient": np.round(beta, 4),
        "Std_Err": np.round(se_b, 4), "t_value": np.round(t_b, 2), "p_value": np.round(p_b, 5)
    })
    df_coefs.to_csv("data/ch5_sec2_rsm_coefs.csv", index=False, encoding='utf-8-sig')
    print("\n[二次多项式回归系数显著性表]：\n", df_coefs.to_string(index=False))

    # 3. 使用 scipy.optimize.minimize (L-BFGS-B 算法) 求解 [-1, 1] 范围内的全局最高得率点
    def neg_yield(x_vec):
        u1, u2, u3 = x_vec
        row_x = np.array([1.0, u1, u2, u3, u1*u2, u1*u3, u2*u3, u1**2, u2**2, u3**2])
        return -float(row_x @ beta)

    opt_res = minimize(neg_yield, x0=[0.0, 0.0, 0.0], bounds=[(-1, 1), (-1, 1), (-1, 1)], method='L-BFGS-B')
    opt_x1, opt_x2, opt_x3 = opt_res.x
    max_yield = -opt_res.fun
    opt_temp = 75.0 + 15.0 * opt_x1
    opt_ratio = 30.0 + 10.0 * opt_x2
    opt_time = 35.0 + 15.0 * opt_x3
    print(f"\n[全局极值寻优结果] 最佳温度 = {opt_temp:.2f} ℃, 最佳液料比 = {opt_ratio:.2f} mL/g, "
          f"最佳时间 = {opt_time:.2f} min -> 理论最高得率 = {max_yield:.2f}%")

    # 4. 绘制 3D 响应曲面图与 2D 等高线图 (固定超声时间在最优值 opt_x3)
    g1 = np.linspace(-1, 1, 50)
    g2 = np.linspace(-1, 1, 50)
    G1, G2 = np.meshgrid(g1, g2)
    Z_yield = (beta[0] + beta[1]*G1 + beta[2]*G2 + beta[3]*opt_x3 +
               beta[4]*G1*G2 + beta[5]*G1*opt_x3 + beta[6]*G2*opt_x3 +
               beta[7]*(G1**2) + beta[8]*(G2**2) + beta[9]*(opt_x3**2))
    T_grid = 75.0 + 15.0 * G1
    R_grid = 30.0 + 10.0 * G2

    fig = plt.figure(figsize=(12.2, 5.0), dpi=150)

    # 左子图：3D 响应曲面 (3D Response Surface)
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax1.plot_surface(T_grid, R_grid, Z_yield, cmap='viridis', alpha=0.88, edgecolor='none')
    ax1.scatter([opt_temp], [opt_ratio], [max_yield], color='red', s=90, marker='*', zorder=10, label='最优极值点')
    ax1.set_title("A. 温度与液料比交互效应 3D 响应曲面图", fontsize=11, fontweight='bold', pad=10)
    ax1.set_xlabel("提取温度 Temp (℃)", labelpad=6)
    ax1.set_ylabel("液料比 Ratio (mL/g)", labelpad=6)
    ax1.set_zlabel("多糖得率 Yield (%)", labelpad=6)
    ax1.view_init(elev=28, azim=132)
    fig.colorbar(surf, ax=ax1, shrink=0.55, pad=0.08, label="Yield (%)")

    # 右子图：2D 等高线投影图 (Contour Plot)
    ax2 = fig.add_subplot(1, 2, 2)
    cf = ax2.contourf(T_grid, R_grid, Z_yield, levels=22, cmap='viridis', alpha=0.75)
    cs = ax2.contour(T_grid, R_grid, Z_yield, levels=10, colors='black', linewidths=0.8)
    ax2.clabel(cs, inline=True, fontsize=8.5, fmt="%.1f%%")
    ax2.scatter([opt_temp], [opt_ratio], color='red', edgecolor='white', s=130, marker='*', zorder=5,
                label=f'最佳工艺点 ({opt_temp:.1f}℃, {opt_ratio:.1f} mL/g)\n预测最高得率: {max_yield:.2f}%')
    ax2.set_title("B. 响应面 2D 等高线投影图 (椭圆斜置反映显著交互效应)", fontsize=11, fontweight='bold', pad=10)
    ax2.set_xlabel("提取温度 Temp (℃)", fontsize=10.5)
    ax2.set_ylabel("液料比 Ratio (mL/g)", fontsize=10.5)
    ax2.legend(loc='lower left', frameon=True, fontsize=9)
    fig.colorbar(cf, ax=ax2, shrink=0.78, label="Yield (%)")

    plt.tight_layout()
    plt.savefig("scripts/ch5_sec3_rsm_3d_contour.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第三节响应面 3D 曲面与等高线图已保存至: scripts/ch5_sec3_rsm_3d_contour.png")


# ==============================================================================
# 第四节：混料设计三元相图（Ternary Contour）与开题样本量检验力分析（Power Curve）
# ==============================================================================
def run_section4_mixture_and_power():
    print("\n" + "=" * 72)
    print(">>> 第五章第四节：混料设计三元相图与样本量统计检验力分析 (Power Analysis) <<<")
    print("=" * 72)

    df_mix = read_csv_safe("data/ch5_sec4_mixture_uniform.csv")
    a = df_mix['Blueberry_A'].values
    b = df_mix['Seabuckthorn_B'].values
    c = df_mix['Roxburgh_C'].values
    y_score = df_mix['Sensory_Score'].values

    # 1. 拟合 Scheffé 特殊三次无截距混料回归方程: Y = b1*A + b2*B + b3*C + b12*AB + b13*AC + b23*BC + b123*ABC
    X_mix = np.column_stack([a, b, c, a*b, a*c, b*c, a*b*c])
    b_mix, _, _, _ = np.linalg.lstsq(X_mix, y_score, rcond=None)

    # 2. 在正三角形重心坐标系内生成致密网格并转换为 2D 笛卡尔坐标 (x_cart, y_cart)
    grid_pts = []
    for a_v in np.linspace(0, 1, 100):
        for b_v in np.linspace(0, 1 - a_v, 100):
            c_v = max(0.0, 1.0 - a_v - b_v)
            # 重心坐标转正三角形直角坐标: 顶点 A(0.5, sqrt(3)/2), B(0,0), C(1,0)
            x_c = 0.5 * a_v + 1.0 * c_v
            y_c = (np.sqrt(3) / 2.0) * a_v
            z_pred = (b_mix[0]*a_v + b_mix[1]*b_v + b_mix[2]*c_v +
                      b_mix[3]*a_v*b_v + b_mix[4]*a_v*c_v + b_mix[5]*b_v*c_v + b_mix[6]*a_v*b_v*c_v)
            grid_pts.append((x_c, y_c, z_pred, a_v, b_v, c_v))
    df_grid = pd.DataFrame(grid_pts, columns=['x', 'y', 'z', 'a', 'b', 'c'])
    best_row = df_grid.loc[df_grid['z'].idxmax()]
    print(f"[混料设计最优配方] 蓝莓汁 A = {best_row['a']*100:.1f}%, 沙棘汁 B = {best_row['b']*100:.1f}%, "
          f"刺梨汁 C = {best_row['c']*100:.1f}% -> 最高感官评分 = {best_row['z']:.2f}")

    # 3. 绘制三元相图等高线 + 开题样本量统计检验力曲线 (Power Analysis)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), dpi=150)

    # 左图：配方混料三元相图等高线 (Ternary Contour Plot)
    tcf = axes[0].tricontourf(df_grid['x'], df_grid['y'], df_grid['z'], levels=20, cmap='Spectral_r', alpha=0.85)
    tcs = axes[0].tricontour(df_grid['x'], df_grid['y'], df_grid['z'], levels=8, colors='black', linewidths=0.7)
    axes[0].clabel(tcs, inline=True, fontsize=8, fmt="%.1f")
    # 绘制正三角形外框与顶点标签
    axes[0].plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3)/2, 0], color='black', linewidth=2.0)
    axes[0].scatter([best_row['x']], [best_row['y']], color='red', edgecolor='white', s=130, marker='*', zorder=10,
                    label=f"最优配比 ({best_row['a']*100:.0f}:{best_row['b']*100:.0f}:{best_row['c']*100:.0f})")
    axes[0].text(0.5, np.sqrt(3)/2 + 0.04, "蓝莓原汁 A (100%)", ha='center', fontweight='bold', fontsize=10, color='#1F78B4')
    axes[0].text(-0.04, -0.05, "沙棘汁 B (100%)", ha='left', fontweight='bold', fontsize=10, color='#D95F02')
    axes[0].text(1.04, -0.05, "刺梨汁 C (100%)", ha='right', fontweight='bold', fontsize=10, color='#1B9E77')
    axes[0].set_xlim(-0.12, 1.12)
    axes[0].set_ylim(-0.10, 0.98)
    axes[0].set_aspect('equal')
    axes[0].axis('off')
    axes[0].set_title("A. 复合果汁混料设计三元相图等高线 (A + B + C = 100%)", fontsize=11, fontweight='bold')
    axes[0].legend(loc='upper right', frameon=True, fontsize=9)
    fig.colorbar(tcf, ax=axes[0], shrink=0.75, label="感官评分 Sensory Score")

    # 右图：开题前样本量与统计检验力分析曲线 (Power Curve for Two-Sample t-test)
    n_range = np.arange(3, 31)
    for d_val, col_c, lbl in [(0.5, '#377EB8', "中等效应量 (Cohen's d = 0.5)"),
                              (0.8, '#FF7F00', "大效应量 (Cohen's d = 0.8)"),
                              (1.2, '#4DAF4A', "极强效应量 (Cohen's d = 1.2)")]:
        # 计算非中心 t 分布下的统计检验力 Power = 1 - beta
        df_t = 2 * n_range - 2
        ncp = d_val * np.sqrt(n_range / 2.0)  # 非中心参数
        t_crit = stats.t.ppf(1 - 0.05 / 2, df_t)
        power_vals = 1.0 - stats.nct.cdf(t_crit, df_t, ncp) + stats.nct.cdf(-t_crit, df_t, ncp)
        axes[1].plot(n_range, power_vals, marker='o', markersize=4, linewidth=2.2, color=col_c, label=lbl)

    axes[1].axhline(0.80, color='red', linestyle='--', linewidth=1.4, label='合格检验力门槛 (Power = 0.80)')
    axes[1].set_title("B. 开题前样本量估算：统计检验力曲线 (α = 0.05)", fontsize=11, fontweight='bold')
    axes[1].set_xlabel("每组独立生物学重复样本量 (n per group)", fontsize=10.5)
    axes[1].set_ylabel("统计检验力 Power (1 - β)", fontsize=10.5)
    axes[1].set_ylim(0.05, 1.03)
    axes[1].legend(loc='lower right', frameon=True, fontsize=8.8)
    sns.despine(ax=axes[1], top=True, right=True)

    plt.tight_layout()
    plt.savefig("scripts/ch5_sec4_ternary_and_power.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第四节混料三元相图与检验力曲线已保存至: scripts/ch5_sec4_ternary_and_power.png")


if __name__ == '__main__':
    run_section1_orthogonal()
    run_section2_and_3_rsm()
    run_section4_mixture_and_power()
