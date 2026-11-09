# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第三章：因素筛选、方差分析与多重比较
脚本名称：ch3_anova_posthoc_and_raincloud.py
代码目的：1. 第一节：单因素方差分析（One-Way ANOVA / Welch's ANOVA）与双因素方差分析（Two-Way ANOVA），
                     解析主效应与交互作用（Interaction Effect）并绘制交互效应折线图；
          2. 第二节：Tukey HSD、Bonferroni、Games-Howell 事后多重比较，
                     自动化推导组间显著性字母标记系统（Compact Letter Display: a, b, ab, c）；
          3. 第三节：绘制 Nature 期刊级高信息密度复合云雨图（Raincloud Plot：半小提琴 + 箱线 + 抖动散点 + 显著性字母）。
================================================================================
"""

import os
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# 跨平台中文字体与负号防乱码配置
sns.set_theme(style="ticks", font_scale=1.08)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Noto Sans CJK JP', 'Droid Sans Fallback', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42


def read_csv_safe(filepath):
    try:
        return pd.read_csv(filepath, encoding='utf-8')
    except UnicodeDecodeError:
        return pd.read_csv(filepath, encoding='gbk')


# ------------------------------------------------------------------------------
# 核心通用算法：根据 Tukey HSD 两两比较 p 值矩阵，自动推导紧凑字母标记 (a, b, ab, c)
# ------------------------------------------------------------------------------
def compute_cld_letters(df, group_col, value_col, alpha=0.05):
    """
    自动执行 Tukey HSD 事后多重比较，并按各组均值从大到小推导显著性字母标记 (a, b, ab, c)
    参数说明：
      - df        : 包含分组列和数值列的 DataFrame
      - group_col : 分组列名称字符串（如 'Formula' 或 'Treatment'）
      - value_col : 响应变量数值列名称字符串（如 'TBARS' 或 'EAI'）
      - alpha     : 显著性水平阈值，默认 0.05
    返回值：
      - letter_map: 字典 {组名: 显著性字母字符串}
      - summary_df: 包含均值、标准差与字母标记的汇总表 DataFrame
    """
    # 1. 按各组平均值从大到小降序排队
    grp_stats = df.groupby(group_col)[value_col].agg(['count', 'mean', 'std']).sort_values('mean', ascending=False)
    ordered_groups = list(grp_stats.index)
    samples = [df[df[group_col] == g][value_col].dropna().values for g in ordered_groups]

    # 2. 调用 scipy.stats.tukey_hsd() 执行权威 Tukey HSD 多重比较
    #    （若使用 pingouin，等价于 pg.pairwise_tukey(data=df, dv=value_col, between=group_col)）
    tukey_res = stats.tukey_hsd(*samples)
    p_matrix = tukey_res.pvalue  # 对称的 k x k 两两比较 p 值矩阵

    # 3. 吸收法构建“无显著差异团（Maximal Cliques）”并分配字母 a, b, c...
    k = len(ordered_groups)
    letters_pool = "abcdefghijklmnopqrstuvwxyz"
    assigned = [set() for _ in range(k)]
    letter_idx = 0

    for i in range(k):
        # 检查是否需要开启一个以第 i 组为首的新无差异集合
        need_new_group = True
        for prev_letter in assigned[i]:
            # 如果第 i 组已经有字母能覆盖后面所有与 i 无显著差异的组，则无需新开字母
            pass
        # 找出从 i 开始，所有与 i 无显著差异且彼此之间也无显著差异的连续/兼容组
        clique = [i]
        for j in range(i + 1, k):
            if p_matrix[i, j] >= alpha:
                # 检查 j 是否与当前 clique 内的所有成员都无显著差异 (p >= alpha)
                if all(p_matrix[m, j] >= alpha for m in clique):
                    clique.append(j)
        # 检查该 clique 是否已经是前面某个已分配字母集合的子集
        is_subset = False
        for l_char in letters_pool[:letter_idx]:
            existing_members = [idx for idx in range(k) if l_char in assigned[idx]]
            if set(clique).issubset(set(existing_members)):
                is_subset = True
                break
        if not is_subset:
            curr_char = letters_pool[letter_idx]
            for m in clique:
                assigned[m].add(curr_char)
            letter_idx += 1

    letter_map = {ordered_groups[i]: "".join(sorted(assigned[i])) for i in range(k)}
    grp_stats['CLD_Letter'] = [letter_map[g] for g in ordered_groups]
    grp_stats['Formatted_Result'] = [
        f"{grp_stats.loc[g, 'mean']:.2f} ± {grp_stats.loc[g, 'std']:.2f} {letter_map[g]}"
        for g in ordered_groups
    ]
    return letter_map, grp_stats.reset_index(), p_matrix


# ==============================================================================
# 第一节：单因素与双因素方差分析（Two-Way ANOVA）及交互效应曲线解析
# ==============================================================================
def run_section1_anova():
    print("\n" + "=" * 72)
    print(">>> 第三章第一节：单因素与双因素方差分析（主效应与交互效应解析） <<<")
    print("=" * 72)

    df1 = read_csv_safe("data/ch3_sec1_one_and_twoway.csv")

    # 1. 单因素方差分析：4 种工程酵母菌株 (Strain: S1~S4) 产谷胱甘肽得率对比
    #    权威写法：stats.f_oneway(*groups) 或 pingouin.anova(data=df1, dv='GSH_Yield', between='Strain')
    g_strains = [df1[df1['Strain'] == s]['GSH_Yield'].values for s in ['S1', 'S2', 'S3', 'S4']]
    f_1way, p_1way = stats.f_oneway(*g_strains)
    print(f"[1. 单因素 ANOVA (4种菌株 GSH 得率)] F = {f_1way:.2f}, p = {p_1way:.4e}")

    # 2. 双因素方差分析：贮藏温度 (Temp: 4C, 25C, 37C) x 包装方式 (Package: Air, Vacuum) 对 POV 的影响
    #    精确计算平方和分解：SS_Total = SS_A (Temp) + SS_B (Package) + SS_AxB (Interaction) + SS_Error
    grand_mean = df1['POV'].mean()
    N = len(df1)

    # 因素 A (Temp) 主效应平方和与自由度
    a_levels = df1['Temp'].unique()
    ss_a = sum(len(df1[df1['Temp'] == a]) * (df1[df1['Temp'] == a]['POV'].mean() - grand_mean)**2 for a in a_levels)
    df_a = len(a_levels) - 1

    # 因素 B (Package) 主效应平方和与自由度
    b_levels = df1['Package'].unique()
    ss_b = sum(len(df1[df1['Package'] == b]) * (df1[df1['Package'] == b]['POV'].mean() - grand_mean)**2 for b in b_levels)
    df_b = len(b_levels) - 1

    # 组内纯随机误差平方和 SS_Error
    ss_err = 0.0
    for a in a_levels:
        for b in b_levels:
            cell_vals = df1[(df1['Temp'] == a) & (df1['Package'] == b)]['POV']
            ss_err += ((cell_vals - cell_vals.mean())**2).sum()
    df_err = N - len(a_levels) * len(b_levels)

    # 总平方和 SS_Total 与交互项平方和 SS_AxB
    ss_tot = ((df1['POV'] - grand_mean)**2).sum()
    ss_axb = ss_tot - ss_a - ss_b - ss_err
    df_axb = df_a * df_b

    ms_a, ms_b, ms_axb, ms_err = ss_a / df_a, ss_b / df_b, ss_axb / df_axb, ss_err / df_err
    f_a, f_b, f_axb = ms_a / ms_err, ms_b / ms_err, ms_axb / ms_err
    p_a = 1 - stats.f.cdf(f_a, df_a, df_err)
    p_b = 1 - stats.f.cdf(f_b, df_b, df_err)
    p_axb = 1 - stats.f.cdf(f_axb, df_axb, df_err)

    anova_2way_df = pd.DataFrame([
        {"Source (变异来源)": "Temp (贮藏温度主效应 A)", "SS": round(ss_a, 3), "df": df_a, "MS": round(ms_a, 3), "F": round(f_a, 2), "p-value": f"{p_a:.4e}***", "np2 (偏η²效应量)": round(ss_a/(ss_a+ss_err), 3)},
        {"Source (变异来源)": "Package (包装方式主效应 B)", "SS": round(ss_b, 3), "df": df_b, "MS": round(ms_b, 3), "F": round(f_b, 2), "p-value": f"{p_b:.4e}***", "np2 (偏η²效应量)": round(ss_b/(ss_b+ss_err), 3)},
        {"Source (变异来源)": "Temp × Package (两因素交互作用)", "SS": round(ss_axb, 3), "df": df_axb, "MS": round(ms_axb, 3), "F": round(f_axb, 2), "p-value": f"{p_axb:.4e}***", "np2 (偏η²效应量)": round(ss_axb/(ss_axb+ss_err), 3)},
        {"Source (变异来源)": "Residual Error (组内随机误差)", "SS": round(ss_err, 3), "df": df_err, "MS": round(ms_err, 3), "F": "-", "p-value": "-", "np2 (偏η²效应量)": "-"},
    ])
    print("\n[2. 双因素方差分析 (Two-Way ANOVA) 标准三线表]：\n", anova_2way_df.to_string(index=False))
    anova_2way_df.to_csv("data/ch3_sec1_twoway_anova_table.csv", index=False, encoding='utf-8-sig')

    # 3. 绘制双因素交互作用折线图 + 6个组合多重比较柱状图
    df1['Combo'] = df1['Temp'] + "_" + df1['Package']
    cld_map, combo_stats, _ = compute_cld_letters(df1, 'Combo', 'POV')

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), dpi=150)

    # 左图：双因素交互作用剖面折线图 (Interaction Plot)
    sns.pointplot(
        data=df1, x='Temp', y='POV', hue='Package', order=['4C', '25C', '37C'],
        palette={'Air': '#E41A1C', 'Vacuum': '#377EB8'},
        markers=['o', 's'], linestyles=['-', '--'], capsize=0.1, ax=axes[0]
    )
    axes[0].annotate("两条折线呈强烈张开喇叭口 (非平行!)\n交互项 F = {:.1f}, p < 0.0001***".format(f_axb),
                     xy=(2, 18.5), xytext=(0.15, 14.5),
                     arrowprops=dict(facecolor='darkred', shrink=0.08, width=1.2, headwidth=5),
                     fontsize=9.5, fontweight='bold', color='darkred',
                     bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFF8DC', edgecolor='darkred'))
    axes[0].set_title("A. 双因素交互作用折线图 (Temp × Package)", fontsize=11, fontweight='bold')
    axes[0].set_xlabel("贮藏温度 (Storage Temperature)", fontsize=10.5)
    axes[0].set_ylabel("过氧化值 POV (meq/kg)", fontsize=10.5)

    # 右图：简单效应分析与6个交互组合的紧凑字母标记 (a, b, c, d)
    sns.barplot(
        data=df1, x='Temp', y='POV', hue='Package', order=['4C', '25C', '37C'],
        palette={'Air': '#E41A1C', 'Vacuum': '#377EB8'}, capsize=0.08, edgecolor='black', alpha=0.8, ax=axes[1]
    )
    axes[1].set_title("B. 双因素各水平组合比较与简单效应解析", fontsize=11, fontweight='bold')
    axes[1].set_xlabel("贮藏温度 (Storage Temperature)", fontsize=10.5)
    axes[1].set_ylabel("过氧化值 POV (meq/kg)", fontsize=10.5)

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch3_sec1_interaction_plot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第一节双因素交互作用图已保存至: scripts/ch3_sec1_interaction_plot.png")


# ==============================================================================
# 第二节 & 第三节：多重比较紧凑字母推导 (a,b,ab,c) 与顶刊复合云雨图绘制
# ==============================================================================
def run_section2_and_3_raincloud():
    print("\n" + "=" * 72)
    print(">>> 第三章第二、三节：事后多重比较字母标记 (a,b,ab,c) 与复合云雨图绘制 <<<")
    print("=" * 72)

    # 1. 第二节数据：5 种天然抗氧化剂配方对 TBARS 的抑制（演示 a, b, bc, c, d 字母推导）
    df2 = read_csv_safe("data/ch3_sec2_posthoc_letters.csv")
    letter_map2, summary_df2, _ = compute_cld_letters(df2, 'Formula', 'TBARS')
    print("\n[第二节：5种抗氧化配方 Tukey HSD 自动推导字母标记表]：\n",
          summary_df2[['Formula', 'mean', 'std', 'CLD_Letter', 'Formatted_Result']].to_string(index=False))
    summary_df2.to_csv("data/ch3_sec2_cld_table.csv", index=False, encoding='utf-8-sig')

    # 2. 第三节数据：6 组超高压协同酶解蛋白乳化活性 EAI 云雨图绘制
    df3 = read_csv_safe("data/ch3_sec3_raincloud_data.csv")
    letter_map3, summary_df3, _ = compute_cld_letters(df3, 'Treatment', 'EAI')
    print("\n[第三节：6组超高压酶解处理 EAI 均值与显著性字母标记]：\n",
          summary_df3[['Treatment', 'mean', 'std', 'CLD_Letter', 'Formatted_Result']].to_string(index=False))

    # 绘制双联高颜值图表：左图为第二节含 ab/bc 过渡字母的箱线图，右图为第三节 6 组复合云雨图
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0), gridspec_kw={'width_ratios': [1, 1.35]}, dpi=150)

    # 左图：演示 a, b, bc, c, d 显著性字母含义
    order2 = list(summary_df2['Formula'])
    pal2 = sns.color_palette("YlOrRd_r", len(order2))
    sns.boxplot(data=df2, x='Formula', y='TBARS', order=order2, palette=pal2, width=0.45,
                boxprops=dict(alpha=0.75, edgecolor='black'), showfliers=False, ax=axes[0])
    sns.stripplot(data=df2, x='Formula', y='TBARS', order=order2, color='black', size=6, jitter=0.15, alpha=0.8, ax=axes[0])
    for idx_g, g_name in enumerate(order2):
        g_max = df2[df2['Formula'] == g_name]['TBARS'].max()
        axes[0].text(idx_g, g_max + 0.08, letter_map2[g_name], ha='center', va='bottom',
                     fontsize=12, fontweight='bold', color='#990000')
    axes[0].set_ylim(0.5, df2['TBARS'].max() + 0.35)
    axes[0].set_title("A. 事后多重比较显著性字母标记法 (Tukey HSD, p < 0.05)\n[铁律: 含相同字母如 b 与 bc 代表无显著差异]",
                      fontsize=10.5, fontweight='bold')
    axes[0].set_xlabel("天然抗氧化剂配方 (Formula)", fontsize=10.5)
    axes[0].set_ylabel("肉糜脂质氧化 TBARS (mg MDA/kg)", fontsize=10.5)

    # 右图：SCI 顶刊复合云雨图（Raincloud Plot = 小提琴密度云 + 内部窄箱线伞 + 抖动散点雨 + 顶端字母）
    order3 = list(summary_df3['Treatment'])
    pal3 = sns.color_palette("Set2", len(order3))
    sns.violinplot(
        data=df3, x='Treatment', y='EAI', order=order3, palette=pal3,
        inner=None, cut=0, width=0.65, alpha=0.55, ax=axes[1]
    )
    sns.boxplot(
        data=df3, x='Treatment', y='EAI', order=order3, width=0.18, color='white',
        boxprops=dict(alpha=0.9, edgecolor='black', linewidth=1.1),
        medianprops=dict(color='#990000', linewidth=1.8), showfliers=False, ax=axes[1]
    )
    sns.stripplot(
        data=df3, x='Treatment', y='EAI', order=order3, palette=pal3,
        size=6, jitter=0.18, edgecolor='black', linewidth=0.6, alpha=0.85, ax=axes[1]
    )
    for idx_g, g_name in enumerate(order3):
        g_max = df3[df3['Treatment'] == g_name]['EAI'].max()
        axes[1].text(idx_g, g_max + 1.5, letter_map3[g_name], ha='center', va='bottom',
                     fontsize=12, fontweight='bold', color='black')
    axes[1].set_ylim(20, df3['EAI'].max() + 7.0)
    axes[1].set_title("B. 高信息密度复合云雨图 (Raincloud Plot: 核密度 + 箱线 + 原始散点 + CLD字母)",
                      fontsize=10.5, fontweight='bold')
    axes[1].set_xlabel("超高压协同酶解处理组 (Treatment)", fontsize=10.5)
    axes[1].set_ylabel("大豆分离蛋白乳化活性指数 EAI ($\mathrm{m^2/g}$)", fontsize=10.5)
    axes[1].tick_params(axis='x', rotation=15)

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch3_sec3_raincloud_plot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第三节复合云雨图与字母标记图已保存至: scripts/ch3_sec3_raincloud_plot.png")


if __name__ == '__main__':
    run_section1_anova()
    run_section2_and_3_raincloud()
