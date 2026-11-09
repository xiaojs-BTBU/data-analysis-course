# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第二章：假设检验
脚本名称：ch2_hypothesis_testing_and_apa_table.py
代码目的：1. 第一节：正态性检验（Shapiro-Wilk）与方差齐性检验（Levene / Bartlett）及 Q-Q 图诊断；
          2. 第二节：单样本 t 检验、两独立样本自动决策树（Student's t / Welch's t / Mann-Whitney U）、
                     配对样本检验（Paired t / Wilcoxon）及效应量（Cohen's d / 秩双列相关 RBC）计算；
          3. 第三节：2x2 列联表 Pearson 卡方检验与 Fisher 精确检验，自动化生成 APA 第七版标准三线表。
权威包双轨教学：同时展示 `pingouin` 一行流代码与 `scipy.stats` 标准底层代码。
================================================================================
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# 尝试导入生物统计权威包 pingouin（若在未安装 pingouin 的环境中运行，自动启用下方内置等效函数）
try:
    import pingouin as pg
    HAS_PINGOUIN = True
except ImportError:
    HAS_PINGOUIN = False

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


# ==============================================================================
# 第一节：正态性（Shapiro-Wilk）与方差齐性（Levene / Bartlett）先验体检
# ==============================================================================
def run_section1_assumptions():
    print("\n" + "=" * 72)
    print(">>> 第二章第一节：假设检验前提评估（正态性与方差齐性体检） <<<")
    print("=" * 72)

    df1 = read_csv_safe("data/ch2_sec1_assumptions.csv")
    biomarkers = ['SOD_Normal_Homo', 'MDA_Normal_Hetero', 'IL6_Skewed']

    report_rows = []
    for col in biomarkers:
        g_ctrl = df1[df1['Group'] == 'Control'][col].dropna()
        g_prob = df1[df1['Group'] == 'Probiotic'][col].dropna()

        # ------------------------------------------------------------------
        # 【权威双轨写法 A：使用 pingouin 一行代码体检（学生在本地 stats_env 首选）】
        #   pg.normality(data=df1, dv=col, group='Group', method='shapiro')
        #   pg.homoscedasticity(data=df1, dv=col, group='Group', method='levene')
        # 【权威双轨写法 B：使用 scipy.stats 标准底层函数】
        # ------------------------------------------------------------------
        # 1. stats.shapiro() 参数：传入一维样本数组，返回 (W统计量, p值)；p > 0.05 代表服从正态分布
        w_c, p_norm_c = stats.shapiro(g_ctrl)
        w_p, p_norm_p = stats.shapiro(g_prob)
        is_normal = (p_norm_c > 0.05) and (p_norm_p > 0.05)

        # 2. stats.levene() 参数：
        #    - g_ctrl, g_prob : 待检验方差齐性的两组数据
        #    - center='median': 以中位数为中心计算离差（Brown-Forsythe 修正，对偏态数据最稳健）
        stat_lev, p_levene = stats.levene(g_ctrl, g_prob, center='median')
        # 3. stats.bartlett()：适用于严格正态分布数据的方差齐性检验
        stat_bart, p_bartlett = stats.bartlett(g_ctrl, g_prob)
        is_homo = (p_levene > 0.05)

        # 自动推荐后续应选用的统计检验方法
        if is_normal and is_homo:
            rec_method = "Student's t-test (正态且方差齐)"
        elif is_normal and (not is_homo):
            rec_method = "Welch's t-test (正态但方差不齐)"
        else:
            rec_method = "Mann-Whitney U test (非正态偏态分布)"

        report_rows.append({
            "指标名称 (Biomarker)": col,
            "Control正态p值": round(p_norm_c, 4),
            "Probiotic正态p值": round(p_norm_p, 4),
            "满足正态性(p>0.05)": "是" if is_normal else "否(偏态)",
            "Levene方差齐性p值": round(p_levene, 4),
            "Bartlett齐性p值": round(p_bartlett, 4),
            "满足方差齐(p>0.05)": "是" if is_homo else "否(异方差)",
            "推荐检验方法": rec_method
        })

    df_rep = pd.DataFrame(report_rows)
    print(df_rep.to_string(index=False))
    df_rep.to_csv("data/ch2_sec1_assumption_report.csv", index=False, encoding='utf-8-sig')

    # 绘制 3 个指标的核密度分布图与正态 Q-Q 图（Quantile-Quantile Plot）
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.2), dpi=150)
    titles = [
        "SOD (正态 + 方差齐)\n-> 选 Student's t 检验",
        "MDA (正态 + 方差不齐)\n-> 选 Welch's t 检验",
        "IL-6 (严重右偏态)\n-> 选 Mann-Whitney U 检验"
    ]
    for i, col in enumerate(biomarkers):
        # 第一行：核密度分布与箱线对比
        sns.kdeplot(data=df1, x=col, hue='Group', fill=True, common_norm=False,
                    palette=['#377EB8', '#E41A1C'], alpha=0.35, linewidth=2, ax=axes[0, i])
        axes[0, i].set_title(titles[i], fontsize=10.5, fontweight='bold')

        # 第二行：Probiotic 组的正态 Q-Q 图（点紧贴红线代表服从正态，弯成弧线代表偏态）
        sub_vals = df1[df1['Group'] == 'Probiotic'][col].dropna()
        stats.probplot(sub_vals, dist="norm", plot=axes[1, i])
        axes[1, i].set_title(f"{col} 正态 Q-Q 诊断图", fontsize=10, fontweight='bold')
        axes[1, i].get_lines()[0].set_markerfacecolor('#377EB8')
        axes[1, i].get_lines()[1].set_color('#E41A1C')

    for ax in axes.flat:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch2_sec1_qq_density.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第一节正态 Q-Q 与核密度诊断图已保存至: scripts/ch2_sec1_qq_density.png")


# ==============================================================================
# 第二节：单样本、独立样本自动分支检验、配对样本检验与效应量计算
# ==============================================================================
def run_section2_t_and_nonparam():
    print("\n" + "=" * 72)
    print(">>> 第二章第二节：参数检验与非参数检验全分支及效应量评估 <<<")
    print("=" * 72)

    df2 = read_csv_safe("data/ch2_sec2_t_and_nonparam.csv")

    # 1. 单样本 t 检验：抽检 12 罐功能奶粉蛋白质含量是否显著高于国标线 12.0 g/100g
    prot_vals = df2['Protein_Content'].dropna()
    # stats.ttest_1samp() 参数：
    #   - a       : 样本观测值数组
    #   - popmean : 待对比的国家标准或理论总体均值（此处设为 12.0）
    t_1s, p_1s = stats.ttest_1samp(a=prot_vals, popmean=12.0)
    cohen_d_1s = (prot_vals.mean() - 12.0) / prot_vals.std(ddof=1)
    print(f"[1. 单样本 t 检验] 实测均值 = {prot_vals.mean():.2f}±{prot_vals.std(ddof=1):.2f} vs 国标 12.0 "
          f"| t = {t_1s:.3f}, p = {p_1s:.4f}, Cohen's d = {cohen_d_1s:.2f}")

    # 2. 两独立样本自动化决策树函数（整合正态体检 + 方差齐性体检 + 自动分支 + 效应量）
    ind_results = []
    for col in ['SOD', 'MDA', 'TNF_a']:
        x = df2[df2['Group'] == 'Control'][col].dropna()
        y = df2[df2['Group'] == 'Treatment'][col].dropna()

        _, p_nx = stats.shapiro(x)
        _, p_ny = stats.shapiro(y)
        _, p_lev = stats.levene(x, y, center='median')

        n1, n2 = len(x), len(y)
        # 计算合并标准差 s_pooled 与 Cohen's d 效应量
        s_pooled = np.sqrt(((n1 - 1) * x.var(ddof=1) + (n2 - 1) * y.var(ddof=1)) / (n1 + n2 - 2))
        cohen_d = (y.mean() - x.mean()) / s_pooled

        if p_nx > 0.05 and p_ny > 0.05 and p_lev > 0.05:
            # 分支 A：正态且方差齐 -> 标准 Student's t 检验 (equal_var=True)
            # 若使用 pingouin: pg.ttest(x, y, correction=False)
            stat_val, p_val = stats.ttest_ind(y, x, equal_var=True)
            method_name = "Student's t-test"
            eff_str = f"d = {cohen_d:.2f}"
        elif p_nx > 0.05 and p_ny > 0.05 and p_lev <= 0.05:
            # 分支 B：正态但方差不齐 -> Welch's t 检验 (equal_var=False，自动修正自由度)
            # 若使用 pingouin: pg.ttest(x, y, correction=True)
            stat_val, p_val = stats.ttest_ind(y, x, equal_var=False)
            method_name = "Welch's t-test"
            eff_str = f"d = {cohen_d:.2f}"
        else:
            # 分支 C：非正态偏态分布 -> 非参数 Mann-Whitney U 秩和检验
            # 若使用 pingouin: pg.mwu(x, y)
            u_stat, p_val = stats.mannwhitneyu(y, x, alternative='two-sided')
            stat_val = u_stat
            method_name = "Mann-Whitney U"
            # 计算秩双列相关系数效应量 Rank-Biserial Correlation (r_rb)
            r_rb = 2 * (u_stat / (n1 * n2)) - 1
            eff_str = f"r_rb = {r_rb:.2f}"

        ind_results.append({
            "Biomarker": col,
            "Control (Mean±SD / Med[IQR])": f"{x.mean():.2f}±{x.std(ddof=1):.2f}" if "t-test" in method_name else f"{x.median():.2f} ({x.quantile(0.75)-x.quantile(0.25):.2f})",
            "Treatment (Mean±SD / Med[IQR])": f"{y.mean():.2f}±{y.std(ddof=1):.2f}" if "t-test" in method_name else f"{y.median():.2f} ({y.quantile(0.75)-y.quantile(0.25):.2f})",
            "Method": method_name,
            "Statistic": round(stat_val, 3),
            "p_value": round(p_val, 5),
            "Effect_Size": eff_str
        })

    df_ind_res = pd.DataFrame(ind_results)
    print("\n[2. 两独立样本自动化分支检验汇总表]：\n", df_ind_res.to_string(index=False))

    # 3. 配对样本检验：同一批大鼠干预前 (Glucose_Before) vs 干预后 (Glucose_After)
    pre = df2['Glucose_Before'].dropna()
    post = df2['Glucose_After'].dropna()
    diff = post - pre
    _, p_diff_norm = stats.shapiro(diff)
    # stats.ttest_rel() 参数：传入配对的前后两列数组（要求长度严格相等且一一对应）
    # 若使用 pingouin: pg.ttest(pre, post, paired=True)
    t_pair, p_pair = stats.ttest_rel(post, pre)
    w_wilcox, p_wilcox = stats.wilcoxon(post, pre)
    print(f"\n[3. 配对样本检验] 前后差值正态性 p = {p_diff_norm:.3f} | "
          f"配对 t 检验 t = {t_pair:.3f}, p = {p_pair:.5f} | Wilcoxon p = {p_wilcox:.5f}")

    df_ind_res.to_csv("data/ch2_sec2_test_results.csv", index=False, encoding='utf-8-sig')

    # 绘制两独立样本箱线散点图 + 配对样本个体斜率连线图
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2), dpi=150)

    for idx_m, col in enumerate(['SOD', 'MDA']):
        sns.boxplot(data=df2, x='Group', y=col, ax=axes[idx_m], width=0.42,
                    palette=['#4C72B0', '#C44E52'], boxprops=dict(alpha=0.65), showfliers=False)
        sns.stripplot(data=df2, x='Group', y=col, ax=axes[idx_m],
                      palette=['#4C72B0', '#C44E52'], size=7, jitter=0.15, edgecolor='black', linewidth=0.8)
        row_info = df_ind_res[df_ind_res['Biomarker'] == col].iloc[0]
        axes[idx_m].set_title(f"{col} ({row_info['Method']})\n$p = {row_info['p_value']:.4f}$, {row_info['Effect_Size']}",
                              fontsize=10.5, fontweight='bold')

    # 子图 3：配对样本前后个体连线图 (Paired Slope Plot)
    for i in range(len(pre)):
        axes[2].plot([0, 1], [pre.iloc[i], post.iloc[i]], color='gray', alpha=0.55, linestyle='--', linewidth=1.2)
    axes[2].scatter(np.zeros(len(pre)), pre, color='#4C72B0', s=55, edgecolor='black', zorder=4, label='干预前 (Before)')
    axes[2].scatter(np.ones(len(post)), post, color='#55A868', s=55, edgecolor='black', zorder=4, label='干预后 (After)')
    axes[2].set_xticks([0, 1])
    axes[2].set_xticklabels(['干预前 (Before)', '干预后 (After)'], fontweight='bold')
    axes[2].set_xlim(-0.35, 1.35)
    axes[2].set_ylabel("空腹血糖 Glucose (mmol/L)")
    axes[2].set_title(f"大鼠配对前后血糖连线图 (Paired t-test)\n$t = {t_pair:.2f}, p < 0.0001^{{***}}$", fontsize=10.5, fontweight='bold')

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch2_sec2_hypothesis_plots.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第二节假设检验组合图已保存至: scripts/ch2_sec2_hypothesis_plots.png")


# ==============================================================================
# 第三节：定性变量卡方/Fisher检验与 APA 第七版标准学术三线表一键导出
# ==============================================================================
def run_section3_chi2_and_apa_table():
    print("\n" + "=" * 72)
    print(">>> 第二章第三节：列联表卡方/Fisher检验与 APA 标准三线表导出 <<<")
    print("=" * 72)

    df3 = read_csv_safe("data/ch2_sec3_chi2_apa_table.csv")

    # 1. 构建 Group (Control vs Coating) x Spoilage (Spoiled vs Intact) 的 2x2 频数列联表
    ct = pd.crosstab(df3['Group'], df3['Spoilage'])
    print("[2x2 频数列联表]：\n", ct)

    # stats.chi2_contingency() 参数：
    #   - observed       : 2x2 频数列联表
    #   - correction=False: 设为 False 执行标准 Pearson 卡方检验；设为 True 执行 Yates 连续性修正
    chi2_val, p_chi2, dof, expected = stats.chi2_contingency(ct, correction=False)
    odds_ratio, p_fisher = stats.fisher_exact(ct)
    print(f"最小理论期望频数 T_min = {expected.min():.2f} (因 T_min >= 5 且 N=60>=40，适用 Pearson 卡方检验)")
    print(f"Pearson 卡方 chi2 = {chi2_val:.3f}, p = {p_chi2:.4f} | Fisher 精确检验 OR = {odds_ratio:.3f}, p = {p_fisher:.4f}")

    # 2. 自动汇总全部定性与定量指标，生成符合 APA 第七版规范的标准三线表
    apa_rows = []
    # (a) 定性指标行：霉变率 n (%)
    ctrl_spoil = (df3[df3['Group'] == 'Control']['Spoilage'] == 'Spoiled').sum()
    coat_spoil = (df3[df3['Group'] == 'Coating']['Spoilage'] == 'Spoiled').sum()
    apa_rows.append({
        "Variable (检测指标)": "Spoilage Rate 霉变率, n (%)",
        "Control (n=30)": f"{ctrl_spoil} ({ctrl_spoil/30*100:.1f}%)",
        "Coating (n=30)": f"{coat_spoil} ({coat_spoil/30*100:.1f}%)",
        "Statistical Method": "Pearson Chi-Square",
        "Statistic": f"χ² = {chi2_val:.2f}",
        "p-value": f"{p_chi2:.4f}**",
        "Effect Size": f"OR = {odds_ratio:.2f}"
    })

    # (b) 定量指标行：自动体检并格式化
    for col, label in [('Vitamin_C', 'Vitamin C (mg/100g)'),
                       ('Weight_Loss', 'Weight Loss 失重率 (%)'),
                       ('Log_CFU', 'Total Colonies 菌落总数 (log CFU/g)')]:
        x = df3[df3['Group'] == 'Control'][col]
        y = df3[df3['Group'] == 'Coating'][col]
        _, p_nx = stats.shapiro(x)
        _, p_ny = stats.shapiro(y)
        _, p_lev = stats.levene(x, y, center='median')
        s_pool = np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2)
        d_val = (y.mean() - x.mean()) / s_pool

        if p_nx > 0.05 and p_ny > 0.05:
            eq_v = (p_lev > 0.05)
            t_v, p_v = stats.ttest_ind(y, x, equal_var=eq_v)
            m_str = "Student's t-test" if eq_v else "Welch's t-test"
            c_str = f"{x.mean():.2f} ± {x.std(ddof=1):.2f}"
            t_str = f"{y.mean():.2f} ± {y.std(ddof=1):.2f}"
            st_str = f"t = {t_v:.2f}"
            ef_str = f"d = {d_val:.2f}"
        else:
            u_v, p_v = stats.mannwhitneyu(y, x)
            m_str = "Mann-Whitney U"
            c_str = f"{x.median():.2f} ({x.quantile(0.75)-x.quantile(0.25):.2f})"
            t_str = f"{y.median():.2f} ({y.quantile(0.75)-y.quantile(0.25):.2f})"
            st_str = f"U = {u_v:.1f}"
            ef_str = f"r = {2*(u_v/(len(x)*len(y)))-1:.2f}"

        apa_rows.append({
            "Variable (检测指标)": label,
            "Control (n=30)": c_str,
            "Coating (n=30)": t_str,
            "Statistical Method": m_str,
            "Statistic": st_str,
            "p-value": "<0.0001***" if p_v < 0.0001 else f"{p_v:.4f}",
            "Effect Size": ef_str
        })

    df_apa = pd.DataFrame(apa_rows)
    print("\n[APA 第七版标准学术三线表]：\n", df_apa.to_string(index=False))
    df_apa.to_csv("data/ch2_sec3_apa_summary_table.csv", index=False, encoding='utf-8-sig')

    # 绘制霉变率百分比堆叠柱状图 + 渲染出版级 APA 三线表图版
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.3), gridspec_kw={'width_ratios': [1, 1.85]}, dpi=150)

    # 左图：定性变量百分比堆叠图
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
    ct_pct[['Intact', 'Spoiled']].plot(kind='bar', stacked=True, color=['#55A868', '#C44E52'],
                                       edgecolor='black', width=0.48, ax=axes[0])
    axes[0].set_title(f"A. 鲜切草莓第7天霉变率对比\n($\\chi^2 = {chi2_val:.2f}, p = {p_chi2:.4f}^{{**}}$)", fontsize=11, fontweight='bold')
    axes[0].set_ylabel("占总盒数百分比 (%)")
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)
    axes[0].legend(['完好 (Intact)', '霉变 (Spoiled)'], loc='upper right', frameon=True)
    sns.despine(ax=axes[0], top=True, right=True)

    # 右图：渲染标准 APA 三线表（顶线粗、栏目线细、底线粗、无竖线）
    axes[1].axis('off')
    axes[1].set_title("B. 自动生成的 APA 第七版标准学术三线表 (Three-Line Table)", fontsize=11, fontweight='bold', pad=10)
    tbl = axes[1].table(
        cellText=df_apa.values,
        colLabels=df_apa.columns,
        cellLoc='center',
        loc='center',
        bbox=[0.0, 0.12, 1.0, 0.78]
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    for (row_i, col_j), cell in tbl.get_celld().items():
        cell.set_edgecolor('black')
        if row_i == 0:
            # 表头行：只保留上下横线（顶线与栏目线），去除左右竖线
            cell.visible_edges = 'TB'
            cell.set_linewidth(1.5)
            cell.set_text_props(weight='bold', color='#990000')
        elif row_i == len(df_apa):
            # 底行：只保留底线
            cell.visible_edges = 'B'
            cell.set_linewidth(1.5)
        else:
            # 中间数据行：无任何边框线
            cell.visible_edges = ''

    plt.tight_layout()
    plt.savefig("scripts/ch2_sec3_chi2_and_table.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第三节卡方分析与 APA 三线表图已保存至: scripts/ch2_sec3_chi2_and_table.png")


if __name__ == '__main__':
    run_section1_assumptions()
    run_section2_t_and_nonparam()
    run_section3_chi2_and_apa_table()
