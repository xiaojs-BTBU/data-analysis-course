# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第一章（第一至三节）
脚本名称：ch1_01_data_wrangling_and_cleaning.py
代码目的：1. 区分技术重复与生物学重复，生成随机化检测顺序并可视化区组（Batch）效应；
          2. 演示 Tidy Data 整洁数据规范：宽表转长表（pd.melt）与长表转宽表（pd.pivot）及分面绘图；
          3. 数据清洗实战：组内中位数稳健插补缺失值（NaN），对比 3σ 原则与箱线图 IQR 离群值判别。
适用对象：食品科学与工程、生物工程及相关专业零编程基础硕士研究生
================================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 跨平台中文字体与负号防乱码配置
sns.set_theme(style="ticks", font_scale=1.08)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Noto Sans CJK JP', 'Droid Sans Fallback', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42


def read_csv_safe(filepath):
    """兼容 UTF-8 与 GBK 编码的安全读取函数"""
    try:
        return pd.read_csv(filepath, encoding='utf-8')
    except UnicodeDecodeError:
        return pd.read_csv(filepath, encoding='gbk')


# ==============================================================================
# 第一节实战：技术重复合并、随机化次序生成与区组效应可视化
# ==============================================================================
def run_section1_principles():
    print("\n" + "=" * 72)
    print(">>> 第一章第一节：试验设计三原则（技术重复合并 + 随机化 + 区组效应） <<<")
    print("=" * 72)

    df_sec1 = read_csv_safe("data/ch1_sec1_principles_demo.csv")
    print("[原始数据预览 (含3列酶标仪技术重复 Tech_Rep1~3)]：\n", df_sec1.head(4))

    # 1. 将同一发酵罐的 3 次技术重复求平均值，得到代表独立生物学重复 (n=18) 的真实观测值
    #    axis=1 参数表示“按行横向求平均值”（即把同一行的 Tech_Rep1, Tech_Rep2, Tech_Rep3 加起来除以 3）
    df_sec1['Lactic_Acid_Mean'] = df_sec1[['Tech_Rep1', 'Tech_Rep2', 'Tech_Rep3']].mean(axis=1).round(2)

    # 2. 为 18 个样品生成完全随机化的仪器测定先后顺序（防止仪器预热或时间漂移带来的系统偏倚）
    np.random.seed(42)  # 固定随机种子，保证每次运行生成的随机顺序完全可复现
    # np.random.permutation() 参数：生成一个从 1 到 n 的随机打乱整数序列
    df_sec1['Random_Test_Order'] = np.random.permutation(np.arange(1, len(df_sec1) + 1))

    # 导出规整后的生物学重复表格（使用 utf-8-sig 编码，保证中文 Excel 双击打开不乱码）
    df_sec1.to_csv("data/ch1_sec1_bio_replicates_out.csv", index=False, encoding='utf-8-sig')
    print("\n[规整后数据预览 (含生物学重复均值与随机上机顺序号)]：\n",
          df_sec1[['Fermenter_ID', 'Batch', 'Group', 'Lactic_Acid_Mean', 'Random_Test_Order']].head(6))

    # 3. 绘制展示“批次区组效应（Batch）”与“处理效应（Group）”的对比图
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), dpi=150)

    # 左图：不区分批次直接混合对比
    sns.boxplot(
        data=df_sec1, x='Group', y='Lactic_Acid_Mean', ax=axes[0],
        palette=['#4C72B0', '#DD8452'], width=0.45, boxprops=dict(alpha=0.75)
    )
    sns.stripplot(
        data=df_sec1, x='Group', y='Lactic_Acid_Mean', hue='Batch', ax=axes[0],
        palette='Set1', size=8, jitter=0.12, edgecolor='black', linewidth=0.8
    )
    axes[0].set_title("A. 按处理组对比 (不同颜色点代表不同原料批次)", fontsize=11, fontweight='bold')
    axes[0].set_xlabel("氮源处理分组 (Group)", fontsize=10.5)
    axes[0].set_ylabel("乳酸产量 Lactic Acid (g/L)", fontsize=10.5)

    # 右图：按区组 (Batch) 分组展示，清晰看见每个批次内部 New_Peptone 均稳定高于 Control
    sns.barplot(
        data=df_sec1, x='Batch', y='Lactic_Acid_Mean', hue='Group', ax=axes[1],
        palette=['#4C72B0', '#DD8452'], capsize=0.08, err_kws={'linewidth': 1.2}, edgecolor='black'
    )
    axes[1].set_title("B. 区组化 (Blocking) 剥离原料批次背景波动", fontsize=11, fontweight='bold')
    axes[1].set_xlabel("原料批次区组 (Batch Block)", fontsize=10.5)
    axes[1].set_ylabel("乳酸产量 Lactic Acid (g/L)", fontsize=10.5)

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch1_sec1_principles_plot.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第一节试验设计原则图表已保存至: scripts/ch1_sec1_principles_plot.png")


# ==============================================================================
# 第二节实战：Tidy Data 宽表转长表（melt）与多指标一键分面绘图
# ==============================================================================
def run_section2_tidy_data():
    print("\n" + "=" * 72)
    print(">>> 第一章第二节：Tidy Data 宽表转长表重塑与多指标分面绘图 <<<")
    print("=" * 72)

    df_wide = read_csv_safe("data/ch1_sec2_tidy_wide.csv")
    print("[原始实验室宽表 (Wide Format) 前 3 行]：\n", df_wide.head(3))

    # 使用 pd.melt() 将宽表重塑为标准 Tidy Data 长表
    # 参数详解：
    #   - frame      : 待转换的宽表 df_wide
    #   - id_vars    : 保持不变的样本标识列与分组列 ['Mouse_ID', 'Group']
    #   - value_vars : 需要被折叠合并的生化指标列 ['SOD', 'MDA', 'CAT', 'GSH_Px']
    #   - var_name   : 新生成的“指标名称列”的列名，设为 'Biomarker'
    #   - value_name : 新生成的“观测数值列”的列名，设为 'Value'
    df_long = pd.melt(
        frame=df_wide,
        id_vars=['Mouse_ID', 'Group'],
        value_vars=['SOD', 'MDA', 'CAT', 'GSH_Px'],
        var_name='Biomarker',
        value_name='Value'
    )
    df_long.to_csv("data/ch1_sec2_tidy_long_out.csv", index=False, encoding='utf-8-sig')
    print("\n[转换后的标准 Tidy 长表 (Long Format) 前 6 行]：\n", df_long.head(6))

    # 演示如何用 pd.pivot() 将长表重新还原回宽表
    df_pivoted = df_long.pivot(index=['Mouse_ID', 'Group'], columns='Biomarker', values='Value').reset_index()
    print("\n[验证长表逆向还原回宽表的前 2 行]：\n", df_pivoted.head(2))

    # 基于 Tidy 长表，调用 sns.catplot 一行代码绘制 4 个指标的分面箱线+散点图
    # 参数详解：
    #   - col='Biomarker': 按 Biomarker 列自动拆分为多个并排子图（分面 Facet）
    #   - sharey=False   : 极其重要！因为 SOD(数值约100) 和 MDA(数值约5) 量纲不同，必须设为 False 让各子图独立纵轴！
    g = sns.catplot(
        data=df_long, x='Group', y='Value', col='Biomarker',
        kind='box', sharey=False, height=3.6, aspect=0.78,
        palette=['#55A868', '#C44E52'], width=0.45, boxprops=dict(alpha=0.7)
    )
    # 在每个分面子图上叠加真实个体散点
    g.map_dataframe(
        sns.stripplot, x='Group', y='Value',
        color='black', alpha=0.7, jitter=0.15, size=6
    )
    g.set_titles("{col_name}", fontweight='bold', size=11)
    g.set_axis_labels("实验分组", "测定值 (独立量纲)")
    plt.subplots_adjust(top=0.82)
    g.fig.suptitle("基于 Tidy Data 长表一键生成的血清 4 项抗氧化指标分面箱线散点图", fontsize=12.5, fontweight='bold')
    g.savefig("scripts/ch1_sec2_tidy_facet.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第二节 Tidy Data 分面图已保存至: scripts/ch1_sec2_tidy_facet.png")


# ==============================================================================
# 第三节实战：缺失值组内中位数插补与 3σ vs IQR 离群值判别对比
# ==============================================================================
def run_section3_data_cleaning():
    print("\n" + "=" * 72)
    print(">>> 第一章第三节：缺失值组内中位数插补与 3σ / IQR 离群点判别 <<<")
    print("=" * 72)

    df_raw = read_csv_safe("data/ch1_sec3_cleaning_raw.csv")
    metrics = ['TNF_alpha', 'IL_6', 'NO_production']

    # 1. 缺失值统计与组内中位数稳健插补
    print("[1] 各指标缺失值 (NaN) 统计：\n", df_raw[metrics].isnull().sum())
    df_imputed = df_raw.copy()
    for col in metrics:
        # groupby('Group')[col].transform('median') 会自动算出每个样本所属组的组内中位数
        group_medians = df_imputed.groupby('Group')[col].transform('median')
        df_imputed[col] = df_imputed[col].fillna(group_medians)

    # 2. 对比 3σ 原则 与 1.5*IQR 箱线图法则对离群点的检出灵敏度
    outlier_records = []
    for col in metrics:
        for grp in df_imputed['Group'].unique():
            sub = df_imputed[df_imputed['Group'] == grp]
            vals = sub[col]

            # (a) 3σ 原则边界: Mean ± 3 * SD
            mu, sd = vals.mean(), vals.std(ddof=1)
            low_3sigma, high_3sigma = mu - 3 * sd, mu + 3 * sd

            # (b) 箱线图 IQR 原则边界: [Q1 - 1.5*IQR, Q3 + 1.5*IQR]
            q1, q3 = vals.quantile(0.25), vals.quantile(0.75)
            iqr = q3 - q1
            low_iqr, high_iqr = q1 - 1.5 * iqr, q3 + 1.5 * iqr

            for idx, row in sub.iterrows():
                v = row[col]
                flag_iqr = (v < low_iqr) or (v > high_iqr)
                flag_3sigma = (v < low_3sigma) or (v > high_3sigma)
                if flag_iqr or flag_3sigma:
                    outlier_records.append({
                        "Well_ID": row['Well_ID'], "Group": grp, "Metric": col, "Value": v,
                        "IQR_Bounds": f"[{low_iqr:.1f}, {high_iqr:.1f}]",
                        "Detected_by_IQR": "是 (Yes)" if flag_iqr else "否 (No)",
                        "3Sigma_Bounds": f"[{low_3sigma:.1f}, {high_3sigma:.1f}]",
                        "Detected_by_3Sigma": "是 (Yes)" if flag_3sigma else "否 (漏检!)"
                    })

    df_outliers = pd.DataFrame(outlier_records)
    print("\n[2] 3σ 原则 vs 箱线图 IQR 准则离群点检出对比表：\n", df_outliers.to_string(index=False))
    df_outliers.to_csv("data/ch1_sec3_outlier_report.csv", index=False, encoding='utf-8-sig')
    df_imputed.to_csv("data/ch1_sec3_cleaned_out.csv", index=False, encoding='utf-8-sig')

    # 3. 绘制离群点清洗诊断可视化图
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), dpi=150)

    # 左图：含极端离群点的 TNF_alpha 分布与 IQR 围栏标示
    sns.boxplot(data=df_imputed, x='Group', y='TNF_alpha', ax=axes[0], width=0.4, palette=['#E41A1C', '#377EB8'], boxprops=dict(alpha=0.6))
    sns.stripplot(data=df_imputed, x='Group', y='TNF_alpha', ax=axes[0], color='black', size=7, jitter=0.12)
    axes[0].annotate("移液故障极端离群点 (W06: 495.0)\n被 IQR 法精准捕获!", xy=(0, 495), xytext=(0.25, 430),
                     arrowprops=dict(facecolor='red', shrink=0.08, width=1.5, headwidth=6),
                     fontsize=9.5, fontweight='bold', color='darkred',
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3CD', edgecolor='red'))
    axes[0].set_title("A. 离群点诊断前：极端异常值严重拉偏坐标轴", fontsize=11, fontweight='bold')

    # 右图：核查实验记录确认为加样气泡后，修正/剔除故障点后的真实分布
    df_fixed = df_imputed[df_imputed['TNF_alpha'] < 400]
    sns.boxplot(data=df_fixed, x='Group', y='TNF_alpha', ax=axes[1], width=0.4, palette=['#E41A1C', '#377EB8'], boxprops=dict(alpha=0.7))
    sns.stripplot(data=df_fixed, x='Group', y='TNF_alpha', ax=axes[1], color='black', size=7, jitter=0.12)
    axes[1].set_title("B. 排查仪器故障剔除异常点后：组间差异清晰呈现", fontsize=11, fontweight='bold')

    for ax in axes:
        sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()
    plt.savefig("scripts/ch1_sec3_outlier_check.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[成功] 第三节数据清洗与离群点诊断图已保存至: scripts/ch1_sec3_outlier_check.png")


if __name__ == '__main__':
    run_section1_principles()
    run_section2_tidy_data()
    run_section3_data_cleaning()
