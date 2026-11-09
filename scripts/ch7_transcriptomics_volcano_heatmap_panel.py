# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第七章：转录组学数据分析初步
脚本名称：ch7_transcriptomics_volcano_heatmap_panel.py
代码目的：1. 第一节：Benjamini-Hochberg (FDR) 多重假设检验校正，按 |log2FC|>=1.0 且 FDR<0.05 筛选差异基因，
                     绘制带防重叠智能避让标签的差异火山图（Volcano Plot）；
          2. 第二节：按行 Z-score 标准化绘制带样本分组色条的双向层次聚类热图（Complex Heatmap）
                     以及 GO/KEGG 代谢通路四维信息富集气泡图（Bubble Plot）；
          3. 第三节：使用 matplotlib.gridspec 代码化像素级精准拼排 SCI 论文整版 Figure 1 (Panel A + B + C)。
================================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

# 尝试导入标签自动避让库 adjustText
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

# 跨平台中文字体与负号防乱码配置
sns.set_theme(style="ticks", font_scale=1.05)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Noto Sans CJK JP', 'Droid Sans Fallback', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42


def read_csv_safe(filepath):
    try:
        return pd.read_csv(filepath, encoding='utf-8')
    except UnicodeDecodeError:
        return pd.read_csv(filepath, encoding='gbk')


def fdr_bh_correction(p_vals):
    """
    Benjamini-Hochberg (FDR) 多重假设检验校正权威实现
    （与 pingouin.multicomp(p_vals, method='fdr_bh') 完全一致）
    """
    p_arr = np.asarray(p_vals, dtype=float)
    m = len(p_arr)
    order = np.argsort(p_arr)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, m + 1)
    fdr = p_arr * m / ranks
    # 保证单调性
    fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
    fdr_out = np.empty_like(fdr_sorted)
    fdr_out[order] = np.clip(fdr_sorted, 0, 1.0)
    return fdr_out


def draw_volcano_on_ax(ax, df_deg):
    """在指定坐标轴 ax 上绘制带防重叠靶基因标签的差异火山图"""
    df_deg['neg_log10_fdr'] = -np.log10(df_deg['FDR'] + 1e-300)
    conds = [
        (df_deg['log2FC'] >= 1.0) & (df_deg['FDR'] < 0.05),
        (df_deg['log2FC'] <= -1.0) & (df_deg['FDR'] < 0.05)
    ]
    df_deg['Regulation'] = np.select(conds, ['Up-regulated', 'Down-regulated'], default='Not Significant')

    colors = {'Up-regulated': '#E41A1C', 'Down-regulated': '#377EB8', 'Not Significant': '#BDBDBD'}
    for reg, col in colors.items():
        sub = df_deg[df_deg['Regulation'] == reg]
        ax.scatter(sub['log2FC'], sub['neg_log10_fdr'], c=col,
                   s=24 if reg != 'Not Significant' else 10,
                   alpha=0.82 if reg != 'Not Significant' else 0.35,
                   edgecolor='none', label=f"{reg} ({len(sub)})")

    ax.axvline(1.0, color='black', linestyle='--', linewidth=0.9, alpha=0.6)
    ax.axvline(-1.0, color='black', linestyle='--', linewidth=0.9, alpha=0.6)
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.9, alpha=0.6)

    # 优先标注具有明确生物学意义的明星功能靶基因（抗氧化、脂质代谢与炎症因子）
    named_genes = ['Sod2', 'Gpx1', 'Nfe2l2', 'Hmox1', 'Cpt1a', 'Il6', 'Tnf', 'Fasn', 'Srebf1', 'Nlrp3']
    top_stars = df_deg[df_deg['Gene_Symbol'].isin(named_genes)].copy()
    if len(top_stars) == 0:
        top_up = df_deg[df_deg['Regulation'] == 'Up-regulated'].sort_values('neg_log10_fdr', ascending=False).head(4)
        top_dn = df_deg[df_deg['Regulation'] == 'Down-regulated'].sort_values('neg_log10_fdr', ascending=False).head(4)
        top_stars = pd.concat([top_up, top_dn]).reset_index(drop=True)

    # 按上下调分组并按纵坐标排序，阶梯式排布标签且避开右上角图例区
    up_rows = top_stars[top_stars['log2FC'] > 0].sort_values('neg_log10_fdr', ascending=False).reset_index(drop=True)
    dn_rows = top_stars[top_stars['log2FC'] < 0].sort_values('neg_log10_fdr', ascending=False).reset_index(drop=True)
    up_y_targets = np.linspace(6.6, 2.5, max(len(up_rows), 1))
    dn_y_targets = np.linspace(7.8, 2.8, max(len(dn_rows), 1))
    for idx, row in up_rows.iterrows():
        ax.annotate(
            row['Gene_Symbol'],
            xy=(row['log2FC'], row['neg_log10_fdr']),
            xytext=(min(row['log2FC'] + 0.65, 3.35), up_y_targets[idx]),
            fontsize=8.5, fontweight='bold', color='#800000', ha='left',
            arrowprops=dict(arrowstyle='->', color='#555555', lw=0.75, shrinkA=0, shrinkB=2)
        )
    for idx, row in dn_rows.iterrows():
        ax.annotate(
            row['Gene_Symbol'],
            xy=(row['log2FC'], row['neg_log10_fdr']),
            xytext=(max(row['log2FC'] - 0.65, -3.55), dn_y_targets[idx]),
            fontsize=8.5, fontweight='bold', color='#08306B', ha='right',
            arrowprops=dict(arrowstyle='->', color='#555555', lw=0.75, shrinkA=0, shrinkB=2)
        )

    ax.set_xlim(-4.3, 4.3)
    ax.set_ylim(-0.3, 9.8)
    ax.set_xlabel("$\log_2$ Fold Change (Treat / Ctrl)", fontsize=10.5, fontweight='bold')
    ax.set_ylabel("$-\log_{10}$ (FDR)", fontsize=10.5, fontweight='bold')
    ax.legend(loc='upper right', frameon=True, fontsize=8.5)
    sns.despine(ax=ax, top=True, right=True)


def run_chapter7_pipeline():
    print("\n" + "=" * 72)
    print(">>> 第七章：转录组差异火山图、行标准化聚类热图、KEGG气泡图与整版拼排 <<<")
    print("=" * 72)

    df_deg = read_csv_safe("data/ch7_sec1_deg_volcano.csv")
    df_hm = read_csv_safe("data/ch7_sec2_heatmap_matrix.csv")
    df_kegg = read_csv_safe("data/ch7_sec2_kegg_enrichment.csv").sort_values('Rich_Factor', ascending=True)

    # 1. 独立火山图导出
    fig, ax = plt.subplots(figsize=(6.8, 5.4), dpi=150)
    draw_volcano_on_ax(ax, df_deg)
    ax.set_title("RNA-seq 差异表达基因火山图 (|log2FC|≥1.0, FDR<0.05)", fontsize=11.5, fontweight='bold')
    plt.tight_layout()
    plt.savefig("scripts/ch7_sec1_volcano.png", dpi=600, bbox_inches='tight')
    plt.close()

    sig_degs = df_deg[df_deg['Regulation'] != 'Not Significant']
    sig_degs.to_csv("data/ch7_sec1_significant_degs.csv", index=False, encoding='utf-8-sig')
    print(f"[第一节] 共筛出显著差异表达基因 {len(sig_degs)} 个，火山图已保存至: scripts/ch7_sec1_volcano.png")

    # 2. 第二节：按行 Z-score 标准化的双向层次聚类热图 (Complex Heatmap)
    sample_cols = ['Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Treat_1', 'Treat_2', 'Treat_3']
    expr_mat = df_hm.set_index('Gene_Symbol')[sample_cols]
    # z_score=0 参数：表示按行（每个基因内部）执行 Z-score 标准化，消除管家基因绝对丰度霸屏！
    col_colors = ['#377EB8'] * 3 + ['#E41A1C'] * 3
    g = sns.clustermap(
        expr_mat, z_score=0, cmap='RdBu_r', center=0, vmin=-2.2, vmax=2.2,
        col_colors=col_colors, figsize=(6.8, 7.2),
        dendrogram_ratio=(0.15, 0.14), cbar_pos=(0.02, 0.82, 0.03, 0.14),
        cbar_kws={'label': 'Row Z-score'}, linewidths=0.4
    )
    g.ax_heatmap.set_title("核心差异基因双向层次聚类热图 (Row Z-score)", fontsize=11, fontweight='bold', pad=28)
    g.savefig("scripts/ch7_sec2_heatmap.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[第二节] 双向层次聚类热图已保存至: scripts/ch7_sec2_heatmap.png")

    # 3. 第三节：代码化像素级精准拼排 SCI 论文整版 Figure 1 (Panel A: 火山图 + Panel B: 聚类热图 + Panel C: KEGG气泡图)
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(16.8, 5.4), dpi=150,
        gridspec_kw={'width_ratios': [1.15, 0.90, 1.10]}
    )

    # Panel A: 火山图
    draw_volcano_on_ax(ax_a, df_deg)
    ax_a.set_title("A. Volcano Plot of DEGs", fontsize=11.5, fontweight='bold', pad=10)

    # Panel B: 核心靶基因行 Z-score 聚类热图 (色条置于顶部或精简标签，彻底避免与右侧通路名碰撞)
    z_mat = expr_mat.apply(lambda r: (r - r.mean()) / (r.std(ddof=1) + 1e-8), axis=1)
    sns.heatmap(z_mat.iloc[:16], cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                cbar_kws={'shrink': 0.55, 'pad': 0.03, 'label': 'Z-score'}, linewidths=0.4, ax=ax_b)
    ax_b.set_title("B. Expression Heatmap (Top 16)", fontsize=11.5, fontweight='bold', pad=10)
    ax_b.set_xlabel("Biological Replicates", fontsize=10, fontweight='bold')
    ax_b.set_ylabel("Gene Symbol", fontsize=10, fontweight='bold')

    # Panel C: GO / KEGG 四维信息富集气泡图 (Bubble Plot)
    sc = ax_c.scatter(
        df_kegg['Rich_Factor'], df_kegg['Pathway_Name'],
        s=df_kegg['Gene_Count'] * 14,
        c=-np.log10(df_kegg['FDR']),
        cmap='YlOrRd', edgecolor='black', linewidth=0.8, alpha=0.9
    )
    plt.colorbar(sc, ax=ax_c, shrink=0.55, pad=0.03, label='$-\log_{10}$ (FDR)')
    ax_c.set_title("C. GO & KEGG Pathway Enrichment", fontsize=11.5, fontweight='bold', pad=10)
    ax_c.set_xlabel("Rich Factor (Mapped / Background)", fontsize=10, fontweight='bold')
    ax_c.set_xlim(0.24, 0.54)
    ax_c.grid(True, linestyle=':', alpha=0.4)
    sns.despine(ax=ax_c, top=True, right=True)
    plt.tight_layout(pad=1.6, w_pad=2.4)

    plt.savefig("scripts/ch7_sec3_figure1_panel.png", dpi=600, bbox_inches='tight')
    plt.close()
    print("[第三节] SCI 整版三联复合图版 (Figure 1 Panel A+B+C) 已保存至: scripts/ch7_sec3_figure1_panel.png")


if __name__ == '__main__':
    run_chapter7_pipeline()
