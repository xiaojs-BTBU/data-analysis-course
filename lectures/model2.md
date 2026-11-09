# 《实验设计与数据处理》讲义
## 第二章 假设检验
### （假设检验原理与正态/方差齐性前提评估、参数与非参数t检验及效应量、定性变量卡方/Fisher检验与APA三线表自动化导出）

* **适用对象**：食品科学与工程、轻工技术与工程、生物工程及相关专业硕士研究生（零编程基础、弱统计基础）
* **核心工具**：Trae IDE + Python 科研统计双轨工具链（权威一行流 `pingouin` + 基础标准库 `scipy.stats`, `tabulate`, `seaborn`）
* **四步闭环标准**：**讲原理（统计逻辑） $\rightarrow$ 写代码（Trae 结构化提示词设计） $\rightarrow$ 做分析看结果（统计指标诊断） $\rightarrow$ 出图和表（SCI 成果交付）**
* **学时分配**：4 学时（对应大纲第 2 次课：理论讲授 2 学时 + 上机实操 2 学时）

---

## 本章导读与整体逻辑架构

在完成第一章的数据清洗与 Tidy Data 规整后，我们正式进入推断统计的核心——**假设检验（Hypothesis Testing）**。

在食品和生物硕士论文中，最常见的审稿人拒稿意见之一就是：“作者未对数据进行正态性和方差齐性检验便直接使用了 Student's $t$ 检验，且仅报告了 $p$ 值，缺乏效应量（Effect Size）与规范的统计三线表”。

为了让零基础同学彻底告别复杂的数学公式，本章全面引入专为生物医学统计设计的 Python 王牌包 **`pingouin`**（只需 1 行代码即可同时输出 $t$ 值、$p$ 值、自由度、效应量 Cohen's $d$、95% 置信区间和统计检验力 Power！），并同步提供 `scipy.stats` 对照写法。本章拆解为三大递进小节，**依次配备 3 套完全独立的合成科研数据集**：

```
                   【第二章：推断统计与假设检验全流程的三阶闭环】
                                           │
  ┌────────────────────────────────────────┼────────────────────────────────────────┐
  │                                        │                                        │
【第一节】假设检验原理与前提评估          【第二节】常用参数与非参数检验及效应量    【第三节】定性变量卡方检验与 APA 三线表
(反证法 + α/β错误 + 正态与方差齐性体检)   (单样本/独立/配对t检验 vs U/Wilcoxon检验) (Pearson卡方 vs Fisher精确检验 + 三线表)
数据: ch2_sec1_assumptions.csv            数据: ch2_sec2_t_and_nonparam.csv        数据: ch2_sec3_chi2_apa_table.csv
产出: 正态Q-Q图 + 齐性检验全景诊断表       产出: 自动分支检验 + 箱线图与配对连线图   产出: 列联表检验 + APA第七版标准三线表
```

---

## 第一节：假设检验原理与前提评估（正态性与方差齐性检验）

### 1. 讲原理：用“法庭判案”与“上岗体检”听懂假设检验前提
* **小概率反证法与两类错误（$\alpha$ 与 $\beta$）**：
  * **原假设（$H_0$，零假设）**：假定实验处理**没有效果**，两组之间的微小差别纯粹是随机抽样碰巧造成的（$\mu_1 = \mu_2$）；
  * **备择假设（$H_1$）**：实验处理**确实有效**，两组总体均值存在真实差异（$\mu_1 \neq \mu_2$）；
  * **$p$ 值的通俗含义**：如果原假设成立（真的没效果），单靠运气抽出眼前这么大差距的概率。当 **$p < 0.05$** 时，说明“靠运气”的概率不到 5%，我们果断推翻 $H_0$，宣布差异显著！
  * **第一类错误（$\alpha$，假阳性/冤枉好人）**：本来没效果却误判为显著（人为卡在 $\alpha = 0.05$）；
  * **第二类错误（$\beta$，假阴性/放跑坏人）**：本来有效果却没测出来。**统计检验力（Power $= 1 - \beta$）**就是抓出真效果的把握度，合格标准为 $\text{Power} \ge 0.80$。
* **做参数检验前的两大“上岗体检”（千万记住：这两项体检 $p > 0.05$ 才是合格！）**：
  1. **正态性检验 —— Shapiro-Wilk 检验**：
     检验数据是否呈左右对称的钟形正态分布。这里 $H_0$ 是“数据服从正态分布”，因此 **$p > 0.05$ 代表满足正态性（体检过关）**；$p \le 0.05$ 代表偏态分布（体检不过关）。
  2. **方差齐性检验 —— Levene 检验（稳健首选）与 Bartlett 检验**：
     检验两组数据的离散波动程度（胖瘦）是否一致。Levene 检验对非正态更稳健，Bartlett 检验对正态数据灵敏。同样，**$p > 0.05$ 代表两组方差齐（体检过关）**；$p \le 0.05$ 代表方差不齐。

### 2. 配套合成数据一：`data/ch2_sec1_assumptions.csv`
* **实验背景**：测定对照组（`Control`）与益生菌发酵提取物组（`Probiotic`）各 15 只小鼠的 3 个不同分布特征的指标：
  1. `SOD_Normal_Homo`（双组正态，且方差齐）；
  2. `MDA_Normal_Hetero`（双组正态，但处理组个体差异大导致方差不齐）；
  3. `IL6_Skewed`（存在重度炎症极值，呈现明显右偏态非正态分布）。

### 3. 写代码：Trae 结构化提示词设计与 `pingouin` + `scipy` 双轨代码
* **Trae 提示词模板**：
  > **【角色】**：生物统计导师。  
  > **【数据输入】**：读取 `data/ch2_sec1_assumptions.csv`。  
  > **【分析逻辑】**：使用 `pingouin.normality()` 和 `pingouin.homoscedasticity()`（同时附上 `scipy.stats.shapiro` 与 `levene`/`bartlett` 对照）对 3 个指标批量进行正态性与方差齐性检验；自动输出每个指标应该选用的后续检验方法名称；绘制 3 个指标的正态 Q-Q 图（Quantile-Quantile Plot）与核密度分布图。  
  > **【输出要求】**：导出先验体检诊断表 `data/ch2_sec1_assumption_report.csv`，保存 600 DPI 诊断图 `scripts/ch2_sec1_qq_density.png`。

* **Python 权威双轨实现片段（`pingouin` 一行流 vs `scipy` 标准写法）**：
  ```python
  import pandas as pd
  import pingouin as pg      # 导入生物统计神器 pingouin，简写为 pg
  from scipy import stats    # 导入基础科学统计库 scipy.stats

  df1 = pd.read_csv("data/ch2_sec1_assumptions.csv")

  # 【方法 A：强烈推荐！用 pingouin 一行代码做完整正态性与方差齐性体检】
  # 1. 按 Group 分组检验 SOD_Normal_Homo 的正态性（直接返回 W 值、p 值和 normal=True/False 布尔判断！）
  norm_res = pg.normality(data=df1, dv='SOD_Normal_Homo', group='Group', method='shapiro')
  print("=== Pingouin 正态性体检结果 ===\n", norm_res)

  # 2. 一行代码做 Levene 方差齐性检验（直接返回 W 值、p 值和 equal_var=True/False 判断！）
  homo_levene = pg.homoscedasticity(data=df1, dv='SOD_Normal_Homo', group='Group', method='levene')
  homo_bartlett = pg.homoscedasticity(data=df1, dv='SOD_Normal_Homo', group='Group', method='bartlett')
  print("=== Pingouin Levene 方差齐性结果 ===\n", homo_levene)

  # 【方法 B：传统 scipy.stats 写法（供对照理解底层）】
  g_ctrl = df1[df1['Group'] == 'Control']['SOD_Normal_Homo']
  g_prob = df1[df1['Group'] == 'Probiotic']['SOD_Normal_Homo']
  w1, p_norm1 = stats.shapiro(g_ctrl)
  w2, p_norm2 = stats.shapiro(g_prob)
  stat_lev, p_lev = stats.levene(g_ctrl, g_prob, center='median')
  ```

---

## 第二节：常用参数检验与非参数检验（单样本、独立样本、配对样本及效应量）

### 1. 讲原理：一张决策树选对所有的 $t$ 检验与非参数检验
根据第一节的“体检结果”和实验设计类型，所有两组对比只有以下 4 条标准路径：

```
                        【两组或单组数据假设检验全选型地图】
                                         │
        ┌────────────────────────────────┼────────────────────────────────┐
        ▼                                ▼                                ▼
【场景1：单样本 vs 标准值】       【场景2：两独立样本对比】         【场景3：同一批对象前后配对】
(如检测奶粉蛋白是否达到国标)     (如对照组12只鼠 vs 给药组12只鼠)   (如15名受试者服药前 vs 服药后)
        │                                │                                │
  单样本 t 检验                    看第一节体检结果                  先算“前后差值 = 后 - 前”
  pg.ttest(x, 12.0)                      │                           检验“差值”是否正态？
                        ┌────────────────┼────────────────┐               │
                        ▼                ▼                ▼         ┌─────┴─────┐
                   双正态+方差齐    双正态+方差不齐   任一组非正态  差值正态    差值非正态
                        │                │                │         │           │
                 【Student's t】   【Welch's t】  【Mann-Whitney U】【配对t检验】【Wilcoxon符号秩】
                 (correction=False)(correction=True) (pg.mwu(x, y)) (paired=True) (pg.wilcoxon(x,y))
```

* **效应量（Effect Size）——审稿人必看的“真实差距标尺”**：
  * 参数 $t$ 检验报告 **Cohen's $d$**：$|d| = 0.2$（小效应）、$0.5$（中等效应）、$\ge 0.8$（大效应）；
  * 非参数检验报告 **秩双列相关系数 $r_{rb}$（Rank-Biserial Correlation）** 或 **CLES（Common Language Effect Size）**——`pingouin` 会自动帮你算好，根本不用手写公式！

### 2. 配套合成数据二：`data/ch2_sec2_t_and_nonparam.csv`
* **实验背景**：综合包含 3 类实验记录：
  1. 12 罐批次功能奶粉的蛋白质含量 `Protein_Content`（对比国标合格线 $12.0\text{ g/100g}$）；
  2. 两独立组小鼠（`Control` vs `Treatment`）的血清 `SOD`（正态齐性）、`MDA`（正态但异方差）、`TNF_a`（严重偏态）；
  3. 12 只大鼠干预前后的空腹血糖配对数据（`Glucose_Before`, `Glucose_After`）。

### 3. 写代码：Trae 结构化提示词设计与 `pingouin` 一行流代码
* **Trae 提示词模板**：
  > **【角色】**：食品营养与生物统计专家。  
  > **【数据输入】**：读取 `data/ch2_sec2_t_and_nonparam.csv`。  
  > **【分析逻辑】**：1. 使用 `pg.ttest(x, 12.0)` 执行单样本 t 检验；2. 编写“先验自检 $\rightarrow$ 自动分支”函数，对 `SOD`、`MDA`、`TNF_a` 自动选用 Student's t (`pg.ttest(correction=False)`)、Welch's t (`pg.ttest(correction=True)`) 或 Mann-Whitney U (`pg.mwu()`)，提取统计量、p 值、效应量 Cohen's d / RBC、95% 置信区间和 Power；3. 对血糖前后数据执行配对检验（`pg.ttest(paired=True)` 与 `pg.wilcoxon()`）；4. 绘制独立样本半透明箱线+散点图与配对样本个体前后斜率连线图。  
  > **【输出要求】**：导出全套检验结果表至 `data/ch2_sec2_test_results.csv`，保存 600 DPI 组合图 `scripts/ch2_sec2_hypothesis_plots.png`。

* **Python 权威双轨实现片段（感受 `pingouin` 一行顶十行的威力！）**：
  ```python
  import pingouin as pg

  # 1. 单样本 t 检验（对比国标 12.0）：一行输出 T, dof, p-val, CI95%, cohen-d, power！
  res_1samp = pg.ttest(df2['Protein_Content'].dropna(), 12.0)

  # 2. 两独立样本 Student's t 检验 (correction=False) 与 Welch's t 检验 (correction=True)
  res_student = pg.ttest(ctrl_sod, treat_sod, paired=False, correction=False)
  res_welch   = pg.ttest(ctrl_mda, treat_mda, paired=False, correction=True)

  # 3. 两独立样本非参数 Mann-Whitney U 检验：一行输出 U-val, p-val, RBC(秩双列效应量), CLES！
  res_mwu = pg.mwu(ctrl_tnf, treat_tnf, alternative='two-sided')

  # 4. 配对样本 t 检验 (paired=True) 与非参数 Wilcoxon 符号秩检验
  res_paired_t = pg.ttest(df2['Glucose_Before'], df2['Glucose_After'], paired=True)
  res_wilcoxon = pg.wilcoxon(df2['Glucose_Before'], df2['Glucose_After'])
  ```

---

## 第三节：定性变量卡方检验与标准学术 APA 三线表自动化导出

### 1. 讲原理：分类频数怎么比？APA 三线表怎么自动生成？
* **列联表分析（Contingency Table）：Pearson 卡方 vs Fisher 精确检验**：
  当因变量是定性分类指标（如：保鲜果实的“未腐烂 vs 腐烂”、实验动物的“存活 vs 死亡”、感官盲测的“喜欢 vs 不喜欢”），我们先统计出 $2 \times 2$ 频数四格表：
  * **选用金标准**：
    1. 总样本量 $N \ge 40$ 且 4 个格子的**理论期望频数 $T$ 全部 $\ge 5$** $\rightarrow$ 用标准 **Pearson 卡方检验**（`pg.chi2_independence()`，不带修正）；
    2. $N \ge 40$ 但有 1 个格子的理论频数在 $1 \le T < 5$ 之间 $\rightarrow$ 用 **Yates 连续性修正卡方**；
    3. 只要总样本量 $N < 40$，或者任一格子的理论频数 $T < 5$ $\rightarrow$ 直接用 `scipy.stats.fisher_exact()` 做 **Fisher 精确检验**，并报告**比值比（Odds Ratio, OR）**！
* **APA 第七版标准学术三线表（Three-Line Table）规范**：
  * 表格只有 3 条横线：**顶线（粗，1.5 pt）、栏目线/表头下线（细，0.75 pt）、底线（粗，1.5 pt）**，绝无任何竖线和内部横格线；
  * 正态分布指标汇报 `Mean ± SD`，偏态分布指标汇报 `Median (IQR)`，定性指标汇报 `n (%)`；
  * 必须同时列出：**指标名、各组描述统计、检验方法名称、统计量具体值（$t / U / \chi^2$）、精确 $p$ 值（如 $p = 0.003$，小于 $0.001$ 写 $<0.001$）、效应量（Cohen's $d$ / OR）**。

### 2. 配套合成数据三：`data/ch2_sec3_chi2_apa_table.csv`
* **实验背景**：天然抑菌涂膜处理鲜切草莓保鲜实验（`Control` 对照组 30 盒 vs `Coating` 涂膜组 30 盒），同时记录了第 7 天的霉变状态（`Spoilage`: Spoiled vs Intact，定性变量）、维生素 C 含量（`Vitamin_C`, 正态定量）、失重率（`Weight_Loss`, 异方差定量）与菌落总数对数值（`Log_CFU`, 偏态定量）。

### 3. 写代码：Trae 结构化提示词设计与 APA 三线表一键导出
* **Trae 提示词模板**：
  > **【角色】**：SCI 期刊统计审稿人与排版专家。  
  > **【数据输入】**：读取 `data/ch2_sec3_chi2_apa_table.csv`。  
  > **【分析逻辑】**：1. 对定性列 `Spoilage` 生成列联表，调用 `pg.chi2_independence()` 与 `scipy.stats.fisher_exact()` 检验霉变率差异，绘制堆叠百分比柱状图；2. 对全部定量与定性指标执行自动化检验，自动组装成一张符合 APA 第七版规范的汇总三线表（同时导出为 Markdown 格式、CSV 格式以及可直接贴入论文的 Word/图片渲染三线表）。  
  > **【输出要求】**：导出三线表 `data/ch2_sec3_apa_summary_table.csv`（编码 `utf-8-sig`），保存 600 DPI 图表 `scripts/ch2_sec3_chi2_and_table.png`。

---

## 四、第二章课后实战作业指引

1. **课后作业数据集**：
   * `data/ch2_hw_sec1_normality.csv`：不同提取工艺下植物蛋白溶解度与粒径分布的正态/方差齐性体检练习；
   * `data/ch2_hw_sec2_group_compare.csv`：益生元干预前后受试者肠道双歧杆菌丰度（配对）及两组血清炎性因子（独立样本）检验练习；
   * `data/ch2_hw_sec3_clinical_table.csv`：功能性食品人群试食试验综合基线与疗效指标 APA 三线表一键生成练习。
2. **思考与简答题**：
   * 为什么在正态性检验（Shapiro-Wilk）和方差齐性检验（Levene）中，$p > 0.05$ 代表满足条件，而在 $t$ 检验中却是 $p < 0.05$ 代表差异显著？请从两者的“原假设 $H_0$ 是什么”来解释。
   * 当你的实验每组只有 6 个平行样品时，发现有一组有一个明显偏高的正常生物个体导致 Shapiro-Wilk 检验 $p = 0.012$，此时你应该选用什么检验方法？在三线表中应该如何表示该组数据的集中趋势和离散程度？
