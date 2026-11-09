# 《实验设计与数据处理》讲义
## 第四章 相关分析、回归分析和拟合
### （Pearson与Spearman相关及显著性星号热图、一元与多元线性回归及VIF共线性/残差诊断、一级/二级反应动力学与S型生长曲线非线性拟合）

* **适用对象**：食品科学与工程、轻工技术与工程、生物工程及相关专业硕士研究生（零编程基础、弱统计基础）
* **核心工具**：Trae IDE + Python 科研分析工具链（`pingouin`, `scipy.optimize`, `statsmodels`, `pandas`, `seaborn`, `matplotlib`）
* **四步闭环标准**：**讲原理（统计逻辑） $\rightarrow$ 写代码（Trae 结构化提示词设计） $\rightarrow$ 做分析看结果（统计指标诊断） $\rightarrow$ 出图和表（SCI 成果交付）**
* **学时分配**：4 学时（对应大纲第 4 次课：理论讲授 2 学时 + 上机实操 2 学时）

---

## 本章导读与整体逻辑架构

在前三章中，我们的自变量都是“分门别类的组名”（如对照组 vs 处理组、不同菌株），主要看组与组之间有没有差异。但科研中还有另一大半核心任务：**当自变量和因变量全都是连续变化的数字时，它们之间有什么定量规律？**

* 比如：发酵乳里的“总酚、黄酮、有机酸、活菌数”与“抗氧化活性、感官评分”之间，**谁跟谁关系最铁**？（**第一节：相关性分析与星号矩阵热图**）
* 比如：能不能建立一个数学公式，输入“总酚、黄酮、多糖含量”，就能**直接算出抗氧化活性**？如果几个指标高度重合怎么办？（**第二节：一元/多元线性回归与 VIF 共线性诊断**）
* 比如：食品保藏过程中维生素 C 的热降解、杀菌时微生物的死亡、发酵罐里菌体的 S 型生长，它们**随时间变化的动力学曲线方程和半衰期**怎么求？（**第三节：非线性回归与科研动力学拟合**）

本章拆解为三大递进小节，**依次配备 3 套完全独立的合成科研数据集**：

```
                   【第四章：从指标关联热图到因果回归及动力学方程的三阶闭环】
                                           │
  ┌────────────────────────────────────────┼────────────────────────────────────────┐
  │                                        │                                        │
【第一节】相关性分析与带星号矩阵热图       【第二节】一元/多元线性回归与共线性诊断  【第三节】非线性回归与科研动力学拟合
(Pearson线性 vs Spearman秩相关 + 伪相关)   (OLS最小二乘 + VIF共线性 + 残差四合一图) (一级/二级降解动力学 + Logistic/Gompertz)
数据: ch4_sec1_correlation.csv             数据: ch4_sec2_regression_vif.csv        数据: ch4_sec3_kinetics_growth.csv
产出: 相关与p值矩阵 + 下三角星号冷暖热图   产出: 回归三线表 + VIF表 + 残差诊断图    产出: 动力学速率k/半衰期 + 95%置信带曲线
```

---

## 第一节：相关性分析（Pearson 线性相关与 Spearman 秩相关及显著性热图）

### 1. 讲原理：两个指标同步涨跌，怎么算相关系数？怎么防“伪相关”？
* **相关性分析的本质（对称关系，不区分谁是因、谁是果）**：
  相关系数（$r$ 或 $\rho$）取值在 $[-1, +1]$ 之间：
  * $r > 0$（正相关，涂暖红色）：你高我也高（如多酚含量越高，自由基清除率越高）；
  * $r < 0$（负相关，涂冷蓝色）：你高我偏低（如过氧化值越高，感官风味评分越低）；
  * $|r|$ 越接近 $1$ 说明关联越紧密（通常 $|r| \ge 0.7$ 为强相关，$0.4 \sim 0.7$ 为中等相关，$< 0.4$ 为弱相关），但**必须同时满足 $p < 0.05$ 才具有统计学意义**！
* **Pearson 相关 vs Spearman 秩相关，到底选哪个？**
  1. **Pearson 线性相关系数（$r$）**：要求两个变量都是**正态分布的连续数值**，且散点图呈**直线趋势**；
  2. **Spearman 秩相关系数（$\rho$）**：先把数值从小到大排成名次（秩次 $1, 2, 3\dots$）再算相关。凡是遇到**偏态分布、有极端离群值、感官评价等级（定序变量），或者呈曲线单调递增/递减**，一律选 Spearman！
* **警惕“伪相关（Spurious Correlation）”**：
  夏天冰淇淋销量与溺水人数呈极显著正相关（$r = 0.89^{***}$），难道吃冰淇淋会导致溺水吗？当然不是！是因为第三个隐藏变量“高温天气”同时推高了两者。所以在论文讨论中，**统计相关必须结合生物化学机理才能推导因果**！

### 2. 配套合成数据一：`data/ch4_sec1_correlation.csv`
* **实验背景**：测定 35 批次传统发酵果酒的 6 项核心理化与感官指标：总酚（`Total_Phenols`）、花青素（`Anthocyanins`）、滴定酸度（`Titratable_Acid`）、残糖（`Residual_Sugar`）、ABTS 抗氧化能力（`ABTS_Capacity`）及感官品评等级分（`Sensory_Rank`）。

### 3. 写代码：Trae 结构化提示词设计与 `pingouin` 一行流代码
* **Trae 提示词模板**：
  > **【角色】**：食品风味与理化关联分析专家。  
  > **【数据输入】**：读取 `data/ch4_sec1_correlation.csv`。  
  > **【分析逻辑】**：1. 使用 `pg.pairwise_corr()` 一行代码批量计算所有指标间的 Pearson 和 Spearman 相关系数、精确 p 值与 95% 置信区间；2. 构建 6x6 相关系数矩阵与显著性星号矩阵（`* p<0.05, ** p<0.01, *** p<0.001`）；3. 绘制 SCI 顶刊风格的下三角相关性热图（使用 `RdBu_r` 冷暖发散色阶，格内同时标注相关系数与星号）。  
  > **【输出要求】**：导出相关分析汇总表 `data/ch4_sec1_corr_table.csv`，保存 600 DPI 热图 `scripts/ch4_sec1_corr_heatmap.png`。

* **Python 权威双轨实现片段**：
  ```python
  import pandas as pd
  import pingouin as pg

  df1 = pd.read_csv("data/ch4_sec1_correlation.csv")
  # 使用 pingouin 一行代码算出所有变量两两之间的相关系数 r、p-val、CI95% 和 power！
  corr_pearson  = pg.pairwise_corr(df1, method='pearson')
  corr_spearman = pg.pairwise_corr(df1, method='spearman')
  print("=== Pingouin 两两相关性分析表 ===\n", corr_pearson[['X', 'Y', 'r', 'p-val', 'CI95%']])
  ```

---

## 第二节：一元与多元线性回归（最小二乘法 OLS、VIF 共线性诊断与残差四合一图）

### 1. 讲原理：从“相关”跨越到“因果预测方程”，怎么防共线性翻车？
* **相关 vs 回归的本质区别**：
  相关分析里 $X$ 和 $Y$ 地位平等；而**回归分析（Regression）**明确规定 **$X$ 是自变量（原因/输入特征），$Y$ 是因变量（结果/预测目标）**，目的是建立定量预测方程：
  * **一元线性回归**：$Y = \beta_0 + \beta_1 X + \epsilon$（一条直线，斜率 $\beta_1$ 代表 $X$ 每增加 1 个单位，$Y$ 平均变化多少）；
  * **多元线性回归**：$Y = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \dots + \beta_p X_p + \epsilon$。
* **怎么求方程系数？普通最小二乘法（OLS, Ordinary Least Squares）**：
  找到一条直线（或平面），让所有真实数据点到这条线的**垂直距离平方和（残差平方和 $\sum (Y_{\text{真实}} - \hat{Y}_{\text{预测}})^2$）最小**！
* **看懂回归表的三个核心指标**：
  1. **决定系数 $R^2$ 与调整后 $R^2_{adj}$**：代表方程能解释因变量百分之多少的变化（越接近 1 越好）。多元回归必看 **$R^2_{adj}$**（因为往方程里乱塞无关变量也会让 $R^2$ 虚胖，而 $R^2_{adj}$ 会惩罚无效变量）；
  2. **模型整体 $F$ 检验 $p$ 值**：$p < 0.05$ 说明方程整体有效；
  3. **各个自变量系数的 $t$ 检验 $p$ 值**：判断哪几个自变量是真正起作用的功臣。
* **多元回归头号陷阱：多重共线性（Multicollinearity）与 VIF 诊断**：
  如果方程里的两个自变量本身就高度相关（比如同时放入“总黄酮”和“芦丁含量”，两者相关性高达 $0.96$），计算机在分配功劳时就会精神分裂——算出来的回归系数标准误暴涨，甚至把本来正作用的物质算出负数！
  * **必做体检：方差膨胀因子（VIF, Variance Inflation Factor）**：
    $$\text{VIF}_j = \frac{1}{1 - R_j^2}$$
    * **$\text{VIF} < 5$**：无共线性，方程健康；
    * **$5 \le \text{VIF} \le 10$**：中度共线性；
    * **$\text{VIF} > 10$**：**严重多重共线性！必须剔除其中一个重复变量，或改用第六章的偏最小二乘（PLS）！**
* **残差诊断四合一图（验证回归前提是否成立）**：
  1. **Residuals vs Fitted（残差 vs 拟合值）**：点应随机均匀分布在 0 上下，若呈弯月形说明该用非线性拟合，若呈喇叭口说明方差不齐；
  2. **Normal Q-Q（残差正态 Q-Q 图）**：点应紧贴 45° 对角线，证明残差服从正态分布；
  3. **Scale-Location（标准化残差平方根图）**：红线越平，方差齐性越好；
  4. **Residuals vs Leverage（残差 vs 杠杆值）**：揪出以一己之力把整条回归线拉歪的强影响异常点（Cook's Distance）。

### 2. 配套合成数据二：`data/ch4_sec2_regression_vif.csv`
* **实验背景**：测定 45 批植物提取物的总酚（`Total_Phenols`）、总黄酮（`Flavonoids`）、儿茶素衍生物（`Catechin_Deriv`，故意设置为与总酚高度共线性 $\text{VIF}>15$）、多糖（`Polysaccharides`）及还原糖（`Reducing_Sugar`），预测抗氧化自由基清除率（`DPPH_Rate`）。

### 3. 写代码：Trae 结构化提示词设计与 Python 实现
* **Trae 提示词模板**：
  > **【角色】**：计量回归建模专家。  
  > **【数据输入】**：读取 `data/ch4_sec2_regression_vif.csv`。  
  > **【分析逻辑】**：1. 使用 `pg.linear_regression()` 和 `statsmodels` 先拟合 `Total_Phenols` 对 `DPPH_Rate` 的一元线性回归，绘制带 95% 置信区间和 95% 预测区间的散点拟合图；2. 拟合全变量多元线性回归，计算每个自变量的 VIF 值，展示 `Catechin_Deriv` 与 `Total_Phenols` 严重共线性（$\text{VIF}>10$）导致的系数失真；3. 剔除冗余共线性变量后重新拟合最优多元回归方程，并绘制标准残差诊断四合一图。  
  > **【输出要求】**：导出共线性整改前后的回归对比三线表 `data/ch4_sec2_regression_table.csv`，保存 600 DPI 残差四合一图 `scripts/ch4_sec2_residual_diagnostics.png`。

---

## 第三节：非线性回归与科研动力学拟合（一级/二级动力学与 S 型生长曲线）

### 1. 讲原理：食品保质期劣变与微生物发酵生长的四大经典非线性方程
在食品科学和生物工程中，营养素降解、酶促反应和微生物生长遵循严格的物理化学或生物学定律，绝不能硬拿直线去套，而要用 `scipy.optimize.curve_fit` 做**非线性最小二乘拟合**：

| 科研经典模型名称 | 数学方程表达式 | 核心动力学参数与物理意义 | 食品与生物典型应用场景 |
| :--- | :--- | :--- | :--- |
| **一级反应动力学**<br>(First-Order Kinetics) | $C_t = C_0 \cdot e^{-k_1 t}$ | $k_1$: 一级速率常数 ($\text{d}^{-1}$)<br>半衰期 $t_{1/2} = \frac{\ln 2}{k_1} \approx \frac{0.693}{k_1}$ | 维生素C、花青素、叶绿素热降解；热杀菌微生物对数失活 |
| **二级反应动力学**<br>(Second-Order Kinetics) | $C_t = \frac{C_0}{1 + C_0 k_2 t}$ | $k_2$: 二级速率常数<br>半衰期 $t_{1/2} = \frac{1}{k_2 C_0}$ | 赖氨酸美拉德褐变损耗、双分子脂质氧化反应 |
| **Logistic S 型方程**<br>(对称 S 型生长曲线) | $N_t = \frac{K}{1 + \left(\frac{K - N_0}{N_0}\right) e^{-r t}}$ | $K$: 最大环境容纳量（最高菌浓）<br>$r$: 最大比生长速率 | 发酵罐细菌/酵母菌经历延滞期$\rightarrow$对数期$\rightarrow$稳定期生长 |
| **Gompertz S 型方程**<br>(非对称 S 型生长曲线) | $N_t = K \cdot \exp\left[-\exp\left(\frac{\mu_m e}{K}(\lambda - t) + 1\right)\right]$ | $\mu_m$: 最大生长速率（拐点斜率）<br>$\lambda$: 延滞期时长 (Lag phase, h) | 食品预测微生物学（精准计算致病菌延滞期 $\lambda$ 与货架期） |

* **拟合优度看什么？**
  看**决定系数 $R^2$**（越接近 $1.0$ 越好，动力学拟合通常要求 $R^2 > 0.95$）与**均方根误差 RMSE**（越小越好）。同时在图上画出 **95% 置信区间带（Confidence Band）**。

### 2. 配套合成数据三：`data/ch4_sec3_kinetics_growth.csv`
* **实验背景**：包含两组经典时间序列实验——① 蓝莓花青素提取液在 70℃ 热处理 0~180 min 内的残留率衰减数据（对比一级动力学 vs 二级动力学拟合优度）；② 植物乳杆菌在 37℃ 发酵 0~24 h 的菌体密度 $\text{OD}_{600}$ 生长数据（对比 Logistic vs Gompertz S 型模型，求延滞期 $\lambda$ 与最大比生长速率）。

### 3. 写代码：Trae 结构化提示词设计与非线性拟合代码
* **Trae 提示词模板**：
  > **【角色】**：食品化学动力学与预测微生物学建模专家。  
  > **【数据输入】**：读取 `data/ch4_sec3_kinetics_growth.csv`。  
  > **【分析逻辑】**：1. 定义一级反应与二级反应动力学函数，调用 `scipy.optimize.curve_fit` 拟合花青素热降解数据，计算速率常数 $k$、参数标准误、半衰期 $t_{1/2}$、$R^2$ 与 RMSE；2. 定义 Logistic 与修正 Gompertz 方程拟合植物乳杆菌生长曲线，求取最大菌浓 $K$、最大生长速率 $\mu_m$ 和延滞期 $\lambda$；3. 基于雅可比矩阵或残差标准差绘制带 95% 置信带的双联动力学拟合大图，图内用文本框标注拟合方程与参数。  
  > **【输出要求】**：导出动力学模型参数对比三线表 `data/ch4_sec3_kinetic_params.csv`，保存 600 DPI 动力学图 `scripts/ch4_sec3_kinetic_curves.png`。

* **Python 权威实现片段**：
  ```python
  import numpy as np
  from scipy.optimize import curve_fit

  # 1. 定义一级降解动力学函数与修正 Gompertz 微生物生长函数
  def first_order_decay(t, C0, k1):
      return C0 * np.exp(-k1 * t)

  def gompertz_growth(t, K, mu_m, lam):
      # K: 渐近最大菌浓, mu_m: 最大比生长速率, lam: 延滞期时长 (Lag time)
      return K * np.exp(-np.exp((mu_m * np.e / K) * (lam - t) + 1.0))

  # 2. 调用 curve_fit 一行代码拟合参数 popt 与协方差矩阵 pcov
  popt, pcov = curve_fit(first_order_decay, t_data, c_data, p0=[100.0, 0.02])
  C0_fit, k1_fit = popt
  se_C0, se_k1 = np.sqrt(np.diag(pcov))  # 对角线开根号即为参数的标准误 (Standard Error)
  t_half = np.log(2) / k1_fit            # 计算半衰期
  ```

---

## 四、第四章课后实战作业指引

1. **课后作业数据集**：
   * `data/ch4_hw_sec1_sensory_corr.csv`：烘焙面包质构剖面分析（TPA：硬度、弹性、咀嚼性）与风味物质的 Pearson/Spearman 相关性练习；
   * `data/ch4_hw_sec2_nir_regr.csv`：多波长近红外吸光度预测奶粉蛋白质含量的多元线性回归与 VIF 共线性剔除练习；
   * `data/ch4_hw_sec3_shelf_life.csv`：鲜榨橙汁在不同储藏温度下维生素 C 降解动力学拟合及金黄色葡萄球菌 Gompertz 生长曲线拟合练习。
2. **思考与简答题**：
   * 为什么在多元线性回归中，即使模型整体的 $R^2$ 高达 $0.92$，但如果其中两个自变量的 $\text{VIF} = 24.5$，你也绝不能直接拿回归系数的大小去解释哪个自变量贡献大？
   * 在预测微生物学中，Gompertz 模型求出的参数 $\lambda$（延滞期）对冷链食品保鲜有什么直接的实际指导意义？
