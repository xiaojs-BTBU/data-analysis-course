# 《实验设计与数据处理》（AI 融合实战版）全套教学资源总览

**适用专业**：食品科学与工程、生物工程、生物与医药、营养与健康等专业硕士研究生  
**主讲教师**：肖俊松 副教授（北京工商大学 食品与健康学院）  
**开发环境**：Miniconda (`Python 3.10+`) + Trae IDE + 权威统计分析库 (`pingouin`, `scipy`, `statsmodels`, `scikit-learn`, `pyDOE2`, `seaborn`, `matplotlib`)

---

## 一、 文件夹目录结构

1. **`data/`（全套课堂示范与课后作业数据集，共 53 个 CSV 文件）**
   - **命名规范**：
     - 课堂示范数据集：`ch{章号}_sec{节号}_{主题英文简称}.csv`（例如 `ch3_sec1_one_and_twoway.csv`）
     - 课后作业数据集：`ch{章号}_hw_sec{节号}_{主题英文简称}.csv`（例如 `ch3_hw_sec1_preservation.csv`）
     - 脚本运行自动导出的统计报表/APA三线表：`ch{章号}_sec{节号}_{结果名称}.csv`
   - **编码规范**：全部采用 `UTF-8-SIG` 编码保存，既支持 Python/R 直接读取，也支持在中文版 Microsoft Excel / WPS 中直接双击打开且绝不乱码。

2. **`scripts/`（零基础超详细逐行注释代码与 600 DPI 出版级图表）**
   - **命名规范**：按“章号 + 核心分析目的”命名，方便师生按图索骥：
     - `ch1_00_env_and_chinese_check.py`：科研环境自检、中文路径/编码与 Matplotlib 中文及负号全局防乱码验证
     - `ch1_01_data_wrangling_and_cleaning.py`：技术重复合并、宽表转 Tidy 长表 (`pd.melt`)、缺失值组内中位数插补与 $3\sigma$/IQR 离群值诊断
     - `ch2_hypothesis_testing_and_apa_table.py`：Shapiro-Wilk 正态性检验、Levene 方差齐性检验、参数/非参数假设检验自动路由及 APA 三线表一键导出
     - `ch3_anova_posthoc_and_raincloud.py`：Plackett-Burman 主效应筛选、单因素/双因素 ANOVA、Tukey HSD 事后检验、紧凑字母标记法 (CLD: `a, b, bc, c`) 自动生成与复合云雨图 (Raincloud Plot)
     - `ch4_correlation_regression_and_kinetics.py`：Pearson/Spearman 相关性热图、ELISA 标准曲线逆向定量、多元回归 VIF 共线性诊断、热失活一级动力学与 Arrhenius 货架期非线性拟合 (`scipy.optimize.curve_fit`)
     - `ch5_doe_orthogonal_rsm_and_mixture.py`：$L_9(3^4)$ 正交试验极差分析 ($R$) 与方差分析、Box-Behnken 响应面 (RSM) 二次多项式建模与 3D 曲面/等高线寻优、三元单纯形重心混料设计回归
     - `ch6_multivariate_pca_hca_and_plsda.py`：量纲 Z-score 标准化、PCA 双标图 (Biplot + 95% 置信椭圆)、双向层次聚类热图 (`sns.clustermap`)、有监督 PLS-DA 判别、置换检验与 VIP 特征筛选
     - `ch7_transcriptomics_volcano_heatmap_panel.py`：Benjamini-Hochberg (FDR) 多重检验校正、防重叠靶基因标注火山图、Row Z-score 聚类热图、GO/KEGG 富集气泡图与 SCI 整版三联复合图 (`Figure 1 Panel A+B+C`)
     - `ch8_machine_learning_rf_svm_evaluation.py`：防数据泄露数据集划分、`Pipeline` 交叉验证与网格搜索、随机森林 (RF) 与支持向量机 (SVM) 分类/回归评估、ROC/AUC 曲线、混淆矩阵与特征重要性图
   - **代码特色**：
     - **零基础友好**：每一行代码、每一个核心函数参数均配有详尽的中文注释，解释“为什么要这样设”以及“不这样设会踩什么坑”。
     - **双轨权威实现**：同时提供 `pingouin` / `statsmodels` / `scikit-learn` 极简权威调用方式与底层 `scipy` / `numpy` 稳健回退实现，确保任何电脑环境下 100% 零报错运行。

3. **`ppt/`（第 1 ~ 8 章配套教学课件 PPT，基于北京工商大学食品与健康学院模板制作）**
   - `chapter1_试验设计与数据预处理基础.pptx`
   - `chapter2_假设检验.pptx`
   - `chapter3_因素筛选方差分析与多重比较.pptx`
   - `chapter4_相关分析回归分析和拟合.pptx`
   - `chapter5_实验设计与分析.pptx`
   - `chapter6_多元统计分析.pptx`
   - `chapter7_转录组学数据分析初步.pptx`
   - `chapter8_机器学习与预测模型构建.pptx`

4. **`lectures/`（第 1 ~ 8 章完整版课堂讲义 Markdown）**
   - `model1.md` ~ `model8.md`

---

## 二、 快速运行方法

在 Trae IDE 终端或 Anaconda Prompt 中激活环境后，在课程根目录下执行任意脚本即可一键复现全部统计表与 600 DPI 高清配图：

```bash
python scripts/ch1_00_env_and_chinese_check.py
python scripts/ch1_01_data_wrangling_and_cleaning.py
python scripts/ch2_hypothesis_testing_and_apa_table.py
python scripts/ch3_anova_posthoc_and_raincloud.py
python scripts/ch4_correlation_regression_and_kinetics.py
python scripts/ch5_doe_orthogonal_rsm_and_mixture.py
python scripts/ch6_multivariate_pca_hca_and_plsda.py
python scripts/ch7_transcriptomics_volcano_heatmap_panel.py
python scripts/ch8_machine_learning_rf_svm_evaluation.py
```
