# 《实验设计与数据处理》讲义
## 第一章 试验设计与数据预处理基础
### （零基础环境搭建与中文避坑、试验设计三原则、变量尺度与Tidy Data重塑、缺失值与离群点清洗）

* **适用对象**：食品科学与工程、轻工技术与工程、生物工程及相关专业硕士研究生（零编程基础、弱统计基础）
* **核心工具**：Trae IDE + Miniconda + Python (`pandas`, `numpy`, `scipy`, `pingouin`, `seaborn`, `matplotlib`) / R
* **四步闭环标准**：**讲原理（统计逻辑） $\rightarrow$ 写代码（Trae 结构化提示词设计） $\rightarrow$ 做分析看结果（统计指标诊断） $\rightarrow$ 出图和表（SCI 成果交付）**
* **学时分配**：4 学时（对应大纲第 1 次课：理论讲授 2 学时 + 上机实操 2 学时）

---

## 本章导读与整体逻辑架构

欢迎来到《实验设计与数据处理》的第一堂课！很多食品和生物专业的同学一听到“写代码处理数据”就心里发怵。其实在 AI 时代，你完全不需要死记硬背复杂的英文编程语法——**Trae IDE（AI 结对编程工具）就是你随叫随到的“免费编程师兄”**，你只需要懂实验设计原则、懂数据结构规范，用中文把需求说清楚，AI 就能自动帮你写出高质量的统计与绘图代码。

本章作为全课程的基石，首先手把手带你装好 **Miniconda、Python、R 与 Trae IDE**，特别是**彻底解决中国学生在 Windows 电脑上最头疼的“四大中文报错难题”**；随后依次通过 3 套独立的食品生物合成数据集，掌握试验设计三原则、Tidy Data 长宽表重塑，以及缺失值与离群点的标准化清洗：

```
                   【第一章：从零基础环境配置到科研数据清洗的四阶闭环】
                                           │
  ┌────────────────────────────────────────┼────────────────────────────────────────┐
  │                                        │                                        │
【第零节】软件安装、Trae配置与中文避坑     【第一节】试验设计三原则与科学逻辑        【第二节】变量尺度与 Tidy Data 规范
(Miniconda + Python/R + Trae + 四大中文坑) (随机/重复/区组 + OFAT单因素轮换局限)    (定类/定序/定量 + 宽表一键转长表)
配套脚本: module1_part0_env_check.py       数据: ch1_sec1_principles_demo.csv       数据: ch1_sec2_tidy_wide.csv
产出: 环境自检通过 + 中文字体测试图        产出: 随机化分组实验工单 + 批次效应图    产出: 标准Tidy长表 + 多指标分面箱线图
  │                                        │                                        │
  └────────────────────────────────────────┴────────────────────────────────────────┘
                                           │
                                           ▼
                      【第三节】数据清洗与质量控制 (缺失与离群点)
                      (组内中位数插补 + 3σ原则 vs 箱线图IQR判别)
                      数据: ch1_sec3_cleaning_raw.csv
                      产出: 离群点核查清单 + 清洗前后分布对比图 + 课后作业
```

---

## 第零节：零基础软件安装、Trae 界面导览与“中文技术大坑”根治指南

### 1. 为什么要装 Miniconda + Python/R + Trae？（打个实验室比方）
* **Python 和 R 语言**：就像实验室里的**超高效全自动分析仪**（Python 擅长统计、机器学习和全能绘图；R 语言在转录组生信包上历史悠久）；
* **Miniconda（环境管家）**：就像实验室里的**独立无菌操作间**。它能把我们这门课需要的 Python 和所有统计包（`pandas`, `pingouin`, `seaborn` 等）装在一个名叫 `stats_env` 的独立房间里，绝不弄乱你电脑原本的系统；
* **Trae IDE（AI 智能操作台）**：字节跳动推出的免费 AI 原生代码编辑器（界面全中文友好）。你在右侧对话框用中文下指令，它在左侧自动帮你写代码、改报错、跑出 SCI 图表！

---

### 2. 第一步：安装 Miniconda、Python 与 R（避开中文路径第 1 坑！）

> ⚠️ **【中国学生头号大坑：Windows 中文用户名与中文安装路径】**  
> 很多同学电脑开机名叫 `张三`，默认安装路径就会变成 `C:\Users\张三\Miniconda3`。很多国际开源统计包底层不支持中文路径，一旦路径里有汉字，装包和画图就会疯狂报 `UnicodeDecodeError` 或找不到文件！  
> **【铁律解决方案】**：安装任何科研软件和存放本课程代码数据，**路径里一律只准用英文字母、数字和下划线**！强烈建议直接在 D 盘或 C 盘根目录新建英文文件夹！

#### （1）下载与安装 Miniconda（5 分钟搞定）
1. 打开清华大学开源软件镜像站（国内极速下载）：访问 `https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/`，下载最新的 `Miniconda3-latest-Windows-x86_64.exe`（苹果电脑选 `MacOSX-arm64.pkg`）。
2. 双击安装，一路点击 `Next` $\rightarrow$ `I Agree` $\rightarrow$ 选择 `Just Me`。
3. **【关键一步】Destination Folder（安装路径）**：千万不要装在带中文的目录！手动修改为纯英文路径，例如：
   `D:\Miniconda3` 或 `C:\Miniconda3`
4. 勾选 **`Add Miniconda3 to my PATH environment variable`**（虽然提示变红，但对零基础同学最省心，勾上它！），点击 `Install` 完成安装。

#### （2）配置国内清华极速镜像源 & 创建课程专属环境
按键盘 `Win + S` 搜索并打开 **Anaconda Prompt (Miniconda3)** 黑窗口，依次复制粘贴以下 3 段命令（每粘贴一段按一次回车）：

```bash
# 第1段：配置清华镜像源（解决国内下载统计包卡死、超时失败问题）
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --set show_channel_urls yes
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple

# 第2段：创建本课程专属的 Python 3.10 独立科研环境（名字叫 stats_env）
conda create -n stats_env python=3.10 -y

# 第3段：激活进入该环境，并一次性装齐本课程全套 8 章权威统计与绘图包！
conda activate stats_env
pip install numpy pandas scipy statsmodels pingouin scikit-learn matplotlib seaborn pyDOE2 adjustText openpyxl tabulate
```
*(注：如果部分章节想同步体验 R 语言，只需在同一窗口运行 `conda install -n stats_env -c conda-forge r-base r-essentials -y` 即可在同一环境内装好 R 语言！)*

---

### 3. 第二步：Trae IDE 界面导览与环境绑定（3 分钟上手）

1. **下载安装**：访问官网 `https://www.trae.com.cn/` 下载安装 Trae IDE，打开后选择**简体中文**界面。
2. **创建纯英文课程工作区文件夹**：
   在电脑 D 盘新建一个纯英文文件夹：`D:\Exp_Design_Course`，在里面再建两个子文件夹：
   * `D:\Exp_Design_Course\data`（专门放老师发的和自己的 CSV 实验数据）
   * `D:\Exp_Design_Course\scripts`（专门放 Python 代码和生成的图表）
   在 Trae 中点击左上角 **“文件 $\rightarrow$ 打开文件夹”**，选中 `D:\Exp_Design_Course`。
3. **认识 Trae 四大核心区域**：
   * **左侧栏（资源管理器）**：像 Windows 文件夹一样，看到你的 `data/` 和 `scripts/`；
   * **中间区（代码编辑与图表预览区）**：点开 `.csv` 表格或 `.py` 代码、`.png` 图片都在这里看；
   * **底部栏（终端 Terminal）**：按快捷键 ``Ctrl + ` ``（数字1左边那个键）呼出，用来看程序运行打印的统计结果；
   * **右侧栏（AI 对话框 Chat / Builder）**：按 `Ctrl + U` 呼出！
     * **Chat 模式**：你问它统计问题、让它解释代码报错；
     * **Builder 模式（最强大！）**：你把讲义里的【Trae 结构化提示词】粘贴进去，它会自动帮你创建 `.py` 代码文件、自动运行、遇到报错还会自己改代码直到跑出图片为止！
4. **绑定我们刚装好的 `stats_env` 解释器**：
   按快捷键 `Ctrl + Shift + P`，输入 `Python: Select Interpreter`（选择解释器），在弹出的列表里选中带有 **`stats_env` (`D:\Miniconda3\envs\stats_env\python.exe`)** 的那一项。从此以后，点右上角的三角形运行按钮 `▷`，就会自动调用装好全套统计包的环境！

---

### 4. 第三步：彻底攻克“中文带来的四大额外技术问题”（必存急救宝典！）

中国研究生用 Python 处理中文实验数据时，95% 的报错都源于以下 4 个中文兼容性问题。请把下面这段**“万能中文防乱码表头”**固定放在每个脚本的最开头：

#### 痛点 1：画出的图片上，中文标题和坐标轴全变成了“口口口”方框豆腐块！
* **原因**：`matplotlib` 和 `seaborn` 是外国人写的，默认字体叫 `DejaVu Sans`，字库里根本没有汉字；而且改用中文字体后，负号 `-` 又会变成方框。
* **根治代码（直接复制到每个脚本开头）**：
  ```python
  import matplotlib.pyplot as plt
  import seaborn as sns

  # 先设置 seaborn 主题，再覆盖字体（顺序千万不能反！否则 seaborn 会把字体洗掉）
  sns.set_theme(style="ticks", font_scale=1.1)
  # Windows 优先用 SimHei(黑体) 或 Microsoft YaHei(微软雅黑)，Mac 优先用 Arial Unicode MS / PingFang SC
  plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS', 'DejaVu Sans']
  plt.rcParams['axes.unicode_minus'] = False  # 解决负号 '-' 显示为方块的问题
  ```

#### 痛点 2：用 Excel 保存的 `.csv` 数据表，Python 读取时报 `UnicodeDecodeError: 'utf-8' codec can't decode`！
* **原因**：国内中文版 Windows 的 Excel 默认把 CSV 存成 `GBK (ANSI)` 编码，而 Python 的 `pd.read_csv()` 默认按国际通用的 `UTF-8` 去读，两边密码本对不上就报错；反之，Python 导出的中文三线表用 Excel 直接双击打开变成乱码。
* **根治代码（写进 Trae 提示词里）**：
  ```python
  import pandas as pd

  # 1. 万能读取函数：先试 utf-8，报错就自动切 gbk，永不报错！
  def read_csv_safe(filepath):
      try:
          return pd.read_csv(filepath, encoding='utf-8')
      except UnicodeDecodeError:
          return pd.read_csv(filepath, encoding='gbk')

  # 2. 导出 CSV 给中文 Excel 打开时，编码一定要写 'utf-8-sig'（带BOM签名），双击打开绝不乱码！
  # df.to_csv("data/统计结果表.csv", index=False, encoding='utf-8-sig')
  ```

#### 痛点 3：表格列名写了中文、括号和斜杠（如 `SOD酶活(U/mL)`），写回归公式时报错！
* **原因**：统计模型公式（如 `statsmodels`）遇到中文括号 `（）`、空格、斜杠 `/` 会把它当成数学除号而崩溃。
* **最佳实践规范**：
  * **表格里的列名（Column Names）尽量用纯英文加下划线**（如 `SOD_activity`, `Group`）；
  * **而在最后画图的坐标轴标签（`xlabel`, `ylabel`）和三线表表头里，再写优雅的中文或标准英文单位**！

#### 痛点 4：自己改代码时，不小心在中文输入法下打了中文逗号 `，`、中文引号 `“”` 或中文括号 `（）`！
* **现象**：报 `SyntaxError: invalid character '，' (U+FF0C)`。
* **识别技巧**：在 Trae 编辑器里，中文字体标点会显得特别宽，且底下会立刻出现**红色波浪线**。选中报错行，直接点 Trae 弹出的 **`AI 修复 (Fix with AI)`**，1 秒钟自动帮你改成英文半角标点！

---

## 第一节：试验设计的基本原则与科学逻辑

### 1. 讲原理：随机、重复、区组三原则与单因素轮换法（OFAT）局限
* **随机化（Randomization）——防顺序偏差**：
  比如测定 20 瓶酸奶的质构硬度，如果按顺序先测完 10 瓶对照组、再测 10 瓶添加果胶组，室温升高会导致后测的酸奶变软。**用随机数打乱测定次序**才能抵消时间趋势。
* **重复（Replication）——认准“独立生物学重复”**：
  * **生物学重复（Biological Replicates，真重复）**：独立配制、独立发酵的 6 个发酵罐（$n=6$）；或者独立饲养的 8 只小鼠（$n=8$）。
  * **技术重复（Technical Replicates，伪重复）**：从同一个发酵罐里舀出 1 管发酵液，在酶标仪上点 3 个孔测吸光度。这 3 个孔只能算平均数代表这 1 罐，**绝不能当成 $n=3$ 去做统计检验**！
* **区组化（Blocking）——消除批次/窝别背景噪音**：
  动物实验分 3 批购买大鼠，或者烘箱上中下三层温度略微不均。把每一层（或每一窝）作为一个“区组（Block）”，每个区组内部都均等安排所有处理组，就能在统计时把区组差异剥离出去。
* **单因素轮换法（OFAT）的局限性**：
  固定时间变温度、再固定温度变时间——完全看不到两个因素之间的**交互协同效应（Interaction）**，且实验次数多、效率低。

### 2. 配套合成数据一：`data/ch1_sec1_principles_demo.csv`
* **实验背景**：考察益生菌发酵产酸实验中，3 个不同原料批次（`Batch`: Batch_1, Batch_2, Batch_3）作为区组，对比对照组（`Control`）与新型氮源组（`New_Peptone`）的乳酸产量（`Lactic_Acid`, g/L），同时包含每个发酵罐的 3 次酶标仪技术重复读数（`Tech_Rep1`, `Tech_Rep2`, `Tech_Rep3`）。

### 3. 写代码：Trae 结构化提示词设计与 Python 实现
* **Trae 提示词模板**：
  > **【角色】**：食品试验设计导师。  
  > **【数据输入】**：读取 `data/ch1_sec1_principles_demo.csv`（使用兼容 utf-8 和 gbk 的安全读取方式）。  
  > **【处理逻辑】**：1. 先将同一发酵罐的 3 列技术重复（`Tech_Rep1~3`）求均值，合并为代表真实生物学重复的列 `Lactic_Acid_Mean`；2. 利用随机数生成器演示如何为 18 个样品生成完全随机化的上机检测顺序号 `Random_Test_Order`；3. 绘制带区组（`Batch`）着色对比的散点箱线图，直观展示批次背景波动与真实处理效应。  
  > **【输出要求】**：导出规整后的生物学重复表至 `data/ch1_sec1_bio_replicates.csv`（编码 `utf-8-sig`），保存 600 DPI 图表至 `scripts/ch1_sec1_principles_plot.png`。

---

## 第二节：变量尺度与科研整洁数据（Tidy Data）规范

### 1. 讲原理：三大变量尺度与长宽表重塑魔术
* **三大变量尺度（决定统计方法选型）**：
  1. **定类尺度（Nominal）**：无大小高低之分（如菌种名称、产地、分组名称）；
  2. **定序尺度（Ordinal）**：有等级顺序但间距不等（如感官评分：差=1、中=2、优=3；腐败等级：轻度、中度、重度）；
  3. **连续定量尺度（Continuous）**：有精确测量单位（如多糖浓度 mg/mL、pH值、色差 $L^*$ 值）。
* **Hadley Wickham 的 Tidy Data（整洁数据）三原则**：
  1. **观测成行**：每一行是一个独立的生物学样品；
  2. **变量成列**：每一列是一个独立的变量；
  3. **长宽表自由重塑**：实验记录本常用**宽表（Wide Format，多列指标并排）**，而高级分面绘图与多指标批量统计必须用 `pd.melt()` 转为**长表（Long Format：`[Sample_ID, Group, Biomarker, Value]`）**。

### 2. 配套合成数据二：`data/ch1_sec2_tidy_wide.csv`
* **实验背景**：蓝莓花青素干预小鼠血清的 4 项抗氧化指标原始宽表（`Mouse_ID`, `Group`, `SOD`, `MDA`, `CAT`, `GSH_Px`）。

### 3. 写代码：Trae 结构化提示词设计与 Python 实现
* **Trae 提示词模板**：
  > **【角色】**：科研数据规整专家。  
  > **【数据输入】**：读取宽表 `data/ch1_sec2_tidy_wide.csv`。  
  > **【操作逻辑】**：1. 使用 `pd.melt()` 将宽表重塑为标准 Tidy 长表；2. 再演示如何用 `pd.pivot()` 将长表还原回宽表；3. 基于 Tidy 长表，调用 `seaborn.catplot` 一行代码绘制 4 项抗氧化指标的分面箱线叠加散点图（设置中文字体防乱码）。  
  > **【输出要求】**：导出 Tidy 长表 `data/ch1_sec2_tidy_long.csv`，保存 600 DPI 分面图 `scripts/ch1_sec2_tidy_facet.png`。

---

## 第三节：数据清洗与质量控制（缺失值插补与 $3\sigma$ / IQR 离群值判别）

### 1. 讲原理：缺失值怎么补？离群值怎么抓？
* **缺失值评估与稳健插补**：
  * 缺失率 $<5\%$ 时，推荐使用**同组内的中位数（Median）**填补（比均值更抗极端值干扰）。
* **离群值（Outliers）两大判别标尺**：
  1. **$3\sigma$ 原则**：落在 $[\mu - 3\sigma,\ \mu + 3\sigma]$ 之外。缺点：小样本下，一个极端大值会把自己所在组的均值 $\mu$ 和标准差 $\sigma$ 同时撑大，导致“自己掩护自己”而漏检！
  2. **箱线图四分位距法（IQR 法，推荐！）**：计算 25% 分位数 $Q_1$ 和 75% 分位数 $Q_3$，$\text{IQR} = Q_3 - Q_1$。凡是落在 $[Q_1 - 1.5\text{IQR},\ Q_3 + 1.5\text{IQR}]$ 之外的点即为离群值，小样本下极其灵敏稳健！

### 2. 配套合成数据三：`data/ch1_sec3_cleaning_raw.csv`
* **实验背景**：24 孔细胞板测定多糖抗炎活性指标（`TNF_alpha`, `IL_6`, `NO_production`），含随机漏测空值（`NaN`）与加样气泡引起的极端离群值。

### 3. 写代码：Trae 结构化提示词设计与 Python 实现
* **Trae 提示词模板**：
  > **【角色】**：生物医药数据质控专家。  
  > **【数据输入】**：读取 `data/ch1_sec3_cleaning_raw.csv`。  
  > **【清洗逻辑】**：1. 统计各列缺失值个数，按 `Group` 分组使用组内中位数插补缺失值；2. 编写对比函数，同时用 $3\sigma$ 原则和 $1.5\times\text{IQR}$ 准则筛查离群值，在终端打印两种方法检出的异常样本编号对比表；3. 绘制清洗前后的箱线散点对比图。  
  > **【输出要求】**：导出清洗后数据 `data/ch1_sec3_cleaned.csv` 与离群点诊断表，保存 600 DPI 对比图 `scripts/ch1_sec3_outlier_check.png`。

---

## 四、第一章课后实战作业指引

1. **环境搭建打卡任务**：在个人电脑完成 Miniconda、`stats_env` 环境及 Trae IDE 的安装配置，运行 `module1_part0_env_check.py` 输出带中文标题的测试图。
2. **课后作业数据集**：
   * `data/ch1_hw_sec1_replicates.csv`：发酵罐技术重复与生物学重复规整练习；
   * `data/ch1_hw_sec2_wide_table.csv`：茶叶 6 种儿茶素组分宽表转 Tidy 长表与分面绘图练习；
   * `data/ch1_hw_sec3_dirty_data.csv`：含缺失值与极端值的果蔬农残检测数据清洗与 IQR 诊断练习。
