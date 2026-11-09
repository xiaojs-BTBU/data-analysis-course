# -*- coding: utf-8 -*-
"""
================================================================================
《实验设计与数据处理》配套代码 —— 第一章（第零节）
脚本名称：ch1_00_env_and_chinese_check.py
代码目的：零基础科研环境自检、中文 CSV 编码（UTF-8 / GBK）安全读写测试、
          Matplotlib / Seaborn 中文字体“方框豆腐块”与负号乱码根治验证。
适用对象：食品科学与工程、生物工程及相关专业零编程基础硕士研究生
================================================================================
"""

# ------------------------------------------------------------------------------
# 1. 导入核心科学计算与绘图模块（测试环境是否安装成功）
# ------------------------------------------------------------------------------
import os                         # Python 内置操作系统模块，用于处理文件路径和创建文件夹
import sys                        # Python 内置系统模块，用于查看当前 Python 解释器版本和路径
import numpy as np                # 科学计算基础库 NumPy，用于数组、矩阵与数学函数运算
import pandas as pd               # 表格处理核心库 Pandas，用于读取和整理 Excel / CSV 表格
import scipy                      # 科学统计计算库 SciPy，提供正态检验、t检验、曲线拟合等核心算法
import matplotlib                 # 基础绘图库 Matplotlib
import matplotlib.pyplot as plt   # Matplotlib 的绘图子模块 pyplot，简写为 plt
import seaborn as sns             # 高级统计可视化库 Seaborn，基于 Matplotlib 封装，一行画美图

# ------------------------------------------------------------------------------
# 2. 【核心避坑代码】设置全局绘图风格与跨平台中文字体防乱码（可复制到每个脚本开头！）
# ------------------------------------------------------------------------------
# sns.set_theme() 参数详解：
#   - style="ticks" : 设定背景为干净的白底带坐标轴刻度线风格（最符合 SCI 期刊要求，无多余网格线）
#   - font_scale=1.1: 将全局字号统一放大 1.1 倍，使导出的图表在 Word 论文中清晰易读
sns.set_theme(style="ticks", font_scale=1.1)

# 重要提醒：必须在调用 sns.set_theme() 【之后】再设置 plt.rcParams 字体，否则会被 Seaborn 覆盖！
# font.sans-serif 参数：按顺序依次寻找电脑上已安装的无衬线中文字体：
#   - 'SimHei'          : Windows 系统自带“黑体”
#   - 'Microsoft YaHei' : Windows 系统自带“微软雅黑”
#   - 'Noto Sans CJK SC': Linux / 思源黑体简体中文
#   - 'Arial Unicode MS': macOS 系统常用全字符集字体
#   - 'DejaVu Sans'     : 英文备用字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Noto Sans CJK JP', 'Droid Sans Fallback', 'Arial Unicode MS', 'DejaVu Sans']

# axes.unicode_minus 参数：
#   - 设为 False：强制使用标准 ASCII 减号作为负号，彻底解决切换中文字体后负数（如 -2.5）变成“方框”的顽疾！
plt.rcParams['axes.unicode_minus'] = False

# pdf.fonttype 与 ps.fonttype 参数：
#   - 设为 42（TrueType 字体嵌入）：保证导出的 PDF 矢量图用 Adobe Illustrator 打开时文字可直接编辑！
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def read_csv_safe(filepath):
    """
    【万能中文 CSV 安全读取函数】
    功能说明：国内中文版 Excel 导出的 CSV 常为 GBK 编码，而国际通用标准为 UTF-8 编码。
              本函数先尝试用 utf-8 读取，若报编码错误则自动切换为 gbk 读取，永不报错！
    参数说明：
      - filepath (str): 要读取的 CSV 文件路径（如 "data/ch1_sec1_principles_demo.csv"）
    返回值：
      - pd.DataFrame : 读取成功的 Pandas 数据表对象
    """
    try:
        # pd.read_csv() 参数详解：
        #   - filepath        : 文件路径字符串
        #   - encoding='utf-8': 指定采用国际标准 UTF-8 编码解码文本
        return pd.read_csv(filepath, encoding='utf-8')
    except UnicodeDecodeError:
        # 若捕获到编码不匹配异常，自动改用中文国标编码 'gbk' 重新读取
        return pd.read_csv(filepath, encoding='gbk')


def main():
    print("=" * 72)
    print(">>> 《实验设计与数据处理》零基础科研环境自检与中文防乱码测试启动 <<<")
    print("=" * 72)

    # 打印当前运行的 Python 版本与核心统计包版本号
    print(f"[1] 当前 Python 解释器路径 : {sys.executable}")
    print(f"[2] Python 版本号          : {sys.version.split()[0]}")
    print(f"[3] NumPy 版本号           : {np.__version__}")
    print(f"[4] Pandas 版本号          : {pd.__version__}")
    print(f"[5] SciPy 版本号           : {scipy.__version__}")
    print(f"[6] Matplotlib 版本号      : {matplotlib.__version__}")
    print(f"[7] Seaborn 版本号         : {sns.__version__}")

    # 检查可选高阶统计包（学生在本地 stats_env 环境中通过 pip 安装后即可显示）
    for opt_pkg in ['pingouin', 'statsmodels', 'sklearn', 'pyDOE2', 'adjustText']:
        try:
            mod = __import__(opt_pkg)
            ver = getattr(mod, '__version__', '已安装')
            print(f"[✓] 扩展包 {opt_pkg:12s} : 状态正常 (版本: {ver})")
        except ImportError:
            print(f"[!] 扩展包 {opt_pkg:12s} : 当前环境未安装，脚本将自动启用内置 SciPy/NumPy 权威算法兜底")

    # 生成一张同时包含“中文标题、中文图例、特殊符号（℃、μg/mL）与负数坐标（-3.0 ~ +3.0）”的验证图
    os.makedirs("scripts", exist_ok=True)
    x = np.linspace(-3.0, 3.0, 200)
    y_norm = (1.0 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x**2)

    # plt.subplots() 参数详解：
    #   - figsize=(7.0, 4.2): 画布宽度 7.0 英寸，高度 4.2 英寸
    #   - dpi=150           : 屏幕预览分辨率
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=150)

    # ax.plot() 参数详解：
    #   - x, y_norm         : 横轴与纵轴数据数组
    #   - color='#1F78B4'   : 线条颜色（十六进制深海蓝）
    #   - linewidth=2.5     : 线条粗细（2.5 磅）
    #   - label=...         : 图例文字说明
    ax.plot(x, y_norm, color='#1F78B4', linewidth=2.5, label='标准正态分布曲线 $N(0, 1)$')

    # ax.fill_between() 参数详解：在曲线下方填充半透明颜色，直观展示 95% 置信区间 [-1.96, +1.96]
    mask = (x >= -1.96) & (x <= 1.96)
    ax.fill_between(x[mask], y_norm[mask], color='#A6CEE3', alpha=0.5, label='95% 置信接受域 ($[-1.96, +1.96]$)')

    # 标注两侧 2.5% 拒绝域（p < 0.05 显著性判别区）
    mask_left = x < -1.96
    mask_right = x > 1.96
    ax.fill_between(x[mask_left], y_norm[mask_left], color='#E31A1C', alpha=0.6, label='双侧拒绝域 ($\\alpha = 0.05$)')
    ax.fill_between(x[mask_right], y_norm[mask_right], color='#E31A1C', alpha=0.6)

    ax.set_title("科研环境自检：中文字体、特殊单位 (℃, μg/mL) 与负号 (-1.96) 显示测试", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("标准化 Z 统计量分值 (含负号测试: -3.0 至 +3.0)", fontsize=10.5)
    ax.set_ylabel("概率密度 (Probability Density)", fontsize=10.5)
    ax.legend(loc='upper right', frameon=True, fontsize=9.5)
    sns.despine(top=True, right=True)  # 移除上方和右侧多余边框线，符合 SCI 极简美学
    plt.tight_layout()                 # 自动紧凑调整边距，防止文字被裁切

    out_fig = "scripts/ch1_sec0_env_test_plot.png"
    # plt.savefig() 参数详解：
    #   - out_fig            : 导出图片路径与文件名
    #   - dpi=600            : 导出分辨率设为 600 DPI（国际主流 SCI 期刊要求的出版级清晰度）
    #   - bbox_inches='tight': 自动裁剪四周多余白边
    plt.savefig(out_fig, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"\n[成功] 中文字体与负号测试图已导出至: {out_fig} (600 DPI)")


if __name__ == '__main__':
    main()
