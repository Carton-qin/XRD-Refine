# XRD Auto-Refine: 自动化 X 射线衍射结构精修与多位点阳离子占位分析系统 (v2.0)

[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-blue.svg)](https://python.org)
[![Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-green.svg)](https://python.org)
[![Portability](https://img.shields.io/badge/Dependencies-Pure%20Python-orange.svg)](https://pypi.org)

`xrd-auto-refine` 是专为**催化分子筛、功能多孔材料、矿物晶体转化及无机粉末衍射**设计的全自动化 X 射线粉末衍射 (XRD) Rietveld 结构精修系统。

针对科研实战中不同样品“换样即崩、初修拟合差、多轮调试、缺乏实验前置沟通”的核心痛点，v2.0 全面重构了精修架构，实现**前置意图与工艺自适应问卷、全谱仪器自诊断、非对称物理峰形、多位点阳离子闭环解析以及出图自检质量门禁**。

---

## 目录结构

```text
xrd-auto-refine/
├── SKILL.md                 # Antigravity Agent 核心技能定义规范 (含 5 问智能对齐表)
├── README.md                # 技能详细使用与跨平台指南
├── requirements.txt         # 纯 Python 环境依赖清单
├── database/                # 预置权威晶体学 CIF 数据库
│   ├── minerals/            # 矿物相 (Kaolinite, Quartz, Calcite, Anatase 等)
│   ├── zeolites/            # 沸石相 (SSZ13_CHA, Faujasite_FAU, Zeolite_A_LTA 等)
│   └── byproducts/          # 副产物与碱相 (LiOH, KOH·H2O, Cs2CO3·3H2O 等)
└── scripts/                 # 高性能纯 Python 物理精修引擎
    ├── engine.py            # 矢量化多相 Rietveld 引擎 (含 Split Pseudo-Voigt 与真实晶格优化)
    ├── cation_siting.py     # 沸石孔道与笼外阳离子多位点落位求解器 (SOF 与电中性反演)
    ├── search.py            # 智能寻峰与 Hanawalt / FOM 交叉相关匹配
    ├── export_origin.py     # OriginLab 专用对齐数据与分层布拉格标记导出器
    └── pipeline.py          # 端到端命令行与自动化调度管线 (含仪器元数据解析与自检)
```

---

## 核心升级特性 (v2.0)

1. **智能前置意图与工艺问卷 (Intelligent Intake Protocol)**:
   - 自动检测并主动对齐：材料体系（如分子筛骨架）、目标补偿阳离子（如 $\text{Na}^+/\text{K}^+$）、前处理状态（水洗高温脱模态/原粉态/离子交换态）、精修核心目标（占位/晶格/定量/作图）及测试模式。
   - 智能继承：若用户提问已说明，自动继承跳过，绝不重复盘问。
2. **纯物理自洽引擎 (Zero Fudge Factors)**:
   - 采用积分强度理论投影求解初始比例因子 $s_{\mathrm{init}}$，第 0 轮彻底消除比例失调。
   - 采用晶胞常数六方/立方公式动态重算布拉格位置，严禁粗暴的全谱线性移位。
   - 默认严格各向同性随机粉末（March-Dollase $r=1.0$），严禁随意加设虚假织构。
3. **低角轴向发散非对称峰形 (Split Pseudo-Voigt)**:
   - 针对低角（$2\theta < 20^\circ$）轴向发散导致的严重非对称拖尾，引入左右独立半高宽与高斯/洛伦兹混合函数，将残差曲线平直压制在零基线附近（中位数残差 $< 20$ counts）。
4. **沸石孔道阳离子多位点落位求解器**:
   - 内置 CHA (SSZ-13)、FAU (NaX/NaY)、LTA (分子筛 A) 等典型骨架的标准 Wyckoff 位点拓扑（如 CHA 的 8-MR 18h、6-MR 6c、超笼 6c、6-MR 外侧 6c）。
   - 联立全谱强度差分与骨架电中性条件（$[\text{AlO}_4]^- = \sum q_i n_i$），自动计算各 Wyckoff 位点 SOF、单胞离子数、元素内百分比、$\text{Na}/\text{K}$ 摩尔比、骨架 $\text{Si}/\text{Al}$ 比与干基氧化物比例。
5. **质量门禁闭环自检 (Quality Gate)**:
   - 自动验证 $R_{\mathrm{Bragg}} < 5.0\%$，验证 Diff 中位数残差 $< 1\%$ 峰高，验证阳离子总量在物理合理区间。

---

## 快速使用

```bash
# 执行端到端精修与阳离子落位
arch -x86_64 /tmp/xrd_env/bin/python3 scripts/pipeline.py \
  --data "sample.uxd" \
  --sample "2-NaK-SSZ-13" \
  --framework "CHA" \
  --cations "Na,K" \
  --state "calcined" \
  --goal "cation_siting" \
  --output "./results"
```

### 交付物包 (5-in-1 Package)
1. `*_Rietveld_Origin.csv`: OriginLab 直接绘图数据（包含观测值、拟合值、残差、背景与单相曲线）。
2. `*_Bragg_Ticks.txt`: 各相垂直对齐的 Bragg 衍射标记。
3. `*_Rietveld_Refinement.png`: 300 DPI 高清双层图件（含完整精修与结构信息卡）。
4. `*_Rietveld_Refinement.pdf`: 矢量 PDF 图件。
5. `*_Refinement_Summary.xlsx`: 4 工作表汇总表（指标、晶胞参数、阳离子多位点落位分布、衍射峰表）。
6. `*_Academic_Draft.md`: SCI 英文论文方法与讨论初稿。
