"""
XRD Auto-Refine: OriginLab Ready Data & Bragg Tick Exporter
Exports publication-standard CSV files and pre-spaced Bragg tick tables for OriginLab.
"""

import os
import pandas as pd
import numpy as np

def export_origin_package(sample_name, two_theta, refinement_result, output_dir, tick_spacing_ratio=0.06):
    """
    Export:
    1. Full profile CSV (2Theta, Y_obs, Y_calc, Y_diff, Y_bg, Y_<phase>)
    2. Staggered Bragg tick text files for each phase (non-overlapping!)
    3. OriginLab step-by-step plotting guide.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    y_obs = refinement_result.get('y_obs', refinement_result['y_calc'] + refinement_result['y_diff'])
    y_calc = refinement_result['y_calc']
    y_diff = refinement_result['y_diff']
    y_bg   = refinement_result['y_bg']
    phase_curves = refinement_result.get('phase_curves', {})
    
    # 1. Build Master CSV
    data_dict = {
        '2Theta (deg)': np.round(two_theta, 4),
        'Y_obs': np.round(y_obs, 2),
        'Y_calc': np.round(y_calc, 2),
        'Y_diff': np.round(y_diff, 2),
        'Y_bg': np.round(y_bg, 2)
    }
    
    for pname, pcurve in phase_curves.items():
        data_dict[f'Y_{pname}'] = np.round(pcurve, 2)
        
    df_profile = pd.DataFrame(data_dict)
    csv_filename = f"{sample_name}_Rietveld_Origin.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    df_profile.to_csv(csv_path, index=False)
    
    # 2. Build Staggered Bragg Tick Files
    max_intensity = float(np.max(y_obs))
    tick_base = max_intensity * 0.90
    spacing = max_intensity * tick_spacing_ratio
    
    refl_dict = refinement_result.get('refl_dict', {})
    tick_files = []
    
    for idx, (pname, refs) in enumerate(refl_dict.items()):
        if not refs:
            continue
        tick_y = tick_base + idx * spacing
        tick_tt = [r['2theta'] for r in refs if two_theta[0] <= r['2theta'] <= two_theta[-1]]
        
        tick_df = pd.DataFrame({
            '2Theta (deg)': np.round(tick_tt, 4),
            'Tick_Level': np.round([tick_y] * len(tick_tt), 2)
        })
        
        tick_fn = f"{sample_name}_{pname}_Bragg_Ticks.txt"
        tick_fp = os.path.join(output_dir, tick_fn)
        tick_df.to_csv(tick_fp, sep='\t', index=False)
        tick_files.append((pname, tick_fp))
        
    # 3. OriginLab Step-by-Step Guide
    guide_content = f"""# OriginLab 精修作图操作指南 ({sample_name})

本文件夹包含已预处理、对齐好的 Origin 专用数据包：
1. **主曲线数据**: `{csv_filename}`
   - 列 A: `2Theta (deg)` (X 轴)
   - 列 B: `Y_obs` (实测数据，黑色散点/圆圈)
   - 列 C: `Y_calc` (拟合全谱，红色实线，线宽 1.5 pt)
   - 列 D: `Y_diff` (残差曲线 $Y_{{obs}} - Y_{{calc}}$，蓝色实线，线宽 1.0 pt)
   - 列 E: `Y_bg` (多项式+无定形背景基线，虚线)
   - 列 F+: 各单相独立计算衍射峰 (`Y_<PhaseName>`)

2. **布拉格衍射位置标记 (Bragg Ticks)**:
"""
    for pname, tf in tick_files:
        guide_content += f"   - `{os.path.basename(tf)}`: {pname} 的特征衍射线位置（已预设分层基线，不会发生重叠）\n"
        
    guide_content += """
## Origin 绘图推荐步骤（2 分钟出图）：
1. **导入数据**:
   - 打开 Origin，直接将主 CSV 文件拖入工作表中。
2. **绘制主谱图与残差上下分层图**:
   - 选中工作表，在菜单栏点击 **Plot -> Multi-Panel/Axis -> Stack...** (或 **Vertical 2 Panels**)。
   - **上层面板 (Panel 1)**:
     - X: `2Theta`
     - Y1: `Y_obs` -> 设为 Scatter (黑点或空心圆，大小 4 pt)
     - Y2: `Y_calc` -> 设为 Line (红色实线，1.5 pt)
     - 将各相的 `Bragg_Ticks.txt` 导入并添加为散点图，符号设为垂线 `|` (Vertical Drop Line)。
   - **下层面板 (Panel 2)**:
     - X: `2Theta`
     - Y: `Y_diff` -> 设为 Line (蓝色实线，1.0 pt)
     - 添加一条 $Y=0$ 的虚线水平基准线。
3. **坐标轴美化**:
   - X 轴范围: 设为实测范围 (如 5° ~ 50°)，主刻度间隔 5° 或 10°。
   - 上面板隐藏 X 轴刻度标签，下层面板保留 X 轴标题 `2θ (degrees)`。
   - 上面板 Y 轴标题: `Intensity (a.u.)` 或 `Intensity (counts)`。
   - 下面板 Y 轴标题: `Diff`。
4. **添加指标文本框**:
   - 在右上角空白处添加文本框注明精修收敛指标与相含量：
     - $R_{wp}$, $χ^2$
     - 各相定量含量 (wt.%)
     - 结晶度 DOC (%)
"""
    guide_path = os.path.join(output_dir, "OriginLab_Plotting_Guide.md")
    with open(guide_path, 'w', encoding='utf-8') as f:
        f.write(guide_content)
        
    return {
        'csv_path': csv_path,
        'tick_files': [tf[1] for tf in tick_files],
        'guide_path': guide_path
    }
