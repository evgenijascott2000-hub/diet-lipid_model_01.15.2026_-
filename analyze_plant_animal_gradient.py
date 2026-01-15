#!/usr/bin/env python3
"""
植物源 vs 动物源饮食模式深度分析脚本

功能:
1. 数据准备: 计算每位受试者的植物源总重量(g)、动物源总重量(g)及植物源占比(%)。
2. 二分类人群分析:
   - K-Means聚类 (k=2): 探索自然分群。
   - 阈值划分 (>50%): 按照植物源占比是否超过50%划分。
   - 分析这两种分群方式下，6种脂质水平的差异。
3. 梯度分析 (Quantiles):
   - 将人群按植物源占比分为5等分 (Quintiles)。
   - 分析随着植物源占比升高，脂质水平的变化趋势 (P for trend)。
4. 非线性关系分析 (RCS):
   - 使用限制性立方样条 (Restricted Cubic Splines) 探索植物源占比与脂质水平的非线性剂量-反应关系。

输出目录: /public3/data/zhangxuanyu/Li_Jinxia_program/results_plant_animal_gradient
"""

import os
import sys
import warnings

# Safe import for numpy: try importing, and if missing attempt to install it automatically.
try:
    import numpy as np
except Exception:
    try:
        import subprocess
        print("numpy not found; attempting to install via pip...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "numpy"])
        import importlib
        np = importlib.import_module("numpy")
        print("numpy installed and imported successfully.")
    except Exception as e:
        raise ImportError("Failed to import or install numpy. Please install numpy manually (e.g., pip install numpy).") from e

import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from patsy import cr  # 用于RCS

warnings.filterwarnings('ignore')

# ============================================================================
# 配置
# ============================================================================
BASE_DIR = Path("/public3/data/zhangxuanyu/Li_Jinxia_program")
OUTPUT_DIR = BASE_DIR / "results_plant_animal_gradient"
EXCEL_PATH = BASE_DIR / "UK biobank 100118_93_classify.xlsx"

# 数据文件
MERGED_DATA_PATH = BASE_DIR / "results_gpt5codex_7lipid-33cate" / "merged_lipid_diet_dataset_with_prs.csv"
RECALL_24H_PATH = BASE_DIR / "UKBB_24hRecallDiet.csv"
BLOOD_DATE_PATH = BASE_DIR / "UKBB_blooddate_process.csv"

# 脂质变量
LIPID_COLUMNS = [
    'triglycerides_mmol_L_log',
    'total_cholesterol_mmol_L',
    'hdl_cholesterol_mmol_L_log',
    'ldl_cholesterol_mmol_L',
    'apolipoprotein_a_g_L_log',
    'apolipoprotein_b_g_L'
]

LIPID_LABELS = {
    'triglycerides_mmol_L_log': 'log(TG)',
    'total_cholesterol_mmol_L': 'TC',
    'hdl_cholesterol_mmol_L_log': 'log(HDL-C)',
    'ldl_cholesterol_mmol_L': 'LDL-C',
    'apolipoprotein_a_g_L_log': 'log(ApoA)',
    'apolipoprotein_b_g_L': 'ApoB'
}

# 协变量
COV_NUMERIC = ['cov_age_at_assessment', 'cov_body_mass_index', 'cov_fasting_time']
COV_CATEGORICAL = ['cov_sex', 'cov_alcohol_status', 'cov_current_smoking']

# PRS变量
PRS_COLUMNS = [
    'prs_triglycerides', 'prs_total_cholesterol', 'prs_hdl_cholesterol',
    'prs_ldl_cholesterol', 'prs_apolipoprotein_a1', 'prs_apolipoprotein_b'
]

COMPLETENESS_THRESHOLD = 0.7
RANDOM_STATE = 42

# ============================================================================
# 工具函数 (复用原脚本逻辑)
# ============================================================================
def setup_output_dirs():
    """创建输出目录结构"""
    dirs = {
        'data': OUTPUT_DIR / '01_data',
        'two_groups': OUTPUT_DIR / '02_two_groups_analysis',
        'quantiles': OUTPUT_DIR / '03_quantile_analysis',
        'rcs': OUTPUT_DIR / '04_rcs_analysis',
        'linear_ratio': OUTPUT_DIR / '05_ratio_linear_trend',
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs

def normalize_var_name(name: str) -> str:
    name = name.lower()
    replacements = {'/': '_', ' ': '_', '-': '_', ',': '', '(': '', ')': '', "'": '', '"': '', '.': ''}
    for old, new in replacements.items():
        name = name.replace(old, new)
    while '__' in name:
        name = name.replace('__', '_')
    return f"diet_{name.strip('_')}_24h"

def find_matching_column(description: str, columns: List[str]) -> Optional[str]:
    for col in columns:
        if description in col and 'Instance 0' in col:
            return col
    for col in columns:
        if description.lower() in col.lower():
            return col
    return None

def load_classification_from_excel() -> Dict[str, List[str]]:
    print("读取变量分类信息...")
    classification = {'animal': [], 'plant': []}
    df_source = pd.read_excel(EXCEL_PATH, sheet_name='animal source vs. plant source ', header=[0,1])
    classification['animal'] = df_source[('animal source', 'Description')].dropna().tolist()
    classification['plant'] = df_source[('plant source', 'Description')].dropna().tolist()
    return classification

def load_and_prepare_data(classification: Dict) -> Tuple[pd.DataFrame, List[str], Dict[str, str]]:
    print("加载并准备数据...")
    main_data = pd.read_csv(MERGED_DATA_PATH)
    recall_data = pd.read_csv(RECALL_24H_PATH)
    
    # 统一eid
    recall_eid = 'Participant ID' if 'Participant ID' in recall_data.columns else 'eid'
    main_eid = 'eid' if 'eid' in main_data.columns else main_data.columns[0]
    recall_data = recall_data.rename(columns={recall_eid: 'eid'})
    main_data = main_data.rename(columns={main_eid: 'eid'})
    
    # 匹配列
    all_vars = classification['animal'] + classification['plant']
    var_mapping = {}
    for desc in all_vars:
        col = find_matching_column(desc, recall_data.columns)
        if col:
            var_mapping[desc] = col
            
    diet_cols = list(var_mapping.values())
    recall_subset = recall_data[['eid'] + diet_cols].copy()
    
    rename_map = {col: normalize_var_name(desc) for desc, col in var_mapping.items()}
    recall_subset = recall_subset.rename(columns=rename_map)
    desc_to_simple = {desc: normalize_var_name(desc) for desc in var_mapping.keys()}
    
    merged = main_data.merge(recall_subset, on='eid', how='inner')
    
    # 完整度筛选
    diet_cols_renamed = list(rename_map.values())
    if diet_cols_renamed:
        merged['diet_completeness'] = merged[diet_cols_renamed].notna().sum(axis=1) / len(diet_cols_renamed)
    else:
        merged['diet_completeness'] = 0
        
    filtered_data = merged[merged['diet_completeness'] >= COMPLETENESS_THRESHOLD].copy()
    print(f"筛选后样本数: {len(filtered_data)}")
    
    return filtered_data, diet_cols_renamed, desc_to_simple

def calculate_weights(data: pd.DataFrame, classification: Dict, desc_to_simple: Dict) -> pd.DataFrame:
    """计算植物源和动物源总重量"""
    print("计算饮食总重量...")
    df = data.copy()
    
    # 获取列名
    animal_cols = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['animal']]
    plant_cols = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['plant']]
    
    # 过滤存在的列
    animal_cols = [c for c in animal_cols if c in df.columns]
    plant_cols = [c for c in plant_cols if c in df.columns]
    
    # 计算总重 (g)
    df['total_animal_g'] = df[animal_cols].sum(axis=1)
    df['total_plant_g'] = df[plant_cols].sum(axis=1)
    df['total_diet_g'] = df['total_animal_g'] + df['total_plant_g']
    
    # 计算植物占比 (%)
    df['plant_ratio_pct'] = (df['total_plant_g'] / df['total_diet_g'].replace(0, np.nan)) * 100
    
    # 移除无效数据
    df = df.dropna(subset=['plant_ratio_pct'])
    print(f"有效计算样本数: {len(df)}")
    return df

# ============================================================================
# 分析模块
# ============================================================================

def run_ols_model(data, target_col, x_cols, cov_cols, custom_formula_term=None):
    """运行OLS回归模型"""
    # 准备数据
    model_data = data[[target_col] + x_cols + cov_cols].dropna()
    if len(model_data) < 50:
        return None
    
    # 构建公式
    # 处理分类变量
    cat_covs = [c for c in cov_cols if c in COV_CATEGORICAL]
    num_covs = [c for c in cov_cols if c not in COV_CATEGORICAL]
    
    if custom_formula_term:
        x_term = custom_formula_term
    else:
        x_term = ' + '.join(x_cols)
        
    formula = f"{target_col} ~ {x_term}"
    if num_covs:
        formula += f" + {' + '.join(num_covs)}"
    if cat_covs:
        formula += f" + {' + '.join([f'C({c})' for c in cat_covs])}"
        
    try:
        model = smf.ols(formula, data=model_data).fit()
        return model
    except Exception as e:
        print(f"  Model failed for {target_col}: {e}")
        return None

def analyze_two_groups(data: pd.DataFrame, dirs: Dict):
    """二分类分析: K-Means 和 >50%"""
    print("\n" + "="*60)
    print("执行二分类分析 (Two-Group Analysis)...")
    print("="*60)
    
    results = []
    
    # 1. K-Means Clustering (k=2)
    X = data[['total_plant_g', 'total_animal_g']].copy()
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    kmeans = KMeans(n_clusters=2, random_state=RANDOM_STATE, n_init=10)
    data['kmeans_group'] = kmeans.fit_predict(X_scaled)
    
    # 确定哪个组是植物为主 (Plant-dominant)
    # 比较两组的 plant_ratio_pct 均值
    grp0_mean = data[data['kmeans_group']==0]['plant_ratio_pct'].mean()
    grp1_mean = data[data['kmeans_group']==1]['plant_ratio_pct'].mean()
    
    plant_dom_label = 1 if grp1_mean > grp0_mean else 0
    data['is_plant_dominant_kmeans'] = (data['kmeans_group'] == plant_dom_label).astype(int)
    
    print(f"K-Means分组: Plant-dominant组 (Label={plant_dom_label}) 均值={max(grp0_mean, grp1_mean):.1f}%")
    
    # 可视化 K-Means
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=data, x='total_plant_g', y='total_animal_g', 
                    hue='is_plant_dominant_kmeans', palette='viridis', alpha=0.6)
    plt.title('K-Means Clustering (k=2) on Diet Weights')
    plt.xlabel('Total Plant Weight (g)')
    plt.ylabel('Total Animal Weight (g)')
    plt.savefig(dirs['two_groups'] / 'kmeans_scatter.png', dpi=300)
    plt.close()
    
    # 2. Threshold Split (>50%)
    data['is_plant_dominant_50pct'] = (data['plant_ratio_pct'] > 50).astype(int)
    print(f">50%分组: 植物为主组样本数 = {data['is_plant_dominant_50pct'].sum()}")
    
    # 运行回归
    covariates = [c for c in COV_NUMERIC + COV_CATEGORICAL + PRS_COLUMNS if c in data.columns]
    
    for lipid in LIPID_COLUMNS:
        if lipid not in data.columns: continue
        
        # K-Means Model
        model_km = run_ols_model(data, lipid, ['is_plant_dominant_kmeans'], covariates)
        if model_km:
            res = {
                'method': 'K-Means',
                'lipid': lipid,
                'beta': model_km.params.get('is_plant_dominant_kmeans', np.nan),
                'p_value': model_km.pvalues.get('is_plant_dominant_kmeans', np.nan),
                'conf_low': model_km.conf_int().loc['is_plant_dominant_kmeans', 0],
                'conf_high': model_km.conf_int().loc['is_plant_dominant_kmeans', 1]
            }
            results.append(res)
            
        # 50% Threshold Model
        model_50 = run_ols_model(data, lipid, ['is_plant_dominant_50pct'], covariates)
        if model_50:
            res = {
                'method': '>50% Cutoff',
                'lipid': lipid,
                'beta': model_50.params.get('is_plant_dominant_50pct', np.nan),
                'p_value': model_50.pvalues.get('is_plant_dominant_50pct', np.nan),
                'conf_low': model_50.conf_int().loc['is_plant_dominant_50pct', 0],
                'conf_high': model_50.conf_int().loc['is_plant_dominant_50pct', 1]
            }
            results.append(res)
            
    # 保存结果
    res_df = pd.DataFrame(results)
    res_df.to_csv(dirs['two_groups'] / 'two_group_associations.csv', index=False)
    
    # 绘制森林图
    if not res_df.empty:
        plt.figure(figsize=(10, 6))
        sns.pointplot(data=res_df, x='beta', y='lipid', hue='method', join=False, dodge=0.5, capsize=0.2)
        plt.axvline(0, color='grey', linestyle='--')
        plt.title('Association of Plant-Dominant Diet with Lipids (Two-Group Models)')
        plt.xlabel('Beta Coefficient (95% CI)')
        plt.tight_layout()
        plt.savefig(dirs['two_groups'] / 'two_group_forest_plot.png', dpi=300)
        plt.close()

def plot_errorbar(x, y, lower, upper, **kwargs):
    """Wrapper for errorbar plot with CI bounds"""
    # Ensure inputs are numpy arrays
    x = np.array(x)
    y = np.array(y)
    l = np.array(lower)
    u = np.array(upper)
    
    # Calculate error deltas (must be positive)
    yerr_lower = y - l
    yerr_upper = u - y
    
    # Handle potential precision issues or bad data
    yerr_lower = np.maximum(0, yerr_lower)
    yerr_upper = np.maximum(0, yerr_upper)
    
    # Remove marker from kwargs if present to avoid conflict with fmt='o'
    if 'marker' in kwargs:
        del kwargs['marker']
        
    plt.errorbar(x, y, yerr=[yerr_lower, yerr_upper], fmt='o', capsize=5, **kwargs)

def format_p(pval: float) -> str:
    """格式化P值: <0.001 显示为 "< 0.001"，否则保留三位小数"""
    if pd.isna(pval):
        return "NA"
    return "< 0.001" if pval < 0.001 else f"{pval:.3f}"


def compute_dense_interval(x: pd.Series, bins: int = 20) -> Tuple[float, float]:
    """找到样本较密集的区间: 以直方图中>=非零中位数频数的连续区间作为稠密区间"""
    x = x.dropna().values
    if len(x) == 0:
        return (np.nan, np.nan)

    counts, edges = np.histogram(x, bins=bins)
    nonzero = counts[counts > 0]
    if len(nonzero) == 0:
        return (np.nan, np.nan)

    threshold = np.median(nonzero)
    dense_mask = counts >= threshold

    # 找到最长的连续True区间
    best_len, best_start, current_len, current_start = 0, None, 0, None
    for i, flag in enumerate(dense_mask):
        if flag:
            if current_start is None:
                current_start = i
                current_len = 1
            else:
                current_len += 1
        else:
            if current_len > best_len:
                best_len, best_start = current_len, current_start
            current_start, current_len = None, 0
    if current_len > best_len:
        best_len, best_start = current_len, current_start

    if best_start is None:
        # 回退为10-90分位区间
        low, high = np.percentile(x, [10, 90])
    else:
        low = edges[best_start]
        high = edges[best_start + best_len]

    return float(low), float(high)


def analyze_quantiles(data: pd.DataFrame, dirs: Dict):
    """梯度分析: 五分位数 (Quintiles) + 十分位数 (Deciles)"""
    print("\n" + "="*60)
    print("执行梯度分析 (Quintiles & Deciles Analysis)...")
    print("="*60)
    
    def _quantile_block(q: int, label_prefix: str):
        """生成 q 分位分析块，label_prefix 用于文件命名"""
        labels = [f"Q{i}" for i in range(1, q+1)]
        colname = f"plant_q{q}"
        data[colname] = pd.qcut(data['plant_ratio_pct'], q=q, labels=labels)

        # 分布图
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=data, x=colname, y='plant_ratio_pct', palette='Greens')
        plt.title(f'Plant Ratio Distribution by {label_prefix}')
        plt.ylabel('Plant Source Ratio (%)')
        plt.savefig(dirs['quantiles'] / f'{label_prefix.lower()}_distribution.png', dpi=300)
        plt.close()

        results = []
        covariates = [c for c in COV_NUMERIC + COV_CATEGORICAL + PRS_COLUMNS if c in data.columns]

        for lipid in LIPID_COLUMNS:
            if lipid not in data.columns:
                continue

            # 参考组 Q1
            formula_term = f"C({colname}, Treatment(reference='Q1'))"
            model = run_ols_model(data, lipid, [colname], covariates, custom_formula_term=formula_term)

            if model:
                # Q2..Qn 系数
                for i in range(2, q+1):
                    term = f"C({colname}, Treatment(reference='Q1'))[T.Q{i}]"
                    if term in model.params:
                        results.append({
                            'lipid': lipid,
                            'quantile': f'Q{i}',
                            'beta': model.params[term],
                            'p_value': model.pvalues[term],
                            'conf_low': model.conf_int().loc[term, 0],
                            'conf_high': model.conf_int().loc[term, 1]
                        })

                # P for trend：将分位当作连续变量
                data[f'{colname}_ord'] = data[colname].cat.codes + 1
                model_trend = run_ols_model(data, lipid, [f'{colname}_ord'], covariates)
                if model_trend:
                    p_trend = model_trend.pvalues[f'{colname}_ord']
                    results.append({
                        'lipid': lipid,
                        'quantile': 'P_trend',
                        'beta': np.nan,
                        'p_value': p_trend,
                        'conf_low': np.nan,
                        'conf_high': np.nan
                    })

        res_df = pd.DataFrame(results)
        res_df.to_csv(dirs['quantiles'] / f'{label_prefix.lower()}_associations.csv', index=False)

        # 绘制趋势图
        plot_data = res_df[res_df['quantile'] != 'P_trend'].copy()
        if plot_data.empty:
            return

        # 添加 Q1=0 参考线
        q1_rows = []
        for lipid in plot_data['lipid'].unique():
            q1_rows.append({'lipid': lipid, 'quantile': 'Q1', 'beta': 0, 'conf_low': 0, 'conf_high': 0})
        plot_data = pd.concat([pd.DataFrame(q1_rows), plot_data], ignore_index=True)

        # 排序
        plot_data['q_idx'] = plot_data['quantile'].apply(lambda x: int(x[1:]))
        plot_data = plot_data.sort_values(['lipid', 'q_idx'])

        g = sns.FacetGrid(plot_data, col="lipid", col_wrap=3, sharey=False, height=4)
        g.map(plot_errorbar, "quantile", "beta", "conf_low", "conf_high")

        # 添加 P-trend 注解（格式化显示）
        for ax, lipid in zip(g.axes.flat, g.col_names):
            p_trend = res_df[(res_df['lipid'] == lipid) & (res_df['quantile'] == 'P_trend')]['p_value'].values
            if len(p_trend) > 0:
                ax.text(0.05, 0.9, f'P-trend = {format_p(p_trend[0])}', transform=ax.transAxes, fontsize=10)
            ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

        plt.savefig(dirs['quantiles'] / f'{label_prefix.lower()}_trend_plot.png', dpi=300)
        plt.close()

    # 五分位
    _quantile_block(q=5, label_prefix='Quintile')
    # 十分位
    _quantile_block(q=10, label_prefix='Decile')

def analyze_rcs(data: pd.DataFrame, dirs: Dict):
    """非线性分析: Restricted Cubic Splines (RCS)"""
    print("\n" + "="*60)
    print("执行RCS分析 (Non-linear Analysis)...")
    print("="*60)
    
    covariates = [c for c in COV_NUMERIC + COV_CATEGORICAL + PRS_COLUMNS if c in data.columns]
    baseline_pct = 50  # 以50%作为基准，输出差值
    # Q1的含义：五分位划分中最低20%植物占比人群，对应的植物占比中位数
    quintiles_tmp = pd.qcut(data['plant_ratio_pct'], q=5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])
    q1_median_pct = float(data.loc[quintiles_tmp == 'Q1', 'plant_ratio_pct'].median())
    
    # 准备绘图用的预测数据
    pred_x = np.linspace(data['plant_ratio_pct'].min(), data['plant_ratio_pct'].max(), 100)
    pred_df = pd.DataFrame({'plant_ratio_pct': pred_x})
    
    # 填充协变量为均值或众数
    for col in covariates:
        if col in COV_NUMERIC:
            pred_df[col] = data[col].mean()
        elif col in COV_CATEGORICAL:
            pred_df[col] = data[col].mode()[0]
    dense_low, dense_high = compute_dense_interval(data['plant_ratio_pct'])

    # 绝对值图
    fig_abs, axes_abs = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes_abs = axes_abs.flatten()
    # 相对50%基准差值图
    fig_diff50, axes_diff50 = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes_diff50 = axes_diff50.flatten()
    # 相对Q1中位数差值图
    fig_diffq1, axes_diffq1 = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes_diffq1 = axes_diffq1.flatten()
    # 绝对值图 + 密度 (新增输出)
    fig_abs_density, axes_abs_density = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes_abs_density = axes_abs_density.flatten()
    # 稠密区间拟合图
    fig_dense_abs, axes_dense_abs = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes_dense_abs = axes_dense_abs.flatten()

    for idx, lipid in enumerate(LIPID_COLUMNS):
        if lipid not in data.columns:
            continue
        ax_abs = axes_abs[idx]
        ax_diff50 = axes_diff50[idx]
        ax_diffq1 = axes_diffq1[idx]
        ax_abs_density = axes_abs_density[idx]
        ax_dense_abs = axes_dense_abs[idx]
        
        # 构建RCS公式 (3 knots)
        # 使用 patsy 的 cr (Natural Cubic Spline) 或 bs (B-Spline)
        # 这里使用 cr, df=3 (3 degrees of freedom, usually 3 knots)
        formula = f"{lipid} ~ cr(plant_ratio_pct, df=3)"
        
        # 添加协变量
        num_covs = [c for c in covariates if c in COV_NUMERIC]
        cat_covs = [c for c in covariates if c in COV_CATEGORICAL]
        
        if num_covs: formula += f" + {' + '.join(num_covs)}"
        if cat_covs: formula += f" + {' + '.join([f'C({c})' for c in cat_covs])}"
        
        try:
            model = smf.ols(formula, data=data).fit()
            # 预测 (绝对值)
            pred_abs = model.get_prediction(pred_df)
            pred_mean_abs = pred_abs.predicted_mean
            pred_ci_abs = pred_abs.conf_int()

            # 基准行构造函数
            def _baseline_pred(baseline_pct_val: float) -> float:
                base_row = {'plant_ratio_pct': baseline_pct_val}
                for col in covariates:
                    if col in COV_NUMERIC:
                        base_row[col] = data[col].mean()
                    elif col in COV_CATEGORICAL:
                        base_row[col] = data[col].mode()[0]
                base_pred_local = model.get_prediction(pd.DataFrame([base_row]))
                return float(np.array(base_pred_local.predicted_mean)[0])

            # 差值: 相对50%
            base_50_val = _baseline_pred(baseline_pct)
            pred_mean_diff50 = pred_mean_abs - base_50_val
            pred_ci_diff50 = pred_ci_abs - base_50_val

            # 差值: 相对Q1中位数
            base_q1_val = _baseline_pred(q1_median_pct)
            pred_mean_diffq1 = pred_mean_abs - base_q1_val
            pred_ci_diffq1 = pred_ci_abs - base_q1_val

            # 绘图 (绝对值)
            ax_abs.plot(pred_df['plant_ratio_pct'], pred_mean_abs, color='blue', label='Predicted')
            ax_abs.fill_between(pred_df['plant_ratio_pct'], pred_ci_abs[:, 0], pred_ci_abs[:, 1], color='blue', alpha=0.1)
            ax_abs.set_xlabel('Plant Source Ratio (%)')
            ax_abs.set_ylabel('Lipid Level')

            # 差值曲线 (相对50%)
            ax_diff50.plot(pred_df['plant_ratio_pct'], pred_mean_diff50, color='purple', label='Diff vs 50%')
            ax_diff50.fill_between(pred_df['plant_ratio_pct'], pred_ci_diff50[:, 0], pred_ci_diff50[:, 1], color='purple', alpha=0.1)
            ax_diff50.axhline(0, color='gray', linestyle='--')
            ax_diff50.set_xlabel('Plant Source Ratio (%)')
            ax_diff50.set_ylabel('Difference from 50%')

            # 差值曲线 (相对Q1中位数)
            ax_diffq1.plot(pred_df['plant_ratio_pct'], pred_mean_diffq1, color='teal', label='Diff vs Q1 median')
            ax_diffq1.fill_between(pred_df['plant_ratio_pct'], pred_ci_diffq1[:, 0], pred_ci_diffq1[:, 1], color='teal', alpha=0.1)
            ax_diffq1.axhline(0, color='gray', linestyle='--')
            ax_diffq1.set_xlabel('Plant Source Ratio (%)')
            ax_diffq1.set_ylabel('Difference from Q1 median')

            # 简单的线性 vs 非线性比较 (LRT)
            lin_formula = formula.replace("cr(plant_ratio_pct, df=3)", "plant_ratio_pct")
            lin_model = smf.ols(lin_formula, data=data).fit()
            lrt_stat = 2 * (model.llf - lin_model.llf)
            p_nonlin = stats.chi2.sf(lrt_stat, df=2)  # df difference approx 2

            ax_abs.set_title(f"{LIPID_LABELS.get(lipid, lipid)}\nP_nonlin = {format_p(p_nonlin)}")
            ax_diff50.set_title(f"{LIPID_LABELS.get(lipid, lipid)} (Δ vs 50%)\nP_nonlin = {format_p(p_nonlin)}")
            ax_diffq1.set_title(f"{LIPID_LABELS.get(lipid, lipid)} (Δ vs Q1 median {q1_median_pct:.1f}%)\nP_nonlin = {format_p(p_nonlin)}")

            # 含人群密度的图: 主轴为预测值，次轴为密度
            ax_abs_density.plot(pred_df['plant_ratio_pct'], pred_mean_abs, color='blue', label='Predicted')
            ax_abs_density.fill_between(pred_df['plant_ratio_pct'], pred_ci_abs[:, 0], pred_ci_abs[:, 1], color='blue', alpha=0.1)
            ax_abs_density.set_xlabel('Plant Source Ratio (%)')
            ax_abs_density.set_ylabel('Lipid Level')
            ax_abs_density.set_title(f"{LIPID_LABELS.get(lipid, lipid)} with Density\nP_nonlin = {format_p(p_nonlin)}")
            ax_abs_density.grid(alpha=0.2)
            axd = ax_abs_density.twinx()
            try:
                kde = stats.gaussian_kde(data['plant_ratio_pct'].dropna())
                density = kde(pred_x)
                axd.fill_between(pred_x, density, color='orange', alpha=0.2, label='Density')
                axd.set_ylabel('Density')
            except Exception as kde_err:
                axd.text(0.5, 0.5, f"Density failed: {kde_err}", ha='center')
            axd.set_ylim(bottom=0)
            axd.legend(loc='upper right')
            ax_abs_density.legend(loc='upper left')

            # 稠密区间拟合 (仅使用密集区间的数据)
            dense_mask = (data['plant_ratio_pct'] >= dense_low) & (data['plant_ratio_pct'] <= dense_high)
            dense_data = data.loc[dense_mask]
            if len(dense_data) >= 50:
                try:
                    dense_pred_x = np.linspace(dense_low, dense_high, 80)
                    dense_pred_df = pd.DataFrame({'plant_ratio_pct': dense_pred_x})
                    for col in covariates:
                        if col in COV_NUMERIC:
                            dense_pred_df[col] = dense_data[col].mean()
                        elif col in COV_CATEGORICAL:
                            dense_pred_df[col] = dense_data[col].mode()[0]

                    dense_model = smf.ols(formula, data=dense_data).fit()
                    dense_pred = dense_model.get_prediction(dense_pred_df)
                    dense_mean = dense_pred.predicted_mean
                    dense_ci = dense_pred.conf_int()
                    ax_dense_abs.plot(dense_pred_df['plant_ratio_pct'], dense_mean, color='darkgreen', label='Predicted (dense)')
                    ax_dense_abs.fill_between(dense_pred_df['plant_ratio_pct'], dense_ci[:, 0], dense_ci[:, 1], color='darkgreen', alpha=0.15)
                    ax_dense_abs.set_xlabel('Plant Source Ratio (%)')
                    ax_dense_abs.set_ylabel('Lipid Level')
                    ax_dense_abs.set_title(f"{LIPID_LABELS.get(lipid, lipid)} (Dense {dense_low:.1f}-{dense_high:.1f}%)")
                    ax_dense_abs.grid(alpha=0.2)
                    ax_dense_abs.legend(loc='upper left')
                except Exception as dense_err:
                    ax_dense_abs.text(0.5, 0.5, f"Dense fit failed: {dense_err}", ha='center')
            else:
                ax_dense_abs.text(0.5, 0.5, "Dense fit skipped (n<50)", ha='center')
            
        except Exception as e:
            print(f"RCS failed for {lipid}: {e}")
            ax_abs.text(0.5, 0.5, "Fit Failed", ha='center')
            ax_diff50.text(0.5, 0.5, "Fit Failed", ha='center')
            ax_diffq1.text(0.5, 0.5, "Fit Failed", ha='center')
            ax_abs_density.text(0.5, 0.5, "Fit Failed", ha='center')
            ax_dense_abs.text(0.5, 0.5, "Fit Failed", ha='center')
            
    fig_abs.savefig(dirs['rcs'] / 'rcs_plots.png', dpi=300)
    plt.close(fig_abs)

    fig_diff50.savefig(dirs['rcs'] / 'rcs_plots_diff_vs50.png', dpi=300)
    plt.close(fig_diff50)

    fig_diffq1.savefig(dirs['rcs'] / 'rcs_plots_diff_vsQ1.png', dpi=300)
    plt.close(fig_diffq1)

    fig_abs_density.savefig(dirs['rcs'] / 'rcs_plots_with_density.png', dpi=300)
    plt.close(fig_abs_density)

    fig_dense_abs.savefig(dirs['rcs'] / 'rcs_plots_dense_range.png', dpi=300)
    plt.close(fig_dense_abs)

    pd.DataFrame({
        'dense_low_pct': [dense_low],
        'dense_high_pct': [dense_high]
    }).to_csv(dirs['rcs'] / 'dense_range_summary.csv', index=False)


def analyze_linear_ratio(data: pd.DataFrame, dirs: Dict):
    """线性趋势分析 (X轴直接使用 plant_ratio_pct, 新增05输出)"""
    print("\n" + "="*60)
    print("执行线性趋势分析 (Plant Ratio as Continuous)...")
    print("="*60)

    covariates = [c for c in COV_NUMERIC + COV_CATEGORICAL + PRS_COLUMNS if c in data.columns]
    results = []

    for lipid in LIPID_COLUMNS:
        if lipid not in data.columns:
            continue
        model = run_ols_model(data, lipid, ['plant_ratio_pct'], covariates)
        if model:
            results.append({
                'lipid': lipid,
                'beta_per_pct': model.params.get('plant_ratio_pct', np.nan),
                'p_value': model.pvalues.get('plant_ratio_pct', np.nan),
                'conf_low': model.conf_int().loc['plant_ratio_pct', 0] if 'plant_ratio_pct' in model.conf_int().index else np.nan,
                'conf_high': model.conf_int().loc['plant_ratio_pct', 1] if 'plant_ratio_pct' in model.conf_int().index else np.nan,
                'n': int(model.nobs)
            })

    res_df = pd.DataFrame(results)
    res_df.to_csv(dirs['linear_ratio'] / 'ratio_linear_associations.csv', index=False)

    # 可视化: 抽样散点 + 回归线
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes = axes.flatten()
    for idx, lipid in enumerate(LIPID_COLUMNS):
        ax = axes[idx]
        if lipid not in data.columns:
            ax.text(0.5, 0.5, "No data", ha='center')
            continue
        plot_df = data[['plant_ratio_pct', lipid]].dropna()
        if len(plot_df) == 0:
            ax.text(0.5, 0.5, "No data", ha='center')
            continue
        plot_sample = plot_df.sample(n=min(5000, len(plot_df)), random_state=RANDOM_STATE)
        sns.regplot(data=plot_sample, x='plant_ratio_pct', y=lipid, lowess=True,
                    scatter_kws={'alpha': 0.25, 's': 12}, line_kws={'color': 'red'}, ax=ax)
        label = LIPID_LABELS.get(lipid, lipid)
        ax.set_title(f"{label}")
        ax.set_xlabel('Plant Source Ratio (%)')
        ax.set_ylabel('Lipid Level')

        # 注释p值
        p_row = res_df[res_df['lipid'] == lipid]
        if not p_row.empty:
            p_val = p_row['p_value'].values[0]
            beta = p_row['beta_per_pct'].values[0]
            ax.text(0.05, 0.93, f"beta={beta:.4f}\nP={format_p(p_val)}", transform=ax.transAxes, fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))

    fig.savefig(dirs['linear_ratio'] / 'ratio_linear_scatter_reg.png', dpi=300)
    plt.close(fig)

def main():
    dirs = setup_output_dirs()
    
    # 1. 加载数据
    classification = load_classification_from_excel()
    data, diet_cols, desc_to_simple = load_and_prepare_data(classification)
    
    # 2. 计算重量和比例
    data = calculate_weights(data, classification, desc_to_simple)
    data.to_csv(dirs['data'] / 'data_with_weights.csv', index=False)
    
    # 3. 二分类分析
    analyze_two_groups(data, dirs)
    
    # 4. 梯度分析
    analyze_quantiles(data, dirs)

    # 5. 连续比例的线性趋势分析 (新增05输出)
    analyze_linear_ratio(data, dirs)
    
    # 6. RCS分析
    analyze_rcs(data, dirs)
    
    print("\n分析全部完成！")

if __name__ == "__main__":
    main()
