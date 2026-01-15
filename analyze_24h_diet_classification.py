#!/usr/bin/env python3
"""
24h饮食变量分类分析脚本

功能:
1. 根据Excel分类表将93个24h饮食变量分为:
   - 食物来源: Animal Source vs Plant Source
   - 加工程度: Extra Processed vs Processed
2. 基于分类后的饮食模式进行聚类分析
3. 分析不同人群的饮食模式占比
4. 对每个人群进行脂质关联分析

输出目录: /public3/data/zhangxuanyu/Li_Jinxia_program/results_24h_93_classify
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, calinski_harabasz_score
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings('ignore')

# ============================================================================
# 配置
# ============================================================================
BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "results_24h_93_classify"
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
    'triglycerides_mmol_L_log': 'log(Triglycerides)',
    'total_cholesterol_mmol_L': 'Total Cholesterol',
    'hdl_cholesterol_mmol_L_log': 'log(HDL Cholesterol)',
    'ldl_cholesterol_mmol_L': 'LDL Cholesterol',
    'apolipoprotein_a_g_L_log': 'log(Apolipoprotein A)',
    'apolipoprotein_b_g_L': 'Apolipoprotein B'
}

# 协变量
COV_NUMERIC = ['cov_age_at_assessment', 'cov_body_mass_index', 'cov_fasting_time']
COV_CATEGORICAL = ['cov_sex', 'cov_alcohol_status', 'cov_current_smoking']

# PRS变量
PRS_COLUMNS = [
    'prs_triglycerides', 'prs_total_cholesterol', 'prs_hdl_cholesterol',
    'prs_ldl_cholesterol', 'prs_apolipoprotein_a1', 'prs_apolipoprotein_b'
]

# 时间匹配参数
TIME_WINDOW_HOURS = 36
COMPLETENESS_THRESHOLD = 0.7

# 聚类参数
N_CLUSTERS_RANGE = range(2, 7)
RANDOM_STATE = 42

# 聚类配色（与展示图保持一致）
CLUSTER_COLORS = {
    0: "#9BBBE1",  # Cluster 0: 淡蓝
    1: "#FEA3A2",  # Cluster 1: 淡粉
    2: "#8E8BFE",  # Cluster 2: 淡紫
}


def get_cluster_color(cluster_id: int) -> str:
    """获取聚类颜色，超出预设时回退到Seaborn调色板。"""
    fallback = sns.color_palette("Set2", n_colors=10)
    return CLUSTER_COLORS.get(cluster_id, fallback[cluster_id % len(fallback)])


# ============================================================================
# 工具函数
# ============================================================================
def setup_output_dirs():
    """创建输出目录结构"""
    dirs = {
        'data_prep': OUTPUT_DIR / '01_data_preparation',
        'classification': OUTPUT_DIR / '02_variable_classification',
        'clustering': OUTPUT_DIR / '03_clustering_analysis',
        'diet_patterns': OUTPUT_DIR / '04_diet_pattern_profiles',
        'lipid_analysis': OUTPUT_DIR / '05_lipid_association',
        'summary': OUTPUT_DIR / '06_summary_reports'
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def normalize_var_name(name: str) -> str:
    """将变量描述转换为标准变量名格式"""
    name = name.lower()
    # 替换特殊字符
    replacements = {
        '/': '_',
        ' ': '_',
        '-': '_',
        ',': '',
        '(': '',
        ')': '',
        "'": '',
        '"': '',
        '.': '',
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    # 移除多余下划线
    while '__' in name:
        name = name.replace('__', '_')
    name = name.strip('_')
    return f"diet_{name}_24h"


def find_matching_column(description: str, columns: List[str]) -> Optional[str]:
    """在列名列表中查找与描述匹配的列(支持Instance格式)"""
    # 精确匹配 (Description | Instance 0 格式)
    for col in columns:
        if description in col and 'Instance 0' in col:
            return col
    # 部分匹配
    for col in columns:
        if description.lower() in col.lower():
            return col
    return None


def load_classification_from_excel() -> Dict[str, Dict[str, List[str]]]:
    """从Excel读取变量分类"""
    print("\n" + "="*60)
    print("读取变量分类信息...")
    print("="*60)
    
    classification = {
        'food_source': {'animal': [], 'plant': []},
        'processing': {'extra_processed': [], 'processed': []}
    }
    
    # 读取Animal vs Plant Source
    df_source = pd.read_excel(EXCEL_PATH, sheet_name='animal source vs. plant source ', header=[0,1])
    animal_vars = df_source[('animal source', 'Description')].dropna().tolist()
    plant_vars = df_source[('plant source', 'Description')].dropna().tolist()
    
    # 存储原始描述名称(用于匹配UKBB数据中的Instance列)
    classification['food_source']['animal'] = animal_vars
    classification['food_source']['plant'] = plant_vars
    
    # 读取Extra Processed vs Processed
    df_process = pd.read_excel(EXCEL_PATH, sheet_name='extra processed vs. processed', header=[0,1])
    extra_vars = df_process[('extra processed', 'Description')].dropna().tolist()
    processed_vars = df_process[('processed', 'Description')].dropna().tolist()
    
    classification['processing']['extra_processed'] = extra_vars
    classification['processing']['processed'] = processed_vars
    
    print(f"\n食物来源分类:")
    print(f"  Animal Source: {len(classification['food_source']['animal'])} 个变量")
    print(f"  Plant Source: {len(classification['food_source']['plant'])} 个变量")
    
    print(f"\n加工程度分类:")
    print(f"  Extra Processed: {len(classification['processing']['extra_processed'])} 个变量")
    print(f"  Processed: {len(classification['processing']['processed'])} 个变量")
    
    return classification


def load_and_prepare_data(classification: Dict) -> Tuple[pd.DataFrame, List[str], Dict[str, str]]:
    """加载并准备数据，进行时间匹配"""
    print("\n" + "="*60)
    print("加载并准备数据...")
    print("="*60)
    
    # 加载主数据 (包含脂质、协变量、PRS)
    print(f"加载主数据: {MERGED_DATA_PATH}")
    main_data = pd.read_csv(MERGED_DATA_PATH)
    print(f"  样本数: {len(main_data)}")
    
    # 加载24h饮食数据
    print(f"加载24h饮食数据: {RECALL_24H_PATH}")
    recall_data = pd.read_csv(RECALL_24H_PATH)
    print(f"  样本数: {len(recall_data)}")
    
    # 加载血液采集日期
    print(f"加载血液日期数据: {BLOOD_DATE_PATH}")
    blood_data = pd.read_csv(BLOOD_DATE_PATH)
    print(f"  样本数: {len(blood_data)}")
    
    # 确定eid列名
    recall_eid_col = 'Participant ID' if 'Participant ID' in recall_data.columns else 'eid'
    main_eid_col = 'eid' if 'eid' in main_data.columns else main_data.columns[0]
    
    # 统一eid列名
    if recall_eid_col != 'eid':
        recall_data = recall_data.rename(columns={recall_eid_col: 'eid'})
    if main_eid_col != 'eid':
        main_data = main_data.rename(columns={main_eid_col: 'eid'})
    
    # 获取所有分类中的变量描述
    all_classified_descriptions = set()
    for cat_type in classification.values():
        for vars_list in cat_type.values():
            all_classified_descriptions.update(vars_list)
    
    print(f"\n分类中的变量描述数: {len(all_classified_descriptions)}")
    
    # 在recall_data中查找匹配的列 (Instance 0)
    var_mapping = {}  # description -> actual column name
    for desc in all_classified_descriptions:
        col = find_matching_column(desc, recall_data.columns)
        if col:
            var_mapping[desc] = col
    
    print(f"成功匹配的变量数: {len(var_mapping)}")
    
    if not var_mapping:
        print("警告: 没有匹配到任何变量！")
        return pd.DataFrame(), [], {}
    
    # 提取需要的列
    diet_cols = list(var_mapping.values())
    recall_subset = recall_data[['eid'] + diet_cols].copy()
    
    # 重命名为简化名称 (便于后续处理)
    rename_map = {col: normalize_var_name(desc) for desc, col in var_mapping.items()}
    recall_subset = recall_subset.rename(columns=rename_map)
    diet_cols_renamed = list(rename_map.values())
    
    # 创建描述到简化名称的映射
    desc_to_simple = {desc: normalize_var_name(desc) for desc in var_mapping.keys()}
    
    # 合并主数据与饮食数据
    merged = main_data.merge(recall_subset, on='eid', how='inner')
    print(f"\n合并后样本数: {len(merged)}")
    
    # 合并血液日期
    if 'blood_sample_date' in blood_data.columns:
        merged = merged.merge(blood_data[['eid', 'blood_sample_date']], on='eid', how='left')
    
    # 计算完整度
    if diet_cols_renamed:
        merged['diet_completeness'] = merged[diet_cols_renamed].notna().sum(axis=1) / len(diet_cols_renamed)
    else:
        merged['diet_completeness'] = 0
    
    # 筛选完整度达标的样本
    valid_mask = merged['diet_completeness'] >= COMPLETENESS_THRESHOLD
    filtered_data = merged[valid_mask].copy()
    
    print(f"\n数据筛选:")
    print(f"  合并后样本数: {len(merged)}")
    print(f"  完整度阈值: {COMPLETENESS_THRESHOLD}")
    print(f"  筛选后样本数: {len(filtered_data)}")
    
    return filtered_data, diet_cols_renamed, desc_to_simple


def compute_classification_scores(data: pd.DataFrame, 
                                   classification: Dict,
                                   diet_cols: List[str],
                                   desc_to_simple: Dict[str, str]) -> pd.DataFrame:
    """计算每个样本的分类得分"""
    print("\n" + "="*60)
    print("计算分类得分...")
    print("="*60)
    
    scores = pd.DataFrame(index=data.index)
    
    # 获取实际存在的变量
    available_vars = set(diet_cols)
    
    # 将分类中的描述转换为简化变量名
    animal_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['food_source']['animal']]
    animal_vars = [v for v in animal_vars if v in available_vars]
    
    plant_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['food_source']['plant']]
    plant_vars = [v for v in plant_vars if v in available_vars]
    
    print(f"  Animal Source变量: {len(animal_vars)}")
    print(f"  Plant Source变量: {len(plant_vars)}")
    
    if animal_vars:
        scores['animal_source_sum'] = data[animal_vars].sum(axis=1)
        scores['animal_source_mean'] = data[animal_vars].mean(axis=1)
    if plant_vars:
        scores['plant_source_sum'] = data[plant_vars].sum(axis=1)
        scores['plant_source_mean'] = data[plant_vars].mean(axis=1)
    
    # 计算比例
    if 'animal_source_sum' in scores.columns and 'plant_source_sum' in scores.columns:
        total = scores['animal_source_sum'] + scores['plant_source_sum']
        scores['animal_ratio'] = scores['animal_source_sum'] / total.replace(0, np.nan)
        scores['plant_ratio'] = scores['plant_source_sum'] / total.replace(0, np.nan)
    
    # 加工程度得分
    extra_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['processing']['extra_processed']]
    extra_vars = [v for v in extra_vars if v in available_vars]
    
    processed_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['processing']['processed']]
    processed_vars = [v for v in processed_vars if v in available_vars]
    
    print(f"  Extra Processed变量: {len(extra_vars)}")
    print(f"  Processed变量: {len(processed_vars)}")
    
    if extra_vars:
        scores['extra_processed_sum'] = data[extra_vars].sum(axis=1)
        scores['extra_processed_mean'] = data[extra_vars].mean(axis=1)
    if processed_vars:
        scores['processed_sum'] = data[processed_vars].sum(axis=1)
        scores['processed_mean'] = data[processed_vars].mean(axis=1)
    
    # 计算比例
    if 'extra_processed_sum' in scores.columns and 'processed_sum' in scores.columns:
        total = scores['extra_processed_sum'] + scores['processed_sum']
        scores['extra_processed_ratio'] = scores['extra_processed_sum'] / total.replace(0, np.nan)
        scores['processed_ratio'] = scores['processed_sum'] / total.replace(0, np.nan)
    
    print(f"计算得分完成，特征数: {len(scores.columns)}")
    print(f"特征列表: {list(scores.columns)}")
    
    return scores


def perform_clustering(scores: pd.DataFrame, dirs: Dict) -> Tuple[pd.Series, Dict]:
    """执行聚类分析"""
    print("\n" + "="*60)
    print("执行聚类分析...")
    print("="*60)
    
    # 准备聚类特征
    cluster_features = ['animal_ratio', 'plant_ratio', 
                        'extra_processed_ratio', 'processed_ratio']
    cluster_features = [f for f in cluster_features if f in scores.columns]
    
    if len(cluster_features) < 2:
        # 如果比例特征不足，使用均值特征
        cluster_features = [c for c in scores.columns if 'mean' in c or 'ratio' in c]
    
    X = scores[cluster_features].dropna()
    print(f"聚类特征: {cluster_features}")
    print(f"有效样本数: {len(X)}")
    
    if len(X) < 100:
        print("警告: 样本数过少，无法进行可靠的聚类分析")
        return pd.Series(dtype=int), {}
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 评估不同聚类数
    results = []
    for n_clusters in N_CLUSTERS_RANGE:
        # KMeans
        kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init=10)
        labels = kmeans.fit_predict(X_scaled)
        
        sil_score = silhouette_score(X_scaled, labels)
        ch_score = calinski_harabasz_score(X_scaled, labels)
        
        results.append({
            'n_clusters': n_clusters,
            'method': 'KMeans',
            'silhouette': sil_score,
            'calinski_harabasz': ch_score,
            'inertia': kmeans.inertia_
        })
        
        # GMM
        gmm = GaussianMixture(n_components=n_clusters, random_state=RANDOM_STATE)
        labels = gmm.fit_predict(X_scaled)
        
        sil_score = silhouette_score(X_scaled, labels)
        ch_score = calinski_harabasz_score(X_scaled, labels)
        
        results.append({
            'n_clusters': n_clusters,
            'method': 'GMM',
            'silhouette': sil_score,
            'calinski_harabasz': ch_score,
            'bic': gmm.bic(X_scaled)
        })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(dirs['clustering'] / 'clustering_evaluation.csv', index=False)
    
    # 选择最优聚类数（基于轮廓系数）
    kmeans_results = results_df[results_df['method'] == 'KMeans']
    best_n = kmeans_results.loc[kmeans_results['silhouette'].idxmax(), 'n_clusters']
    best_n = int(best_n)
    
    print(f"\n最优聚类数: {best_n} (基于轮廓系数)")
    
    # 执行最终聚类
    final_kmeans = KMeans(n_clusters=best_n, random_state=RANDOM_STATE, n_init=10)
    cluster_labels = pd.Series(final_kmeans.fit_predict(X_scaled), index=X.index)
    
    # 绘制评估图
    metrics_to_plot = []
    if 'silhouette' in results_df.columns:
        metrics_to_plot.append(('silhouette', 'Silhouette Score', 'higher is better'))
    if 'calinski_harabasz' in results_df.columns:
        metrics_to_plot.append(('calinski_harabasz', 'Calinski-Harabasz Index', 'higher is better'))
    if 'inertia' in results_df.columns:
        metrics_to_plot.append(('inertia', 'Inertia (KMeans)', 'lower is better'))
    if 'bic' in results_df.columns:
        metrics_to_plot.append(('bic', 'BIC (GMM)', 'lower is better'))

    n_metrics = len(metrics_to_plot)
    n_cols = 2 if n_metrics > 1 else 1
    n_rows = int(np.ceil(n_metrics / n_cols)) if n_metrics else 1
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
    axes = np.array(axes).reshape(-1)

    for ax_idx, (metric, title, note) in enumerate(metrics_to_plot):
        ax = axes[ax_idx]
        for method in results_df['method'].unique():
            subset = results_df[(results_df['method'] == method) & results_df[metric].notna()]
            if subset.empty:
                continue
            ax.plot(subset['n_clusters'], subset[metric], marker='o', label=method)
        ax.set_xlabel('Number of Clusters')
        ax.set_ylabel(title)
        suffix = f" ({note})" if note else ''
        ax.set_title(f"{title}{suffix}")
        ax.axvline(best_n, color='red', linestyle='--', alpha=0.5)
        ax.legend()

    # 隐藏未使用的子图
    for ax in axes[n_metrics:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(dirs['clustering'] / 'clustering_evaluation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 聚类结果可视化 - 食物来源 vs 加工程度
    plot_df = X.copy()
    ratio_cols = [c for c in cluster_features if 'ratio' in c]
    plot_df[ratio_cols] = plot_df[ratio_cols] * 100

    if len(cluster_features) >= 4:
        fig, ax = plt.subplots(figsize=(10, 8))
        colors = [get_cluster_color(i) for i in range(best_n)]
        for i in range(best_n):
            mask = cluster_labels == i
            ax.scatter(plot_df.loc[mask, cluster_features[0]], 
                      plot_df.loc[mask, cluster_features[2]], 
                      c=[colors[i]], label=f'Cluster {i}', alpha=0.6, s=20)
        ax.set_xlabel('Animal Source Ratio (%)')
        ax.set_ylabel('Extra Processed Ratio (%)')
        ax.set_title(f'Clustering: Food Source vs Processing Level (n={best_n} clusters)')
        ax.set_xticks(np.arange(0, 101, 20))
        ax.set_yticks(np.arange(0, 101, 20))
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.legend(title='Cluster', loc='best')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(dirs['clustering'] / 'cluster_scatter_food_source.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 聚类结果可视化 - 植物来源 vs 最小加工
    if len(cluster_features) >= 4:
        fig, ax = plt.subplots(figsize=(10, 8))
        colors = [get_cluster_color(i) for i in range(best_n)]
        for i in range(best_n):
            mask = cluster_labels == i
            ax.scatter(plot_df.loc[mask, cluster_features[1]], 
                      plot_df.loc[mask, cluster_features[3]], 
                      c=[colors[i]], label=f'Cluster {i}', alpha=0.6, s=20)
        ax.set_xlabel('Plant Source Ratio (%)')
        ax.set_ylabel('Processed (Minimally) Ratio (%)')
        ax.set_title(f'Clustering: Plant Source vs Processing Level (n={best_n} clusters)')
        ax.set_xticks(np.arange(0, 101, 20))
        ax.set_yticks(np.arange(0, 101, 20))
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.legend(title='Cluster', loc='best')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(dirs['clustering'] / 'cluster_scatter_processing.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    cluster_info = {
        'n_clusters': best_n,
        'features': cluster_features,
        'scaler': scaler,
        'model': final_kmeans
    }
    
    return cluster_labels, cluster_info


def analyze_diet_patterns(data: pd.DataFrame, 
                          cluster_labels: pd.Series,
                          classification: Dict,
                          diet_cols: List[str],
                          dirs: Dict,
                          desc_to_simple: Dict[str, str]):
    """分析每个聚类的饮食模式特征"""
    print("\n" + "="*60)
    print("分析各聚类的饮食模式...")
    print("="*60)
    
    # 合并聚类标签
    data_with_clusters = data.loc[cluster_labels.index].copy()
    data_with_clusters['cluster'] = cluster_labels
    
    n_clusters = cluster_labels.nunique()
    
    # 计算每个聚类的饮食变量均值
    cluster_profiles = []
    for cluster_id in range(n_clusters):
        cluster_data = data_with_clusters[data_with_clusters['cluster'] == cluster_id]
        
        profile = {'cluster': cluster_id, 'n_samples': len(cluster_data)}
        
        # 计算各分类的占比
        available_vars = set(diet_cols)
        
        # Animal vs Plant
        animal_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['food_source']['animal']]
        animal_vars = [v for v in animal_vars if v in available_vars]
        
        plant_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['food_source']['plant']]
        plant_vars = [v for v in plant_vars if v in available_vars]
        
        if animal_vars:
            animal_total = cluster_data[animal_vars].sum().sum()
        else:
            animal_total = 0
        if plant_vars:
            plant_total = cluster_data[plant_vars].sum().sum()
        else:
            plant_total = 0
            
        total_source = animal_total + plant_total
        profile['animal_pct'] = (animal_total / total_source * 100) if total_source > 0 else 0
        profile['plant_pct'] = (plant_total / total_source * 100) if total_source > 0 else 0
        
        # Extra Processed vs Processed
        extra_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['processing']['extra_processed']]
        extra_vars = [v for v in extra_vars if v in available_vars]
        
        proc_vars = [desc_to_simple.get(v, normalize_var_name(v)) for v in classification['processing']['processed']]
        proc_vars = [v for v in proc_vars if v in available_vars]
        
        if extra_vars:
            extra_total = cluster_data[extra_vars].sum().sum()
        else:
            extra_total = 0
        if proc_vars:
            proc_total = cluster_data[proc_vars].sum().sum()
        else:
            proc_total = 0
            
        total_proc = extra_total + proc_total
        profile['extra_processed_pct'] = (extra_total / total_proc * 100) if total_proc > 0 else 0
        profile['processed_pct'] = (proc_total / total_proc * 100) if total_proc > 0 else 0
        
        cluster_profiles.append(profile)
    
    profiles_df = pd.DataFrame(cluster_profiles)
    profiles_df.to_csv(dirs['diet_patterns'] / 'cluster_diet_profiles.csv', index=False)
    
    print("\n各聚类饮食模式特征:")
    print(profiles_df.to_string(index=False))
    
    # 可视化
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 食物来源占比
    x = np.arange(n_clusters)
    width = 0.35
    
    axes[0].bar(x - width/2, profiles_df['animal_pct'], width, label='Animal Source', color='#E74C3C')
    axes[0].bar(x + width/2, profiles_df['plant_pct'], width, label='Plant Source', color='#27AE60')
    axes[0].set_xlabel('Cluster')
    axes[0].set_ylabel('Percentage (%)')
    axes[0].set_title('Food Source Distribution by Cluster')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f'Cluster {i}' for i in range(n_clusters)])
    axes[0].legend()
    
    # 加工程度占比
    axes[1].bar(x - width/2, profiles_df['extra_processed_pct'], width, 
                label='Extra Processed', color='#9B59B6')
    axes[1].bar(x + width/2, profiles_df['processed_pct'], width, 
                label='Processed', color='#3498DB')
    axes[1].set_xlabel('Cluster')
    axes[1].set_ylabel('Percentage (%)')
    axes[1].set_title('Processing Level Distribution by Cluster')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f'Cluster {i}' for i in range(n_clusters)])
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(dirs['diet_patterns'] / 'cluster_diet_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 为每个聚类命名
    cluster_names = {}
    for _, row in profiles_df.iterrows():
        cluster_id = int(row['cluster'])
        
        # 基于占比命名
        source_type = 'Animal-based' if row['animal_pct'] > row['plant_pct'] else 'Plant-based'
        process_type = 'Ultra-processed' if row['extra_processed_pct'] > row['processed_pct'] else 'Minimally-processed'
        
        cluster_names[cluster_id] = f"{source_type}, {process_type}"
    
    print("\n聚类命名:")
    for cid, name in cluster_names.items():
        print(f"  Cluster {cid}: {name}")
    
    return profiles_df, cluster_names


def analyze_lipid_by_cluster(data: pd.DataFrame,
                              cluster_labels: pd.Series,
                              cluster_names: Dict[int, str],
                              dirs: Dict):
    """分析各聚类与脂质的关联"""
    print("\n" + "="*60)
    print("分析各聚类与脂质的关联...")
    print("="*60)
    
    # 合并数据
    analysis_data = data.loc[cluster_labels.index].copy()
    analysis_data['cluster'] = cluster_labels
    analysis_data['cluster_name'] = analysis_data['cluster'].map(cluster_names)
    
    results = []
    
    for lipid in LIPID_COLUMNS:
        if lipid not in analysis_data.columns:
            print(f"  跳过 {lipid}: 列不存在")
            continue
            
        lipid_label = LIPID_LABELS.get(lipid, lipid)
        
        # 获取各聚类的脂质水平
        cluster_stats = analysis_data.groupby('cluster')[lipid].agg(['mean', 'std', 'count'])
        
        # ANOVA检验
        groups = [analysis_data[analysis_data['cluster'] == c][lipid].dropna() 
                  for c in sorted(analysis_data['cluster'].unique())]
        
        if all(len(g) > 10 for g in groups):
            f_stat, p_value = stats.f_oneway(*groups)
        else:
            f_stat, p_value = np.nan, np.nan
        
        # OLS回归（控制协变量）
        available_cov = [c for c in COV_NUMERIC + COV_CATEGORICAL if c in analysis_data.columns]
        available_prs = [c for c in PRS_COLUMNS if c in analysis_data.columns]
        
        model_vars = ['cluster'] + available_cov + available_prs
        subset = analysis_data[[lipid] + model_vars].copy()
        
        # 确保数据类型正确
        subset[lipid] = pd.to_numeric(subset[lipid], errors='coerce')
        for col in available_prs:
            if col in subset.columns:
                subset[col] = pd.to_numeric(subset[col], errors='coerce')
        for col in available_cov:
            if col in subset.columns and col in COV_NUMERIC:
                subset[col] = pd.to_numeric(subset[col], errors='coerce')
        
        # 删除所有包含NaN的行
        subset = subset.dropna()
        
        if len(subset) > 100:
            # 为每个聚类分别计算与脂质的关联（不使用参考组）
            try:
                # 对每个聚类单独进行回归分析
                for cluster_id in sorted(analysis_data['cluster'].unique()):
                    cluster_subset = subset.copy()
                    
                    # 创建该聚类的二元指示变量
                    cluster_subset['is_target_cluster'] = (cluster_subset['cluster'] == cluster_id).astype(int)
                    
                    # 准备协变量（不包括cluster变量）
                    cov_cols = [c for c in available_cov if c in cluster_subset.columns]
                    prs_cols = [c for c in available_prs if c in cluster_subset.columns]
                    
                    # 创建设计矩阵
                    X_vars = ['is_target_cluster'] + cov_cols + prs_cols
                    X = cluster_subset[X_vars].copy()
                    
                    # 对分类协变量进行虚拟编码
                    cat_cov = [c for c in cov_cols if c in COV_CATEGORICAL]
                    if cat_cov:
                        X = pd.get_dummies(X, columns=cat_cov, drop_first=True)
                    
                    # 删除可能产生的NaN和Inf
                    X = X.replace([np.inf, -np.inf], np.nan).dropna(axis=1, how='any')
                    X = X.astype(float)
                    
                    # 对齐y
                    valid_idx = X.index
                    y = cluster_subset.loc[valid_idx, lipid].astype(float)
                    
                    # 删除NaN
                    valid_mask = ~(X.isna().any(axis=1) | y.isna())
                    X = X[valid_mask]
                    y = y[valid_mask]
                    
                    if len(X) < 100:
                        continue
                    
                    X = sm.add_constant(X)
                    
                    model = sm.OLS(y, X).fit()
                    
                    # 提取该聚类的效应
                    if 'is_target_cluster' in model.params.index:
                        results.append({
                            'lipid': lipid,
                            'lipid_label': lipid_label,
                            'cluster': cluster_id,
                            'cluster_name': cluster_names.get(cluster_id, f'Cluster {cluster_id}'),
                            'beta': model.params['is_target_cluster'],
                            'std_error': model.bse['is_target_cluster'],
                            't_value': model.tvalues['is_target_cluster'],
                            'p_value': model.pvalues['is_target_cluster'],
                            'n_obs': len(y),
                            'n_in_cluster': int(X['is_target_cluster'].sum()) if 'is_target_cluster' in X.columns else 0,
                            'r_squared': model.rsquared,
                            'anova_f': f_stat,
                            'anova_p': p_value
                        })
            except Exception as e:
                print(f"  回归分析失败 ({lipid}): {e}")
    
    if results:
        results_df = pd.DataFrame(results)
        
        # FDR校正
        from statsmodels.stats.multitest import multipletests
        _, fdr_pvals, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        results_df['p_fdr'] = fdr_pvals
        results_df['fdr_significant'] = results_df['p_fdr'] < 0.05
        
        results_df.to_csv(dirs['lipid_analysis'] / 'cluster_lipid_associations.csv', index=False)
        
        print("\n显著的聚类-脂质关联 (FDR < 0.05):")
        sig_results = results_df[results_df['fdr_significant']]
        if not sig_results.empty:
            for _, row in sig_results.iterrows():
                direction = "↑" if row['beta'] > 0 else "↓"
                print(f"  {row['lipid_label']}: Cluster {row['cluster']} ({row['cluster_name']}) "
                      f"{direction} β={row['beta']:.4f}, p_FDR={row['p_fdr']:.4f}")

        # 以 cluster 0 为参照的聚类-脂质关联（统一参考组）
        ref_results = []
        cluster_ref = 0
        clusters_sorted = sorted(analysis_data['cluster'].dropna().unique())
        if cluster_ref not in clusters_sorted:
            print("\n提示: 数据中不存在 cluster 0，跳过参考组回归输出。")
        else:
            for lipid in LIPID_COLUMNS:
                if lipid not in analysis_data.columns:
                    continue

                subset = analysis_data[[lipid, 'cluster'] + available_cov + available_prs].copy()
                subset[lipid] = pd.to_numeric(subset[lipid], errors='coerce')
                for col in available_prs:
                    if col in subset.columns:
                        subset[col] = pd.to_numeric(subset[col], errors='coerce')
                for col in available_cov:
                    if col in subset.columns and col in COV_NUMERIC:
                        subset[col] = pd.to_numeric(subset[col], errors='coerce')

                subset = subset.dropna()
                if len(subset) < 100:
                    continue

                # 设置cluster为有序分类，确保0为基准
                cat_order = [cluster_ref] + [c for c in clusters_sorted if c != cluster_ref]
                subset['cluster'] = pd.Categorical(subset['cluster'], categories=cat_order, ordered=True)
                cluster_dummies = pd.get_dummies(subset['cluster'], prefix='cluster', drop_first=True)

                X = pd.concat([cluster_dummies, subset[available_cov + available_prs]], axis=1)
                cat_cov = [c for c in available_cov if c in COV_CATEGORICAL and c in X.columns]
                if cat_cov:
                    X = pd.get_dummies(X, columns=cat_cov, drop_first=True)

                X = X.apply(pd.to_numeric, errors='coerce')
                X = X.replace([np.inf, -np.inf], np.nan).dropna()
                y = subset.loc[X.index, lipid].astype(float)

                # 对齐并过滤缺失
                valid_mask = ~(X.isna().any(axis=1) | y.isna())
                X = X[valid_mask]
                y = y[valid_mask]
                if len(X) < 100:
                    continue

                X = X.astype(float)
                y = y.astype(float)

                X = sm.add_constant(X)
                model = sm.OLS(y, X).fit()

                # 记录每个非基准聚类的系数
                for cname in cluster_dummies.columns:
                    if cname in model.params.index:
                        try:
                            cluster_id = int(cname.split('_')[1])
                        except Exception:
                            cluster_id = cname
                        ref_results.append({
                            'ref_cluster': cluster_ref,
                            'cluster': cluster_id,
                            'lipid': lipid,
                            'lipid_label': LIPID_LABELS.get(lipid, lipid),
                            'beta': model.params[cname],
                            'std_error': model.bse.get(cname, np.nan),
                            't_value': model.tvalues.get(cname, np.nan),
                            'p_value': model.pvalues.get(cname, np.nan),
                            'n_obs': len(y)
                        })

            if ref_results:
                ref_df = pd.DataFrame(ref_results)
                _, fdr_ref, _, _ = multipletests(ref_df['p_value'], method='fdr_bh')
                ref_df['p_fdr'] = fdr_ref
                ref_df['fdr_significant'] = ref_df['p_fdr'] < 0.05
                ref_df.to_csv(dirs['lipid_analysis'] / 'cluster_lipid_associations_ref_cluster0.csv', index=False)

                print("\n以cluster 0为参照的显著关联 (FDR < 0.05):")
                sig_ref = ref_df[ref_df['fdr_significant']]
                if not sig_ref.empty:
                    for _, row in sig_ref.iterrows():
                        direction = "↑" if row['beta'] > 0 else "↓"
                        print(f"  {row['lipid_label']}: Cluster {row['cluster']} 相对 Cluster {row['ref_cluster']} "
                              f"{direction} β={row['beta']:.4f}, p_FDR={row['p_fdr']:.4f}")

                # 绘制参考组可视化
                plot_lipid_cluster_associations_ref(ref_df, cluster_names, dirs, ref_cluster=cluster_ref)
    else:
        print("  无显著关联")
        
        # 可视化
        plot_lipid_cluster_associations(results_df, cluster_names, dirs)
        
        return results_df
    
    return pd.DataFrame()


def plot_lipid_cluster_associations(results_df: pd.DataFrame, 
                                     cluster_names: Dict,
                                     dirs: Dict):
    """绘制脂质-聚类关联图"""
    
    # 热图 - 只显示实际有数据的聚类
    pivot_data = results_df.pivot_table(
        values='beta', 
        index='lipid_label', 
        columns='cluster',  # 使用cluster数字而非cluster_name
        aggfunc='first'
    )
    
    # 重命名列为聚类名称
    pivot_data.columns = [f'Cluster {int(c)}\n{cluster_names.get(int(c), "")}' 
                          for c in pivot_data.columns]
    
    fig, ax = plt.subplots(figsize=(max(12, len(pivot_data.columns)*3), 8))
    sns.heatmap(pivot_data, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                ax=ax, cbar_kws={'label': 'Beta Coefficient'})
    ax.set_title('Lipid Levels by Diet Pattern Cluster\n(Adjusted for Covariates)')
    ax.set_xlabel('Diet Pattern Cluster')
    ax.set_ylabel('Lipid')
    
    plt.tight_layout()
    plt.savefig(dirs['lipid_analysis'] / 'lipid_cluster_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 森林图
    fig, axes = plt.subplots(len(LIPID_LABELS), 1, figsize=(10, 4*len(LIPID_LABELS)))
    if len(LIPID_LABELS) == 1:
        axes = [axes]
    
    for ax, (lipid, label) in zip(axes, LIPID_LABELS.items()):
        subset = results_df[results_df['lipid'] == lipid].sort_values('cluster')
        
        if subset.empty:
            ax.set_visible(False)
            continue
        
        y_pos = np.arange(len(subset))
        ax.errorbar(subset['beta'], y_pos, 
                    xerr=1.96*subset['std_error'],
                    fmt='o', capsize=5, capthick=2,
                    color='steelblue', ecolor='gray')
        ax.axvline(0, color='red', linestyle='--', alpha=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"Cluster {int(c)}" for c in subset['cluster']])
        ax.set_xlabel('Beta Coefficient (95% CI)')
        ax.set_title(label)
    
    plt.tight_layout()
    plt.savefig(dirs['lipid_analysis'] / 'lipid_cluster_forest.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_lipid_cluster_associations_ref(ref_df: pd.DataFrame,
                                        cluster_names: Dict,
                                        dirs: Dict,
                                        ref_cluster: int = 0):
    """绘制以参考组为基准的脂质-聚类关联可视化（热图 + lollipop）。"""
    if ref_df.empty:
        print("  参考组可视化: 无数据可绘制")
        return

    # 只保留非参考组聚类
    ref_df = ref_df[ref_df['cluster'] != ref_cluster].copy()
    if ref_df.empty:
        print("  参考组可视化: 仅有参考组，跳过")
        return

    ref_df['annot'] = ref_df.apply(
        lambda r: f"{r['beta']:.3f}{'*' if r.get('fdr_significant', False) else ''}", axis=1
    )

    # 热图数据
    pivot = ref_df.pivot_table(
        values='beta', index='lipid_label', columns='cluster', aggfunc='first'
    )
    annot = ref_df.pivot_table(
        values='annot', index='lipid_label', columns='cluster', aggfunc='first'
    )

    # 按字母顺序排序Y轴
    sorted_index = sorted(pivot.index, key=lambda x: x.lower())
    pivot = pivot.reindex(sorted_index)
    annot = annot.reindex(index=sorted_index, columns=pivot.columns)

    col_labels = [f"Cluster {int(c)}\n{cluster_names.get(int(c), '')}" for c in pivot.columns]

    fig, ax = plt.subplots(figsize=(max(8, len(col_labels)*2.6), 7))
    sns.heatmap(
        pivot,
        annot=annot,
        fmt="",
        cmap='RdBu_r',
        center=0,
        cbar_kws={'label': 'β Coefficient'},
        ax=ax,
        linewidths=0.5,
        linecolor='white'
    )
    ax.set_title(f"β Coefficient Heatmap (vs Cluster {ref_cluster})\n* = FDR < 0.05", fontsize=14, fontweight='bold')
    ax.set_xlabel('Diet Pattern Cluster', fontsize=12, fontweight='bold')
    ax.set_ylabel('Lipid', fontsize=12, fontweight='bold')
    ax.set_xticklabels(col_labels, rotation=0)
    plt.tight_layout()
    plt.savefig(dirs['lipid_analysis'] / 'lipid_cluster_heatmap_ref.png', dpi=300, bbox_inches='tight')
    plt.savefig(dirs['lipid_analysis'] / 'lipid_cluster_heatmap_ref.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # Lollipop 图
    fig, ax = plt.subplots(figsize=(10, 5 + 0.6*len(pivot.index)))
    colors = []
    for c in pivot.columns:
        try:
            cid = int(c)
        except Exception:
            cid = c
        colors.append(get_cluster_color(cid))
    xlims = []
    for idx, c in enumerate(pivot.columns):
        sub = ref_df[ref_df['cluster'] == c].sort_values('lipid_label', key=lambda s: s.str.lower())
        if sub.empty:
            continue
        y_pos = np.arange(len(sub)) + idx*0.1  # 微偏移避免重合
        ax.hlines(y=y_pos, xmin=0, xmax=sub['beta'], colors=colors[idx], linestyles='-', alpha=0.6)
        ax.plot(sub['beta'], y_pos, 'o', color=colors[idx], label=f"Cluster {int(c)} vs {ref_cluster}")
        for x, y, beta, sig in zip(sub['beta'], y_pos, sub['beta'], sub['fdr_significant']):
            txt = f"{beta:.3f}{'*' if sig else ''}"
            ax.text(x + (0.002 if x>=0 else -0.002), y, txt, va='center', ha='left' if x>=0 else 'right', fontsize=9)
        xlims.extend(sub['beta'].tolist())

    ax.axvline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel('β Coefficient (Cluster vs Reference)', fontsize=12, fontweight='bold')
    ax.set_title(f'Lipid Associations by Cluster (Reference = Cluster {ref_cluster})\n* = FDR < 0.05', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9)
    if xlims:
        max_abs = max(abs(min(xlims)), abs(max(xlims)))
        ax.set_xlim(-max_abs*1.2, max_abs*1.2)
    plt.tight_layout()
    plt.savefig(dirs['lipid_analysis'] / 'lipid_cluster_lollipop_ref.png', dpi=300, bbox_inches='tight')
    plt.savefig(dirs['lipid_analysis'] / 'lipid_cluster_lollipop_ref.pdf', dpi=300, bbox_inches='tight')
    plt.close()


def generate_summary_report(profiles_df: pd.DataFrame,
                            cluster_names: Dict,
                            lipid_results: pd.DataFrame,
                            dirs: Dict):
    """生成汇总报告"""
    print("\n" + "="*60)
    print("生成汇总报告...")
    print("="*60)
    
    report_path = dirs['summary'] / 'analysis_summary.txt'
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("24h饮食变量分类分析报告\n")
        f.write("="*80 + "\n\n")
        
        f.write("1. 分析概述\n")
        f.write("-"*40 + "\n")
        f.write("本分析将93个24h饮食变量按照以下两个维度进行分类:\n")
        f.write("  - 食物来源: Animal Source vs Plant Source\n")
        f.write("  - 加工程度: Extra Processed vs Processed\n\n")
        
        f.write("2. 聚类结果\n")
        f.write("-"*40 + "\n")
        f.write(f"最优聚类数: {len(cluster_names)}\n\n")
        
        f.write("聚类特征:\n")
        for cid, name in sorted(cluster_names.items()):
            profile = profiles_df[profiles_df['cluster'] == cid].iloc[0]
            f.write(f"\n  Cluster {cid}: {name}\n")
            f.write(f"    样本数: {int(profile['n_samples'])}\n")
            f.write(f"    Animal Source: {profile['animal_pct']:.1f}%\n")
            f.write(f"    Plant Source: {profile['plant_pct']:.1f}%\n")
            f.write(f"    Extra Processed: {profile['extra_processed_pct']:.1f}%\n")
            f.write(f"    Processed: {profile['processed_pct']:.1f}%\n")
        
        f.write("\n3. 脂质关联分析\n")
        f.write("-"*40 + "\n")
        
        if not lipid_results.empty:
            sig_results = lipid_results[lipid_results['fdr_significant']]
            f.write(f"显著关联数 (FDR < 0.05): {len(sig_results)}\n\n")
            
            if not sig_results.empty:
                f.write("显著结果:\n")
                for _, row in sig_results.iterrows():
                    direction = "增加" if row['beta'] > 0 else "降低"
                    f.write(f"  - {row['lipid_label']}: {row['cluster_name']} "
                            f"{direction} (β={row['beta']:.4f}, p_FDR={row['p_fdr']:.4f})\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("分析完成\n")
    
    print(f"报告已保存: {report_path}")


# ============================================================================
# 主函数
# ============================================================================
def main():
    """主函数"""
    print("\n" + "="*80)
    print("  24h饮食变量分类分析")
    print("="*80)
    
    # 创建输出目录
    dirs = setup_output_dirs()
    
    # 读取分类信息
    classification = load_classification_from_excel()
    
    # 保存分类信息
    classification_summary = []
    for cat_type, categories in classification.items():
        for cat_name, vars_list in categories.items():
            for var in vars_list:
                classification_summary.append({
                    'category_type': cat_type,
                    'category': cat_name,
                    'variable': var
                })
    pd.DataFrame(classification_summary).to_csv(
        dirs['classification'] / 'variable_classification.csv', index=False
    )
    
    # 加载数据
    data, diet_cols, desc_to_simple = load_and_prepare_data(classification)
    
    if len(data) < 100:
        print("错误: 有效样本数过少，无法进行分析")
        return
    
    # 计算分类得分
    scores = compute_classification_scores(data, classification, diet_cols, desc_to_simple)
    scores.to_csv(dirs['data_prep'] / 'classification_scores.csv')
    
    # 聚类分析
    cluster_labels, cluster_info = perform_clustering(scores, dirs)
    
    if cluster_labels.empty:
        print("聚类分析失败")
        return
    
    # 保存聚类标签
    cluster_labels.to_csv(dirs['clustering'] / 'cluster_labels.csv')
    
    # 分析饮食模式
    profiles_df, cluster_names = analyze_diet_patterns(
        data, cluster_labels, classification, diet_cols, dirs, desc_to_simple
    )
    
    # 分析脂质关联
    lipid_results = analyze_lipid_by_cluster(
        data, cluster_labels, cluster_names, dirs
    )
    
    # 生成汇总报告
    generate_summary_report(profiles_df, cluster_names, lipid_results, dirs)
    
    print("\n" + "="*80)
    print(f"分析完成！结果保存在: {OUTPUT_DIR}")
    print("="*80)


if __name__ == "__main__":
    main()
