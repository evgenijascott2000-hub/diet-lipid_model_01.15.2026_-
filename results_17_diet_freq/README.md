# 17个FFQ饮食频率变量分析

## 分析配置

### 数据来源
- **输入数据**: `results_gpt5codex_7lipid-33cate/merged_lipid_diet_dataset_with_prs.csv`
- **变量来源**: FFQ (Food Frequency Questionnaire) 饮食频率问卷
- **样本数量**: 72,811 (完整PRS数据)

### 饮食变量 (17个FFQ daily变量)
1. diet_beef_intake_daily
2. diet_lamb_mutton_intake_daily
3. diet_pork_intake_daily
4. diet_poultry_intake_daily
5. diet_processed_meat_intake_daily
6. diet_oily_fish_intake_daily
7. diet_non_oily_fish_intake_daily
8. diet_cheese_intake_daily
9. diet_cooked_vegetable_intake_daily
10. diet_salad_raw_vegetable_intake_daily
11. diet_fresh_fruit_intake_daily
12. diet_dried_fruit_intake_daily
13. diet_tea_intake_daily
14. diet_coffee_intake_daily
15. diet_water_intake_daily (cups/day)
16. diet_bread_intake_daily
17. diet_cereal_intake_daily

### 脂质指标 (6个)
1. Triglycerides (log)
2. Total cholesterol
3. HDL cholesterol (log)
4. LDL cholesterol
5. Apolipoprotein A (log)
6. Apolipoprotein B

### 分析方法
- ✅ **OLS回归**: 完整模型 + 饮食变量 + 协变量
- ✅ **LASSO正则化**: 变量选择
- ✅ **PCA分析**: 饮食模式降维 (已启用)
- ✅ **PC-脂质相关性**: 协变量调整后的偏相关
- ❌ **Bootstrap**: 已禁用 (节省时间)

### 协变量 (12个)
**数值型 (9个)**:
- age
- bmi
- waist_circumference
- physical_activity_met_minutes_per_week
- alcohol_intake_frequency
- smoking_pack_years
- dietary_fiber_g_per_day
- total_energy_kcal_per_day
- 6个脂质PRS

**分类型 (3个)**:
- sex
- ethnicity
- education_level

## 关键特点
1. **无需时间匹配**: FFQ变量不依赖24h recall,无需与抽血时间匹配
2. **完整数据集**: 使用全部72,811个有完整PRS的样本
3. **标准化单位**: 所有频率变量转换为每日摄入量
4. **PCA启用**: 可以看到饮食模式与脂质的关系

## 输出目录结构
```
results_17_diet_freq/
├── linear_models/
│   ├── 01_data_preparation/
│   ├── 02_pca_analysis/              # PCA碎石图、loadings热图
│   ├── 03_ols_regression/            # OLS系数、图表
│   ├── 04_lasso_regularization/      # LASSO系数、路径图
│   ├── 05_bootstrap_selection/       # (已禁用)
│   ├── 06_model_diagnostics/         # 模型拟合指标
│   ├── 07_pc_lipid_correlations/     # PC与脂质的相关性
│   └── 08_sensitivity_analysis/      # 敏感性分析
└── non-linear_models/                # (如果生成)
```

## 运行时间
- 预计: ~10-15分钟 (无Bootstrap)
- 进程ID: 82259
- 日志文件: `run_17_diet_freq.log`

## 与24h分析的对比
| 特征 | 17个FFQ变量 | 93个24h变量 |
|------|-------------|-------------|
| 数据类型 | 频率问卷 | 24小时膳食回顾 |
| 变量数 | 17 | 93 |
| 样本数 | 72,811 | 10,943 |
| 时间匹配 | 不需要 | ≤36小时 |
| 完整度要求 | - | ≥80% |
| 数据粒度 | 食物组 | 具体食物 |
