"""
Analysis 2 Method Comparison: Pairwise comparisons between methods on Original group
- Compares different methods on the same isolates (paired data)
- Uses McNemar's test for paired categorical data
- Evaluates success rates and correct identification rates
- MALDI evaluated at cutoffs 2.0, 1.9, and 1.8
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
import os

# Set up paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_file = os.path.join(base_dir, 'data', 'AI_database.csv')
output_dir = os.path.join(base_dir, 'results')

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load data
df = pd.read_csv(data_file)

print("="*80)
print("ANALYSIS 2: METHOD COMPARISON ON ORIGINAL GROUP")
print("Pairwise comparisons between methods on the same isolates")
print("="*80)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def wilson_score_interval(successes, total, confidence=0.95):
    """Calculate Wilson Score confidence interval for proportion"""
    if total == 0:
        return (0, 0)
    z = stats.norm.ppf((1 + confidence) / 2)
    p = successes / total
    denominator = 1 + (z**2 / total)
    centre_adjusted_probability = (p + (z**2 / (2 * total))) / denominator
    adjusted_standard_deviation = np.sqrt((p * (1 - p) + z**2 / (4 * total)) / total) / denominator
    lower_bound = centre_adjusted_probability - z * adjusted_standard_deviation
    upper_bound = centre_adjusted_probability + z * adjusted_standard_deviation
    return (max(0, lower_bound), min(1, upper_bound))

def holm_bonferroni_correction(p_values):
    """Apply Holm-Bonferroni correction for multiple comparisons"""
    rejected, p_corrected, _, _ = multipletests(p_values, method='holm')
    return p_corrected

def mcnemar_test_paired(data1, data2):
    """Perform McNemar's test on paired binary data"""
    # Create contingency table
    both_success = ((data1 == 1) & (data2 == 1)).sum()
    method1_only = ((data1 == 1) & (data2 == 0)).sum()
    method2_only = ((data1 == 0) & (data2 == 1)).sum()
    both_fail = ((data1 == 0) & (data2 == 0)).sum()
    
    # Create contingency table for McNemar's test
    contingency_table = np.array([[both_success, method1_only],
                                   [method2_only, both_fail]])
    
    # Perform McNemar's test
    try:
        result = mcnemar(contingency_table, exact=False, correction=True)
        return result.statistic, result.pvalue
    except:
        # If test fails, return None
        return None, None

# ==============================================================================
# FILTER TO ORIGINAL GROUP ONLY
# ==============================================================================

df_original = df[df['Plate'] == 'Original'].copy()
print(f"\nTotal Original samples: {len(df_original)}")
print(f"Unique isolates: {df_original['Isolate_ID'].nunique()}")

# ==============================================================================
# DEFINE METHODS AND CUTOFFS
# ==============================================================================

CUTOFF_N_BIOCHEM = 0.85
CUTOFF_P_BIOCHEM = 0.85
CUTOFF_VITEK = 85
CUTOFF_MALDI = [2.0, 1.9, 1.8]

# Define method configurations
method_configs = {
    'N_biochem': {
        'score_col': 'N_biochem_score',
        'id_col': 'N_biochem_ID',
        'cutoff': CUTOFF_N_BIOCHEM
    },
    'P_biochem': {
        'score_col': 'P_biochem_score',
        'id_col': 'P_biochem_ID',
        'cutoff': CUTOFF_P_BIOCHEM
    },
    'Vitek': {
        'score_col': 'Vitek_score',
        'id_col': 'Vitek_ID',
        'cutoff': CUTOFF_VITEK
    }
}

# Add MALDI methods for each cutoff
for cutoff in CUTOFF_MALDI:
    method_configs[f'MALDI_{cutoff}'] = {
        'score_col': 'MALDI_score',
        'id_col': 'MALDI_ID',
        'cutoff': cutoff
    }

# ==============================================================================
# CALCULATE SUCCESS AND CORRECT RATES FOR EACH METHOD
# ==============================================================================

method_results = {}
method_data = {}

for method_name, method_info in method_configs.items():
    # Calculate success (passes cutoff)
    df_original[f'{method_name}_success'] = (
        df_original[method_info['score_col']] >= method_info['cutoff']
    ).astype(int)
    
    # Calculate correct (successful AND matches Final_ID)
    df_original[f'{method_name}_correct'] = (
        (df_original[f'{method_name}_success'] == 1) &
        (df_original[method_info['id_col']] == df_original['Final_ID'])
    ).astype(int)
    
    # Store data
    method_data[method_name] = {
        'success': df_original[f'{method_name}_success'].values,
        'correct': df_original[f'{method_name}_correct'].values
    }
    
    # Calculate rates
    total = len(df_original)
    success_count = df_original[f'{method_name}_success'].sum()
    correct_count = df_original[f'{method_name}_correct'].sum()
    
    success_rate = success_count / total if total > 0 else 0
    correct_rate = correct_count / success_count if success_count > 0 else 0
    
    ci_success = wilson_score_interval(success_count, total)
    ci_correct = wilson_score_interval(correct_count, success_count) if success_count > 0 else (0, 0)
    
    method_results[method_name] = {
        'total': total,
        'success_count': success_count,
        'success_rate': success_rate,
        'success_ci_lower': ci_success[0],
        'success_ci_upper': ci_success[1],
        'correct_count': correct_count,
        'correct_rate': correct_rate,
        'correct_ci_lower': ci_correct[0],
        'correct_ci_upper': ci_correct[1]
    }

# ==============================================================================
# PAIRWISE COMPARISONS
# ==============================================================================

comparison_results = []
method_names = list(method_configs.keys())

# Compare all pairs of methods
for i, method1 in enumerate(method_names):
    for method2 in method_names[i+1:]:
        # Success rate comparison
        success1 = method_data[method1]['success']
        success2 = method_data[method2]['success']
        
        stat_success, p_success = mcnemar_test_paired(success1, success2)
        
        # Correct rate comparison (only among those who succeeded in both methods)
        # We need to compare correct rates among samples where both methods succeeded
        both_successful = (success1 == 1) & (success2 == 1)
        
        if both_successful.sum() > 0:
            correct1_both = method_data[method1]['correct'][both_successful]
            correct2_both = method_data[method2]['correct'][both_successful]
            
            stat_correct, p_correct = mcnemar_test_paired(correct1_both, correct2_both)
        else:
            stat_correct, p_correct = None, None
        
        comparison_results.append({
            'Method1': method1,
            'Method2': method2,
            'Metric': 'Success Rate',
            'Test': 'McNemar',
            'Statistic': stat_success,
            'p_value': p_success,
            'Method1_Rate': method_results[method1]['success_rate'],
            'Method2_Rate': method_results[method2]['success_rate'],
            'Method1_Count': f"{method_results[method1]['success_count']}/{method_results[method1]['total']}",
            'Method2_Count': f"{method_results[method2]['success_count']}/{method_results[method2]['total']}"
        })
        
        if stat_correct is not None:
            comparison_results.append({
                'Method1': method1,
                'Method2': method2,
                'Metric': 'Correct Rate',
                'Test': 'McNemar',
                'Statistic': stat_correct,
                'p_value': p_correct,
                'Method1_Rate': method_results[method1]['correct_rate'],
                'Method2_Rate': method_results[method2]['correct_rate'],
                'Method1_Count': f"{method_results[method1]['correct_count']}/{method_results[method1]['success_count']}",
                'Method2_Count': f"{method_results[method2]['correct_count']}/{method_results[method2]['success_count']}"
            })

# Apply multiple comparison correction
comparison_df = pd.DataFrame(comparison_results)
if len(comparison_df) > 0:
    # Separate by metric for correction
    for metric in ['Success Rate', 'Correct Rate']:
        mask = comparison_df['Metric'] == metric
        if mask.sum() > 0:
            p_values = comparison_df.loc[mask, 'p_value'].values
            p_values = [p if p is not None and not np.isnan(p) else 1.0 for p in p_values]
            p_corrected = holm_bonferroni_correction(p_values)
            comparison_df.loc[mask, 'p_value_corrected'] = p_corrected

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

# Save method results
results_df = pd.DataFrame([{
    'Method': method,
    **results
} for method, results in method_results.items()])

results_df.to_csv(os.path.join(output_dir, 'analysis2_method_comparison_original_results.csv'), index=False)

# Save comparison results
comparison_df.to_csv(os.path.join(output_dir, 'analysis2_method_comparison_original_pairwise.csv'), index=False)

print("\n" + "="*80)
print("RESULTS SAVED")
print("="*80)
print(f"\nMethod results saved to:")
print(f"  - {os.path.join(output_dir, 'analysis2_method_comparison_original_results.csv')}")
print(f"\nPairwise comparisons saved to:")
print(f"  - {os.path.join(output_dir, 'analysis2_method_comparison_original_pairwise.csv')}")

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

print("\n" + "="*80)
print("SUMMARY OF RESULTS")
print("="*80)

print("\nSuccess Rates (passing cutoff):")
print("-" * 80)
for method, results in sorted(method_results.items(), key=lambda x: x[1]['success_rate'], reverse=True):
    print(f"{method:15s}: {results['success_rate']:.1%} ({results['success_count']}/{results['total']}) "
          f"[95% CI: {results['success_ci_lower']:.1%} - {results['success_ci_upper']:.1%}]")

print("\nCorrect Rates (among successful identifications):")
print("-" * 80)
for method, results in sorted(method_results.items(), key=lambda x: x[1]['correct_rate'], reverse=True):
    if results['success_count'] > 0:
        print(f"{method:15s}: {results['correct_rate']:.1%} ({results['correct_count']}/{results['success_count']}) "
              f"[95% CI: {results['correct_ci_lower']:.1%} - {results['correct_ci_upper']:.1%}]")
    else:
        print(f"{method:15s}: N/A (no successful identifications)")

print("\n" + "="*80)
print("SIGNIFICANT PAIRWISE DIFFERENCES (p < 0.05 after correction)")
print("="*80)

significant_comparisons = comparison_df[
    (comparison_df['p_value_corrected'] < 0.05) & 
    (comparison_df['p_value_corrected'].notna())
]

if len(significant_comparisons) > 0:
    for _, row in significant_comparisons.iterrows():
        print(f"\n{row['Metric']}: {row['Method1']} vs {row['Method2']}")
        print(f"  {row['Method1']}: {row['Method1_Rate']:.1%} ({row['Method1_Count']})")
        print(f"  {row['Method2']}: {row['Method2_Rate']:.1%} ({row['Method2_Count']})")
        print(f"  p-value (corrected): {row['p_value_corrected']:.6f}")
else:
    print("\nNo significant differences found after multiple comparison correction.")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)

