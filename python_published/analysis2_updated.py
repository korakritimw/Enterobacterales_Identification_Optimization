"""
Analysis 2 (Updated): Main analysis with Positive control group
- Success and correct identification rates for all methods
- Pairwise comparisons for each method across groups
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
print("ANALYSIS 2 (UPDATED): MAIN ANALYSIS WITH POSITIVE CONTROL")
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

# ==============================================================================
# ANALYSIS 2: MAIN ANALYSIS
# ==============================================================================

# Define cutoffs
CUTOFF_N_BIOCHEM = 0.85
CUTOFF_P_BIOCHEM = 0.85
CUTOFF_VITEK = 85
CUTOFF_MALDI = [2.0, 1.9, 1.8]  # Evaluate these cutoffs for MALDI

methods = {
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

results_analysis2 = []
comparison_results = []

# All groups including Positive
all_groups = ['Original', 'Subculture', 'Positive']

# Calculate rates for each method and group
for method_name, method_info in methods.items():
    for plate_type in all_groups:
        subset = df[df['Plate'] == plate_type].copy()
        
        # Calculate success rate (passing cutoff)
        subset['passes_cutoff'] = subset[method_info['score_col']] >= method_info['cutoff']
        success_count = subset['passes_cutoff'].sum()
        total_count = len(subset)
        success_rate = success_count / total_count if total_count > 0 else 0
        ci_success = wilson_score_interval(success_count, total_count)
        
        # Calculate correct identification rate (among those passing cutoff)
        successful = subset[subset['passes_cutoff']]
        if len(successful) > 0:
            # Check if ID matches Final_ID
            correct = (successful[method_info['id_col']] == successful['Final_ID']).sum()
            correct_rate = correct / len(successful) if len(successful) > 0 else 0
            ci_correct = wilson_score_interval(correct, len(successful))
        else:
            correct = 0
            correct_rate = 0
            ci_correct = (0, 0)
        
        results_analysis2.append({
            'Method': method_name,
            'Group': plate_type,
            'Total_Entries': total_count,
            'Success_Count': success_count,
            'Success_Rate': success_rate,
            'Success_CI_Lower': ci_success[0],
            'Success_CI_Upper': ci_success[1],
            'Correct_Count': correct,
            'Correct_Rate': correct_rate,
            'Correct_CI_Lower': ci_correct[0],
            'Correct_CI_Upper': ci_correct[1]
        })

# MALDI with multiple cutoffs
for cutoff in CUTOFF_MALDI:
    for plate_type in all_groups:
        subset = df[df['Plate'] == plate_type].copy()
        
        subset['passes_cutoff'] = subset['MALDI_score'] >= cutoff
        success_count = subset['passes_cutoff'].sum()
        total_count = len(subset)
        success_rate = success_count / total_count if total_count > 0 else 0
        ci_success = wilson_score_interval(success_count, total_count)
        
        successful = subset[subset['passes_cutoff']]
        if len(successful) > 0:
            correct = (successful['MALDI_ID'] == successful['Final_ID']).sum()
            correct_rate = correct / len(successful) if len(successful) > 0 else 0
            ci_correct = wilson_score_interval(correct, len(successful))
        else:
            correct = 0
            correct_rate = 0
            ci_correct = (0, 0)
        
        results_analysis2.append({
            'Method': f'MALDI_{cutoff}',
            'Group': plate_type,
            'Total_Entries': total_count,
            'Success_Count': success_count,
            'Success_Rate': success_rate,
            'Success_CI_Lower': ci_success[0],
            'Success_CI_Upper': ci_success[1],
            'Correct_Count': correct,
            'Correct_Rate': correct_rate,
            'Correct_CI_Lower': ci_correct[0],
            'Correct_CI_Upper': ci_correct[1]
        })

results_df2 = pd.DataFrame(results_analysis2)
results_df2.to_csv(os.path.join(output_dir, 'analysis2_updated_results.csv'), index=False)

print("\n--- Pairwise Comparisons for Each Method ---")

# ==============================================================================
# PAIRWISE COMPARISONS FOR EACH METHOD
# ==============================================================================

# For each method, compare groups pairwise
for method_name, method_info in methods.items():
    print(f"\n{method_name}:")
    
    # Get data for all groups
    group_data = {}
    for group in all_groups:
        group_data[group] = df[df['Plate'] == group].copy()
    
    # Success rate comparisons
    # Original vs Subculture (paired if same Isolate_IDs)
    if 'Isolate_ID' in df.columns:
        # Check if we can do paired comparison
        orig_ids = set(group_data['Original']['Isolate_ID'].unique())
        sub_ids = set(group_data['Subculture']['Isolate_ID'].unique())
        
        if orig_ids == sub_ids:
            # Paired comparison using McNemar
            merged = pd.merge(
                group_data['Original'][['Isolate_ID', method_info['score_col']]],
                group_data['Subculture'][['Isolate_ID', method_info['score_col']]],
                on='Isolate_ID',
                suffixes=('_orig', '_sub')
            )
            merged['orig_pass'] = merged[f"{method_info['score_col']}_orig"] >= method_info['cutoff']
            merged['sub_pass'] = merged[f"{method_info['score_col']}_sub"] >= method_info['cutoff']
            
            table = pd.crosstab(merged['orig_pass'], merged['sub_pass'])
            if table.shape == (2, 2):
                result = mcnemar(table, exact=False, correction=True)
                comparison_results.append({
                    'Method': method_name,
                    'Comparison': 'Original vs Subculture',
                    'Metric': 'Success Rate',
                    'Test': 'McNemar',
                    'Statistic': result.statistic,
                    'p_value': result.pvalue
                })
        else:
            # Unpaired comparison
            orig_pass = (group_data['Original'][method_info['score_col']] >= method_info['cutoff']).sum()
            orig_total = len(group_data['Original'])
            sub_pass = (group_data['Subculture'][method_info['score_col']] >= method_info['cutoff']).sum()
            sub_total = len(group_data['Subculture'])
            
            contingency = [[orig_pass, orig_total - orig_pass],
                          [sub_pass, sub_total - sub_pass]]
            
            if min([min(row) for row in contingency]) < 5:
                oddsratio, p_value = stats.fisher_exact(contingency)
                comparison_results.append({
                    'Method': method_name,
                    'Comparison': 'Original vs Subculture',
                    'Metric': 'Success Rate',
                    'Test': "Fisher's Exact",
                    'Statistic': oddsratio,
                    'p_value': p_value
                })
            else:
                chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
                comparison_results.append({
                    'Method': method_name,
                    'Comparison': 'Original vs Subculture',
                    'Metric': 'Success Rate',
                    'Test': 'Chi-squared (Yates)',
                    'Statistic': chi2,
                    'p_value': p_value
                })
        
        # Original vs Positive (unpaired)
        orig_pass = (group_data['Original'][method_info['score_col']] >= method_info['cutoff']).sum()
        orig_total = len(group_data['Original'])
        pos_pass = (group_data['Positive'][method_info['score_col']] >= method_info['cutoff']).sum()
        pos_total = len(group_data['Positive'])
        
        contingency = [[orig_pass, orig_total - orig_pass],
                      [pos_pass, pos_total - pos_pass]]
        
        if min([min(row) for row in contingency]) < 5:
            oddsratio, p_value = stats.fisher_exact(contingency)
            comparison_results.append({
                'Method': method_name,
                'Comparison': 'Original vs Positive',
                'Metric': 'Success Rate',
                'Test': "Fisher's Exact",
                'Statistic': oddsratio,
                'p_value': p_value
            })
        else:
            chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
            comparison_results.append({
                'Method': method_name,
                'Comparison': 'Original vs Positive',
                'Metric': 'Success Rate',
                'Test': 'Chi-squared (Yates)',
                'Statistic': chi2,
                'p_value': p_value
            })
        
        # Subculture vs Positive (unpaired)
        sub_pass = (group_data['Subculture'][method_info['score_col']] >= method_info['cutoff']).sum()
        sub_total = len(group_data['Subculture'])
        pos_pass = (group_data['Positive'][method_info['score_col']] >= method_info['cutoff']).sum()
        pos_total = len(group_data['Positive'])
        
        contingency = [[sub_pass, sub_total - sub_pass],
                      [pos_pass, pos_total - pos_pass]]
        
        if min([min(row) for row in contingency]) < 5:
            oddsratio, p_value = stats.fisher_exact(contingency)
            comparison_results.append({
                'Method': method_name,
                'Comparison': 'Subculture vs Positive',
                'Metric': 'Success Rate',
                'Test': "Fisher's Exact",
                'Statistic': oddsratio,
                'p_value': p_value
            })
        else:
            chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
            comparison_results.append({
                'Method': method_name,
                'Comparison': 'Subculture vs Positive',
                'Metric': 'Success Rate',
                'Test': 'Chi-squared (Yates)',
                'Statistic': chi2,
                'p_value': p_value
            })
        
        # Correct identification rate comparisons (only for those passing cutoff)
        for group1, group2 in [('Original', 'Subculture'), ('Original', 'Positive'), ('Subculture', 'Positive')]:
            g1_success = group_data[group1][group_data[group1][method_info['score_col']] >= method_info['cutoff']].copy()
            g2_success = group_data[group2][group_data[group2][method_info['score_col']] >= method_info['cutoff']].copy()
            
            if len(g1_success) > 0 and len(g2_success) > 0:
                g1_success['correct'] = (g1_success[method_info['id_col']] == g1_success['Final_ID'])
                g2_success['correct'] = (g2_success[method_info['id_col']] == g2_success['Final_ID'])
                
                g1_correct = g1_success['correct'].sum()
                g1_total = len(g1_success)
                g2_correct = g2_success['correct'].sum()
                g2_total = len(g2_success)
                
                contingency = [[g1_correct, g1_total - g1_correct],
                              [g2_correct, g2_total - g2_correct]]
                
                if min([min(row) for row in contingency]) < 5:
                    oddsratio, p_value = stats.fisher_exact(contingency)
                    comparison_results.append({
                        'Method': method_name,
                        'Comparison': f'{group1} vs {group2}',
                        'Metric': 'Correct Rate',
                        'Test': "Fisher's Exact",
                        'Statistic': oddsratio,
                        'p_value': p_value
                    })
                else:
                    chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
                    comparison_results.append({
                        'Method': method_name,
                        'Comparison': f'{group1} vs {group2}',
                        'Metric': 'Correct Rate',
                        'Test': 'Chi-squared (Yates)',
                        'Statistic': chi2,
                        'p_value': p_value
                    })

# MALDI comparisons for each cutoff
for cutoff in CUTOFF_MALDI:
    method_name = f'MALDI_{cutoff}'
    print(f"\n{method_name}:")
    
    group_data = {}
    for group in all_groups:
        group_data[group] = df[df['Plate'] == group].copy()
    
    # Success rate comparisons
    for group1, group2 in [('Original', 'Subculture'), ('Original', 'Positive'), ('Subculture', 'Positive')]:
        g1_pass = (group_data[group1]['MALDI_score'] >= cutoff).sum()
        g1_total = len(group_data[group1])
        g2_pass = (group_data[group2]['MALDI_score'] >= cutoff).sum()
        g2_total = len(group_data[group2])
        
        contingency = [[g1_pass, g1_total - g1_pass],
                      [g2_pass, g2_total - g2_pass]]
        
        if min([min(row) for row in contingency]) < 5:
            oddsratio, p_value = stats.fisher_exact(contingency)
            comparison_results.append({
                'Method': method_name,
                'Comparison': f'{group1} vs {group2}',
                'Metric': 'Success Rate',
                'Test': "Fisher's Exact",
                'Statistic': oddsratio,
                'p_value': p_value
            })
        else:
            chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
            comparison_results.append({
                'Method': method_name,
                'Comparison': f'{group1} vs {group2}',
                'Metric': 'Success Rate',
                'Test': 'Chi-squared (Yates)',
                'Statistic': chi2,
                'p_value': p_value
            })
        
        # Correct rate comparisons
        g1_success = group_data[group1][group_data[group1]['MALDI_score'] >= cutoff].copy()
        g2_success = group_data[group2][group_data[group2]['MALDI_score'] >= cutoff].copy()
        
        if len(g1_success) > 0 and len(g2_success) > 0:
            g1_correct = (g1_success['MALDI_ID'] == g1_success['Final_ID']).sum()
            g1_total = len(g1_success)
            g2_correct = (g2_success['MALDI_ID'] == g2_success['Final_ID']).sum()
            g2_total = len(g2_success)
            
            contingency = [[g1_correct, g1_total - g1_correct],
                          [g2_correct, g2_total - g2_correct]]
            
            if min([min(row) for row in contingency]) < 5:
                oddsratio, p_value = stats.fisher_exact(contingency)
                comparison_results.append({
                    'Method': method_name,
                    'Comparison': f'{group1} vs {group2}',
                    'Metric': 'Correct Rate',
                    'Test': "Fisher's Exact",
                    'Statistic': oddsratio,
                    'p_value': p_value
                })
            else:
                chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
                comparison_results.append({
                    'Method': method_name,
                    'Comparison': f'{group1} vs {group2}',
                    'Metric': 'Correct Rate',
                    'Test': 'Chi-squared (Yates)',
                    'Statistic': chi2,
                    'p_value': p_value
                })

# Apply Holm-Bonferroni correction
if comparison_results:
    p_values = [r['p_value'] for r in comparison_results]
    p_corrected = holm_bonferroni_correction(p_values)
    for i, result in enumerate(comparison_results):
        result['p_value_corrected'] = p_corrected[i]

comparison_df = pd.DataFrame(comparison_results)
comparison_df.to_csv(os.path.join(output_dir, 'analysis2_updated_pairwise_comparisons.csv'), index=False)

print("\n" + "="*80)
print("RESULTS SAVED")
print("="*80)
print(f"\nAnalysis 2 results saved to:")
print(f"  - {os.path.join(output_dir, 'analysis2_updated_results.csv')}")
print(f"  - {os.path.join(output_dir, 'analysis2_updated_pairwise_comparisons.csv')}")

# Display summary
print("\n" + "="*80)
print("SUMMARY OF RESULTS")
print("="*80)

print("\nSuccess Rates by Method and Group:")
summary_df = results_df2.pivot(index='Method', columns='Group', values='Success_Rate')
print(summary_df.round(3))

print("\nCorrect Rates by Method and Group:")
correct_df = results_df2.pivot(index='Method', columns='Group', values='Correct_Rate')
print(correct_df.round(3))

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)

