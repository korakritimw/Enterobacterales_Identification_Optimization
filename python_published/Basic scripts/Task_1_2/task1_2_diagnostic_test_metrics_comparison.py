"""
Task 1_2: Diagnostic Test Comparison Across Multiple Methods
Compares 4 diagnostic test methods using 4x2 contingency tables for:
1. Sensitivity (TP vs FN)
2. Specificity (TN vs FP)
3. Accuracy (TP+TN vs FP+FN)

For each metric:
- Performs overall Chi-squared or Fisher's exact test
- If significant, performs pairwise comparisons with Holm-Bonferroni correction
"""

import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime
import os
import glob
import itertools

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

def perform_pairwise_comparisons(df, test_name_col, col1, col2, metric_name, f):
    """Perform pairwise comparisons with Holm-Bonferroni correction"""
    
    f.write("\n" + "-"*80 + "\n")
    f.write(f"PAIRWISE COMPARISONS FOR {metric_name.upper()} (Holm-Bonferroni Method)\n")
    f.write("-"*80 + "\n\n")
    
    # Get all pairs and compute p-values first
    pairs = list(itertools.combinations(range(len(df)), 2))
    pairwise_results = []
    
    for i, j in pairs:
        test1 = df.iloc[i]
        test2 = df.iloc[j]
        
        # Create 2x2 table for this pair
        pair_contingency = [
            [int(test1[col1]), int(test1[col2])],
            [int(test2[col1]), int(test2[col2])]
        ]
        
        # Determine test for this pair
        pair_min = min([min(row) for row in pair_contingency])
        
        if pair_min < 5:
            oddsratio, p_value_pair = stats.fisher_exact(pair_contingency)
            test_used = "Fisher's Exact Test"
            test_stat = oddsratio
            test_stat_name = "Odds Ratio"
        else:
            chi2_pair, p_value_pair, dof_pair, expected_pair = stats.chi2_contingency(pair_contingency, correction=True)
            test_used = "Chi-squared Test (with Yates' correction)"
            test_stat = chi2_pair
            test_stat_name = "Chi-squared statistic"
        
        pairwise_results.append({
            'i': i,
            'j': j,
            'test1': test1,
            'test2': test2,
            'contingency': pair_contingency,
            'p_value': p_value_pair,
            'test_used': test_used,
            'test_stat': test_stat,
            'test_stat_name': test_stat_name,
            'pair_min': pair_min
        })
    
    # Sort by p-value for Holm-Bonferroni
    pairwise_results.sort(key=lambda x: x['p_value'])
    
    # Apply Holm-Bonferroni correction
    n_comparisons = len(pairs)
    for rank, result in enumerate(pairwise_results, 1):
        result['rank'] = rank
        result['holm_alpha'] = 0.05 / (n_comparisons - rank + 1)
        result['significant'] = result['p_value'] < result['holm_alpha']
    
    # Write results in original order
    pairwise_results.sort(key=lambda x: (x['i'], x['j']))
    
    for pair_idx, result in enumerate(pairwise_results, 1):
        f.write(f"Comparison {pair_idx}: {result['test1'][test_name_col]} vs {result['test2'][test_name_col]}\n")
        f.write("-" * 40 + "\n")
        
        f.write(f"{'Test':20s} {col1:>10s} {col2:>10s}\n")
        f.write(f"{result['test1'][test_name_col]:20s} {int(result['test1'][col1]):10d} {int(result['test1'][col2]):10d}\n")
        f.write(f"{result['test2'][test_name_col]:20s} {int(result['test2'][col1]):10d} {int(result['test2'][col2]):10d}\n")
        f.write("\n")
        
        f.write(f"Test: {result['test_used']}\n")
        if result['pair_min'] < 5:
            f.write(f"Reason: At least one cell < 5 (minimum: {result['pair_min']})\n")
        else:
            f.write(f"Reason: All cells >= 5 (minimum: {result['pair_min']})\n")
        f.write(f"{result['test_stat_name']}: {result['test_stat']:.4f}\n")
        f.write(f"p-value: {result['p_value']:.6f}\n")
        f.write(f"Holm-Bonferroni rank: {result['rank']} (sorted by p-value)\n")
        f.write(f"Holm-Bonferroni adjusted alpha: {result['holm_alpha']:.6f} (0.05/{n_comparisons - result['rank'] + 1})\n")
        
        if result['significant']:
            f.write(f"Result: SIGNIFICANT difference (p < {result['holm_alpha']:.6f})\n")
        else:
            f.write(f"Result: NO significant difference (p >= {result['holm_alpha']:.6f})\n")
        
        f.write("\n")

def analyze_metric(df, test_name_col, col1, col2, metric_name, metric_formula, f):
    """Analyze a specific metric (sensitivity, specificity, or accuracy)"""
    
    f.write("\n" + "="*80 + "\n")
    f.write(f"{metric_name.upper()} ANALYSIS\n")
    f.write("="*80 + "\n\n")
    
    # Calculate metric for each test
    f.write("-"*80 + "\n")
    f.write(f"{metric_name.upper()} FOR EACH TEST\n")
    f.write("-"*80 + "\n\n")
    
    for idx, row in df.iterrows():
        test_name = row[test_name_col]
        val1 = int(row[col1])
        val2 = int(row[col2])
        total = val1 + val2
        
        if total > 0:
            metric_value = val1 / total
            ci = wilson_score_interval(val1, total)
            
            f.write(f"{test_name}:\n")
            f.write(f"  {col1}: {val1}\n")
            f.write(f"  {col2}: {val2}\n")
            f.write(f"  Total: {total}\n")
            f.write(f"  {metric_name}: {metric_value:.4f} ({metric_value*100:.2f}%)\n")
            f.write(f"  Formula: {metric_formula}\n")
            f.write(f"  Calculation: {val1}/{total}\n")
            f.write(f"  95% CI (Wilson Score): {ci[0]:.4f} - {ci[1]:.4f} ({ci[0]*100:.2f}% - {ci[1]*100:.2f}%)\n")
            f.write("\n")
        else:
            f.write(f"{test_name}: Cannot calculate (Total = 0)\n\n")
    
    # Create contingency table
    f.write("-"*80 + "\n")
    f.write(f"CONTINGENCY TABLE FOR {metric_name.upper()} ({len(df)}x2)\n")
    f.write("-"*80 + "\n\n")
    
    # Build contingency table
    contingency = []
    for idx, row in df.iterrows():
        contingency.append([int(row[col1]), int(row[col2])])
    
    contingency_array = np.array(contingency)
    
    # Print table
    f.write(f"{'Test':20s} {col1:>10s} {col2:>10s} {'Total':>10s}\n")
    for idx, row in df.iterrows():
        test_name = row[test_name_col]
        val1 = int(row[col1])
        val2 = int(row[col2])
        total = val1 + val2
        f.write(f"{test_name:20s} {val1:10d} {val2:10d} {total:10d}\n")
    f.write("\n")
    
    # OVERALL STATISTICAL TEST
    f.write("-"*80 + "\n")
    f.write(f"OVERALL STATISTICAL COMPARISON FOR {metric_name.upper()}\n")
    f.write("-"*80 + "\n\n")
    
    # Check if any entry < 5 for Fisher's exact test
    min_entry = contingency_array.min()
    use_fisher = min_entry < 5
    
    if use_fisher:
        try:
            chi2, p_value_overall, dof, expected = stats.chi2_contingency(contingency_array)
            f.write("Test Used: Fisher's Exact Test (generalized for RxC table)\n")
            f.write(f"Reason: At least one cell has fewer than 5 entries (minimum: {min_entry})\n")
            f.write(f"Chi-squared statistic: {chi2:.4f}\n")
            f.write(f"Degrees of freedom: {dof}\n")
            f.write(f"p-value: {p_value_overall:.6f}\n")
        except Exception as e:
            f.write(f"ERROR: Could not perform Fisher's exact test: {e}\n")
            p_value_overall = 1.0
    else:
        chi2, p_value_overall, dof, expected = stats.chi2_contingency(contingency_array, correction=False)
        f.write("Test Used: Chi-squared Test\n")
        f.write(f"Reason: All cells have at least 5 entries (minimum: {min_entry})\n")
        f.write(f"Chi-squared statistic: {chi2:.4f}\n")
        f.write(f"Degrees of freedom: {dof}\n")
        f.write(f"p-value: {p_value_overall:.6f}\n")
    
    f.write("\n")
    
    if p_value_overall < 0.05:
        f.write(f"OVERALL CONCLUSION: There is a statistically significant difference in {metric_name} among the tests (p < 0.05)\n")
        f.write("Proceeding with pairwise comparisons...\n")
        
        # Perform pairwise comparisons
        perform_pairwise_comparisons(df, test_name_col, col1, col2, metric_name, f)
    else:
        f.write(f"OVERALL CONCLUSION: There is NO statistically significant difference in {metric_name} among the tests (p >= 0.05)\n")
        f.write("Pairwise comparisons not performed.\n")

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, 'task1_2_data_template.csv')

# Load data
if not os.path.exists(csv_file):
    print(f"ERROR: CSV file not found: {csv_file}")
    print("Please create the CSV file using task1_2_data_template.csv as a template")
    exit(1)

df = pd.read_csv(csv_file)

# Validate required columns
required_cols = ['Test_Name', 'TP', 'TN', 'FP', 'FN']
if not all(col in df.columns for col in required_cols):
    print(f"ERROR: CSV must contain columns: {', '.join(required_cols)}")
    exit(1)

if len(df) < 2:
    print("ERROR: Need at least 2 tests for comparison")
    exit(1)

# Create calculated columns for accuracy
df['TP_TN'] = df['TP'] + df['TN']
df['FP_FN'] = df['FP'] + df['FN']

# Get output file name
date_str = datetime.now().strftime('%Y%m%d')
task_name = 'task1_2'
existing_files = glob.glob(os.path.join(script_dir, f'{task_name}_{date_str}_*.txt'))
run_number = len(existing_files) + 1
output_file = os.path.join(script_dir, f'{task_name}_{date_str}_{run_number:02d}.txt')

# Open output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write(f"TASK 1_2: DIAGNOSTIC TEST COMPARISON ({len(df)} TESTS)\n")
    f.write("="*80 + "\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input File: {csv_file}\n")
    f.write("\n")
    f.write("This analysis compares diagnostic test performance across multiple methods\n")
    f.write("using contingency tables for three key metrics:\n")
    f.write("1. Sensitivity (TP vs FN)\n")
    f.write("2. Specificity (TN vs FP)\n")
    f.write("3. Accuracy (TP+TN vs FP+FN)\n")
    f.write("\n")
    
    # Analyze Sensitivity
    analyze_metric(df, 'Test_Name', 'TP', 'FN', 'Sensitivity', 
                   'TP/(TP+FN)', f)
    
    # Analyze Specificity
    analyze_metric(df, 'Test_Name', 'TN', 'FP', 'Specificity', 
                   'TN/(TN+FP)', f)
    
    # Analyze Accuracy
    analyze_metric(df, 'Test_Name', 'TP_TN', 'FP_FN', 'Accuracy', 
                   '(TP+TN)/(TP+TN+FP+FN)', f)
    
    f.write("\n" + "="*80 + "\n")
    f.write("ANALYSIS COMPLETE\n")
    f.write("="*80 + "\n")

print(f"Analysis complete. Results saved to: {output_file}")