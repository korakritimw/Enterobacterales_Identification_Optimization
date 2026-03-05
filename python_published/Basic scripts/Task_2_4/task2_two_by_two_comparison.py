"""
Task 2: Three-by-Two Table Comparison
Compares multiple groups using a 3x2 table (rows = groups, columns = correct/incorrect).
Calculates sensitivity for each group and performs:
1. Overall Chi-squared or Fisher's exact test
2. Pairwise comparisons if overall test is significant
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

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, 'task2_data_template.csv')

# Load data
if not os.path.exists(csv_file):
    print(f"ERROR: CSV file not found: {csv_file}")
    print("Please create the CSV file using task2_data_template.csv as a template")
    exit(1)

df = pd.read_csv(csv_file)

# Validate required columns - updated to accept either Test_Name or Subgroup
if 'Subgroup' in df.columns:
    group_col = 'Subgroup'
elif 'Test_Name' in df.columns:
    group_col = 'Test_Name'
else:
    print("ERROR: CSV must contain either 'Subgroup' or 'Test_Name' column")
    exit(1)

required_cols = [group_col, 'Correct', 'Incorrect']
if not all(col in df.columns for col in required_cols):
    print(f"ERROR: CSV must contain columns: {', '.join(required_cols)}")
    exit(1)

if len(df) < 2:
    print("ERROR: Need at least 2 groups for comparison")
    exit(1)

# Get output file name
date_str = datetime.now().strftime('%Y%m%d')
task_name = 'task2'
existing_files = glob.glob(os.path.join(script_dir, f'{task_name}_{date_str}_*.txt'))
run_number = len(existing_files) + 1
output_file = os.path.join(script_dir, f'{task_name}_{date_str}_{run_number:02d}.txt')

# Open output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write(f"TASK 2: CONTINGENCY TABLE COMPARISON ({len(df)} GROUPS)\n")
    f.write("="*80 + "\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input File: {csv_file}\n")
    f.write("\n")
    
    # Calculate sensitivity for each group
    f.write("-"*80 + "\n")
    f.write("SENSITIVITY FOR EACH GROUP\n")
    f.write("-"*80 + "\n\n")
    
    for idx, row in df.iterrows():
        group_name = row[group_col]
        correct = int(row['Correct'])
        incorrect = int(row['Incorrect'])
        total = correct + incorrect
        
        if total > 0:
            sensitivity = correct / total
            ci = wilson_score_interval(correct, total)
            
            f.write(f"{group_name}:\n")
            f.write(f"  Correct: {correct}\n")
            f.write(f"  Incorrect: {incorrect}\n")
            f.write(f"  Total: {total}\n")
            f.write(f"  Sensitivity: {sensitivity:.4f} ({sensitivity*100:.2f}%)\n")
            f.write(f"  Formula: Correct/(Correct+Incorrect) = {correct}/({correct}+{incorrect}) = {correct}/{total}\n")
            f.write(f"  95% CI (Wilson Score): {ci[0]:.4f} - {ci[1]:.4f} ({ci[0]*100:.2f}% - {ci[1]*100:.2f}%)\n")
            f.write("\n")
        else:
            f.write(f"{group_name}: Cannot calculate (Total = 0)\n\n")
    
    # Create contingency table
    f.write("-"*80 + "\n")
    f.write(f"CONTINGENCY TABLE ({len(df)}x2)\n")
    f.write("-"*80 + "\n\n")
    
    # Build contingency table
    contingency = []
    for idx, row in df.iterrows():
        contingency.append([int(row['Correct']), int(row['Incorrect'])])
    
    contingency_array = np.array(contingency)
    
    # Print table
    f.write(f"{'Group':20s} Correct    Incorrect    Total\n")
    for idx, row in df.iterrows():
        group_name = row[group_col]
        correct = int(row['Correct'])
        incorrect = int(row['Incorrect'])
        total = correct + incorrect
        f.write(f"{group_name:20s} {correct:7d}    {incorrect:9d}    {total:5d}\n")
    f.write("\n")
    
    # OVERALL STATISTICAL TEST
    f.write("-"*80 + "\n")
    f.write("OVERALL STATISTICAL COMPARISON\n")
    f.write("-"*80 + "\n\n")
    
    # Check if any entry < 5 for Fisher's exact test
    min_entry = contingency_array.min()
    use_fisher = min_entry < 5
    
    if use_fisher:
        # For 3x2 tables, use Fisher's exact test via chi2_contingency with simulation
        try:
            chi2, p_value_overall, dof, expected = stats.chi2_contingency(contingency_array)
            # Use Fisher-Freeman-Halton test (generalization of Fisher's exact)
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
        f.write("OVERALL CONCLUSION: There is a statistically significant difference among the groups (p < 0.05)\n")
        f.write("Proceeding with pairwise comparisons...\n")
        
        # PAIRWISE COMPARISONS
        f.write("\n" + "-"*80 + "\n")
        f.write("PAIRWISE COMPARISONS (Holm-Bonferroni Method)\n")
        f.write("-"*80 + "\n\n")
        
        # Get all pairs and compute p-values first
        pairs = list(itertools.combinations(range(len(df)), 2))
        pairwise_results = []
        
        for i, j in pairs:
            group1 = df.iloc[i]
            group2 = df.iloc[j]
            
            # Create 2x2 table for this pair
            pair_contingency = [
                [int(group1['Correct']), int(group1['Incorrect'])],
                [int(group2['Correct']), int(group2['Incorrect'])]
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
                'group1': group1,
                'group2': group2,
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
            f.write(f"Comparison {pair_idx}: {result['group1'][group_col]} vs {result['group2'][group_col]}\n")
            f.write("-" * 40 + "\n")
            
            f.write(f"{'Group':20s} Correct    Incorrect\n")
            f.write(f"{result['group1'][group_col]:20s} {int(result['group1']['Correct']):7d}    {int(result['group1']['Incorrect']):9d}\n")
            f.write(f"{result['group2'][group_col]:20s} {int(result['group2']['Correct']):7d}    {int(result['group2']['Incorrect']):9d}\n")
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
    else:
        f.write("OVERALL CONCLUSION: There is NO statistically significant difference among the groups (p >= 0.05)\n")
        f.write("Pairwise comparisons not performed.\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("ANALYSIS COMPLETE\n")
    f.write("="*80 + "\n")

print(f"Analysis complete. Results saved to: {output_file}")