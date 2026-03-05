"""
Task 5: Two-Group Comparison (Study vs Control)
Compares two groups of organisms using a 2x2 table (rows = groups, columns = yes/no).
Calculates proportion for each group and performs statistical comparison.
"""

import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime
import os
import glob

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
csv_file = os.path.join(script_dir, 'task5_data_template.csv')

# Load data
if not os.path.exists(csv_file):
    print(f"ERROR: CSV file not found: {csv_file}")
    print("Please create the CSV file using task5_data_template.csv as a template")
    exit(1)

df = pd.read_csv(csv_file)

# Validate required columns
required_cols = ['Group', 'Yes', 'No']
if not all(col in df.columns for col in required_cols):
    print(f"ERROR: CSV must contain columns: {', '.join(required_cols)}")
    print(f"Found columns: {', '.join(df.columns)}")
    exit(1)

if len(df) != 2:
    print("ERROR: Need exactly 2 groups for comparison (Study and Control)")
    exit(1)

# Get output file name
date_str = datetime.now().strftime('%Y%m%d')
task_name = 'task5'
existing_files = glob.glob(os.path.join(script_dir, f'{task_name}_{date_str}_*.txt'))
run_number = len(existing_files) + 1
output_file = os.path.join(script_dir, f'{task_name}_{date_str}_{run_number:02d}.txt')

# Open output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write("TASK 5: TWO-GROUP COMPARISON (STUDY VS CONTROL)\n")
    f.write("="*80 + "\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input File: {csv_file}\n")
    f.write("\n")
    
    # Calculate proportion for each group
    f.write("-"*80 + "\n")
    f.write("PROPORTION FOR EACH GROUP\n")
    f.write("-"*80 + "\n\n")
    
    for idx, row in df.iterrows():
        group_name = row['Group']
        yes_count = int(row['Yes'])
        no_count = int(row['No'])
        total = yes_count + no_count
        
        if total > 0:
            proportion = yes_count / total
            ci = wilson_score_interval(yes_count, total)
            
            f.write(f"{group_name}:\n")
            f.write(f"  Yes: {yes_count}\n")
            f.write(f"  No: {no_count}\n")
            f.write(f"  Total: {total}\n")
            f.write(f"  Proportion (Yes): {proportion:.4f} ({proportion*100:.2f}%)\n")
            f.write(f"  Formula: Yes/(Yes+No) = {yes_count}/({yes_count}+{no_count}) = {yes_count}/{total}\n")
            f.write(f"  95% CI (Wilson Score): {ci[0]:.4f} - {ci[1]:.4f} ({ci[0]*100:.2f}% - {ci[1]*100:.2f}%)\n")
            f.write("\n")
        else:
            f.write(f"{group_name}: Cannot calculate (Total = 0)\n\n")
    
    # Create 2x2 contingency table
    f.write("-"*80 + "\n")
    f.write("TWO-BY-TWO CONTINGENCY TABLE\n")
    f.write("-"*80 + "\n\n")
    
    group1 = df.iloc[0]
    group2 = df.iloc[1]
    
    contingency = [
        [int(group1['Yes']), int(group1['No'])],
        [int(group2['Yes']), int(group2['No'])]
    ]
    
    f.write(f"{'Group':15s} {'Yes':>10s} {'No':>10s} {'Total':>10s}\n")
    f.write(f"{group1['Group']:15s} {int(group1['Yes']):10d} {int(group1['No']):10d} {int(group1['Yes'])+int(group1['No']):10d}\n")
    f.write(f"{group2['Group']:15s} {int(group2['Yes']):10d} {int(group2['No']):10d} {int(group2['Yes'])+int(group2['No']):10d}\n")
    f.write("\n")
    
    # Perform statistical test
    f.write("-"*80 + "\n")
    f.write("STATISTICAL COMPARISON\n")
    f.write("-"*80 + "\n\n")
    
    # Check if any entry < 5
    min_entry = min([min(row) for row in contingency])
    use_fisher = min_entry < 5
    
    if use_fisher:
        oddsratio, p_value = stats.fisher_exact(contingency)
        f.write("Test Used: Fisher's Exact Test\n")
        f.write(f"Reason: At least one cell has fewer than 5 entries (minimum: {min_entry})\n")
        f.write(f"Odds Ratio: {oddsratio:.4f}\n")
        f.write(f"p-value: {p_value:.6f}\n")
        
        # Calculate confidence interval for odds ratio
        if oddsratio > 0:
            # Log odds ratio CI
            log_or = np.log(oddsratio)
            # Standard error approximation
            se_log_or = np.sqrt(1/contingency[0][0] + 1/contingency[0][1] + 
                               1/contingency[1][0] + 1/contingency[1][1])
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)
            f.write(f"95% CI for Odds Ratio: {ci_lower:.4f} - {ci_upper:.4f}\n")
    else:
        chi2, p_value, dof, expected = stats.chi2_contingency(contingency, correction=True)
        f.write("Test Used: Chi-squared Test (with Yates' correction)\n")
        f.write(f"Reason: All cells have at least 5 entries (minimum: {min_entry})\n")
        f.write(f"Chi-squared statistic: {chi2:.4f}\n")
        f.write(f"Degrees of freedom: {dof}\n")
        f.write(f"p-value: {p_value:.6f}\n")
        
        # Calculate odds ratio
        oddsratio = (contingency[0][0] * contingency[1][1]) / (contingency[0][1] * contingency[1][0])
        f.write(f"Odds Ratio: {oddsratio:.4f}\n")
        
        # Calculate confidence interval for odds ratio
        if oddsratio > 0:
            log_or = np.log(oddsratio)
            se_log_or = np.sqrt(1/contingency[0][0] + 1/contingency[0][1] + 
                               1/contingency[1][0] + 1/contingency[1][1])
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)
            f.write(f"95% CI for Odds Ratio: {ci_lower:.4f} - {ci_upper:.4f}\n")
    
    f.write("\n")
    
    # Interpretation
    f.write("-"*80 + "\n")
    f.write("INTERPRETATION\n")
    f.write("-"*80 + "\n\n")
    
    if p_value < 0.05:
        f.write("CONCLUSION: There is a statistically significant difference between the two groups (p < 0.05)\n")
        if oddsratio > 1:
            f.write(f"The {group1['Group']} group has {oddsratio:.2f} times the odds of 'Yes' compared to the {group2['Group']} group.\n")
        elif oddsratio < 1:
            f.write(f"The {group2['Group']} group has {1/oddsratio:.2f} times the odds of 'Yes' compared to the {group1['Group']} group.\n")
    else:
        f.write("CONCLUSION: There is NO statistically significant difference between the two groups (p >= 0.05)\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("ANALYSIS COMPLETE\n")
    f.write("="*80 + "\n")

print(f"Analysis complete. Results saved to: {output_file}")