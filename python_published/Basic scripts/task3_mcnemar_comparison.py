"""
Task 3: McNemar Test Comparison
Same as Task 2, but uses McNemar test and associated exact test for paired data.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
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
csv_file = os.path.join(script_dir, 'task3_data_template.csv')

# Load data
if not os.path.exists(csv_file):
    print(f"ERROR: CSV file not found: {csv_file}")
    print("Please create the CSV file using task3_data_template.csv as a template")
    exit(1)

df = pd.read_csv(csv_file)

# Validate required columns - McNemar needs paired data structure
# Expected format: Test1_Correct_Test2_Correct, Test1_Correct_Test2_Incorrect, 
#                  Test1_Incorrect_Test2_Correct, Test1_Incorrect_Test2_Incorrect
required_cols = ['Test1_Correct_Test2_Correct', 'Test1_Correct_Test2_Incorrect',
                 'Test1_Incorrect_Test2_Correct', 'Test1_Incorrect_Test2_Incorrect']
if not all(col in df.columns for col in required_cols):
    print(f"ERROR: CSV must contain columns: {', '.join(required_cols)}")
    print("These represent the 2x2 table for paired data:")
    print("  Test1_Correct_Test2_Correct: Both tests correct")
    print("  Test1_Correct_Test2_Incorrect: Test1 correct, Test2 incorrect")
    print("  Test1_Incorrect_Test2_Correct: Test1 incorrect, Test2 correct")
    print("  Test1_Incorrect_Test2_Incorrect: Both tests incorrect")
    exit(1)

# Get output file name
date_str = datetime.now().strftime('%Y%m%d')
task_name = 'task3'
existing_files = glob.glob(os.path.join(script_dir, f'{task_name}_{date_str}_*.txt'))
run_number = len(existing_files) + 1
output_file = os.path.join(script_dir, f'{task_name}_{date_str}_{run_number:02d}.txt')

# Open output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write("TASK 3: MCNEMAR TEST COMPARISON\n")
    f.write("="*80 + "\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input File: {csv_file}\n")
    f.write("\n")
    
    # Get the first row (should contain the 2x2 table)
    row = df.iloc[0]
    
    both_correct = int(row['Test1_Correct_Test2_Correct'])
    test1_only = int(row['Test1_Correct_Test2_Incorrect'])
    test2_only = int(row['Test1_Incorrect_Test2_Correct'])
    both_incorrect = int(row['Test1_Incorrect_Test2_Incorrect'])
    
    # Calculate sensitivity for each test
    f.write("-"*80 + "\n")
    f.write("SENSITIVITY FOR EACH TEST\n")
    f.write("-"*80 + "\n\n")
    
    test1_correct = both_correct + test1_only
    test1_incorrect = test2_only + both_incorrect
    test1_total = test1_correct + test1_incorrect
    
    test2_correct = both_correct + test2_only
    test2_incorrect = test1_only + both_incorrect
    test2_total = test2_correct + test2_incorrect
    
    total_pairs = both_correct + test1_only + test2_only + both_incorrect
    
    if test1_total > 0:
        test1_sensitivity = test1_correct / test1_total
        ci_test1 = wilson_score_interval(test1_correct, test1_total)
        f.write("Test 1:\n")
        f.write(f"  Correct: {test1_correct}\n")
        f.write(f"  Incorrect: {test1_incorrect}\n")
        f.write(f"  Total: {test1_total}\n")
        f.write(f"  Sensitivity: {test1_sensitivity:.4f} ({test1_sensitivity*100:.2f}%)\n")
        f.write(f"  Formula: Correct/(Correct+Incorrect) = {test1_correct}/({test1_correct}+{test1_incorrect}) = {test1_correct}/{test1_total}\n")
        f.write(f"  95% CI (Wilson Score): {ci_test1[0]:.4f} - {ci_test1[1]:.4f} ({ci_test1[0]*100:.2f}% - {ci_test1[1]*100:.2f}%)\n\n")
    
    if test2_total > 0:
        test2_sensitivity = test2_correct / test2_total
        ci_test2 = wilson_score_interval(test2_correct, test2_total)
        f.write("Test 2:\n")
        f.write(f"  Correct: {test2_correct}\n")
        f.write(f"  Incorrect: {test2_incorrect}\n")
        f.write(f"  Total: {test2_total}\n")
        f.write(f"  Sensitivity: {test2_sensitivity:.4f} ({test2_sensitivity*100:.2f}%)\n")
        f.write(f"  Formula: Correct/(Correct+Incorrect) = {test2_correct}/({test2_correct}+{test2_incorrect}) = {test2_correct}/{test2_total}\n")
        f.write(f"  95% CI (Wilson Score): {ci_test2[0]:.4f} - {ci_test2[1]:.4f} ({ci_test2[0]*100:.2f}% - {ci_test2[1]*100:.2f}%)\n\n")
    
    # Create 2x2 contingency table for McNemar
    f.write("-"*80 + "\n")
    f.write("MCNEMAR TEST CONTINGENCY TABLE\n")
    f.write("(Paired data format)\n")
    f.write("-"*80 + "\n\n")
    
    contingency = np.array([
        [both_correct, test1_only],
        [test2_only, both_incorrect]
    ])
    
    f.write(f"                Test2_Correct    Test2_Incorrect    Total\n")
    f.write(f"Test1_Correct   {both_correct:8d}         {test1_only:8d}      {test1_correct:5d}\n")
    f.write(f"Test1_Incorrect {test2_only:8d}         {both_incorrect:8d}      {test1_incorrect:5d}\n")
    f.write(f"Total           {test2_correct:8d}         {test2_incorrect:8d}      {total_pairs:5d}\n")
    f.write("\n")
    f.write("NOTE: This table shows paired comparisons (same subjects tested with both methods).\n\n")
    
    # Perform McNemar test
    f.write("-"*80 + "\n")
    f.write("MCNEMAR TEST\n")
    f.write("-"*80 + "\n\n")
    
    # Check if exact test needed (if discordant pairs < 25)
    discordant = contingency[0, 1] + contingency[1, 0]
    use_exact = discordant < 25
    
    try:
        if use_exact:
            result = mcnemar(contingency, exact=True)
            f.write("Test Used: McNemar's Exact Test\n")
            f.write(f"Reason: Number of discordant pairs ({discordant}) < 25\n")
        else:
            result = mcnemar(contingency, exact=False, correction=True)
            f.write("Test Used: McNemar's Test (with continuity correction)\n")
            f.write(f"Reason: Number of discordant pairs ({discordant}) >= 25\n")
        
        f.write(f"Chi-squared statistic: {result.statistic:.4f}\n")
        f.write(f"p-value: {result.pvalue:.6f}\n")
        f.write("\n")
        
        if result.pvalue < 0.05:
            f.write("CONCLUSION: There is a statistically significant difference between the two tests (p < 0.05)\n")
        else:
            f.write("CONCLUSION: There is NO statistically significant difference between the two tests (p >= 0.05)\n")
    except Exception as e:
        f.write(f"ERROR performing McNemar test: {str(e)}\n")
        f.write("NOTE: McNemar test requires paired data structure.\n")
        f.write("Please ensure your data represents paired comparisons.\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("ANALYSIS COMPLETE\n")
    f.write("="*80 + "\n")

print(f"Analysis complete. Results saved to: {output_file}")

