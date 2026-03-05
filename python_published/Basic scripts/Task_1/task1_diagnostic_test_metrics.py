"""
Task 1: Diagnostic Test Performance Metrics
Calculates sensitivity, specificity, and accuracy with 95% confidence intervals
for multiple diagnostic tests.
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

def exact_binomial_test(successes, total, p0=0.5):
    """Perform exact binomial test"""
    result = stats.binomtest(successes, total, p0, alternative='two-sided')
    return result.pvalue

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, 'task1_data_template.csv')

# Load data
if not os.path.exists(csv_file):
    print(f"ERROR: CSV file not found: {csv_file}")
    print("Please create the CSV file using task1_data_template.csv as a template")
    exit(1)

df = pd.read_csv(csv_file)

# Validate required columns
required_cols = ['Test_Name', 'TP', 'TN', 'FP', 'FN']
if not all(col in df.columns for col in required_cols):
    print(f"ERROR: CSV must contain columns: {', '.join(required_cols)}")
    exit(1)

# Get output file name
date_str = datetime.now().strftime('%Y%m%d')
task_name = 'task1'
existing_files = glob.glob(os.path.join(script_dir, f'{task_name}_{date_str}_*.txt'))
run_number = len(existing_files) + 1
output_file = os.path.join(script_dir, f'{task_name}_{date_str}_{run_number:02d}.txt')

# Open output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write("TASK 1: DIAGNOSTIC TEST PERFORMANCE METRICS\n")
    f.write("="*80 + "\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input File: {csv_file}\n")
    f.write("\n")
    
    results = []
    
    for idx, row in df.iterrows():
        test_name = row['Test_Name']
        tp = int(row['TP'])
        tn = int(row['TN'])
        fp = int(row['FP'])
        fn = int(row['FN'])
        
        f.write("-"*80 + "\n")
        f.write(f"Test: {test_name}\n")
        f.write("-"*80 + "\n")
        f.write(f"True Positive (TP):  {tp}\n")
        f.write(f"True Negative (TN):  {tn}\n")
        f.write(f"False Positive (FP): {fp}\n")
        f.write(f"False Negative (FN): {fn}\n")
        f.write("\n")
        
        # Calculate Sensitivity
        sensitivity_total = tp + fn
        if sensitivity_total > 0:
            sensitivity = tp / sensitivity_total
            ci_sens = wilson_score_interval(tp, sensitivity_total)
            
            # Check if exact test needed
            use_exact_sens = (tp < 5) or (fn < 5)
            
            f.write("SENSITIVITY:\n")
            f.write(f"  Value: {sensitivity:.4f} ({sensitivity*100:.2f}%)\n")
            f.write(f"  Formula: TP/(TP+FN) = {tp}/({tp}+{fn}) = {tp}/{sensitivity_total}\n")
            f.write(f"  95% CI (Wilson Score): {ci_sens[0]:.4f} - {ci_sens[1]:.4f} ({ci_sens[0]*100:.2f}% - {ci_sens[1]*100:.2f}%)\n")
            if use_exact_sens:
                p_exact = exact_binomial_test(tp, sensitivity_total)
                f.write(f"  [Note: Exact binomial test used (p-value: {p_exact:.6f}) due to small sample size]\n")
            f.write("\n")
        else:
            sensitivity = np.nan
            ci_sens = (np.nan, np.nan)
            f.write("SENSITIVITY: Cannot calculate (TP + FN = 0)\n\n")
        
        # Calculate Specificity
        specificity_total = tn + fp
        if specificity_total > 0:
            specificity = tn / specificity_total
            ci_spec = wilson_score_interval(tn, specificity_total)
            
            # Check if exact test needed
            use_exact_spec = (tn < 5) or (fp < 5)
            
            f.write("SPECIFICITY:\n")
            f.write(f"  Value: {specificity:.4f} ({specificity*100:.2f}%)\n")
            f.write(f"  Formula: TN/(TN+FP) = {tn}/({tn}+{fp}) = {tn}/{specificity_total}\n")
            f.write(f"  95% CI (Wilson Score): {ci_spec[0]:.4f} - {ci_spec[1]:.4f} ({ci_spec[0]*100:.2f}% - {ci_spec[1]*100:.2f}%)\n")
            if use_exact_spec:
                p_exact = exact_binomial_test(tn, specificity_total)
                f.write(f"  [Note: Exact binomial test used (p-value: {p_exact:.6f}) due to small sample size]\n")
            f.write("\n")
        else:
            specificity = np.nan
            ci_spec = (np.nan, np.nan)
            f.write("SPECIFICITY: Cannot calculate (TN + FP = 0)\n\n")
        
        # Calculate Accuracy
        accuracy_total = tp + tn + fp + fn
        if accuracy_total > 0:
            accuracy = (tp + tn) / accuracy_total
            ci_acc = wilson_score_interval(tp + tn, accuracy_total)
            
            # Check if exact test needed
            use_exact_acc = (tp < 5) or (tn < 5) or (fp < 5) or (fn < 5)
            
            f.write("ACCURACY:\n")
            f.write(f"  Value: {accuracy:.4f} ({accuracy*100:.2f}%)\n")
            f.write(f"  Formula: (TP+TN)/(TP+TN+FP+FN) = ({tp}+{tn})/({tp}+{tn}+{fp}+{fn}) = {tp+tn}/{accuracy_total}\n")
            f.write(f"  95% CI (Wilson Score): {ci_acc[0]:.4f} - {ci_acc[1]:.4f} ({ci_acc[0]*100:.2f}% - {ci_acc[1]*100:.2f}%)\n")
            if use_exact_acc:
                p_exact = exact_binomial_test(tp + tn, accuracy_total)
                f.write(f"  [Note: Exact binomial test used (p-value: {p_exact:.6f}) due to small sample size]\n")
            f.write("\n")
        else:
            accuracy = np.nan
            ci_acc = (np.nan, np.nan)
            f.write("ACCURACY: Cannot calculate (Total = 0)\n\n")
        
        results.append({
            'Test_Name': test_name,
            'Sensitivity': sensitivity,
            'Sensitivity_CI_Lower': ci_sens[0],
            'Sensitivity_CI_Upper': ci_sens[1],
            'Specificity': specificity,
            'Specificity_CI_Lower': ci_spec[0],
            'Specificity_CI_Upper': ci_spec[1],
            'Accuracy': accuracy,
            'Accuracy_CI_Lower': ci_acc[0],
            'Accuracy_CI_Upper': ci_acc[1]
        })
    
    # Summary table
    f.write("\n" + "="*80 + "\n")
    f.write("SUMMARY TABLE\n")
    f.write("="*80 + "\n\n")
    
    results_df = pd.DataFrame(results)
    f.write(results_df.to_string(index=False))
    f.write("\n\n")
    
    f.write("="*80 + "\n")
    f.write("ANALYSIS COMPLETE\n")
    f.write("="*80 + "\n")

print(f"Analysis complete. Results saved to: {output_file}")

