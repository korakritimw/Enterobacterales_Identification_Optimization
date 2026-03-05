"""
Task 4: Numerical Data Comparison
Compares two columns of numerical data using appropriate statistical tests.
Performs Shapiro-Wilk test for normality, then uses appropriate descriptive
statistics and comparison tests.
"""

import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime
import os
import glob

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, 'task4_data_template.csv')

# Load data
if not os.path.exists(csv_file):
    print(f"ERROR: CSV file not found: {csv_file}")
    print("Please create the CSV file using task4_data_template.csv as a template")
    exit(1)

df = pd.read_csv(csv_file)

# Validate required columns
if len(df.columns) < 2:
    print("ERROR: CSV must contain at least 2 columns of numerical data")
    exit(1)

col1_name = df.columns[0]
col2_name = df.columns[1]

# Get output file name
date_str = datetime.now().strftime('%Y%m%d')
task_name = 'task4'
existing_files = glob.glob(os.path.join(script_dir, f'{task_name}_{date_str}_*.txt'))
run_number = len(existing_files) + 1
output_file = os.path.join(script_dir, f'{task_name}_{date_str}_{run_number:02d}.txt')

# Open output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write("TASK 4: NUMERICAL DATA COMPARISON\n")
    f.write("="*80 + "\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input File: {csv_file}\n")
    f.write("\n")
    
    # Extract data (remove missing values)
    col1_data = df[col1_name].dropna().values
    col2_data = df[col2_name].dropna().values
    
    f.write(f"Column 1: {col1_name} (n = {len(col1_data)})\n")
    f.write(f"Column 2: {col2_name} (n = {len(col2_data)})\n")
    f.write("\n")
    
    results = {}
    
    # Analyze each column
    for col_name, col_data in [(col1_name, col1_data), (col2_name, col2_data)]:
        f.write("-"*80 + "\n")
        f.write(f"ANALYSIS: {col_name}\n")
        f.write("-"*80 + "\n\n")
        
        if len(col_data) == 0:
            f.write("ERROR: No valid data points\n\n")
            results[col_name] = {'normal': False, 'n': 0}
            continue
        
        # Shapiro-Wilk test for normality
        if len(col_data) >= 3 and len(col_data) <= 5000:
            shapiro_stat, shapiro_p = stats.shapiro(col_data)
            is_normal = shapiro_p >= 0.05
        else:
            # For very small (<3) or very large (>5000) samples, use alternative
            if len(col_data) < 3:
                f.write("NOTE: Sample size too small for Shapiro-Wilk test (n < 3)\n")
                f.write("Assuming non-normal distribution\n")
                is_normal = False
                shapiro_stat, shapiro_p = np.nan, np.nan
            else:
                # For large samples, use Anderson-Darling or assume normal if symmetric
                f.write("NOTE: Sample size too large for Shapiro-Wilk test (n > 5000)\n")
                f.write("Using alternative normality assessment\n")
                # Check skewness and kurtosis
                skew = stats.skew(col_data)
                kurt = stats.kurtosis(col_data)
                is_normal = (abs(skew) < 2) and (abs(kurt) < 2)
                shapiro_stat, shapiro_p = np.nan, np.nan
        
        f.write(f"Shapiro-Wilk Test:\n")
        if not np.isnan(shapiro_stat):
            f.write(f"  Statistic: {shapiro_stat:.4f}\n")
            f.write(f"  p-value: {shapiro_p:.6f}\n")
        f.write(f"  Conclusion: Data is {'NORMALLY DISTRIBUTED' if is_normal else 'NOT NORMALLY DISTRIBUTED'} (p {'>=' if is_normal else '<'} 0.05)\n")
        f.write("\n")
        
        # Calculate descriptive statistics
        if is_normal:
            mean_val = np.mean(col_data)
            std_val = np.std(col_data, ddof=1)
            n = len(col_data)
            se = std_val / np.sqrt(n)
            t_critical = stats.t.ppf(0.975, n - 1)
            ci_lower = mean_val - t_critical * se
            ci_upper = mean_val + t_critical * se
            
            f.write("Descriptive Statistics (Normal Distribution):\n")
            f.write(f"  Mean: {mean_val:.4f}\n")
            f.write(f"  Standard Deviation: {std_val:.4f}\n")
            f.write(f"  95% Confidence Interval: {ci_lower:.4f} - {ci_upper:.4f}\n")
            f.write(f"  Formula: Mean ± t(0.975, df={n-1}) × SE\n")
            f.write(f"           {mean_val:.4f} ± {t_critical:.4f} × {se:.4f}\n")
            
            results[col_name] = {
                'normal': True,
                'mean': mean_val,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'n': n
            }
        else:
            median_val = np.median(col_data)
            min_val = np.min(col_data)
            max_val = np.max(col_data)
            
            f.write("Descriptive Statistics (Non-Normal Distribution):\n")
            f.write(f"  Median: {median_val:.4f}\n")
            f.write(f"  Range: {min_val:.4f} - {max_val:.4f}\n")
            f.write(f"  Minimum: {min_val:.4f}\n")
            f.write(f"  Maximum: {max_val:.4f}\n")
            
            results[col_name] = {
                'normal': False,
                'median': median_val,
                'min': min_val,
                'max': max_val,
                'n': len(col_data)
            }
        
        f.write("\n")
    
    # Perform comparison test
    f.write("-"*80 + "\n")
    f.write("STATISTICAL COMPARISON\n")
    f.write("-"*80 + "\n\n")
    
    col1_info = results[col1_name]
    col2_info = results[col2_name]
    
    # Determine which test to use
    use_mann_whitney = (
        not col1_info['normal'] or 
        not col2_info['normal'] or 
        col1_info['n'] < 30 or 
        col2_info['n'] < 30
    )
    
    if use_mann_whitney:
        # Mann-Whitney U test
        u_stat, p_value = stats.mannwhitneyu(col1_data, col2_data, alternative='two-sided')
        
        # Calculate z-score using normal approximation
        n1, n2 = len(col1_data), len(col2_data)
        mean_u = n1 * n2 / 2
        std_u = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        if std_u > 0:
            z_score = (u_stat - mean_u) / std_u
        else:
            z_score = 0
        
        f.write("Test Used: Mann-Whitney U Test\n")
        f.write("Reason: ")
        reasons = []
        if not col1_info['normal']:
            reasons.append(f"{col1_name} is not normally distributed")
        if not col2_info['normal']:
            reasons.append(f"{col2_name} is not normally distributed")
        if col1_info['n'] < 30:
            reasons.append(f"{col1_name} has fewer than 30 entries (n={col1_info['n']})")
        if col2_info['n'] < 30:
            reasons.append(f"{col2_name} has fewer than 30 entries (n={col2_info['n']})")
        f.write("; ".join(reasons) + "\n")
        f.write(f"\nMann-Whitney U statistic: {u_stat:.4f}\n")
        f.write(f"z-score: {z_score:.4f}\n")
        f.write(f"p-value: {p_value:.6f}\n")
    else:
        # Student's t-test
        t_stat, p_value = stats.ttest_ind(col1_data, col2_data)
        
        f.write("Test Used: Student's t-test (Independent Samples)\n")
        f.write("Reason: Both columns are normally distributed AND both have at least 30 entries\n")
        f.write(f"  {col1_name}: normally distributed, n = {col1_info['n']}\n")
        f.write(f"  {col2_name}: normally distributed, n = {col2_info['n']}\n")
        f.write(f"\nt-statistic: {t_stat:.4f}\n")
        f.write(f"p-value: {p_value:.6f}\n")
    
    f.write("\n")
    
    if p_value < 0.05:
        f.write("CONCLUSION: There is a statistically significant difference between the two columns (p < 0.05)\n")
    else:
        f.write("CONCLUSION: There is NO statistically significant difference between the two columns (p >= 0.05)\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("ANALYSIS COMPLETE\n")
    f.write("="*80 + "\n")

print(f"Analysis complete. Results saved to: {output_file}")

