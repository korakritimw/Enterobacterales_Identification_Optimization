"""
Cost-Effectiveness Analysis for Bacterial Identification Methods (Updated)
================================================================

This script performs comprehensive cost-effectiveness analysis comparing:
- P_Biochem (Positive Biochemical)
- Vitek2 (Automated system)
- MALDI_2.0 (MALDI-TOF MS with cutoff 2.0)
- MALDI_1.9 (MALDI-TOF MS with cutoff 1.9)
- MALDI_1.8 (MALDI-TOF MS with cutoff 1.8)

Analysis includes:
1. Direct cost comparison
2. Expected cost including error consequences
3. Cost per correct identification
4. Dominance analysis
5. Incremental cost-effectiveness ratios (ICERs)
6. Cost per QALY analysis
7. Net Monetary Benefit (NMB)
8. Sensitivity analysis
9. Budget impact analysis

All costs in USD.
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

# Set up paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cea_dir = os.path.join(base_dir, 'data', 'CEA')
output_dir = os.path.join(base_dir, 'results', 'CEA_Analysis_Updated')
os.makedirs(output_dir, exist_ok=True)

print("="*80)
print("COST-EFFECTIVENESS ANALYSIS (UPDATED)")
print("Bacterial Identification Methods")
print(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
print()

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================
print("PART 1: Loading Data")
print("-"*80)

# Load methods data
methods_file = os.path.join(cea_dir, 'Methods_KI.csv')
methods_df = pd.read_csv(methods_file)
print("Methods data loaded:")
print(methods_df)
print()

# Load risk and cost data
risk_file = os.path.join(cea_dir, 'RIsk_and_cost_KI.csv')
risk_df = pd.read_csv(risk_file)
# Remove empty rows
risk_df = risk_df[risk_df['Risk_type'].notna() & (risk_df['Risk_type'] != 'Notes:')]
risk_df = risk_df.dropna(how='all')  # Remove empty rows

# Convert numeric columns to proper types
numeric_columns = ['Chance', 'QALY_loss', 'Monetized_QALY_loss_lower', 
                   'Monetized_QALY_loss_upper', 'Fixed_cost']
for col in numeric_columns:
    risk_df[col] = pd.to_numeric(risk_df[col], errors='coerce')

print("Risk and cost data loaded:")
print(risk_df[['Risk_type', 'Description', 'Chance', 'QALY_loss', 
               'Monetized_QALY_loss_lower', 'Monetized_QALY_loss_upper', 'Fixed_cost']])
print()

# =============================================================================
# PART 2: CALCULATE EXPECTED COSTS OF ERRORS
# =============================================================================
print("PART 2: Calculating Expected Costs of Errors")
print("-"*80)

# Calculate expected cost and QALY loss for each error type
# We'll use the midpoint of lower and upper bounds for base case

# Very Major (VM) errors
vm_complications = risk_df[(risk_df['Risk_type'] == 'VM') & (risk_df['Description'] == 'Complication')].iloc[0]
vm_death = risk_df[(risk_df['Risk_type'] == 'VM') & (risk_df['Description'] == 'Death')].iloc[0]
vm_baseline = risk_df[(risk_df['Risk_type'] == 'VM') & (risk_df['Description'] == 'Baseline')].iloc[0]

# Expected cost per VM error (midpoint)
vm_cost_mid = (
    vm_complications['Chance'] * (
        vm_complications['Fixed_cost'] + 
        (vm_complications['Monetized_QALY_loss_lower'] + vm_complications['Monetized_QALY_loss_upper']) / 2
    ) +
    vm_death['Chance'] * (
        vm_death['Fixed_cost'] + 
        (vm_death['Monetized_QALY_loss_lower'] + vm_death['Monetized_QALY_loss_upper']) / 2
    ) +
    vm_baseline['Chance'] * 0
)

# Expected cost per VM error (lower bound - best case)
vm_cost_lower = (
    vm_complications['Chance'] * (vm_complications['Fixed_cost'] + vm_complications['Monetized_QALY_loss_lower']) +
    vm_death['Chance'] * (vm_death['Fixed_cost'] + vm_death['Monetized_QALY_loss_lower']) +
    vm_baseline['Chance'] * 0
)

# Expected cost per VM error (upper bound - worst case)
vm_cost_upper = (
    vm_complications['Chance'] * (vm_complications['Fixed_cost'] + vm_complications['Monetized_QALY_loss_upper']) +
    vm_death['Chance'] * (vm_death['Fixed_cost'] + vm_death['Monetized_QALY_loss_upper']) +
    vm_baseline['Chance'] * 0
)

# Expected QALY loss per VM error
vm_qaly_loss = (
    vm_complications['Chance'] * vm_complications['QALY_loss'] +
    vm_death['Chance'] * vm_death['QALY_loss'] +
    vm_baseline['Chance'] * 0
)

# Major (M) errors
m_aki = risk_df[(risk_df['Risk_type'] == 'M') & (risk_df['Description'] == 'AKI')].iloc[0]
m_baseline = risk_df[(risk_df['Risk_type'] == 'M') & (risk_df['Description'] == 'Baseline')].iloc[0]

# Expected cost per M error (midpoint)
m_cost_mid = (
    m_aki['Chance'] * (
        m_aki['Fixed_cost'] + 
        (m_aki['Monetized_QALY_loss_lower'] + m_aki['Monetized_QALY_loss_upper']) / 2
    ) +
    m_baseline['Chance'] * 0
)

# Expected cost per M error (lower bound)
m_cost_lower = (
    m_aki['Chance'] * (m_aki['Fixed_cost'] + m_aki['Monetized_QALY_loss_lower']) +
    m_baseline['Chance'] * 0
)

# Expected cost per M error (upper bound)
m_cost_upper = (
    m_aki['Chance'] * (m_aki['Fixed_cost'] + m_aki['Monetized_QALY_loss_upper']) +
    m_baseline['Chance'] * 0
)

# Expected QALY loss per M error
m_qaly_loss = (
    m_aki['Chance'] * m_aki['QALY_loss'] +
    m_baseline['Chance'] * 0
)

print(f"Expected cost per Very Major error:")
print(f"  Base case (midpoint): ${vm_cost_mid:,.2f}")
print(f"  Lower bound (best case): ${vm_cost_lower:,.2f}")
print(f"  Upper bound (worst case): ${vm_cost_upper:,.2f}")
print(f"  Expected QALY loss: {vm_qaly_loss:.4f}")
print()

print(f"Expected cost per Major error:")
print(f"  Base case (midpoint): ${m_cost_mid:,.2f}")
print(f"  Lower bound (best case): ${m_cost_lower:,.2f}")
print(f"  Upper bound (worst case): ${m_cost_upper:,.2f}")
print(f"  Expected QALY loss: {m_qaly_loss:.4f}")
print()

# =============================================================================
# PART 3: CALCULATE TOTAL EXPECTED COSTS FOR EACH METHOD
# =============================================================================
print("PART 3: Total Expected Costs per Test (Base Case)")
print("-"*80)

results = []

for idx, row in methods_df.iterrows():
    method = row['Methods']
    unit_cost = row['Cost_per_test_USD']
    vm_prob = row['VM_chance']
    m_prob = row['M_chance']
    correct_prob = row['Correct']
    turnaround = row['Turnaround_time_hours']
    
    # Expected cost including error consequences
    expected_cost_mid = unit_cost + (vm_prob * vm_cost_mid) + (m_prob * m_cost_mid)
    expected_cost_lower = unit_cost + (vm_prob * vm_cost_lower) + (m_prob * m_cost_lower)
    expected_cost_upper = unit_cost + (vm_prob * vm_cost_upper) + (m_prob * m_cost_upper)
    
    # Expected QALY loss per test
    expected_qaly_loss = (vm_prob * vm_qaly_loss) + (m_prob * m_qaly_loss)
    
    # Cost per correct identification
    cost_per_correct = expected_cost_mid / correct_prob
    
    results.append({
        'Method': method,
        'Unit_Cost_USD': unit_cost,
        'VM_Error_Rate': vm_prob,
        'Major_Error_Rate': m_prob,
        'Correct_Rate': correct_prob,
        'Expected_Error_Cost_Mid': (vm_prob * vm_cost_mid) + (m_prob * m_cost_mid),
        'Total_Expected_Cost_Mid': expected_cost_mid,
        'Total_Expected_Cost_Lower': expected_cost_lower,
        'Total_Expected_Cost_Upper': expected_cost_upper,
        'Expected_QALY_Loss': expected_qaly_loss,
        'Cost_per_Correct_ID': cost_per_correct,
        'Turnaround_Hours': turnaround
    })

results_df = pd.DataFrame(results)

print(results_df[['Method', 'Unit_Cost_USD', 'Total_Expected_Cost_Mid', 'Cost_per_Correct_ID']].to_string(index=False))
print()

# =============================================================================
# PART 4: DOMINANCE ANALYSIS
# =============================================================================
print("PART 4: Dominance Analysis")
print("-"*80)

# Sort by cost
results_df_sorted = results_df.sort_values('Total_Expected_Cost_Mid').copy()
results_df_sorted['Dominated'] = False
results_df_sorted['Dominance_Type'] = 'On frontier'

# Check for dominance
for i in range(len(results_df_sorted)):
    current_method = results_df_sorted.iloc[i]
    
    for j in range(len(results_df_sorted)):
        if i != j:
            compare_method = results_df_sorted.iloc[j]
            
            # Check if current method is dominated
            if (compare_method['Total_Expected_Cost_Mid'] <= current_method['Total_Expected_Cost_Mid'] and
                compare_method['Correct_Rate'] > current_method['Correct_Rate']):
                results_df_sorted.at[results_df_sorted.index[i], 'Dominated'] = True
                results_df_sorted.at[results_df_sorted.index[i], 'Dominance_Type'] = 'Dominated'
                break

# Check for dominant methods
for i in range(len(results_df_sorted)):
    current_method = results_df_sorted.iloc[i]
    is_dominant = False
    
    for j in range(len(results_df_sorted)):
        if i != j:
            compare_method = results_df_sorted.iloc[j]
            
            if (current_method['Total_Expected_Cost_Mid'] < compare_method['Total_Expected_Cost_Mid'] and
                current_method['Correct_Rate'] >= compare_method['Correct_Rate']):
                is_dominant = True
                break
    
    if is_dominant and not results_df_sorted.iloc[i]['Dominated']:
        results_df_sorted.at[results_df_sorted.index[i], 'Dominance_Type'] = 'Dominant'

print(results_df_sorted[['Method', 'Total_Expected_Cost_Mid', 'Correct_Rate', 'Dominance_Type']].to_string(index=False))
print()

# =============================================================================
# PART 5: INCREMENTAL COST-EFFECTIVENESS RATIOS (ICERs)
# =============================================================================
print("PART 5: Incremental Cost-Effectiveness Ratios (ICERs)")
print("-"*80)

# Remove dominated methods for ICER calculation
frontier_df = results_df_sorted[~results_df_sorted['Dominated']].copy()
frontier_df = frontier_df.sort_values('Total_Expected_Cost_Mid').reset_index(drop=True)

print("Methods on efficiency frontier:")
print(frontier_df[['Method', 'Total_Expected_Cost_Mid', 'Correct_Rate']].to_string(index=False))
print()

# Calculate ICERs
icers = []

for i in range(1, len(frontier_df)):
    method_current = frontier_df.iloc[i]['Method']
    method_previous = frontier_df.iloc[i-1]['Method']
    
    cost_diff = frontier_df.iloc[i]['Total_Expected_Cost_Mid'] - frontier_df.iloc[i-1]['Total_Expected_Cost_Mid']
    effect_diff = frontier_df.iloc[i]['Correct_Rate'] - frontier_df.iloc[i-1]['Correct_Rate']
    
    if effect_diff > 0:
        icer = cost_diff / effect_diff
    else:
        icer = np.inf
    
    icers.append({
        'Comparison': f"{method_current} vs {method_previous}",
        'Cost_Difference_USD': cost_diff,
        'Effectiveness_Difference': effect_diff,
        'ICER_per_Additional_Correct_ID': icer
    })

icers_df = pd.DataFrame(icers)
print("ICERs (Cost per additional correct identification):")
print(icers_df.to_string(index=False))
print()

# =============================================================================
# PART 6: COST PER QALY ANALYSIS
# =============================================================================
print("PART 6: Cost per QALY Analysis")
print("-"*80)

# Calculate ICERs in terms of QALY
qaly_icers = []

for i in range(1, len(frontier_df)):
    method_current = frontier_df.iloc[i]['Method']
    method_previous = frontier_df.iloc[i-1]['Method']
    
    cost_diff = frontier_df.iloc[i]['Total_Expected_Cost_Mid'] - frontier_df.iloc[i-1]['Total_Expected_Cost_Mid']
    qaly_diff = frontier_df.iloc[i-1]['Expected_QALY_Loss'] - frontier_df.iloc[i]['Expected_QALY_Loss']
    
    if qaly_diff > 0:
        icer_qaly = cost_diff / qaly_diff
    else:
        icer_qaly = np.inf
    
    qaly_icers.append({
        'Comparison': f"{method_current} vs {method_previous}",
        'Cost_Difference_USD': cost_diff,
        'QALY_Gained': qaly_diff,
        'ICER_per_QALY_USD': icer_qaly
    })

qaly_icers_df = pd.DataFrame(qaly_icers)
print("ICERs (Cost per QALY gained):")
print(qaly_icers_df.to_string(index=False))
print()

# Interpretation against standard thresholds
print("Interpretation (using standard willingness-to-pay thresholds):")
print("  $50,000/QALY threshold (lower bound, USA):")
for idx, row in qaly_icers_df.iterrows():
    if row['ICER_per_QALY_USD'] < 50000:
        print(f"    [OK] {row['Comparison']}: Cost-effective (ICER = ${row['ICER_per_QALY_USD']:,.2f}/QALY)")
    elif row['ICER_per_QALY_USD'] < np.inf:
        print(f"    [NO] {row['Comparison']}: NOT cost-effective (ICER = ${row['ICER_per_QALY_USD']:,.2f}/QALY)")
print()
print("  $100,000/QALY threshold (upper bound, USA):")
for idx, row in qaly_icers_df.iterrows():
    if row['ICER_per_QALY_USD'] < 100000:
        print(f"    [OK] {row['Comparison']}: Cost-effective (ICER = ${row['ICER_per_QALY_USD']:,.2f}/QALY)")
    elif row['ICER_per_QALY_USD'] < np.inf:
        print(f"    [NO] {row['Comparison']}: NOT cost-effective (ICER = ${row['ICER_per_QALY_USD']:,.2f}/QALY)")
print()

# =============================================================================
# PART 7: NET MONETARY BENEFIT (NMB)
# =============================================================================
print("PART 7: Net Monetary Benefit Analysis")
print("-"*80)

# Use least expensive method as reference
reference_method = results_df.loc[results_df['Total_Expected_Cost_Mid'].idxmin()]
print(f"Reference method (least costly): {reference_method['Method']}")
print()

# Calculate NMB at different WTP thresholds
wtp_thresholds = [50000, 100000, 150000]

nmb_all_results = []

for wtp in wtp_thresholds:
    print(f"Net Monetary Benefit at WTP = ${wtp:,}/QALY:")
    print("-" * 70)
    
    nmb_results = []
    
    for idx, row in results_df.iterrows():
        qaly_gained = reference_method['Expected_QALY_Loss'] - row['Expected_QALY_Loss']
        cost_diff = row['Total_Expected_Cost_Mid'] - reference_method['Total_Expected_Cost_Mid']
        nmb = (qaly_gained * wtp) - cost_diff
        
        nmb_results.append({
            'Method': row['Method'],
            'NMB_USD': nmb,
            'QALY_Gained': qaly_gained,
            'Cost_Difference': cost_diff
        })
        
        nmb_all_results.append({
            'WTP_Threshold': wtp,
            'Method': row['Method'],
            'NMB_USD': nmb,
            'QALY_Gained': qaly_gained,
            'Cost_Difference': cost_diff
        })
    
    nmb_df = pd.DataFrame(nmb_results)
    nmb_df = nmb_df.sort_values('NMB_USD', ascending=False)
    
    print(nmb_df.to_string(index=False))
    print()
    print(f"  -> Best value: {nmb_df.iloc[0]['Method']} (highest NMB)")
    print()

nmb_all_df = pd.DataFrame(nmb_all_results)

# =============================================================================
# PART 8: SENSITIVITY ANALYSIS
# =============================================================================
print("PART 8: Sensitivity Analysis")
print("-"*80)

print("One-Way Sensitivity Analysis: Varying Error Costs")
print()

# Best case (lower bound)
print("BEST CASE (Lower bound error costs):")
results_best = results_df[['Method', 'Total_Expected_Cost_Lower', 'Correct_Rate']].copy()
results_best = results_best.sort_values('Total_Expected_Cost_Lower')
print(results_best.to_string(index=False))
print(f"  -> Least costly: {results_best.iloc[0]['Method']} (${results_best.iloc[0]['Total_Expected_Cost_Lower']:.2f}/test)")
print()

# Worst case (upper bound)
print("WORST CASE (Upper bound error costs):")
results_worst = results_df[['Method', 'Total_Expected_Cost_Upper', 'Correct_Rate']].copy()
results_worst = results_worst.sort_values('Total_Expected_Cost_Upper')
print(results_worst.to_string(index=False))
print(f"  -> Least costly: {results_worst.iloc[0]['Method']} (${results_worst.iloc[0]['Total_Expected_Cost_Upper']:.2f}/test)")
print()

# Check if ranking changes
best_ranking = list(results_best['Method'])
worst_ranking = list(results_worst['Method'])
base_ranking = list(results_df.sort_values('Total_Expected_Cost_Mid')['Method'])

print("Sensitivity of Rankings:")
if best_ranking == base_ranking == worst_ranking:
    print("  [OK] Rankings remain CONSISTENT across all scenarios")
    print("  -> Conclusions are ROBUST to uncertainty in error costs")
else:
    print("  [WARNING] Rankings CHANGE across scenarios:")
    print(f"    Best case order: {' > '.join(best_ranking)}")
    print(f"    Base case order: {' > '.join(base_ranking)}")
    print(f"    Worst case order: {' > '.join(worst_ranking)}")
print()

# =============================================================================
# PART 9: BUDGET IMPACT ANALYSIS
# =============================================================================
print("PART 9: Budget Impact Analysis")
print("-"*80)

annual_volume = methods_df['Expected_number_of_test'].iloc[0]
print(f"Annual testing volume: {annual_volume:,} tests")
print()

print("Annual Budget Requirements (Base Case):")
budget_results = []
for idx, row in results_df.iterrows():
    annual_cost = row['Total_Expected_Cost_Mid'] * annual_volume
    budget_results.append({
        'Method': row['Method'],
        'Annual_Cost_USD': annual_cost
    })
    print(f"  {row['Method']:15s}: ${annual_cost:,.2f}")

budget_df = pd.DataFrame(budget_results)
print()

# Calculate savings compared to each method
print("Potential Annual Savings (compared to each method):")
savings_results = []
for i in range(len(results_df)):
    reference = results_df.iloc[i]
    
    for j in range(len(results_df)):
        if i != j:
            alternative = results_df.iloc[j]
            savings_per_test = reference['Total_Expected_Cost_Mid'] - alternative['Total_Expected_Cost_Mid']
            annual_savings = savings_per_test * annual_volume
            
            savings_results.append({
                'Reference_Method': reference['Method'],
                'Alternative_Method': alternative['Method'],
                'Savings_per_Test_USD': savings_per_test,
                'Annual_Savings_USD': annual_savings
            })
            
            if annual_savings > 0:
                print(f"  {reference['Method']:15s} -> {alternative['Method']:15s}: Save ${annual_savings:,.2f}/year (${savings_per_test:.2f}/test)")
            else:
                print(f"  {reference['Method']:15s} -> {alternative['Method']:15s}: Cost ${abs(annual_savings):,.2f}/year more (${abs(savings_per_test):.2f}/test)")

savings_df = pd.DataFrame(savings_results)
print()

# =============================================================================
# PART 10: COMPREHENSIVE SUMMARY
# =============================================================================
print("\n")
print("="*80)
print("COMPREHENSIVE SUMMARY")
print("="*80)
print()

# Find best methods by different criteria
best_unit_cost = results_df.loc[results_df['Unit_Cost_USD'].idxmin()]
best_accuracy = results_df.loc[results_df['Correct_Rate'].idxmax()]
best_expected_cost = results_df.loc[results_df['Total_Expected_Cost_Mid'].idxmin()]
best_cost_per_correct = results_df.loc[results_df['Cost_per_Correct_ID'].idxmin()]
fastest_turnaround = results_df.loc[results_df['Turnaround_Hours'].idxmin()]

print("BEST METHODS BY CRITERIA:")
print(f"  Lowest unit cost:                {best_unit_cost['Method']} (${best_unit_cost['Unit_Cost_USD']:.2f}/test)")
print(f"  Highest accuracy:                {best_accuracy['Method']} ({best_accuracy['Correct_Rate']*100:.2f}% correct)")
print(f"  Lowest expected cost (w/ errors): {best_expected_cost['Method']} (${best_expected_cost['Total_Expected_Cost_Mid']:.2f}/test)")
print(f"  Best cost per correct ID:         {best_cost_per_correct['Method']} (${best_cost_per_correct['Cost_per_Correct_ID']:.2f})")
print(f"  Fastest turnaround:               {fastest_turnaround['Method']} ({fastest_turnaround['Turnaround_Hours']:.0f} hours)")
print()

print("KEY FINDINGS:")
print()

# Finding 1: Direct cost vs. expected cost
print("1. Impact of Error Costs:")
for idx, row in results_df.iterrows():
    error_cost_contribution = ((row['Total_Expected_Cost_Mid'] - row['Unit_Cost_USD']) / row['Total_Expected_Cost_Mid']) * 100
    print(f"   {row['Method']:15s}: Error costs add ${row['Total_Expected_Cost_Mid'] - row['Unit_Cost_USD']:.2f} " +
          f"({error_cost_contribution:.1f}% of total expected cost)")
print()

# Finding 2: Cost range
cost_range = results_df['Total_Expected_Cost_Mid'].max() - results_df['Total_Expected_Cost_Mid'].min()
print(f"2. Cost Range: ${cost_range:.2f} between least and most expensive methods")
print(f"   ({(cost_range/results_df['Total_Expected_Cost_Mid'].min())*100:.1f}% difference)")
print()

# Finding 3: Accuracy range
acc_range = results_df['Correct_Rate'].max() - results_df['Correct_Rate'].min()
print(f"3. Accuracy Range: {acc_range*100:.2f}% between least and most accurate methods")
print(f"   ({(acc_range/results_df['Correct_Rate'].min())*100:.1f}% relative difference)")
print()

# Finding 4: Dominance
dominated_methods = results_df_sorted[results_df_sorted['Dominated']]['Method'].tolist()
if dominated_methods:
    print(f"4. Dominated Methods (should NOT be used): {', '.join(dominated_methods)}")
else:
    print("4. No methods are strictly dominated")
print()

# Finding 5: Overall recommendation
print("5. OVERALL RECOMMENDATION:")
print(f"   -> Most cost-effective method: {best_expected_cost['Method']}")
print(f"     - Expected cost: ${best_expected_cost['Total_Expected_Cost_Mid']:.2f}/test")
print(f"     - Accuracy: {best_expected_cost['Correct_Rate']*100:.2f}%")
print(f"     - Turnaround: {best_expected_cost['Turnaround_Hours']:.0f} hours")
print(f"     - Annual cost for {annual_volume:,} tests: ${best_expected_cost['Total_Expected_Cost_Mid']*annual_volume:,.2f}")
print()

# =============================================================================
# SAVE RESULTS
# =============================================================================
print("="*80)
print("SAVING RESULTS")
print("="*80)
print()

# Save all results to CSV
results_df.to_csv(os.path.join(output_dir, "detailed_results.csv"), index=False)
print(f"[OK] Detailed results saved to: {os.path.join(output_dir, 'detailed_results.csv')}")

results_df_sorted.to_csv(os.path.join(output_dir, "dominance_analysis.csv"), index=False)
print(f"[OK] Dominance analysis saved to: {os.path.join(output_dir, 'dominance_analysis.csv')}")

icers_df.to_csv(os.path.join(output_dir, "icers_correct_id.csv"), index=False)
print(f"[OK] ICERs (correct ID) saved to: {os.path.join(output_dir, 'icers_correct_id.csv')}")

qaly_icers_df.to_csv(os.path.join(output_dir, "icers_qaly.csv"), index=False)
print(f"[OK] ICERs (QALY) saved to: {os.path.join(output_dir, 'icers_qaly.csv')}")

nmb_all_df.to_csv(os.path.join(output_dir, "net_monetary_benefit.csv"), index=False)
print(f"[OK] Net Monetary Benefit saved to: {os.path.join(output_dir, 'net_monetary_benefit.csv')}")

budget_df.to_csv(os.path.join(output_dir, "budget_impact.csv"), index=False)
print(f"[OK] Budget impact saved to: {os.path.join(output_dir, 'budget_impact.csv')}")

savings_df.to_csv(os.path.join(output_dir, "potential_savings.csv"), index=False)
print(f"[OK] Potential savings saved to: {os.path.join(output_dir, 'potential_savings.csv')}")

# Create summary table for publication
summary_table = results_df[[
    'Method',
    'Unit_Cost_USD',
    'Correct_Rate',
    'VM_Error_Rate',
    'Major_Error_Rate',
    'Total_Expected_Cost_Mid',
    'Cost_per_Correct_ID',
    'Expected_QALY_Loss',
    'Turnaround_Hours'
]].copy()

summary_table.columns = [
    'Method',
    'Unit Cost (USD)',
    'Accuracy (%)',
    'VM Error Rate (%)',
    'Major Error Rate (%)',
    'Expected Cost (USD)',
    'Cost per Correct ID (USD)',
    'Expected QALY Loss',
    'Turnaround (hours)'
]

# Format percentages
summary_table['Accuracy (%)'] = (summary_table['Accuracy (%)'] * 100).round(2)
summary_table['VM Error Rate (%)'] = (summary_table['VM Error Rate (%)'] * 100).round(4)
summary_table['Major Error Rate (%)'] = (summary_table['Major Error Rate (%)'] * 100).round(4)

summary_table.to_csv(os.path.join(output_dir, "summary_table_for_publication.csv"), index=False)
print(f"[OK] Publication-ready summary saved to: {os.path.join(output_dir, 'summary_table_for_publication.csv')}")

print()
print("="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print()
print(f"All results saved to: {output_dir}/")
print()

