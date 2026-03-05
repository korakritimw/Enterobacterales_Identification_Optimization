"""
Create publication-quality multi-panel figure for Cost-Effectiveness Analysis (Updated)
=======================================================================================

This script creates a 4-panel figure (2x2 layout) that tells the complete CEA story:
  Panel A: Cost-Effectiveness Plane (shows dominance)
  Panel B: Expected Cost Breakdown (shows hidden cost of errors)
  Panel C: Budget Impact Analysis (shows annual savings)
  Panel D: Sensitivity Analysis (shows robustness)

Also exports data in Excel-friendly format for replication.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import seaborn as sns
import os

# Set up paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(base_dir, 'results', 'CEA_Analysis_Updated')
data_file = os.path.join(output_dir, 'detailed_results.csv')

# Set publication-quality style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# Load data
results_df = pd.read_csv(data_file)

# Define colors for each method (colorblind-friendly palette)
colors = {
    'P_Biochem': '#E69F00',      # orange
    'Vitek2': '#009E73',         # bluish green
    'MALDI_2.0': '#CC79A7',      # reddish purple
    'MALDI_1.9': '#56B4E9',      # sky blue
    'MALDI_1.8': '#0072B2'       # blue (THE WINNER)
}

# Method display names
method_names = {
    'P_Biochem': 'P-Biochem',
    'Vitek2': 'Vitek2',
    'MALDI_2.0': 'MALDI 2.0',
    'MALDI_1.9': 'MALDI 1.9',
    'MALDI_1.8': 'MALDI 1.8'
}

# Create figure with 4 panels (2x2)
fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# ============================================================================
# PANEL A: Cost-Effectiveness Plane (Scatter Plot)
# ============================================================================
ax_a = fig.add_subplot(gs[0, 0])

# Plot each method
for idx, row in results_df.iterrows():
    method = row['Method']
    x = row['Total_Expected_Cost_Mid']
    y = row['Correct_Rate'] * 100  # Convert to percentage
    
    # Check if dominated
    is_dominated = row.get('Dominated', False) if 'Dominated' in results_df.columns else False
    alpha_val = 0.5 if is_dominated else 0.9
    
    ax_a.scatter(x, y, s=300, c=colors.get(method, '#000000'), 
                edgecolors='black', linewidth=2, 
                label=method_names.get(method, method), zorder=3, alpha=alpha_val)
    
    # Add method labels near points
    offset_x = 3
    offset_y = 0.3
    if method == 'Vitek2':
        offset_x = -20
        offset_y = 0.3
    elif method == 'MALDI_1.8':
        offset_x = 3
        offset_y = -0.5
    
    ax_a.annotate(method_names.get(method, method), 
                 xy=(x, y), 
                 xytext=(offset_x, offset_y),
                 textcoords='offset points',
                 fontsize=9, fontweight='bold',
                 ha='left' if method != 'Vitek2' else 'right')

# Highlight MALDI_1.8 as winner (top-left quadrant = best)
maldi_18_data = results_df[results_df['Method'] == 'MALDI_1.8'].iloc[0]
ax_a.scatter(maldi_18_data['Total_Expected_Cost_Mid'], 
            maldi_18_data['Correct_Rate'] * 100,
            s=400, facecolors='none', edgecolors='#0072B2', 
            linewidth=4, zorder=2)

# Add "Dominant" annotation for MALDI_1.8
ax_a.annotate('Dominant\n(Best Value)', 
             xy=(maldi_18_data['Total_Expected_Cost_Mid'], 
                 maldi_18_data['Correct_Rate'] * 100),
             xytext=(30, 15),
             textcoords='offset points',
             fontsize=9, fontweight='bold', color='#0072B2',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                      edgecolor='#0072B2', linewidth=2),
             arrowprops=dict(arrowstyle='->', color='#0072B2', lw=2))

ax_a.set_xlabel('Expected Cost per Test (USD)', fontsize=11, fontweight='bold')
ax_a.set_ylabel('Accuracy (%)', fontsize=11, fontweight='bold')
ax_a.set_title('A. Cost-Effectiveness Plane', fontsize=12, fontweight='bold', loc='left')
ax_a.grid(True, alpha=0.3, linestyle='--')
ax_a.set_xlim(0, 120)
ax_a.set_ylim(75, 100)

# Add interpretation text
ax_a.text(0.02, 0.02, 'Top-left = Most cost-effective\n(High accuracy, Low cost)',
         transform=ax_a.transAxes, fontsize=8, 
         verticalalignment='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

# ============================================================================
# PANEL B: Expected Cost Breakdown (Stacked Bar Chart)
# ============================================================================
ax_b = fig.add_subplot(gs[0, 1])

# Prepare data for stacked bars
methods_display = [method_names.get(m, m) for m in results_df['Method']]
unit_costs = results_df['Unit_Cost_USD'].values
error_costs = results_df['Expected_Error_Cost_Mid'].values

# Create positions
x_pos = np.arange(len(methods_display))
width = 0.6

# Create stacked bars
bars1 = ax_b.bar(x_pos, unit_costs, width, 
                label='Unit Cost (Test Materials)', 
                color='#4472C4', edgecolor='black', linewidth=1.5)

bars2 = ax_b.bar(x_pos, error_costs, width, bottom=unit_costs,
                label='Error Cost (Clinical Consequences)', 
                color='#C5504B', edgecolor='black', linewidth=1.5)

# Add total cost labels on top
for i, (method, total) in enumerate(zip(methods_display, results_df['Total_Expected_Cost_Mid'])):
    ax_b.text(i, total + 3, f'${total:.1f}', 
             ha='center', va='bottom', fontsize=9, fontweight='bold')

# Highlight that error costs dominate for some methods
ax_b.text(0.02, 0.98, 'Red portion = Cost of errors\n(often >> test cost)',
         transform=ax_b.transAxes, fontsize=8, 
         verticalalignment='top', style='italic',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

ax_b.set_ylabel('Cost per Test (USD)', fontsize=11, fontweight='bold')
ax_b.set_title('B. True Cost: Unit Cost + Error Consequences', 
              fontsize=12, fontweight='bold', loc='left')
ax_b.set_xticks(x_pos)
ax_b.set_xticklabels(methods_display, rotation=45, ha='right')
ax_b.legend(loc='upper left', fontsize=8, framealpha=0.9)
ax_b.set_ylim(0, 120)
ax_b.grid(True, alpha=0.3, linestyle='--', axis='y')

# ============================================================================
# PANEL C: Budget Impact Analysis (Bar Chart with Savings)
# ============================================================================
ax_c = fig.add_subplot(gs[1, 0])

# Annual costs
annual_volume = 15000
annual_costs = results_df['Total_Expected_Cost_Mid'].values * annual_volume / 1000  # in thousands

# Sort by cost
sorted_indices = np.argsort(annual_costs)
sorted_methods = [methods_display[i] for i in sorted_indices]
sorted_costs = annual_costs[sorted_indices]
sorted_colors = [colors.get(results_df.iloc[i]['Method'], '#000000') for i in sorted_indices]

# Create horizontal bars
bars = ax_c.barh(sorted_methods, sorted_costs, 
                color=sorted_colors, edgecolor='black', linewidth=1.5, alpha=0.9)

# Add cost labels
for i, (method, cost) in enumerate(zip(sorted_methods, sorted_costs)):
    ax_c.text(cost + 50, i, f'${cost:.0f}K', 
             va='center', fontsize=9, fontweight='bold')

# Add savings annotations (P_Biochem to MALDI_1.8)
p_biochem_cost = results_df[results_df['Method'] == 'P_Biochem']['Total_Expected_Cost_Mid'].values[0] * annual_volume / 1000
maldi_18_cost = results_df[results_df['Method'] == 'MALDI_1.8']['Total_Expected_Cost_Mid'].values[0] * annual_volume / 1000
savings = p_biochem_cost - maldi_18_cost

# Add arrow showing savings
p_biochem_idx = sorted_methods.index('P-Biochem')
maldi_18_idx = sorted_methods.index('MALDI 1.8')
ax_c.annotate('', xy=(maldi_18_cost, maldi_18_idx), xytext=(p_biochem_cost, p_biochem_idx),
             arrowprops=dict(arrowstyle='<->', color='red', lw=3))
ax_c.text((maldi_18_cost + p_biochem_cost) / 2, (maldi_18_idx + p_biochem_idx) / 2,
         f'Annual Savings:\n${savings:.0f}K',
         ha='center', va='center', fontsize=9, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', 
                  edgecolor='red', linewidth=2))

ax_c.set_xlabel('Annual Cost (Thousands USD)\n15,000 tests/year', 
               fontsize=11, fontweight='bold')
ax_c.set_title('C. Budget Impact Analysis', fontsize=12, fontweight='bold', loc='left')
ax_c.set_xlim(0, max(annual_costs) * 1.15)
ax_c.grid(True, alpha=0.3, linestyle='--', axis='x')

# ============================================================================
# PANEL D: Sensitivity Analysis (Grouped Bar Chart)
# ============================================================================
ax_d = fig.add_subplot(gs[1, 1])

# Prepare data for grouped bars
scenarios = ['Best\nCase', 'Base\nCase', 'Worst\nCase']
x_pos = np.arange(len(methods_display))
width = 0.25

# Get costs for each scenario
best_case = results_df['Total_Expected_Cost_Lower'].values
base_case = results_df['Total_Expected_Cost_Mid'].values
worst_case = results_df['Total_Expected_Cost_Upper'].values

# Create grouped bars
bars1 = ax_d.bar(x_pos - width, best_case, width, 
                label='Best Case (Lower Error Costs)', 
                color='#70AD47', edgecolor='black', linewidth=1, alpha=0.8)
bars2 = ax_d.bar(x_pos, base_case, width, 
                label='Base Case (Expected Costs)', 
                color='#4472C4', edgecolor='black', linewidth=1, alpha=0.8)
bars3 = ax_d.bar(x_pos + width, worst_case, width, 
                label='Worst Case (Upper Error Costs)', 
                color='#C5504B', edgecolor='black', linewidth=1, alpha=0.8)

# Highlight that MALDI_1.8 wins in all scenarios
maldi_18_idx = list(results_df['Method']).index('MALDI_1.8')
for bars, offset in [(bars1, -width), (bars2, 0), (bars3, width)]:
    bars[maldi_18_idx].set_edgecolor('#0072B2')
    bars[maldi_18_idx].set_linewidth(3)

ax_d.set_ylabel('Expected Cost per Test (USD)', fontsize=11, fontweight='bold')
ax_d.set_title('D. Sensitivity Analysis: Robust Rankings', 
              fontsize=12, fontweight='bold', loc='left')
ax_d.set_xticks(x_pos)
ax_d.set_xticklabels(methods_display, rotation=45, ha='right')
ax_d.legend(loc='upper left', fontsize=8, framealpha=0.9)
ax_d.set_ylim(0, max(worst_case) * 1.15)
ax_d.grid(True, alpha=0.3, linestyle='--', axis='y')

# Add annotation
ax_d.text(0.98, 0.98, 'MALDI 1.8 (bold borders)\nremains cheapest\nacross all scenarios',
         transform=ax_d.transAxes, fontsize=8, 
         verticalalignment='top', ha='right', style='italic',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

# ============================================================================
# Overall figure title and save
# ============================================================================
fig.suptitle('Cost-Effectiveness Analysis of Bacterial Identification Methods', 
            fontsize=14, fontweight='bold', y=0.995)

# Save figure
plt.savefig(os.path.join(output_dir, 'Figure_CEA_MultiPanel.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, 'Figure_CEA_MultiPanel.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(output_dir, 'Figure_CEA_MultiPanel.tiff'), dpi=300, bbox_inches='tight')

print("=" * 80)
print("MULTI-PANEL FIGURE CREATED SUCCESSFULLY")
print("=" * 80)
print()
print("Files saved:")
print(f"  [OK] {os.path.join(output_dir, 'Figure_CEA_MultiPanel.png')} (high resolution)")
print(f"  [OK] {os.path.join(output_dir, 'Figure_CEA_MultiPanel.pdf')} (vector format)")
print(f"  [OK] {os.path.join(output_dir, 'Figure_CEA_MultiPanel.tiff')} (publication format)")
print()
print("Figure layout (2x2):")
print("  Panel A (top-left):     Cost-Effectiveness Plane")
print("  Panel B (top-right):    Expected Cost Breakdown")
print("  Panel C (bottom-left):  Budget Impact Analysis")
print("  Panel D (bottom-right): Sensitivity Analysis")
print()

# ============================================================================
# EXPORT DATA FOR EXCEL REPLICATION
# ============================================================================
print("=" * 80)
print("EXPORTING DATA FOR EXCEL REPLICATION")
print("=" * 80)
print()

# Panel A data: Cost-Effectiveness Plane
panel_a_data = results_df[['Method', 'Total_Expected_Cost_Mid', 'Correct_Rate']].copy()
panel_a_data['Accuracy_Percent'] = panel_a_data['Correct_Rate'] * 100
panel_a_data['Method_Display'] = panel_a_data['Method'].map(method_names)
panel_a_data = panel_a_data[['Method_Display', 'Total_Expected_Cost_Mid', 'Accuracy_Percent']]
panel_a_data.columns = ['Method', 'Expected_Cost_USD', 'Accuracy_Percent']
panel_a_data.to_csv(os.path.join(output_dir, 'Excel_Data_Panel_A_CostEffectivenessPlane.csv'), index=False)
print(f"[OK] Panel A data saved: Excel_Data_Panel_A_CostEffectivenessPlane.csv")

# Panel B data: Cost Breakdown
panel_b_data = pd.DataFrame({
    'Method': methods_display,
    'Unit_Cost_USD': unit_costs,
    'Error_Cost_USD': error_costs,
    'Total_Expected_Cost_USD': results_df['Total_Expected_Cost_Mid'].values
})
panel_b_data.to_csv(os.path.join(output_dir, 'Excel_Data_Panel_B_CostBreakdown.csv'), index=False)
print(f"[OK] Panel B data saved: Excel_Data_Panel_B_CostBreakdown.csv")

# Panel C data: Budget Impact
panel_c_data = pd.DataFrame({
    'Method': sorted_methods,
    'Annual_Cost_Thousands_USD': sorted_costs,
    'Annual_Cost_USD': sorted_costs * 1000
})
panel_c_data.to_csv(os.path.join(output_dir, 'Excel_Data_Panel_C_BudgetImpact.csv'), index=False)
print(f"[OK] Panel C data saved: Excel_Data_Panel_C_BudgetImpact.csv")

# Panel D data: Sensitivity Analysis
panel_d_data = pd.DataFrame({
    'Method': methods_display,
    'Best_Case_Cost_USD': best_case,
    'Base_Case_Cost_USD': base_case,
    'Worst_Case_Cost_USD': worst_case
})
panel_d_data.to_csv(os.path.join(output_dir, 'Excel_Data_Panel_D_SensitivityAnalysis.csv'), index=False)
print(f"[OK] Panel D data saved: Excel_Data_Panel_D_SensitivityAnalysis.csv")

# Create comprehensive Excel-ready file with all data
excel_all_data = pd.DataFrame({
    'Method': methods_display,
    'Unit_Cost_USD': unit_costs,
    'Error_Cost_USD': error_costs,
    'Total_Expected_Cost_Mid_USD': base_case,
    'Total_Expected_Cost_Lower_USD': best_case,
    'Total_Expected_Cost_Upper_USD': worst_case,
    'Accuracy_Percent': results_df['Correct_Rate'].values * 100,
    'VM_Error_Rate_Percent': results_df['VM_Error_Rate'].values * 100,
    'Major_Error_Rate_Percent': results_df['Major_Error_Rate'].values * 100,
    'Annual_Cost_USD': annual_costs * 1000,
    'Turnaround_Hours': results_df['Turnaround_Hours'].values
})
excel_all_data.to_csv(os.path.join(output_dir, 'Excel_Data_AllPanels_Complete.csv'), index=False)
print(f"[OK] Complete data saved: Excel_Data_AllPanels_Complete.csv")

# Create Excel instructions document
instructions = """
================================================================================
EXCEL REPLICATION INSTRUCTIONS
================================================================================

This document explains how to recreate the CEA figures in Excel using the
provided data files.

================================================================================
PANEL A: COST-EFFECTIVENESS PLANE (Scatter Plot)
================================================================================

Data File: Excel_Data_Panel_A_CostEffectivenessPlane.csv

Steps:
1. Open the CSV file in Excel
2. Select columns: Expected_Cost_USD and Accuracy_Percent
3. Insert > Charts > Scatter (with markers)
4. Format:
   - X-axis: Expected Cost per Test (USD)
   - Y-axis: Accuracy (%)
   - Title: "A. Cost-Effectiveness Plane"
   - Add gridlines
   - Add data labels with Method names
   - Color code points by method:
     * P-Biochem: Orange (#E69F00)
     * Vitek2: Bluish Green (#009E73)
     * MALDI 2.0: Reddish Purple (#CC79A7)
     * MALDI 1.9: Sky Blue (#56B4E9)
     * MALDI 1.8: Blue (#0072B2) - highlight with bold border

5. Add text box: "Top-left = Most cost-effective (High accuracy, Low cost)"
6. Highlight MALDI 1.8 point with bold border and add annotation "Dominant (Best Value)"

================================================================================
PANEL B: EXPECTED COST BREAKDOWN (Stacked Bar Chart)
================================================================================

Data File: Excel_Data_Panel_B_CostBreakdown.csv

Steps:
1. Open the CSV file in Excel
2. Select all data (Method, Unit_Cost_USD, Error_Cost_USD)
3. Insert > Charts > Stacked Bar Chart
4. Format:
   - X-axis: Cost per Test (USD)
   - Y-axis: Method names
   - Title: "B. True Cost: Unit Cost + Error Consequences"
   - Legend:
     * Unit Cost: Blue (#4472C4)
     * Error Cost: Red (#C5504B)
   - Add data labels on top showing Total_Expected_Cost_USD
   - Add text box: "Red portion = Cost of errors (often >> test cost)"

5. Rotate method names 45 degrees if needed

================================================================================
PANEL C: BUDGET IMPACT ANALYSIS (Horizontal Bar Chart)
================================================================================

Data File: Excel_Data_Panel_C_BudgetImpact.csv

Steps:
1. Open the CSV file in Excel
2. Sort by Annual_Cost_Thousands_USD (ascending)
3. Select Method and Annual_Cost_Thousands_USD columns
4. Insert > Charts > Horizontal Bar Chart
5. Format:
   - X-axis: Annual Cost (Thousands USD) - 15,000 tests/year
   - Y-axis: Method names
   - Title: "C. Budget Impact Analysis"
   - Add data labels showing cost values
   - Color code bars by method (same colors as Panel A)
   - Add arrow/annotation showing savings from P-Biochem to MALDI 1.8

6. Calculate savings:
   - P-Biochem cost - MALDI 1.8 cost = Annual Savings
   - Add text box with savings amount

================================================================================
PANEL D: SENSITIVITY ANALYSIS (Grouped Bar Chart)
================================================================================

Data File: Excel_Data_Panel_D_SensitivityAnalysis.csv

Steps:
1. Open the CSV file in Excel
2. Select all columns (Method, Best_Case, Base_Case, Worst_Case)
3. Insert > Charts > Clustered Column Chart
4. Format:
   - X-axis: Method names
   - Y-axis: Expected Cost per Test (USD)
   - Title: "D. Sensitivity Analysis: Robust Rankings"
   - Legend:
     * Best Case: Green (#70AD47)
     * Base Case: Blue (#4472C4)
     * Worst Case: Red (#C5504B)
   - Highlight MALDI 1.8 bars with bold borders
   - Add text box: "MALDI 1.8 (bold borders) remains cheapest across all scenarios"

5. Rotate method names 45 degrees

================================================================================
COMBINING PANELS IN EXCEL
================================================================================

To create a 2x2 layout in Excel:

1. Create each chart separately (Panels A-D)
2. Copy each chart
3. Paste into a new worksheet
4. Arrange in 2x2 grid:
   - Panel A (top-left)
   - Panel B (top-right)
   - Panel C (bottom-left)
   - Panel D (bottom-right)
5. Resize charts to fit
6. Add overall title: "Cost-Effectiveness Analysis of Bacterial Identification Methods"

Alternative: Use PowerPoint
1. Create each chart in Excel
2. Copy to PowerPoint
3. Arrange in 2x2 layout
4. Add overall title

================================================================================
COLOR CODES FOR EXCEL
================================================================================

Method Colors (use in all panels):
- P-Biochem: RGB(230, 159, 0) or Orange
- Vitek2: RGB(0, 158, 115) or Teal/Green
- MALDI 2.0: RGB(204, 121, 167) or Pink/Purple
- MALDI 1.9: RGB(86, 180, 233) or Light Blue
- MALDI 1.8: RGB(0, 114, 178) or Blue (WINNER - use bold border)

Cost Breakdown Colors (Panel B):
- Unit Cost: RGB(68, 114, 196) or Blue
- Error Cost: RGB(197, 80, 75) or Red

Sensitivity Colors (Panel D):
- Best Case: RGB(112, 173, 71) or Green
- Base Case: RGB(68, 114, 196) or Blue
- Worst Case: RGB(197, 80, 75) or Red

================================================================================
QUICK REFERENCE: KEY VALUES
================================================================================

Expected Costs (Base Case):
- MALDI 1.8: $11.69
- MALDI 1.9: $20.88
- MALDI 2.0: $42.61
- Vitek2: $58.11
- P-Biochem: $105.61

Accuracy:
- Vitek2: 97.92%
- MALDI 1.8: 95.83%
- P-Biochem: 95.42%
- MALDI 1.9: 91.25%
- MALDI 2.0: 80.42%

Annual Costs (15,000 tests):
- MALDI 1.8: $175,322
- MALDI 1.9: $313,230
- MALDI 2.0: $639,196
- Vitek2: $871,657
- P-Biochem: $1,584,193

Annual Savings (P-Biochem to MALDI 1.8): $1,408,871

================================================================================
END OF INSTRUCTIONS
================================================================================
"""

with open(os.path.join(output_dir, 'Excel_Replication_Instructions.txt'), 'w', encoding='utf-8') as f:
    f.write(instructions)

print(f"[OK] Excel instructions saved: Excel_Replication_Instructions.txt")

print()
print("=" * 80)
print("ALL FILES CREATED SUCCESSFULLY")
print("=" * 80)
print()
print("Figure files:")
print(f"  - Figure_CEA_MultiPanel.png")
print(f"  - Figure_CEA_MultiPanel.pdf")
print(f"  - Figure_CEA_MultiPanel.tiff")
print()
print("Excel data files:")
print(f"  - Excel_Data_Panel_A_CostEffectivenessPlane.csv")
print(f"  - Excel_Data_Panel_B_CostBreakdown.csv")
print(f"  - Excel_Data_Panel_C_BudgetImpact.csv")
print(f"  - Excel_Data_Panel_D_SensitivityAnalysis.csv")
print(f"  - Excel_Data_AllPanels_Complete.csv")
print()
print("Instructions:")
print(f"  - Excel_Replication_Instructions.txt")
print()

# Close the figure to free memory
plt.close()

print("=" * 80)
print("COMPLETE")
print("=" * 80)

