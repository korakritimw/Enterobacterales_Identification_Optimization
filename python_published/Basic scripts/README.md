# Basic Scripts

This directory contains modular statistical analysis scripts used throughout the study. Each script is designed to be reusable for different analyses by providing data through CSV templates.

## Script Overview

### Diagnostic Test Metrics

#### `task1_diagnostic_test_metrics.py`
**Purpose:** Calculate sensitivity, specificity, accuracy with 95% Wilson score confidence intervals

**Input:** CSV with columns: `Test_Name`, `TP`, `TN`, `FP`, `FN`

**Output:** Text file with diagnostic metrics for each test

**Used for:** Table S7 (Benchmark performance data)

**Example usage:**
```bash
cd "Basic scripts/Task_1"
python task1_diagnostic_test_metrics.py
```

---

#### `task1_2_diagnostic_test_metrics_comparison.py`
**Purpose:** Compare diagnostic metrics between two groups using proportion tests

**Input:** CSV with columns: `Test_Name`, `Group`, `TP`, `TN`, `FP`, `FN`

**Output:** Text file with pairwise comparisons and p-values

**Used for:** Comparing Original vs Subculture vs Positive groups

---

### Two-by-Two Table Comparisons

#### `task2_two_by_two_comparison.py`
**Purpose:** Chi-square and Fisher's exact tests for 2×2 contingency tables

**Input:** CSV with columns: `Group`, `Category_A`, `Category_B`

**Output:** Text file with chi-square statistic, p-value, effect size (Cramér's V)

**Used for:**
- Table S1: Species distribution comparison
- Table S2: Specimen type distribution
- Resistance prevalence comparisons

**Example usage:**
```bash
cd "Basic scripts/Task_2_1"
python task2_two_by_two_comparison.py
```

**Subdirectories:**
- `Task_2_1/`: Species distribution (study vs benchmark)
- `Task_2_2/`: Specimen type distribution
- `Task_2_3/`: ESBL E. coli prevalence
- `Task_2_4/`: ESBL K. pneumoniae prevalence
- `Task_2_5/`: CRE prevalence

---

### McNemar's Test for Paired Data

#### `task3_mcnemar_comparison.py`
**Purpose:** McNemar's test for paired binary outcomes (Original vs Subculture)

**Input:** CSV with columns: `Isolate_ID`, `Original_Result`, `Subculture_Result`

**Output:** Text file with McNemar statistic, p-value, odds ratio

**Used for:** Comparing identification success rates between Original and Subculture plates

---

### Numerical Comparisons

#### `task4_numerical_comparison.py`
**Purpose:** Mann-Whitney U test for comparing continuous variables between groups

**Input:** CSV with columns: `Group`, `Value`

**Output:** Text file with U statistic, p-value, median values

**Used for:** Comparing MALDI-TOF scores between groups (Figure 2)

---

### Two-Group Comparisons with Multiple Tests

#### `task5_two_group_comparison.py`
**Purpose:** Compare proportions between two groups with Holm-Bonferroni correction

**Input:** CSV with columns: `Group`, `Success`, `Total`

**Output:** Text file with proportion test results and corrected p-values

**Used for:** Resistance prevalence comparisons across multiple organisms

**Subdirectories:**
- `Task_5_1/`: ESBL E. coli comparison
- `Task_5_2/`: ESBL K. pneumoniae comparison
- `Task_5_3/`: CRE E. coli comparison
- `Task_5_4/`: CRE K. pneumoniae comparison

---

## Data Templates

Each script directory contains a `*_data_template.csv` file showing the required format. Copy and populate these templates with your data before running the scripts.

## Common Statistical Methods

All scripts use:
- **Wilson score intervals** for proportion confidence intervals (more accurate than normal approximation)
- **Holm-Bonferroni correction** for multiple comparisons
- **Exact tests** when cell counts are small (Fisher's exact, exact binomial)

## Output Format

All scripts generate timestamped text files with:
- Input data summary
- Statistical test results
- Interpretation notes
- Timestamp and script version

## Dependencies

See `requirements.txt` in repository root:
```
pandas>=1.5.3
numpy>=1.23.5
scipy>=1.10.0
statsmodels>=0.14.0
```

## Notes

- **Subdirectories (Task_X_Y):** These contain the same scripts with different input data for specific analyses. The canonical script versions are in the main Basic scripts/ directory.
- **Output files:** Named with timestamp (e.g., `task1_20251115_01.txt`) to track analysis versions
- **Reproducibility:** All random seeds are set where applicable
