# Enterobacterales Identification Optimization

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

Python analysis scripts for the manuscript **"Optimizing species identification after antimicrobial pretreatment in Enterobacterales: a prospective cohort study"** accepted in *BMC Microbiology*.

## Study Overview

This prospective cohort study (May-November 2025, Siriraj Hospital, Thailand) evaluated the impact of antimicrobial pretreatment on bacterial identification accuracy for 240 Enterobacterales isolates showing visible pretreatment effects on culture plates.

**Key Findings:**
- MALDI-TOF MS with reduced cutoff (1.8) achieved 93.8% success rate and 99.6% accuracy
- Vitek2 maintained high accuracy comparable to benchmark
- Optimized MALDI-TOF MS was most cost-effective method
- Pretreatment effect associated with lower antimicrobial resistance prevalence

**Methods Compared:**
1. Conventional phenotypic (biochemical tests with/without empirical priors)
2. Vitek2 automated system
3. MALDI-TOF MS (cutoffs: 2.0, 1.9, 1.8)

## Citation

If you use this code, please cite:

```bibtex
@article{kijsinthopchai2026enterobacterales,
  title={Optimizing species identification after antimicrobial pretreatment in Enterobacterales: a prospective cohort study},
  author={Kijsinthopchai, Usa and Song, Yizhe and Yongyod, Samaporn and Wensanthia, Thidarat and Disthaporn, Pensiri and Pimnon, Siriporn and Imwattana, Korakrit},
  journal={BMC Microbiology},
  year={2026},
  note={Accepted for publication}
}
```

## Repository Structure

```
.
├── python_published/
│   ├── analysis2_updated.py              # Main analysis script
│   ├── analysis2_method_comparison_original.py  # Original version
│   ├── Basic scripts/                    # Modular statistical tests
│   │   ├── README.md                     # Detailed script documentation
│   │   ├── task1_diagnostic_test_metrics.py
│   │   ├── task2_two_by_two_comparison.py
│   │   ├── task3_mcnemar_comparison.py
│   │   ├── task4_numerical_comparison.py
│   │   ├── task5_two_group_comparison.py
│   │   └── Task_*/                       # Analysis-specific subdirectories
│   ├── CEA analysis/                     # Cost-effectiveness analysis
│   │   ├── cost_effectiveness_analysis_updated.py
│   │   ├── create_cea_figure_updated.py
│   │   ├── Methods_KI.csv                # Method performance data
│   │   └── RIsk_and_cost_KI.csv          # Cost and QALY parameters
│   ├── data/                             # Data directory (not included)
│   └── results/                          # Output directory
├── notebooks/                            # Jupyter notebooks (see below)
├── docs/
│   └── DATA_STRUCTURE.md                 # Complete data dictionary
├── requirements.txt                      # Python dependencies
├── LICENSE                               # MIT License
├── CITATION.cff                          # Citation metadata
└── README.md                             # This file
```

## Installation

### Prerequisites
- Python 3.9 or higher
- pip package manager

### Setup

1. **Clone the repository:**
```bash
git clone https://github.com/sscien/Enterobacterales_Identification_Optimization.git
cd Enterobacterales_Identification_Optimization
```

2. **Create virtual environment (recommended):**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies:**
```bash
pip install -r requirements.txt
```

## Data Structure

Due to patient privacy considerations, the actual dataset is **not included** in this repository. However, complete data structure documentation is provided in [`docs/DATA_STRUCTURE.md`](docs/DATA_STRUCTURE.md).

### Required Data Files

To run the analysis scripts, you need:

1. **`data/AI_database.csv`** - Main dataset with 33 columns:
   - Isolate identification (ID, plate type, specimen type, final ID)
   - Biochemical test results (14 tests)
   - Identification results from 3 methods (conventional, Vitek2, MALDI-TOF MS)
   - Antimicrobial susceptibility data

2. **`data/CEA/Methods_KI.csv`** - Method performance metrics
3. **`data/CEA/RIsk_and_cost_KI.csv`** - Cost-effectiveness parameters

See [`docs/DATA_STRUCTURE.md`](docs/DATA_STRUCTURE.md) for complete column descriptions, valid values, and example data rows.

### Data Access

Researchers interested in accessing the dataset should contact:

**Korakrit Imwattana**
Department of Microbiology, Faculty of Medicine Siriraj Hospital
Mahidol University, Bangkok, Thailand
Email: korakrit.imw@mahidol.ac.th

## Usage

### Main Analysis

The primary analysis script compares identification methods across Original, Subculture, and Positive (benchmark) groups:

```bash
cd python_published
python analysis2_updated.py
```

**Outputs:**
- `results/analysis2_success_rates.txt` - Success rates for each method
- `results/analysis2_accuracy.txt` - Accuracy among successful identifications
- `results/analysis2_comparisons.txt` - Pairwise statistical comparisons
- `results/analysis2_maldi_cutoff_optimization.txt` - ROC analysis for MALDI cutoffs

**Key analyses performed:**
- Diagnostic test metrics (sensitivity, specificity, accuracy)
- McNemar tests for paired comparisons
- Holm-Bonferroni correction for multiple comparisons
- ROC curve analysis for MALDI-TOF cutoff optimization

### Cost-Effectiveness Analysis

```bash
cd python_published/CEA\ analysis
python cost_effectiveness_analysis_updated.py
```

**Outputs:**
- `results/CEA_Analysis_Updated/cea_summary.txt` - Cost-effectiveness metrics
- `results/CEA_Analysis_Updated/dominance_analysis.txt` - Method dominance
- `results/CEA_Analysis_Updated/icer_analysis.txt` - Incremental cost-effectiveness ratios
- `results/CEA_Analysis_Updated/sensitivity_analysis.txt` - Parameter sensitivity

**Generate CEA figures:**
```bash
python create_cea_figure_updated.py
```

**Outputs:**
- `results/CEA_Analysis_Updated/cea_plane.png` - Cost-effectiveness plane (Figure 4)
- `results/CEA_Analysis_Updated/budget_impact.png` - Budget impact analysis (Figure 5)

### Modular Statistical Tests

Individual statistical tests can be run using scripts in `Basic scripts/`. See [`python_published/Basic scripts/README.md`](python_published/Basic%20scripts/README.md) for detailed documentation.

**Example - Diagnostic test metrics:**
```bash
cd python_published/Basic\ scripts/Task_1
python task1_diagnostic_test_metrics.py
```

## Reproducing Manuscript Figures

### Figure 2A: MALDI-TOF Score Distributions

Box plots comparing MALDI-TOF log scores across groups:

```python
# See notebooks/figure2A_maldi_scores.ipynb
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('data/AI_database.csv')
fig, ax = plt.subplots(figsize=(8, 6))
sns.boxplot(data=df, x='Plate', y='MALDI_score', ax=ax)
ax.set_ylabel('MALDI-TOF Log Score')
ax.axhline(y=2.0, color='r', linestyle='--', label='Standard cutoff')
ax.axhline(y=1.8, color='g', linestyle='--', label='Optimized cutoff')
plt.legend()
plt.savefig('results/figure2A.png', dpi=300, bbox_inches='tight')
```

### Figure 2B: ROC Curves

ROC analysis for MALDI-TOF cutoff optimization:

```python
# See notebooks/figure2B_roc_curves.ipynb
from sklearn.metrics import roc_curve, auc

# Calculate ROC for Original plate
fpr, tpr, thresholds = roc_curve(y_true, y_scores)
roc_auc = auc(fpr, tpr)

plt.plot(fpr, tpr, label=f'AUC = {roc_auc:.3f}')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve - MALDI-TOF Cutoff Optimization')
```

## Jupyter Notebooks

Interactive notebooks demonstrating the analysis workflow with embedded figures:

1. **`notebooks/figure2A_maldi_scores.ipynb`** - Figure 2A: MALDI-TOF score distributions (box plots)
2. **`notebooks/figure2B_roc_curves.ipynb`** - Figure 2B: ROC curve analysis for cutoff optimization

## Statistical Methods

### Key Statistical Approaches

- **Wilson score intervals:** Confidence intervals for proportions (more accurate than normal approximation)
- **McNemar's test:** Paired comparisons between Original and Subculture plates
- **Holm-Bonferroni correction:** Multiple comparison adjustment
- **Mann-Whitney U test:** Non-parametric comparison of MALDI scores
- **ROC analysis:** Optimal cutoff determination
- **Bayesian cost-effectiveness:** Monte Carlo simulation with 10,000 iterations

### Software Versions

- Python: 3.13.9
- pandas: 1.5.3+
- numpy: 1.23.5+
- scipy: 1.10.0+
- statsmodels: 0.14.0+
- matplotlib: 3.7.0+
- seaborn: 0.12.2+

## Results Summary

### Identification Success Rates (Original Plate)

| Method | Cutoff | Success Rate | 95% CI |
|--------|--------|--------------|--------|
| P_Biochem | ≥0.85 | 97.1% | 94.2-98.7% |
| Vitek2 | ≥85% | 97.5% | 94.7-99.0% |
| MALDI (2.0) | ≥2.0 | 75.0% | 69.1-80.3% |
| MALDI (1.9) | ≥1.9 | 86.7% | 81.7-90.7% |
| MALDI (1.8) | ≥1.8 | 93.8% | 89.9-96.4% |

### Accuracy (Among Successful Identifications)

All methods achieved >99% accuracy among isolates passing their respective cutoffs.

### Cost-Effectiveness

- **Most cost-effective:** MALDI-TOF MS with cutoff 1.8
- **Least cost-effective:** Conventional phenotypic method
- **ICER:** MALDI 1.8 dominates other methods (lower cost, higher effectiveness)

## Contributing

This repository contains analysis code for a published study. For questions or issues:

1. Open a GitHub issue
2. Contact the corresponding author (see above)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Department of Microbiology, Faculty of Medicine Siriraj Hospital, Mahidol University
- All study participants and clinical staff who contributed specimens
- Funding sources: [To be added]

## Contact

**Corresponding Author:**
Korakrit Imwattana
Department of Microbiology
Faculty of Medicine Siriraj Hospital
Mahidol University, Bangkok, Thailand
Email: korakrit.imw@mahidol.ac.th

---

**Last Updated:** March 2026
