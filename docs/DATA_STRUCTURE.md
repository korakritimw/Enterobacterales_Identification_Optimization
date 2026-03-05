# Data Structure Documentation

## Overview

This document describes the data structure used in the Enterobacterales identification optimization study. Due to patient privacy considerations, the actual dataset is not included in this repository. However, this documentation provides complete information about the data format to understand the analysis scripts.

## Main Dataset: AI_database.csv

### File Structure
- **Format:** CSV (comma-separated values)
- **Encoding:** UTF-8 with BOM
- **Total isolates:** 240 Enterobacterales with antimicrobial pretreatment effect
- **Rows per isolate:** 2 (Original plate + Subculture plate)
- **Total rows:** 480 + 150 benchmark isolates = 630 rows
- **Columns:** 33

### Column Descriptions

#### Identification Columns
| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| `Isolate_ID` | Integer | Unique isolate identifier | 1-240 |
| `Plate` | String | Culture plate type | "Original", "Subculture", "Positive" |
| `Specimen type` | String | Clinical specimen source | "Urine", "Body fluid", "Sputum", "Swab", "Pus", "Tissue" |
| `Final_ID` | String | Reference identification (consensus) | Species name (e.g., "Escherichia coli") |

#### Biochemical Test Results (Conventional Phenotypic Method)
| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| `TSI_KA` | Binary | Triple Sugar Iron: Alkaline slant | 0 (negative), 1 (positive) |
| `TSI_AA` | Binary | Triple Sugar Iron: Acid slant | 0 (negative), 1 (positive) |
| `TSI_Gas` | Binary | Triple Sugar Iron: Gas production | 0 (negative), 1 (positive) |
| `TSI_H2S` | Binary | Triple Sugar Iron: H2S production | 0 (negative), 1 (positive) |
| `Lactose` | Binary | Lactose fermentation (MacConkey) | 0 (negative), 1 (positive) |
| `Indole` | Binary | Indole production | 0 (negative), 1 (positive) |
| `LDM` | Binary | Lysine deaminase | 0 (negative), 1 (positive) |
| `LDC` | Binary | Lysine decarboxylase | 0 (negative), 1 (positive) |
| `Motility` | Binary | Motility test | 0 (negative), 1 (positive) |
| `Urease` | Binary | Urease production | 0 (negative), 1 (positive) |
| `Citrate` | Binary | Citrate utilization | 0 (negative), 1 (positive) |
| `VP` | Binary | Voges-Proskauer test | 0 (negative), 1 (positive) |
| `Malonate` | Binary | Malonate utilization | 0 (negative), 1 (positive) |
| `Mannitol` | Binary | Mannitol fermentation | 0 (negative), 1 (positive) |

#### Conventional Phenotypic Identification Results
| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| `Profile_ID` | Scientific notation | Numeric profile from biochemical tests | e.g., 1.10E+12 |
| `N_biochem_ID` | String | ID without priors (uniform probability) | Species name or "Unidentified" |
| `N_biochem_score` | Float | Posterior probability (no priors) | 0.0 to 1.0 |
| `P_biochem_ID` | String | ID with empirical priors | Species name or "Unidentified" |
| `P_biochem_score` | Float | Posterior probability (with priors) | 0.0 to 1.0 |

#### Vitek2 Automated System Results
| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| `Vitek_ID` | String | Vitek2 identification | Species name or "Unidentified" |
| `Vitek_Bionumber` | String | Vitek2 bionumber code | 16-character alphanumeric |
| `Vitek_score` | Integer | Vitek2 confidence percentage | 0-99 (%) |

#### MALDI-TOF MS Results
| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| `MALDI_ID` | String | MALDI-TOF MS identification | Species name or "Unidentified" |
| `MALDI_score` | Float | MALDI-TOF MS log score | 0.0 to 3.0 |

#### Antimicrobial Susceptibility Testing (AST) Results
| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| `AST_performed` | String | Whether AST was performed | "YES", "NO" |
| `Predicted_AmpC` | String | AmpC β-lactamase producer | "YES" or empty |
| `Predicted_ESBL` | String | ESBL producer | "YES" or empty |
| `Confirmed_CRE` | String | Carbapenem-resistant Enterobacterales | "YES" or empty |

### Data Relationships

#### Paired Observations
Each isolate has two rows:
- **Original:** Results from primary culture plate (with pretreatment effect)
- **Subculture:** Results from overnight subculture (reference standard)

The `Isolate_ID` links these paired observations.

#### Benchmark Group
An additional 150 isolates without pretreatment effect:
- `Plate` = "Positive"
- Used to establish baseline performance of identification methods

### Identification Method Cutoffs

#### Success Criteria (Passing Cutoff)
- **N_biochem:** `N_biochem_score >= 0.85`
- **P_biochem:** `P_biochem_score >= 0.85`
- **Vitek2:** `Vitek_score >= 85`
- **MALDI-TOF MS:** `MALDI_score >= 2.0` (standard), also evaluated at 1.9 and 1.8

#### Correct Identification
Among isolates passing cutoff, identification is correct if:
- Method ID matches `Final_ID` (reference identification)

### Example Data Rows

```csv
Isolate_ID,Plate,Specimen type,Final_ID,TSI_KA,TSI_AA,TSI_Gas,TSI_H2S,Lactose,Indole,LDM,LDC,Motility,Urease,Citrate,VP,Malonate,Mannitol,Profile_ID,N_biochem_ID,N_biochem_score,P_biochem_ID,P_biochem_score,Vitek_ID,Vitek_Bionumber,Vitek_score,MALDI_ID,MALDI_score,AST_performed,Predicted_AmpC,Predicted_ESBL,Confirmed_CRE
1,Original,Body fluid,Escherichia coli,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1.10E+12,Escherichia coli,0.999670006,Escherichia coli,0.999999821,Unidentified,B0635710554526611,0,Escherichia coli,1.93,YES,,,
1,Subculture,Body fluid,Escherichia coli,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1.10E+12,Escherichia coli,0.999670006,Escherichia coli,0.999999821,Escherichia coli,B0405610044006611,99,Escherichia coli,2.17,YES,,,
```

### Species Distribution

**Study Group (n=240):**
- Escherichia coli: 146 (60.8%)
- Klebsiella pneumoniae: 63 (26.3%)
- Proteus mirabilis: 11 (4.6%)
- Enterobacter cloacae complex: 8 (3.3%)
- Others: 12 (5.0%)

**Benchmark Group (n=150):**
- Similar distribution to study group (see Table S1 in manuscript)

## Cost-Effectiveness Analysis Data

### Methods_KI.csv

Contains performance metrics for each identification method:

| Column | Description |
|--------|-------------|
| `Method` | Method name (P_Biochem, Vitek, MALDI_2.0, MALDI_1.9, MALDI_1.8) |
| `Direct_cost` | Direct cost per test (USD) |
| `Success_rate` | Proportion passing cutoff |
| `Accuracy` | Proportion correct among successful |
| `Major_error_rate` | Rate of major errors |
| `Very_major_error_rate` | Rate of very major errors |

### RIsk_and_cost_KI.csv

Contains cost and QALY parameters for error consequences:

| Column | Description |
|--------|-------------|
| `Risk_type` | Error type (VM = very major, M = major) |
| `Description` | Outcome description |
| `Chance` | Probability of outcome |
| `QALY_loss` | Quality-adjusted life years lost |
| `Monetized_QALY_loss_lower` | Lower bound of monetized QALY (USD) |
| `Monetized_QALY_loss_upper` | Upper bound of monetized QALY (USD) |
| `Fixed_cost` | Direct treatment cost (USD) |

## Data Access

The actual dataset contains patient specimens and is not publicly available due to privacy regulations. Researchers interested in accessing the data should contact the corresponding author:

**Korakrit Imwattana**
Department of Microbiology
Faculty of Medicine Siriraj Hospital
Mahidol University, Bangkok, Thailand
Email: korakrit.imw@mahidol.ac.th

## Data Generation

Data was collected prospectively at Siriraj Hospital between May-November 2025 following the protocol described in the manuscript. All specimens were processed according to standard clinical microbiology protocols with institutional review board approval.
