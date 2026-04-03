# EOCancerRiskAging
Code repository for article "Generational Shifts in Early-Onset Cancer Risk: The Role of Systemic and Organ-Specific Aging"

# Biological age and early-onset cancer: publication-ready code layout

This folder reorganizes the analysis into four publication-facing files:

- `01_data_cohort_cleaning.Rmd`
  - raw data import
  - variable harmonization
  - cohort inclusion/exclusion
  - biological age derivation
  - creation of manuscript-aligned analysis-ready datasets
- `02_main_analysis.Rmd`
  - birth cohort trends
  - PhenoAge models
  - KDM models
  - metabolomic aging models
  - organ-specific aging models
- `03_sensitivity_analysis.Rmd`
  - follow-up restriction
  - EO definition `<50`
  - additional adjustment models
  - repeated-measure stability
- `R/functions.R`
  - reusable helper functions
  - generic enough to adapt to external datasets

## Design principles used

1. **One-way data flow**: raw/cohort cleaning first, all downstream analysis reads frozen analysis-ready files.
2. **Manuscript alignment**: exclusions and analysis families mirror the current manuscript.
3. **Simple adaptation**: core clock and survival helpers accept user-supplied column names.
4. **Minimal duplication**: repetitive Cox/FDR/table logic moved into `R/functions.R`.

## Suggested repository structure

```text
project/
  data_raw/                 # not shared
  data_derived/             # analysis-ready files created in 01
  results/
    tables/
    figures/
    model_objects/
  R/
    functions.R
  01_data_cohort_cleaning.Rmd
  02_main_analysis.Rmd
  03_sensitivity_analysis.Rmd
  README.md
```

## Important notes

- The scripts are intentionally written as **publication templates** and **clean analysis scaffolds**.
- You will still need to replace local path objects and, in a few places, exact variable names if they differ from the latest working dataset.
- The cohort definitions and analysis families were organized to match the current manuscript provided in this conversation.
