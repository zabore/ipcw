# Package index

## Single-event IPCW

Functions for IPCW analysis of single-event survival data, including
data simulation, weight estimation, weighted Cox regression, and
weighted Kaplan-Meier curves.

<!-- end list -->

  - `sim_data_SE()` : Simulate single-event survival data with
    informative censoring
  - `get_ipcw_wgt()` : Compute IPCW weights for single-event survival
    data
  - `get_ipcw_cox_fit()` : Fit an IPCW-weighted Cox proportional hazards
    model
  - `get_ipcw_km_prob_x()` : Estimate IPCW Kaplan-Meier survival
    probabilities by binary covariate
  - `get_cox_fit()` : Fit a standard (unweighted) Cox proportional
    hazards model
  - `get_boot_var()` : Compute bootstrap variance of the log hazard
    ratio
  - `get_boot_pci()` : Compute bootstrap percentile confidence interval
    for the log hazard ratio

## Competing-risks IPCW

Functions for IPCW analysis of competing risks data, including weight
estimation, cumulative incidence estimation, and Fine-Gray regression.

<!-- end list -->

  - `sim_data_CR()` : Simulate competing risks survival data
  - `wide_to_long_CR()` : Convert wide competing risks data to long
    (counting-process) format
  - `add_ipcw_weights()` : Add IPCW weights to competing risks
    long-format data
  - `cuminc_naive()` : Naive (unweighted) cumulative incidence estimate
  - `cuminc_waverage()` : Weighted-average (non-parametric) IPCW
    cumulative incidence estimate
  - `cuminc_ipcw()` : Cox model IPCW cumulative incidence estimate
  - `fg_split()` : Prepare long-format data for Fine-Gray weighted
    regression
  - `add_fg_weights()` : Add Fine-Gray IPCW weights to Fine-Gray split
    data
  - `fg_naive()` : Naive Fine-Gray sub-distribution hazard regression
  - `fg_weighted()` : IPCW-weighted Fine-Gray sub-distribution hazard
    regression

## Data

Example datasets included with the package.

<!-- end list -->

  - `single_example_dat` : Single-event example dataset
  - `single_example_ipcw_dat` : Single-event example dataset in IPCW
    long format
