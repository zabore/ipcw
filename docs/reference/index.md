# Package index

## Single-event IPCW

Functions for IPCW analysis of single-event survival data, including
data simulation, weight estimation, weighted Cox regression, and
weighted Kaplan-Meier curves.

- [`sim_data_se()`](https://www.emilyzabor.com/ipcw/reference/sim_data_se.md)
  : Simulate single-event survival data with informative censoring
- [`get_ipcw_wgt_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_wgt_se.md)
  : Compute IPCW weights for single-event survival data
- [`get_ipcw_cox_fit_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_cox_fit_se.md)
  : Fit an IPCW-weighted Cox proportional hazards model
- [`get_ipcw_km_prob_x_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_km_prob_x_se.md)
  : Estimate IPCW Kaplan-Meier survival probabilities by covariate
- [`get_cox_fit_se()`](https://www.emilyzabor.com/ipcw/reference/get_cox_fit_se.md)
  : Fit a standard (unweighted) Cox proportional hazards model
- [`get_boot_var_se()`](https://www.emilyzabor.com/ipcw/reference/get_boot_var_se.md)
  : Compute bootstrap variance of the log hazard ratio
- [`get_boot_pci_se()`](https://www.emilyzabor.com/ipcw/reference/get_boot_pci_se.md)
  : Compute bootstrap percentile confidence interval for the log hazard
  ratio
- [`get_ipcw_boot_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot_se.md)
  : Bootstrap IPCW weighted data for single-event survival analysis
- [`plot_ipcw_km_boot_ci_se()`](https://www.emilyzabor.com/ipcw/reference/plot_ipcw_km_boot_ci_se.md)
  : Plot IPCW Kaplan-Meier curves with bootstrap percentile confidence
  intervals
- [`cox_boot_lhr`](https://www.emilyzabor.com/ipcw/reference/cox_boot_lhr.md)
  : Bootstrap log hazard ratios from IPCW Cox models

## Competing-risks IPCW

Functions for IPCW analysis of competing risks data, including weight
estimation, cumulative incidence estimation, and Fine-Gray regression.

- [`sim_data_cr()`](https://www.emilyzabor.com/ipcw/reference/sim_data_cr.md)
  : Simulate competing risks survival data
- [`wide_to_long_cr()`](https://www.emilyzabor.com/ipcw/reference/wide_to_long_cr.md)
  : Convert wide competing risks data to long (counting-process) format
- [`add_ipcw_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights_cr.md)
  : Add IPCW weights to competing risks long-format data
- [`cuminc_naive_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_naive_cr.md)
  : Naive (unweighted) cumulative incidence estimate
- [`cuminc_waverage_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_waverage_cr.md)
  : Weighted-average (non-parametric) IPCW cumulative incidence estimate
- [`cuminc_ipcw_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_ipcw_cr.md)
  : Cox model IPCW cumulative incidence estimate
- [`fg_split_cr()`](https://www.emilyzabor.com/ipcw/reference/fg_split_cr.md)
  : Prepare long-format data for Fine-Gray regression
- [`add_fg_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_fg_weights_cr.md)
  : Add Fine-Gray weights to Fine-Gray split data
- [`fg_naive_cr()`](https://www.emilyzabor.com/ipcw/reference/fg_naive_cr.md)
  : Naive Fine-Gray sub-distribution hazard regression
- [`fg_weighted_cr()`](https://www.emilyzabor.com/ipcw/reference/fg_weighted_cr.md)
  : FG-weighted Fine-Gray sub-distribution hazard regression
- [`get_ipcw_boot_cr()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot_cr.md)
  : Bootstrap IPCW weighted data for competing risks survival analysis
- [`get_boot_pci_cr()`](https://www.emilyzabor.com/ipcw/reference/get_boot_pci_cr.md)
  : Compute bootstrap percentile confidence intervals for competing
  risks estimates
- [`plot_ipcw_cuminc_boot_ci_cr()`](https://www.emilyzabor.com/ipcw/reference/plot_ipcw_cuminc_boot_ci_cr.md)
  : Plot IPCW cumulative incidence with bootstrap percentile confidence
  intervals
- [`fg_strat_boot`](https://www.emilyzabor.com/ipcw/reference/fg_strat_boot.md)
  : Bootstrap Fine-Gray coefficients from stratified FG weights
