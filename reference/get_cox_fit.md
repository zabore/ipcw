# Fit a standard (unweighted) Cox proportional hazards model

Convenience wrapper around
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) using
counting-process (`tstart`, `tstop`) time variables. Returns log-hazard
ratio estimates with standard errors and the exponentiated hazard ratio
with 95% confidence intervals.

## Usage

``` r
get_cox_fit(data, covariate = "x")
```

## Arguments

- data:

  A data frame in long (counting-process) format containing columns
  `tstart`, `tstop`, `delta`, and the covariate named by `covariate`.

- covariate:

  Character string. Name of the predictor covariate column. Default is
  `"x"`.

## Value

A data frame with one row per term containing:

- term:

  Covariate name.

- log_hr:

  Log hazard ratio estimate.

- log_hr_se:

  Standard error of the log hazard ratio.

- hr:

  Hazard ratio estimate.

- hr_ci_low:

  Lower 95% confidence limit.

- hr_ci_high:

  Upper 95% confidence limit.

## Examples

``` r
set.seed(20240429)
dat <- sim_data_se(n = 500)
dat_long <- get_ipcw_wgt(dat)
get_cox_fit(dat_long)
#> # A tibble: 1 × 6
#>   term  log_hr log_hr_se    hr hr_ci_low hr_ci_high
#>   <chr>  <dbl>     <dbl> <dbl>     <dbl>      <dbl>
#> 1 x      -1.85     0.147 0.157     0.117      0.209
```
