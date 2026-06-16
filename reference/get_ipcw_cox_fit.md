# Fit an IPCW-weighted Cox proportional hazards model

Fits a weighted Cox model with a robust sandwich variance and returns
log-hazard ratio estimates together with both model-based and robust
standard errors, plus the exponentiated hazard ratio and 95% confidence
interval based on the robust standard error.

## Usage

``` r
get_ipcw_cox_fit(data, covariate = "x", weight = "wgt")
```

## Arguments

- data:

  A data frame in long (counting-process) format, as returned by
  [`get_ipcw_wgt()`](https://zabore.github.io/ipcw/reference/get_ipcw_wgt.md).
  Must contain columns `tstart`, `tstop`, `delta`, `id`, the covariate
  named by `covariate`, and the weight column named by `weight`.

- covariate:

  Character string. Name of the predictor covariate column. Default is
  `"x"`.

- weight:

  Character string. Name of the weight column. Default is `"wgt"`.

## Value

A data frame with one row per term containing:

- term:

  Covariate name.

- log_hr:

  Log hazard ratio estimate.

- log_hr_se:

  Model-based standard error of the log hazard ratio.

- log_hr_rob_se:

  Robust (sandwich) standard error of the log hazard ratio.

- hr:

  Hazard ratio estimate.

- hr_ci_low:

  Lower 95% confidence limit (based on robust SE).

- hr_ci_high:

  Upper 95% confidence limit (based on robust SE).

## Examples

``` r
data(single_example_ipcw_dat)
get_ipcw_cox_fit(single_example_ipcw_dat, weight = "wgt")
#> # A tibble: 1 × 7
#>   term  log_hr log_hr_se log_hr_rob_se    hr hr_ci_low hr_ci_high
#>   <chr>  <dbl>     <dbl>         <dbl> <dbl>     <dbl>      <dbl>
#> 1 x      -1.21     0.103         0.162 0.298     0.217      0.410
```
