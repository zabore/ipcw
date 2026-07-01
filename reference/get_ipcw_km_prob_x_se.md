# Estimate IPCW Kaplan-Meier survival probabilities by covariate

Fits a weighted Kaplan-Meier estimator stratified by the specified
covariate and returns survival probabilities evaluated at a
pre-specified set of times.

## Usage

``` r
get_ipcw_km_prob_x_se(
  data,
  covariate = "x",
  weight_var = "wgt",
  event_var = "delta",
  pre_times = seq(0, 50, 1)
)
```

## Arguments

- data:

  A data frame in long (counting-process) format with an IPCW weight
  column, as returned by
  [`get_ipcw_wgt_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_wgt_se.md).
  Must contain time columns named `tstart` and `tstop` for the counting
  process time intervals.

- covariate:

  Character string. Name of the stratification covariate column. Default
  is `"x"`.

- weight_var:

  Character string. Name of the IPCW weight column. Default is `"wgt"`.

- event_var:

  Character string. Name of the event indicator column (1 = event, 0 =
  censored). Default is `"delta"`.

- pre_times:

  Numeric vector of times at which to evaluate survival probabilities.
  Defaults to `seq(0, 50, 1)`. Choose the second number to cover the
  entire range of times observed in `data`

## Value

A tibble with columns `time`, `surv` (survival probability), and a
column named after `covariate` containing the stratum labels.

## Examples

``` r
set.seed(20240429)
dat <- sim_data_se(n = 500)
dat_long <- get_ipcw_wgt_se(dat)
get_ipcw_km_prob_x_se(dat_long, pre_times = seq(0, 2429, 1))
#> # A tibble: 3,724 × 3
#>     time  surv x    
#>    <dbl> <dbl> <fct>
#>  1     0     1 x=0  
#>  2     1     1 x=0  
#>  3     2     1 x=0  
#>  4     3     1 x=0  
#>  5     4     1 x=0  
#>  6     5     1 x=0  
#>  7     6     1 x=0  
#>  8     7     1 x=0  
#>  9     8     1 x=0  
#> 10     9     1 x=0  
#> # ℹ 3,714 more rows
```
