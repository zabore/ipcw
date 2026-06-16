# Estimate IPCW Kaplan-Meier survival probabilities by covariate

Fits a weighted Kaplan-Meier estimator stratified by the specified
covariate and returns survival probabilities evaluated at a
pre-specified set of times.

## Usage

``` r
get_ipcw_km_prob_x(
  data,
  covariate = "x",
  weight_var = "wgt",
  pre_times = seq(0, 50, 1)
)
```

## Arguments

- data:

  A data frame in long (counting-process) format with an IPCW weight
  column, as returned by
  [`get_ipcw_wgt()`](https://zabore.github.io/ipcw/reference/get_ipcw_wgt.md).

- covariate:

  Character string. Name of the stratification covariate column. Default
  is `"x"`.

- weight_var:

  Character string. Name of the IPCW weight column. Default is `"wgt"`.

- pre_times:

  Numeric vector of times at which to evaluate survival probabilities.
  Defaults to `seq(0, 50, 1)`.

## Value

A tibble with columns `time`, `surv` (survival probability), and a
column named after `covariate` containing the stratum labels.

## Examples

``` r
data(single_example_ipcw_dat)
get_ipcw_km_prob_x(single_example_ipcw_dat, pre_times = seq(0, 500, 10))
#> # A tibble: 102 × 3
#>     time  surv x    
#>    <dbl> <dbl> <fct>
#>  1     0 1     x=0  
#>  2    10 1     x=0  
#>  3    20 0.995 x=0  
#>  4    30 0.995 x=0  
#>  5    40 0.989 x=0  
#>  6    50 0.989 x=0  
#>  7    60 0.973 x=0  
#>  8    70 0.965 x=0  
#>  9    80 0.956 x=0  
#> 10    90 0.947 x=0  
#> # ℹ 92 more rows
```
