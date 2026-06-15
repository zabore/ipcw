# Estimate IPCW Kaplan-Meier survival probabilities by binary covariate

Fits a weighted Kaplan-Meier estimator stratified by the binary
covariate `x` and returns survival probabilities evaluated at a
pre-specified set of times.

## Usage

``` r
get_ipcw_km_prob_x(data, pre_times = seq(0, 50, 1))
```

## Arguments

  - data:
    
    A data frame in long (counting-process) format with an IPCW weight
    column `wgt`, as returned by `get_ipcw_wgt()`.

  - pre\_times:
    
    Numeric vector of times at which to evaluate survival probabilities.
    Defaults to `seq(0, 50, 1)`.

## Value

A tibble with columns `time`, `surv` (survival probability), and `x`
(stratum).

## Examples

``` r
data(single_example_ipcw_dat)
get_ipcw_km_prob_x(single_example_ipcw_dat, pre_times = seq(0, 500, 10))
```
