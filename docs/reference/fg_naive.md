# Naive Fine-Gray sub-distribution hazard regression

Fits a Fine-Gray model using the standard `survival::finegray()`
approach, without any IPCW adjustment. Serves as a comparison to
`fg_weighted()`.

## Usage

``` r
fg_naive(dat)
```

## Arguments

  - dat:
    
    A wide-format competing risks data frame with columns `t`, `delta`,
    and `z1`.

## Value

A matrix with one row per term and two columns: the log sub-distribution
hazard ratio and its standard error.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 200, censoring = "independent")
fg_naive(dat)
```
