# Naive Fine-Gray sub-distribution hazard regression

Fits a Fine-Gray model using the standard
[`survival::finegray()`](https://rdrr.io/pkg/survival/man/finegray.html)
approach, without any IPCW adjustment. Serves as a comparison to
[`fg_weighted_cr()`](https://www.emilyzabor.com/ipcw/reference/fg_weighted_cr.md).

## Usage

``` r
fg_naive_cr(dat, time_var = "t", event_var = "delta", covariate = "z1")
```

## Arguments

- dat:

  A wide-format competing risks data frame containing the columns
  specified by `time_var`, `event_var`, and `covariate`.

- time_var:

  Character string. Name of the event/censoring time column. Default is
  `"t"`.

- event_var:

  Character string. Name of the event indicator column. Default is
  `"delta"`.

- covariate:

  Character string. Name of the covariate column. Default is `"z1"`.

## Value

A matrix with one row per term and two columns: the log sub-distribution
hazard ratio and its standard error.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 200, censoring = "independent")
fg_naive_cr(dat)
#>          coef  se(coef)
#> z11 0.9452312 0.3624174
#> z12 1.3652479 0.3799441
#> z13 1.7795374 0.3446775
```
