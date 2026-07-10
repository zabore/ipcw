# FG-weighted Fine-Gray sub-distribution hazard regression

Fits a Fine-Gray model on the Fine-Gray split dataset, weighting by the
probability of remaining uncensored after the competing event
(`p_notcens_after_death`), and uses a robust sandwich variance via
`cluster(id)`.

## Usage

``` r
fg_weighted_cr(
  data_long_fg,
  covariate = "z1",
  extend = TRUE,
  event1_level = "event_1"
)
```

## Arguments

- data_long_fg:

  A data frame in Fine-Gray format with weights, as returned by
  [`add_fg_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_fg_weights_cr.md).

- covariate:

  Character string. Name of the covariate column. Default is `"z1"`.

- extend:

  Logical. If `FALSE`, data are truncated at the minimum of the
  stratum-specific maximum follow-up times before fitting. Default is
  `TRUE`.

- event1_level:

  Character string. Factor level in the `delta` column representing the
  primary event. Default is `"event_1"`.

## Value

A matrix with one row per term and two columns: the log sub-distribution
hazard ratio and its robust standard error.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 300, censoring = "baseline")
dat_long    <- wide_to_long_cr(dat)
dat_long_fg <- fg_split_cr(dat_long)
dat_long_fg <- add_fg_weights_cr(dat_long_fg, strat = "no")
fg_weighted_cr(dat_long_fg)
#>          coef robust se
#> z11 0.9135669 0.3975405
#> z12 1.3004739 0.3955882
#> z13 1.4468114 0.3928682
```
