# IPCW-weighted Fine-Gray sub-distribution hazard regression

Fits a Fine-Gray model on the Fine-Gray split dataset, weighting by the
inverse of the probability of remaining uncensored after event_2
(`p_notcens_after_death`), and uses a robust sandwich variance via
`cluster(id)`.

## Usage

``` r
fg_weighted(data_long_fg, extend = TRUE)
```

## Arguments

- data_long_fg:

  A data frame in Fine-Gray format with weights, as returned by
  [`add_fg_weights()`](https://zabore.github.io/ipcw/reference/add_fg_weights.md).

- extend:

  Logical. If `FALSE`, data are truncated at the minimum of the
  stratum-specific maximum follow-up times before fitting. Default is
  `TRUE`.

## Value

A matrix with one row per term and two columns: the log sub-distribution
hazard ratio and its robust standard error.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 200, censoring = "baseline")
dat_long    <- wide_to_long_CR(dat)
dat_long_fg <- fg_split(dat_long)
dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")
#> Warning: Loglik converged before variable  3 ; beta may be infinite. 
fg_weighted(dat_long_fg)
#>         coef robust se
#> z11 1.364480 0.5617613
#> z12 1.798696 0.5626535
#> z13 2.222293 0.5323173
```
