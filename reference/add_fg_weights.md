# Add Fine-Gray IPCW weights to Fine-Gray split data

Computes the probability of remaining uncensored after a type-2 event
and appends it as column `p_notcens_after_death`.

## Usage

``` r
add_fg_weights(data_long_fg, covariate = "z1", strat = "no")
```

## Arguments

- data_long_fg:

  A data frame in Fine-Gray format, as returned by
  [`fg_split()`](https://www.emilyzabor.com/ipcw/reference/fg_split.md).

- covariate:

  Character string. Name of the covariate column. Default is `"z1"`.

- strat:

  Character. Passed to
  [`add_ipcw_weights()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights.md).
  `"no"` (default) uses a Cox model; `"yes"` uses stratum-specific KM
  estimates.

## Value

`data_long_fg` with an additional column `p_notcens_after_death`.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 100, censoring = "baseline")
dat_long    <- wide_to_long_cr(dat)
dat_long_fg <- fg_split(dat_long)
dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")
summary(dat_long_fg$p_notcens_after_death)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0002632 1.0000000 1.0000000 0.9525524 1.0000000 1.0000000 
```
