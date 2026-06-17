# Prepare long-format data for Fine-Gray weighted regression

After a type-2 event, subjects are artificially re-entered into the risk
set (as in the Fine-Gray sub-distribution hazard model). This function
appends those additional rows, split at every censoring time.

## Usage

``` r
fg_split(data_long, covariate = "z1", event2_level = "event_2")
```

## Arguments

- data_long:

  A long-format data frame, as returned by
  [`wide_to_long_cr()`](https://www.emilyzabor.com/ipcw/reference/wide_to_long_CR.md).

- covariate:

  Character string. Name of the covariate column. Default is `"z1"`.

- event2_level:

  Character string. Factor level in the `delta` column representing the
  competing event. Default is `"event_2"`.

## Value

A data frame in Fine-Gray format, with additional rows for subjects who
experienced the competing event, sorted by `id` and `tstart`.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 100, censoring = "baseline")
dat_long <- wide_to_long_cr(dat)
dat_long_fg <- fg_split(dat_long)
nrow(dat_long_fg) > nrow(dat_long)
#> [1] TRUE
```
