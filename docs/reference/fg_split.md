# Prepare long-format data for Fine-Gray weighted regression

After a type-2 event (event\_2), subjects are artificially re-entered
into the risk set (as in the Fine-Gray sub-distribution hazard model).
This function appends those additional rows, split at every censoring
time.

## Usage

``` r
fg_split(data_long)
```

## Arguments

  - data\_long:
    
    A long-format data frame, as returned by `wide_to_long_CR()`.

## Value

A data frame in Fine-Gray format, with additional rows for subjects who
experienced event\_2, sorted by `id` and `tstart`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 100, censoring = "baseline")
dat_long <- wide_to_long_CR(dat)
dat_long_fg <- fg_split(dat_long)
nrow(dat_long_fg) > nrow(dat_long)
```
