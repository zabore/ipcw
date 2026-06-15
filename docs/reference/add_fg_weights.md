# Add Fine-Gray IPCW weights to Fine-Gray split data

Computes the probability of remaining uncensored after a type-2 event
and appends it as column `p_notcens_after_death`.

## Usage

``` r
add_fg_weights(data_long_fg, strat = "no")
```

## Arguments

  - data\_long\_fg:
    
    A data frame in Fine-Gray format, as returned by `fg_split()`.

  - strat:
    
    Character. Passed to `add_ipcw_weights()`. `"no"` (default) uses a
    Cox model; `"yes"` uses stratum-specific KM estimates.

## Value

`data_long_fg` with an additional column `p_notcens_after_death`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 100, censoring = "baseline")
dat_long    <- wide_to_long_CR(dat)
dat_long_fg <- fg_split(dat_long)
dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")
summary(dat_long_fg$p_notcens_after_death)
```
