# Add IPCW weights to competing risks long-format data

Estimates the probability of remaining uncensored and appends it as
column `p_notcens` to the dataset. Supports two estimation strategies: a
Cox proportional hazards model (`strat = "no"`) or non-parametric (KM)
estimates within each level of the covariate (`strat = "yes"`).

## Usage

``` r
add_ipcw_weights_cr(
  data_long,
  covariate = "z1",
  strat = "no",
  new_data = NULL,
  by.start = TRUE
)
```

## Arguments

- data_long:

  A data frame in long format, as returned by
  [`wide_to_long_cr()`](https://www.emilyzabor.com/ipcw/reference/wide_to_long_CR.md).
  Must contain columns `tstart`, `tstop`, `censor`, and the covariate
  named by `covariate`.

- covariate:

  Character string. Name of the covariate column used to model the
  censoring distribution. Default is `"z1"`.

- strat:

  Character. `"no"` (default) fits a single Cox model for the censoring
  distribution using `covariate`. `"yes"` estimates the censoring
  distribution non-parametrically within each stratum of `covariate`.

- new_data:

  Optional data frame to which weights are applied. If `NULL` (default),
  weights are computed for `data_long` itself.

- by.start:

  Logical. If `TRUE` (default), the weight is P(not censored by
  `tstart`). If `FALSE`, the weight is P(not censored during the
  interval from `tstart` to `tstop`).

## Value

`new_data` (or `data_long` if `new_data` is `NULL`) with an additional
column `p_notcens`.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 100, censoring = "baseline")
dat_long <- wide_to_long_cr(dat)
dat_long <- add_ipcw_weights_cr(dat_long, strat = "no")
summary(dat_long$p_notcens)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0547  0.7641  0.9126  0.8413  0.9757  1.0000 
```
