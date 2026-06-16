# Cox model IPCW cumulative incidence estimate

Estimates the marginal cumulative incidence of event type 1 using a
weighted Aalen-Johansen estimator, where the weights are the IPCW
weights stored in the `p_notcens` column (as added by
[`add_ipcw_weights()`](https://zabore.github.io/ipcw/reference/add_ipcw_weights.md)).

## Usage

``` r
cuminc_ipcw(
  data_long,
  esttimes = seq(from = 0, to = 10, length.out = 100),
  extend = TRUE,
  covariate = "z1"
)
```

## Arguments

- data_long:

  A long-format data frame with IPCW weights, as returned by
  [`add_ipcw_weights()`](https://zabore.github.io/ipcw/reference/add_ipcw_weights.md).
  Must contain columns `tstart`, `tstop`, `delta`, `id`, and
  `p_notcens`.

- esttimes:

  Numeric vector of times at which to return estimates. Defaults to 100
  equally spaced points from 0 to 10.

- extend:

  Logical. If `FALSE`, estimates beyond the minimum of the
  stratum-specific maximum follow-up times are set to `NA`. Default is
  `TRUE`.

- covariate:

  Character string. Name of the covariate column, used only when
  `extend = FALSE` to compute the truncation time. Default is `"z1"`.

## Value

A numeric vector of IPCW cumulative incidence estimates at `esttimes`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 200, censoring = "baseline")
dat_long <- wide_to_long_CR(dat)
dat_long <- add_ipcw_weights(dat_long, strat = "no")
#> Warning: Loglik converged before variable  3 ; beta may be infinite. 
cuminc_ipcw(dat_long, esttimes = seq(0, 5, 0.5))
#>  [1] 0.0000000 0.2525084 0.3656652 0.4349979 0.4489364 0.4768644 0.4768644
#>  [8] 0.4768644 0.4847719 0.4847719 0.4847719
```
