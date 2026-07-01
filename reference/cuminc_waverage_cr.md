# Weighted-average (non-parametric) IPCW cumulative incidence estimate

Estimates the marginal cumulative incidence of event type 1 by computing
stratum-specific Aalen-Johansen estimates and then combining them with
sample-proportion weights. This is the non-parametric IPCW approach.

## Usage

``` r
cuminc_waverage_cr(
  dat,
  esttimes = seq(from = 0, to = 10, length.out = 100),
  time_var = "t",
  event_var = "delta",
  covariate = "z1"
)
```

## Arguments

- dat:

  A wide-format competing risks data frame containing the columns
  specified by `time_var`, `event_var`, and `covariate`.

- esttimes:

  Numeric vector of times at which to return estimates. Defaults to 100
  equally spaced points from 0 to 10.

- time_var:

  Character string. Name of the event/censoring time column. Default is
  `"t"`.

- event_var:

  Character string. Name of the event indicator column. Default is
  `"delta"`.

- covariate:

  Character string. Name of the stratification covariate column. Default
  is `"z1"`.

## Value

A numeric vector of weighted-average cumulative incidence estimates at
`esttimes`. Values beyond the minimum of the stratum-specific maximum
follow-up times are set to `NA`.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 200, censoring = "baseline")
cuminc_waverage_cr(dat, esttimes = seq(0, 5, 0.5))
#>  [1] 0.0000000 0.2549215 0.3706694 0.4349190 0.4462671 0.4689078 0.4689078
#>  [8] 0.4689078        NA        NA        NA
```
