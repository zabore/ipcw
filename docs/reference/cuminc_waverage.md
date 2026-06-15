# Weighted-average (non-parametric) IPCW cumulative incidence estimate

Estimates the marginal cumulative incidence of event type 1 by computing
stratum-specific Aalen-Johansen estimates (stratified by `z1`) and then
combining them with sample-proportion weights. This is the
non-parametric IPCW approach.

## Usage

``` r
cuminc_waverage(dat, esttimes = seq(from = 0, to = 10, length.out = 100))
```

## Arguments

  - dat:
    
    A wide-format competing risks data frame with columns `t`, `delta`,
    and `z1`.

  - esttimes:
    
    Numeric vector of times at which to return estimates. Defaults to
    100 equally spaced points from 0 to 10.

## Value

A numeric vector of weighted-average cumulative incidence estimates at
`esttimes`. Values beyond the minimum of the stratum-specific maximum
follow-up times are set to `NA`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 200, censoring = "baseline")
cuminc_waverage(dat, esttimes = seq(0, 5, 0.5))
```
