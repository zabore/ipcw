# Naive (unweighted) cumulative incidence estimate

Estimates the cumulative incidence of event type 1 using the standard
Aalen-Johansen estimator without any IPCW adjustment. Serves as a
comparison method when censoring is independent of covariates.

## Usage

``` r
cuminc_naive(dat, esttimes)
```

## Arguments

  - dat:
    
    A wide-format competing risks data frame with columns `t` and
    `delta` (factor with levels `"censor"`, `"event_1"`, `"event_2"`).

  - esttimes:
    
    Numeric vector of times at which to return estimates.

## Value

A numeric vector of cumulative incidence estimates at `esttimes`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 200, censoring = "independent")
cuminc_naive(dat, esttimes = seq(0, 5, 0.5))
```
