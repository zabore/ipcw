# Naive (unweighted) cumulative incidence estimate

Estimates the cumulative incidence of event type 1 using the standard
Aalen-Johansen estimator without any IPCW adjustment. Serves as a
comparison method when censoring is independent of covariates.

## Usage

``` r
cuminc_naive(dat, esttimes, time_var = "t", event_var = "delta")
```

## Arguments

- dat:

  A wide-format competing risks data frame containing the columns
  specified by `time_var` and `event_var`.

- esttimes:

  Numeric vector of times at which to return estimates.

- time_var:

  Character string. Name of the event/censoring time column. Default is
  `"t"`.

- event_var:

  Character string. Name of the event indicator column. Default is
  `"delta"`.

## Value

A numeric vector of cumulative incidence estimates at `esttimes`.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 200, censoring = "independent")
cuminc_naive(dat, esttimes = seq(0, 5, 0.5))
#>  [1] 0.0000000 0.2705481 0.3935556 0.4512375 0.4692767 0.4945862 0.5009775
#>  [8] 0.5079067 0.5155229 0.5155229 0.5155229
```
