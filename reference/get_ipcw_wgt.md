# Compute IPCW weights for single-event survival data

Fits a Cox model for the censoring distribution and returns the original
dataset in counting-process (long) format with an unstabilized IPCW
weight column appended.

## Usage

``` r
get_ipcw_wgt(data, time_var = "t", event_var = "delta", cens_cov = "W2")
```

## Arguments

- data:

  A data frame containing the follow-up time, event indicator, and
  censoring covariate columns.

- time_var:

  Character string. Name of the observed follow-up time column. Default
  is `"t"`.

- event_var:

  Character string. Name of the event indicator column (1 = event, 0 =
  censored). Default is `"delta"`.

- cens_cov:

  Character string. Name of the covariate column used in the Cox model
  for the censoring distribution. Default is `"W2"`.

## Value

A data frame in long (counting-process) format with columns `tstart`,
`tstop`, `delta` (event indicator), `id`, `wgt` (the unstabilized IPCW
weight), and all original columns.

## Examples

``` r
set.seed(20240429)
dat <- sim_data_se(n = 500)
dat_long <- get_ipcw_wgt(dat)
head(dat_long)
#>          S x        W2 id    tstart     tstop delta censor      wgt
#> 1 426.7239 0 0.5326369  1 0.0000000 0.6883977     0      0 1.000000
#> 2 426.7239 0 0.5326369  1 0.6883977 1.2823328     0      0 1.000648
#> 3 426.7239 0 0.5326369  1 1.2823328 2.7738294     0      0 1.001298
#> 4 426.7239 0 0.5326369  1 2.7738294 5.1543230     0      0 1.001951
#> 5 426.7239 0 0.5326369  1 5.1543230 5.7907172     0      0 1.002608
#> 6 426.7239 0 0.5326369  1 5.7907172 6.1787837     0      0 1.003269
```
