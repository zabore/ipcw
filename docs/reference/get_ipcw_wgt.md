# Compute IPCW weights for single-event survival data

Fits a Cox model for the censoring distribution and returns the original
dataset in counting-process (long) format with an unstabilized IPCW
weight column appended. The dataset must contain columns named `t`
(event/censoring time), `delta` (event indicator, 1 = event), and `W2`
(the covariate used to model the censoring distribution).

## Usage

``` r
get_ipcw_wgt(data)
```

## Arguments

- data:

  A data frame with columns `t`, `delta`, and `W2`, as produced by the
  data-generation code in the package vignette.

## Value

A data frame in long (counting-process) format with columns `tstart`,
`tstop`, `delta`, `id`, `wgt` (the unstabilized IPCW weight), and all
original columns.

## Examples

``` r
data(single_example_dat)
dat_long <- get_ipcw_wgt(single_example_dat)
head(dat_long)
#>          S x        W2 id    tstart     tstop delta censor      wgt
#> 1 426.7239 0 0.5326369  1 0.0000000 0.6883977     0      0 1.000000
#> 2 426.7239 0 0.5326369  1 0.6883977 1.2823328     0      0 1.000648
#> 3 426.7239 0 0.5326369  1 1.2823328 2.7738294     0      0 1.001298
#> 4 426.7239 0 0.5326369  1 2.7738294 5.1543230     0      0 1.001951
#> 5 426.7239 0 0.5326369  1 5.1543230 5.7907172     0      0 1.002608
#> 6 426.7239 0 0.5326369  1 5.7907172 6.1787837     0      0 1.003269
```
