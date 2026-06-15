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
```
