# Convert wide competing risks data to long (counting-process) format

Splits the wide dataset at every censoring time, creating one or more
rows per subject. The resulting long dataset is suitable for fitting Cox
models for the censoring distribution and for use with
`add_ipcw_weights()`.

## Usage

``` r
wide_to_long_CR(dat)
```

## Arguments

  - dat:
    
    A data frame with columns `t` (event/censoring time), `delta`
    (factor with levels `"censor"`, `"event_1"`, `"event_2"`), and `z1`
    (covariate). Typically the output of `sim_data_CR()`.

## Value

A data frame in long format with columns `id`, `z1`, `delta`, `censor`,
`event2_time`, `tstart`, and `tstop`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 100, censoring = "baseline")
dat_long <- wide_to_long_CR(dat)
head(dat_long)
```
