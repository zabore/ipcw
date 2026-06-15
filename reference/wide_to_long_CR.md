# Convert wide competing risks data to long (counting-process) format

Splits the wide dataset at every censoring time, creating one or more
rows per subject. The resulting long dataset is suitable for fitting Cox
models for the censoring distribution and for use with
[`add_ipcw_weights()`](https://zabore.github.io/ipcw/reference/add_ipcw_weights.md).

## Usage

``` r
wide_to_long_CR(dat)
```

## Arguments

- dat:

  A data frame with columns `t` (event/censoring time), `delta` (factor
  with levels `"censor"`, `"event_1"`, `"event_2"`), and `z1`
  (covariate). Typically the output of
  [`sim_data_CR()`](https://zabore.github.io/ipcw/reference/sim_data_CR.md).

## Value

A data frame in long format with columns `id`, `z1`, `delta`, `censor`,
`event2_time`, `tstart`, and `tstop`.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 100, censoring = "baseline")
dat_long <- wide_to_long_CR(dat)
head(dat_long)
#>   z1 id censor event2_time       tstart        tstop  delta
#> 1  0  1      0          NA 0.0000000000 0.0008870787 censor
#> 2  0  1      0          NA 0.0008870787 0.0236331037 censor
#> 3  0  1      0          NA 0.0236331037 0.0295500713 censor
#> 4  0  1      0          NA 0.0295500713 0.0724101320 censor
#> 5  0  1      0          NA 0.0724101320 0.0812450008 censor
#> 6  0  1      0          NA 0.0812450008 0.1824410379 censor
```
