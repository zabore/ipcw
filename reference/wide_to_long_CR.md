# Convert wide competing risks data to long (counting-process) format

Splits the wide dataset at every censoring time, creating one or more
rows per subject. The resulting long dataset is suitable for fitting Cox
models for the censoring distribution and for use with
[`add_ipcw_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights_cr.md).

## Usage

``` r
wide_to_long_cr(
  dat,
  time_var = "t",
  event_var = "delta",
  covariate = "z1",
  cens_level = "censor",
  event2_level = "event_2"
)
```

## Arguments

- dat:

  A wide-format competing risks data frame. Must contain the columns
  specified by `time_var`, `event_var`, and `covariate`.

- time_var:

  Character string. Name of the event/censoring time column. Default is
  `"t"`.

- event_var:

  Character string. Name of the event indicator column (factor with
  levels for censoring and the two event types). Default is `"delta"`.

- covariate:

  Character string. Name of the covariate column. Default is `"z1"`.

- cens_level:

  Character string. Factor level in `event_var` representing censoring.
  Default is `"censor"`.

- event2_level:

  Character string. Factor level in `event_var` representing the
  competing event (event type 2). Default is `"event_2"`.

## Value

A data frame in long format with columns `id`, `delta`, `censor`,
`event2_time`, `tstart`, `tstop`, and the covariate column.

## Examples

``` r
set.seed(42)
dat <- sim_data_cr(n = 100, censoring = "baseline")
dat_long <- wide_to_long_cr(dat)
head(dat_long)
#>   z1 id censor event2_time       tstart        tstop  delta
#> 1  0  1      0          NA 0.0000000000 0.0008870787 censor
#> 2  0  1      0          NA 0.0008870787 0.0236331037 censor
#> 3  0  1      0          NA 0.0236331037 0.0295500713 censor
#> 4  0  1      0          NA 0.0295500713 0.0724101320 censor
#> 5  0  1      0          NA 0.0724101320 0.0812450008 censor
#> 6  0  1      0          NA 0.0812450008 0.1824410379 censor
```
