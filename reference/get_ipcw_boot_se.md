# Bootstrap IPCW weighted data for single-event survival analysis

Draws `B` bootstrap samples from the original wide-format data, fits
IPCW weights to each sample via
[`get_ipcw_wgt_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_wgt_se.md),
and returns the results as a list of long-format data frames. The list
can be passed to
[`plot_ipcw_km_boot_ci_se()`](https://www.emilyzabor.com/ipcw/reference/plot_ipcw_km_boot_ci_se.md)
for plotting or used directly to bootstrap other quantities (e.g. Cox
hazard ratios via
[`get_ipcw_cox_fit_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_cox_fit_se.md)).

## Usage

``` r
get_ipcw_boot_se(
  data,
  B = 500,
  time_var = "t",
  event_var = "delta",
  cens_cov = "W2",
  seed = NULL
)
```

## Arguments

- data:

  A wide-format data frame containing the columns specified by
  `time_var`, `event_var`, and `cens_cov`.

- B:

  Integer. Number of bootstrap samples. Default is `500`.

- time_var:

  Character string. Name of the observed follow-up time column. Default
  is `"t"`.

- event_var:

  Character string. Name of the event indicator column (1 = event, 0 =
  censored). Default is `"delta"`.

- cens_cov:

  Character string. Name of the covariate column used to model the
  censoring distribution. Default is `"W2"`.

- seed:

  Optional integer seed for reproducibility. Default is `NULL`.

## Value

A list of `B` data frames, each in long (counting-process) format as
returned by
[`get_ipcw_wgt_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_wgt_se.md).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
dat <- sim_data_se(n = 200)
boot_list <- get_ipcw_boot_se(dat, B = 50)
} # }
```
