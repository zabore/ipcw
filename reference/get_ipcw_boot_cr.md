# Bootstrap IPCW weighted data for competing risks survival analysis

Draws `B` bootstrap samples from the original wide-format competing
risks data, converts each sample to long (counting-process) format via
[`wide_to_long_cr()`](https://www.emilyzabor.com/ipcw/reference/wide_to_long_CR.md),
and appends IPCW weights via
[`add_ipcw_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights_cr.md).
Returns the results as a list of long-format data frames that can be
used directly with
[`cuminc_ipcw_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_ipcw_cr.md),
or further processed with
[`fg_split_cr()`](https://www.emilyzabor.com/ipcw/reference/fg_split_cr.md)
and
[`add_fg_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_fg_weights_cr.md)
for a Fine-Gray bootstrap.

## Usage

``` r
get_ipcw_boot_cr(
  data,
  B = 500,
  time_var = "t",
  event_var = "delta",
  covariate = "z1",
  cens_level = "censor",
  event2_level = "event_2",
  strat = "no",
  seed = NULL
)
```

## Arguments

- data:

  A wide-format competing risks data frame containing the columns
  specified by `time_var`, `event_var`, and `covariate`.

- B:

  Integer. Number of bootstrap samples. Default is `500`.

- time_var:

  Character string. Name of the event/censoring time column. Default is
  `"t"`.

- event_var:

  Character string. Name of the event indicator column. Default is
  `"delta"`.

- covariate:

  Character string. Name of the covariate column. Default is `"z1"`.

- cens_level:

  Character string. Factor level in `event_var` representing censoring.
  Default is `"censor"`.

- event2_level:

  Character string. Factor level in `event_var` representing the
  competing event (event type 2). Default is `"event_2"`.

- strat:

  Character. Passed to
  [`add_ipcw_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights_cr.md).
  `"no"` (default) fits a single Cox model for the censoring
  distribution; `"yes"` estimates the censoring distribution
  non-parametrically within each stratum of `covariate`.

- seed:

  Optional integer seed for reproducibility. Default is `NULL`.

## Value

A list of `B` data frames, each in long (counting-process) format with
an IPCW weight column `p_notcens`, as returned by
[`add_ipcw_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights_cr.md).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
dat <- sim_data_cr(n = 200, censoring = "baseline")
boot_list <- get_ipcw_boot_cr(dat, B = 50)
} # }
```
