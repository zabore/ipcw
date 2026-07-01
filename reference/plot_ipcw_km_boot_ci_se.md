# Plot IPCW Kaplan-Meier curves with bootstrap percentile confidence intervals

Combines bootstrap survival curves from
[`get_ipcw_boot_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot_se.md)
with point estimates from the original IPCW-weighted data to produce a
Kaplan-Meier plot with 95% percentile confidence intervals.

## Usage

``` r
plot_ipcw_km_boot_ci_se(
  boot_data,
  orig_data,
  pre_times = seq(0, 50, 1),
  covariate = "x",
  weight_var = "wgt",
  event_var = "delta"
)
```

## Arguments

- boot_data:

  A list of long-format data frames as returned by
  [`get_ipcw_boot_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot_se.md).

- orig_data:

  A single long-format data frame for the original (non-bootstrapped)
  dataset, as returned by
  [`get_ipcw_wgt_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_wgt_se.md).
  Used for the point estimates.

- pre_times:

  Numeric vector of times at which to evaluate survival probabilities.
  Default is `seq(0, 50, 1)`.

- covariate:

  Character string. Name of the stratification covariate column. Default
  is `"x"`.

- weight_var:

  Character string. Name of the IPCW weight column. Default is `"wgt"`.

- event_var:

  Character string. Name of the event indicator column in the
  long-format data. Default is `"delta"`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. Add layers or themes to customise.

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run due to computational time. 
# Use a smaller n and/or smaller B for testing.
set.seed(1)
dat      <- sim_data_se(n = 500)
dat_long <- get_ipcw_wgt_se(dat)
boot_list <- get_ipcw_boot_se(dat, B = 500)
plot_ipcw_km_boot_ci_se(boot_list, dat_long, pre_times = seq(0, 500, 10))
} # }
```
