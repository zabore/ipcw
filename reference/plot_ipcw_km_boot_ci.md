# Plot IPCW Kaplan-Meier curves with bootstrap percentile confidence intervals

Combines bootstrap survival curves from
[`get_ipcw_boot()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot.md)
with point estimates from the original IPCW-weighted data to produce a
Kaplan-Meier plot with 95% percentile confidence intervals.

## Usage

``` r
plot_ipcw_km_boot_ci(
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
  [`get_ipcw_boot()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot.md).

- orig_data:

  A single long-format data frame for the original (non-bootstrapped)
  dataset, as returned by
  [`get_ipcw_wgt()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_wgt.md).
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
set.seed(1)
dat      <- sim_data_se(n = 200)
dat_long <- get_ipcw_wgt(dat)
boot_list <- get_ipcw_boot(dat, B = 50)
plot_ipcw_km_boot_ci(boot_list, dat_long, pre_times = seq(0, 500, 10))
} # }
```
