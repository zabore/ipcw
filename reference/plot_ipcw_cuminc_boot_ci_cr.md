# Plot IPCW cumulative incidence with bootstrap percentile confidence intervals

Combines bootstrap cumulative incidence curves from
[`get_ipcw_boot_cr()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot_cr.md)
with the point estimate from the original IPCW-weighted data to produce
a cumulative incidence plot with 95% percentile confidence intervals,
mirroring
[`plot_ipcw_km_boot_ci_se()`](https://www.emilyzabor.com/ipcw/reference/plot_ipcw_km_boot_ci_se.md)
for the single-event case.

## Usage

``` r
plot_ipcw_cuminc_boot_ci_cr(
  boot_data,
  orig_data,
  esttimes = seq(from = 0, to = 10, length.out = 100),
  extend = TRUE,
  covariate = "z1"
)
```

## Arguments

- boot_data:

  A list of long-format data frames with IPCW weights, as returned by
  [`get_ipcw_boot_cr()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_boot_cr.md).

- orig_data:

  A single long-format data frame with IPCW weights for the original
  (non-bootstrapped) dataset, as returned by
  [`add_ipcw_weights_cr()`](https://www.emilyzabor.com/ipcw/reference/add_ipcw_weights_cr.md).
  Used for the point estimate.

- esttimes:

  Numeric vector of times at which to evaluate cumulative incidence.
  Defaults to 100 equally spaced points from 0 to 10.

- extend:

  Logical. Passed to
  [`cuminc_ipcw_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_ipcw_cr.md).
  Default is `TRUE`.

- covariate:

  Character string. Name of the covariate column, used only when
  `extend = FALSE` to compute the truncation time. Default is `"z1"`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. Add layers or themes to customise.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
dat <- sim_data_cr(n = 200, censoring = "baseline")
dat_long  <- add_ipcw_weights_cr(wide_to_long_cr(dat), strat = "no")
boot_list <- get_ipcw_boot_cr(dat, B = 50)
plot_ipcw_cuminc_boot_ci_cr(boot_list, dat_long, esttimes = seq(0, 5, 0.1))
} # }
```
