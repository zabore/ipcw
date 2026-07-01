# Compute bootstrap percentile confidence intervals for competing risks estimates

Computes the 2.5th and 97.5th percentiles of the bootstrap distribution
for each column of a matrix of bootstrap estimates. Intended for
competing risks quantities (e.g. cumulative incidence at a set of time
points, or Fine-Gray regression terms), where each bootstrap replicate
produces more than one estimate, unlike the single scalar log hazard
ratio handled by
[`get_boot_pci_se()`](https://www.emilyzabor.com/ipcw/reference/get_boot_pci_se.md).

## Usage

``` r
get_boot_pci_cr(boot_mat)
```

## Arguments

- boot_mat:

  A numeric matrix of bootstrap estimates, with one row per bootstrap
  replicate and one column per estimand (e.g. time point or regression
  term), as returned by row-binding the results of multiple calls to
  [`cuminc_waverage_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_waverage_cr.md),
  [`cuminc_ipcw_cr()`](https://www.emilyzabor.com/ipcw/reference/cuminc_ipcw_cr.md),
  or
  [`fg_weighted_cr()`](https://www.emilyzabor.com/ipcw/reference/fg_weighted_cr.md).

## Value

A 2-row numeric matrix with rows `"lower"` and `"upper"` giving the
2.5th and 97.5th percentiles of the bootstrap distribution for each
column of `boot_mat`. Columns with any `NA` bootstrap estimate return
`NA` for both bounds.

## Examples

``` r
set.seed(1)
boot_mat <- matrix(rnorm(500 * 3), nrow = 500, ncol = 3)
get_boot_pci_cr(boot_mat)
#>            [,1]      [,2]      [,3]
#> lower -1.995029 -2.242602 -1.840346
#> upper  2.030194  2.003819  2.046006
```
