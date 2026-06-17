# Compute bootstrap percentile confidence interval for the log hazard ratio

Compute bootstrap percentile confidence interval for the log hazard
ratio

## Usage

``` r
get_boot_pci(data)
```

## Arguments

- data:

  A data frame with a column `log_hr` containing bootstrap log hazard
  ratio estimates, as returned by multiple calls to
  [`get_ipcw_cox_fit()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_cox_fit.md)
  or
  [`get_cox_fit()`](https://www.emilyzabor.com/ipcw/reference/get_cox_fit.md).

## Value

A named numeric vector of length 2 giving the 2.5th and 97.5th
percentiles of the bootstrap distribution.

## Examples

``` r
boot_results <- data.frame(log_hr = rnorm(500, mean = -0.5, sd = 0.2))
get_boot_pci(boot_results)
#>       2.5%      97.5% 
#> -0.9221499 -0.1397219 
```
