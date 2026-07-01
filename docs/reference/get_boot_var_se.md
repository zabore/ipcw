# Compute bootstrap variance of the log hazard ratio

Compute bootstrap variance of the log hazard ratio

## Usage

``` r
get_boot_var_se(data, B)
```

## Arguments

- data:

  A data frame with a column `log_hr` containing bootstrap log hazard
  ratio estimates, as returned by multiple calls to
  [`get_ipcw_cox_fit_se()`](https://www.emilyzabor.com/ipcw/reference/get_ipcw_cox_fit_se.md)
  or
  [`get_cox_fit_se()`](https://www.emilyzabor.com/ipcw/reference/get_cox_fit_se.md).

- B:

  Integer. Number of bootstrap samples used.

## Value

A single numeric value: the bootstrap variance estimate.

## Examples

``` r
# Toy example with 10 simulated bootstrap log HRs
boot_results <- data.frame(log_hr = rnorm(10, mean = -0.5, sd = 0.2))
get_boot_var_se(boot_results, B = 10)
#> [1] 0.05129527
```
