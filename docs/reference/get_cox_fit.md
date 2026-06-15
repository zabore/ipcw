# Fit a standard (unweighted) Cox proportional hazards model

Convenience wrapper around `survival::coxph()` for the covariate `x`
using counting-process (`tstart`, `tstop`) time variables. Returns
log-hazard ratio estimates with standard errors and the exponentiated
hazard ratio with 95% confidence intervals.

## Usage

``` r
get_cox_fit(data)
```

## Arguments

  - data:
    
    A data frame in long (counting-process) format containing columns
    `tstart`, `tstop`, `delta`, and `x`.

## Value

A data frame with one row per term containing:

  - term:
    
    Covariate name.

  - log\_hr:
    
    Log hazard ratio estimate.

  - log\_hr\_se:
    
    Standard error of the log hazard ratio.

  - hr:
    
    Hazard ratio estimate.

  - hr\_ci\_low:
    
    Lower 95% confidence limit.

  - hr\_ci\_high:
    
    Upper 95% confidence limit.

## Examples

``` r
data(single_example_ipcw_dat)
get_cox_fit(single_example_ipcw_dat)
```
