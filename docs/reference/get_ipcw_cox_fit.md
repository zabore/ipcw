# Fit an IPCW-weighted Cox proportional hazards model

Fits a weighted Cox model with a robust sandwich variance and returns
log-hazard ratio estimates together with both model-based and robust
standard errors, plus the exponentiated hazard ratio and 95% confidence
interval based on the robust standard error.

## Usage

``` r
get_ipcw_cox_fit(data, weight)
```

## Arguments

  - data:
    
    A data frame in long (counting-process) format, as returned by
    `get_ipcw_wgt()`. Must contain columns `tstart`, `tstop`, `delta`,
    `id`, and the weight column named by `weight`.

  - weight:
    
    A character string giving the name of the weight column in `data`.

## Value

A data frame with one row per term containing:

  - term:
    
    Covariate name.

  - log\_hr:
    
    Log hazard ratio estimate.

  - log\_hr\_se:
    
    Model-based standard error of the log hazard ratio.

  - log\_hr\_rob\_se:
    
    Robust (sandwich) standard error of the log hazard ratio.

  - hr:
    
    Hazard ratio estimate.

  - hr\_ci\_low:
    
    Lower 95% confidence limit (based on robust SE).

  - hr\_ci\_high:
    
    Upper 95% confidence limit (based on robust SE).

## Examples

``` r
data(single_example_ipcw_dat)
get_ipcw_cox_fit(single_example_ipcw_dat, weight = "wgt")
```
