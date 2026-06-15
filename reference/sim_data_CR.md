# Simulate competing risks survival data

Generates a dataset with two competing events under a sub-distribution
hazard model parameterized by a four-level baseline covariate `z1`.
Optionally adds independent or covariate-dependent administrative
censoring.

## Usage

``` r
sim_data_CR(
  n = 100,
  censoring = "none",
  beta1 = log(1.5),
  beta2 = log(2.25),
  beta3 = log(3.4),
  p = 0.3
)
```

## Arguments

- n:

  Integer. Number of subjects to simulate. Default is 100.

- censoring:

  Character string specifying the censoring mechanism. One of `"none"`,
  `"independent"`, or `"baseline"` (censoring rate depends on `z1`).
  Default is `"none"`.

- beta1:

  Numeric. Log sub-distribution hazard ratio for `z1 = 1` vs `z1 = 0`
  for event type 1. Default is `log(1.5)`.

- beta2:

  Numeric. Log sub-distribution hazard ratio for `z1 = 2` vs `z1 = 0`
  for event type 1. Default is `log(2.25)`.

- beta3:

  Numeric. Log sub-distribution hazard ratio for `z1 = 3` vs `z1 = 0`
  for event type 1. Default is `log(3.4)`.

- p:

  Numeric in (0, 1). Baseline probability of event type 1 in the
  reference group (`z1 = 0`). Default is `0.3`.

## Value

A data frame with columns:

- z1:

  Four-level factor covariate (levels 0-3).

- delta:

  Factor event indicator with levels `"censor"`, `"event_1"`,
  `"event_2"`.

- t:

  Event or censoring time.

## Examples

``` r
set.seed(42)
dat <- sim_data_CR(n = 200, censoring = "baseline")
table(dat$delta)
#> 
#>  censor event_1 event_2 
#>      74      85      41 
```
