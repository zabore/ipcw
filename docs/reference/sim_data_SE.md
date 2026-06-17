# Simulate single-event survival data with informative censoring

Generates a dataset where event times follow a Weibull distribution
parameterized by a binary treatment covariate `x`, and censoring depends
on an ancillary biomarker `W2` through a shared gamma frailty, inducing
informative censoring using Clayton's copula. This is the
data-generating mechanism used in the single-event guided example.

## Usage

``` r
sim_data_se(
  n = 500,
  alpha = 0.05,
  x_prop = 0.5,
  a = 2,
  sigma = 500,
  beta = log(0.25),
  lambda = 0.01,
  phi = -5
)
```

## Arguments

- n:

  Integer. Number of subjects to simulate. Default is 500.

- alpha:

  Numeric. Shape parameter for the gamma frailty distribution, which has
  a fixed rate of 1. Here, `alpha` denotes the level of dependence
  between the survival time and the covariate that drives censoring. A
  smaller `alpha` leads to a stronger dependence. This parameter is also
  related to Kendall's correlation coefficient, \\\tau\\, such that
  \\\tau = (1/\alpha)/(1/\alpha + 2)\\. Default is `0.05`.

- x_prop:

  Numeric in (0, 1). Probability of treatment (`x = 1`). Default is
  `0.5`.

- a:

  Numeric. Weibull shape parameter for the event time distribution.
  Default is `2`.

- sigma:

  Numeric. Weibull scale parameter for the event time distribution.
  Default is `500`.

- beta:

  Numeric. Log hazard ratio for the treatment effect on event time.
  Default is `log(0.25)`.

- lambda:

  Numeric. Baseline censoring rate when W2=0. Default is `0.01`.

- phi:

  Numeric. The log hazard ratio for the association between a one-unit
  increase in `W2` and censoring time. Negative values mean higher `W2`
  leads to a lower censoring rate (i.e., longer follow-up for high-`W2`
  subjects). Default is `-5`.

## Value

A data frame with columns:

- S:

  True (latent) event time.

- t:

  Observed time (minimum of `S` and the censoring time).

- delta:

  Event indicator: 1 = event occurred, 0 = censored.

- x:

  Binary treatment indicator (0 or 1).

- W2:

  Ancillary biomarker covariate that drives informative censoring.

## Examples

``` r
set.seed(20240429)
dat <- sim_data_se(n = 500)
table(dat$delta)
#> 
#>   0   1 
#> 198 302 
summary(dat$t)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#>    0.6884  172.7074  441.1759  552.9262  803.2563 2429.5429 
```
