# Simulate single-event survival data with informative censoring

Generates a dataset where event times follow a Weibull distribution
parameterized by a binary treatment covariate `x`, and censoring depends
on an ancillary biomarker `W2` through a shared gamma frailty, inducing
informative censoring. This is the data-generating mechanism used in the
single-event guided example.

## Usage

``` r
sim_data_SE(
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
    
    Numeric. Shape parameter for the gamma frailty distribution, also
    used as the Pareto exponent for `W1` and `W2`. Default is `0.05`.

  - x\_prop:
    
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
    
    Numeric. Baseline censoring rate. Default is `0.01`.

  - phi:
    
    Numeric. Effect of `W2` on the censoring rate (log scale). Negative
    values mean higher `W2` leads to a lower censoring rate (i.e.,
    longer follow-up for high-`W2` subjects). Default is `-5`.

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
dat <- sim_data_SE(n = 500)
table(dat$delta)
summary(dat$t)
```
