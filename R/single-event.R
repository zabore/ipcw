#' Simulate single-event survival data with informative censoring
#'
#' Generates a dataset where event times follow a Weibull distribution
#' parameterized by a binary treatment covariate `x`, and censoring depends
#' on an ancillary biomarker `W2` through a shared gamma frailty, inducing
#' informative censoring. This is the data-generating mechanism used in the
#' single-event guided example.
#'
#' @param n Integer. Number of subjects to simulate. Default is 500.
#' @param alpha Numeric. Shape parameter for the gamma frailty distribution,
#'   also used as the Pareto exponent for `W1` and `W2`. Default is `0.05`.
#' @param x_prop Numeric in (0, 1). Probability of treatment (`x = 1`).
#'   Default is `0.5`.
#' @param a Numeric. Weibull shape parameter for the event time distribution.
#'   Default is `2`.
#' @param sigma Numeric. Weibull scale parameter for the event time
#'   distribution. Default is `500`.
#' @param beta Numeric. Log hazard ratio for the treatment effect on event
#'   time. Default is `log(0.25)`.
#' @param lambda Numeric. Baseline censoring rate. Default is `0.01`.
#' @param phi Numeric. Effect of `W2` on the censoring rate (log scale).
#'   Negative values mean higher `W2` leads to a lower censoring rate
#'   (i.e., longer follow-up for high-`W2` subjects). Default is `-5`.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{S}{True (latent) event time.}
#'     \item{t}{Observed time (minimum of `S` and the censoring time).}
#'     \item{delta}{Event indicator: 1 = event occurred, 0 = censored.}
#'     \item{x}{Binary treatment indicator (0 or 1).}
#'     \item{W2}{Ancillary biomarker covariate that drives informative
#'       censoring.}
#'   }
#'
#' @examples
#' set.seed(20240429)
#' dat <- sim_data_SE(n = 500)
#' table(dat$delta)
#' summary(dat$t)
#'
#' @importFrom stats rbinom rexp rgamma
#' @export
sim_data_SE <- function(n = 500, alpha = 0.05, x_prop = 0.5,
                        a = 2, sigma = 500, beta = log(0.25),
                        lambda = 0.01, phi = -5) {
  Y1 <- rexp(n, rate = 1)
  Y2 <- rexp(n, rate = 1)
  Z  <- rgamma(n, shape = alpha, rate = 1)

  W1 <- (1 + Y1 / Z)^(-alpha)
  W2 <- (1 + Y2 / Z)^(-alpha)

  x <- rbinom(n, 1, prob = x_prop)

  S <- sigma * (-log(1 - W1) / exp(beta * x))^(1 / a)
  C <- rexp(n, rate = lambda * exp(phi * W2))

  t     <- pmin(S, C)
  delta <- as.integer(S <= C)

  data.frame(S = S, t = t, delta = delta, x = x, W2 = W2)
}


#' Compute IPCW weights for single-event survival data
#'
#' Fits a Cox model for the censoring distribution and returns the original
#' dataset in counting-process (long) format with an unstabilized IPCW weight
#' column appended. The dataset must contain columns named `t` (event/censoring
#' time), `delta` (event indicator, 1 = event), and `W2` (the covariate used to
#' model the censoring distribution).
#'
#' @param data A data frame with columns `t`, `delta`, and `W2`, as produced by
#'   the data-generation code in the package vignette.
#'
#' @return A data frame in long (counting-process) format with columns
#'   `tstart`, `tstop`, `delta`, `id`, `wgt` (the unstabilized IPCW weight),
#'   and all original columns.
#'
#' @examples
#' data(single_example_dat)
#' dat_long <- get_ipcw_wgt(single_example_dat)
#' head(dat_long)
#'
#' @importFrom survival survSplit coxph Surv
#' @importFrom dplyr mutate arrange select rename full_join case_when row_number
#' @importFrom tibble add_column
#' @importFrom stats predict
#' @export
get_ipcw_wgt <- function(data) {

  times <- sort(unique(data$t[data$delta == 0]))

  dat_prep <-
    data |>
    mutate(
      censor = 1 - delta,
      tstart = 0,
      id = row_number()
    )

  dat_long1 <-
    survSplit(
      dat_prep,
      cut = times,
      end = "t",
      start = "tstart",
      event = "delta"
    ) |>
    arrange(id, t)

  dat_long2 <-
    survSplit(
      dat_prep,
      cut = times,
      end = "t",
      start = "tstart",
      event = "censor"
    ) |>
    arrange(id, t)

  dat_long0 <-
    dat_long1 |>
    select(-censor) |>
    add_column(censor = dat_long2$censor) |>
    rename(tstop = t)

  cens_mod <- coxph(Surv(tstart, tstop, censor) ~ W2,
                    data = dat_long0, timefix = FALSE)

  dat_long3 <-
    dat_long0 |>
    mutate(
      tstop = tstart,
      tstart = 0,
      inv_wgt = case_when(tstop == 0 ~ 1)
    )

  dat_long3$inv_wgt[is.na(dat_long3$inv_wgt)] <-
    exp(-predict(cens_mod,
                 newdata = dat_long3[is.na(dat_long3$inv_wgt), ],
                 type = "expected"))

  dat_long3$wgt <- 1 / dat_long3$inv_wgt

  dat_long0 |>
    full_join(
      dat_long3 |> select(id, tstop, wgt),
      by = c("id", "tstart" = "tstop")
    )
}


#' Fit an IPCW-weighted Cox proportional hazards model
#'
#' Fits a weighted Cox model with a robust sandwich variance and returns
#' log-hazard ratio estimates together with both model-based and robust
#' standard errors, plus the exponentiated hazard ratio and 95% confidence
#' interval based on the robust standard error.
#'
#' @param data A data frame in long (counting-process) format, as returned by
#'   [get_ipcw_wgt()]. Must contain columns `tstart`, `tstop`, `delta`, `id`,
#'   and the weight column named by `weight`.
#' @param weight A character string giving the name of the weight column in
#'   `data`.
#'
#' @return A data frame with one row per term containing:
#'   \describe{
#'     \item{term}{Covariate name.}
#'     \item{log_hr}{Log hazard ratio estimate.}
#'     \item{log_hr_se}{Model-based standard error of the log hazard ratio.}
#'     \item{log_hr_rob_se}{Robust (sandwich) standard error of the log hazard ratio.}
#'     \item{hr}{Hazard ratio estimate.}
#'     \item{hr_ci_low}{Lower 95% confidence limit (based on robust SE).}
#'     \item{hr_ci_high}{Upper 95% confidence limit (based on robust SE).}
#'   }
#'
#' @examples
#' data(single_example_ipcw_dat)
#' get_ipcw_cox_fit(single_example_ipcw_dat, weight = "wgt")
#'
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom dplyr full_join rename
#' @export
get_ipcw_cox_fit <- function(data, weight) {
  ipcw_cox_fit <- coxph(Surv(tstart, tstop, delta) ~ x + cluster(id),
                        data = data,
                        weights = data[[weight]],
                        timefix = FALSE)
  full_join(
    tidy(ipcw_cox_fit)[, 1:4] |>
      rename(log_hr = estimate, log_hr_se = std.error, log_hr_rob_se = robust.se),
    tidy(ipcw_cox_fit, exponentiate = TRUE, conf.int = TRUE)[, c(1, 2, 7, 8)] |>
      rename(hr = estimate, hr_ci_low = conf.low, hr_ci_high = conf.high),
    by = "term"
  )
}


#' Estimate IPCW Kaplan-Meier survival probabilities by binary covariate
#'
#' Fits a weighted Kaplan-Meier estimator stratified by the binary covariate `x`
#' and returns survival probabilities evaluated at a pre-specified set of times.
#'
#' @param data A data frame in long (counting-process) format with an IPCW
#'   weight column `wgt`, as returned by [get_ipcw_wgt()].
#' @param pre_times Numeric vector of times at which to evaluate survival
#'   probabilities. Defaults to `seq(0, 50, 1)`.
#'
#' @return A tibble with columns `time`, `surv` (survival probability), and
#'   `x` (stratum).
#'
#' @examples
#' data(single_example_ipcw_dat)
#' get_ipcw_km_prob_x(single_example_ipcw_dat, pre_times = seq(0, 500, 10))
#'
#' @importFrom survival survfit Surv
#' @importFrom tibble tibble
#' @export
get_ipcw_km_prob_x <- function(data, pre_times = seq(0, 50, 1)) {
  ipcw_km_surv_fit <- survfit(Surv(tstart, tstop, delta) ~ x, data = data,
                              weights = wgt, timefix = FALSE)
  tibble(
    time = summary(ipcw_km_surv_fit, times = pre_times)$time,
    surv = summary(ipcw_km_surv_fit, times = pre_times)$surv,
    x    = summary(ipcw_km_surv_fit, times = pre_times)$strata
  )
}


#' Fit a standard (unweighted) Cox proportional hazards model
#'
#' Convenience wrapper around [survival::coxph()] for the covariate `x` using
#' counting-process (`tstart`, `tstop`) time variables. Returns log-hazard
#' ratio estimates with standard errors and the exponentiated hazard ratio with
#' 95% confidence intervals.
#'
#' @param data A data frame in long (counting-process) format containing
#'   columns `tstart`, `tstop`, `delta`, and `x`.
#'
#' @return A data frame with one row per term containing:
#'   \describe{
#'     \item{term}{Covariate name.}
#'     \item{log_hr}{Log hazard ratio estimate.}
#'     \item{log_hr_se}{Standard error of the log hazard ratio.}
#'     \item{hr}{Hazard ratio estimate.}
#'     \item{hr_ci_low}{Lower 95% confidence limit.}
#'     \item{hr_ci_high}{Upper 95% confidence limit.}
#'   }
#'
#' @examples
#' data(single_example_ipcw_dat)
#' get_cox_fit(single_example_ipcw_dat)
#'
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom dplyr full_join rename
#' @export
get_cox_fit <- function(data) {
  cox_fit <- coxph(Surv(tstart, tstop, delta) ~ x, data = data,
                   timefix = FALSE)
  full_join(
    tidy(cox_fit)[, 1:3] |> rename(log_hr = estimate, log_hr_se = std.error),
    tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE)[, c(1, 2, 6, 7)] |>
      rename(hr = estimate, hr_ci_low = conf.low, hr_ci_high = conf.high),
    by = "term"
  )
}


#' Compute bootstrap variance of the log hazard ratio
#'
#' @param data A data frame with a column `log_hr` containing bootstrap
#'   log hazard ratio estimates, as returned by multiple calls to
#'   [get_ipcw_cox_fit()] or [get_cox_fit()].
#' @param B Integer. Number of bootstrap samples used.
#'
#' @return A single numeric value: the bootstrap variance estimate.
#'
#' @examples
#' # Toy example with 10 simulated bootstrap log HRs
#' boot_results <- data.frame(log_hr = rnorm(10, mean = -0.5, sd = 0.2))
#' get_boot_var(boot_results, B = 10)
#'
#' @export
get_boot_var <- function(data, B) {
  sum((data$log_hr - mean(data$log_hr))^2) / B
}


#' Compute bootstrap percentile confidence interval for the log hazard ratio
#'
#' @param data A data frame with a column `log_hr` containing bootstrap
#'   log hazard ratio estimates, as returned by multiple calls to
#'   [get_ipcw_cox_fit()] or [get_cox_fit()].
#'
#' @return A named numeric vector of length 2 giving the 2.5th and 97.5th
#'   percentiles of the bootstrap distribution.
#'
#' @examples
#' boot_results <- data.frame(log_hr = rnorm(500, mean = -0.5, sd = 0.2))
#' get_boot_pci(boot_results)
#'
#' @importFrom stats quantile
#' @export
get_boot_pci <- function(data) {
  quantile(data$log_hr, c(0.025, 0.975))
}
