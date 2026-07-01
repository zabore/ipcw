#' Simulate single-event survival data with informative censoring
#'
#' Generates a dataset where event times follow a Weibull distribution
#' parameterized by a binary treatment covariate `x`, and censoring depends
#' on an ancillary biomarker `W2` through a shared gamma frailty, inducing
#' informative censoring using Clayton's copula. 
#' This is the data-generating mechanism used in the
#' single-event guided example.
#'
#' @param n Integer. Number of subjects to simulate. Default is 500.
#' @param alpha Numeric. Shape parameter for the gamma frailty distribution, 
#' which has a fixed rate of 1. Here, `alpha` denotes the level of dependence 
#' between the survival time and the covariate that drives censoring. 
#' A smaller `alpha` leads to a stronger dependence. 
#' This parameter is also related to Kendall's correlation coefficient, 
#' \eqn{\tau}, such that \eqn{\tau = (1/\alpha)/(1/\alpha + 2)}.
#' Default is `0.05`.
#' @param x_prop Numeric in (0, 1). Probability of treatment (`x = 1`).
#'   Default is `0.5`.
#' @param a Numeric. Weibull shape parameter for the event time distribution.
#'   Default is `2`.
#' @param sigma Numeric. Weibull scale parameter for the event time
#'   distribution. Default is `500`.
#' @param beta Numeric. Log hazard ratio for the treatment effect on event
#'   time. Default is `log(0.25)`.
#' @param lambda Numeric. Baseline censoring rate when W2=0. Default is `0.01`.
#' @param phi Numeric. The log hazard ratio for the association between a 
#' one-unit increase in `W2` and censoring time. Negative values mean higher 
#' `W2` leads to a lower censoring rate 
#' (i.e., longer follow-up for high-`W2` subjects). Default is `-5`.
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
#' dat <- sim_data_se(n = 500)
#' table(dat$delta)
#' summary(dat$t)
#'
#' @importFrom stats rbinom rexp rgamma
#' @export
sim_data_se <- function(n = 500, alpha = 0.05, x_prop = 0.5,
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
#' column appended.
#'
#' @param data A data frame containing the follow-up time, event indicator, and
#'   censoring covariate columns.
#' @param time_var Character string. Name of the observed follow-up time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column
#'   (1 = event, 0 = censored). Default is `"delta"`.
#' @param cens_cov Character string. Name of the covariate column used in the
#'   Cox model for the censoring distribution. Default is `"W2"`.
#'
#' @return A data frame in long (counting-process) format with columns
#'   `tstart`, `tstop`, `delta` (event indicator), `id`, `wgt` (the
#'   unstabilized IPCW weight), and all original columns.
#'
#' @examples
#' set.seed(20240429)
#' dat <- sim_data_se(n = 500)
#' dat_long <- get_ipcw_wgt_se(dat)
#' head(dat_long)
#'
#' @importFrom survival survSplit coxph Surv
#' @importFrom dplyr mutate arrange select full_join case_when row_number
#' @importFrom tibble add_column
#' @importFrom stats predict
#' @export
get_ipcw_wgt_se <- function(data, time_var = "t", event_var = "delta",
                          cens_cov = "W2") {
  times <- sort(unique(data[[time_var]][data[[event_var]] == 0]))

  dat_prep <- 
    data |>
    mutate(
      censor = 1L - .data[[event_var]],
      tstart = 0,
      id     = row_number()
    )

  dat_long1 <- 
    survSplit(
      dat_prep, 
      cut = times, 
      end = time_var,
      start = "tstart", 
      event = event_var) |>
    arrange(id, .data[[time_var]])

  dat_long2 <- 
    survSplit(
      dat_prep, 
      cut = times, 
      end = time_var,
      start = "tstart", 
      event = "censor") |>
    arrange(id, .data[[time_var]])

  dat_long0 <- 
    dat_long1 |>
    select(-censor) |>
    add_column(censor = dat_long2$censor)

  names(dat_long0)[names(dat_long0) == time_var] <- "tstop"
  
  if (event_var != "delta")
    names(dat_long0)[names(dat_long0) == event_var] <- "delta"

  cens_formula <- as.formula(paste("Surv(tstart, tstop, censor) ~", cens_cov))
  
  cens_mod <- coxph(cens_formula, data = dat_long0, timefix = FALSE)

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
      dat_long3[, c("id", "tstop", "wgt")],
      by = c("id", "tstart" = "tstop")
    )
}


#' Estimate IPCW Kaplan-Meier survival probabilities by covariate
#'
#' Fits a weighted Kaplan-Meier estimator stratified by the specified covariate
#' and returns survival probabilities evaluated at a pre-specified set of times.
#'
#' @param data A data frame in long (counting-process) format with an IPCW
#'   weight column, as returned by [get_ipcw_wgt_se()]. Must contain time columns
#'   named `tstart` and `tstop` for the counting process time intervals.
#' @param covariate Character string. Name of the stratification covariate
#'   column. Default is `"x"`.
#' @param weight_var Character string. Name of the IPCW weight column. Default
#'   is `"wgt"`.
#' @param event_var Character string. Name of the event indicator column
#'   (1 = event, 0 = censored). Default is `"delta"`.
#' @param pre_times Numeric vector of times at which to evaluate survival
#'   probabilities. Defaults to `seq(0, 50, 1)`. Choose the second number
#'   to cover the entire range of times observed in `data`
#'
#' @return A tibble with columns `time`, `surv` (survival probability), and a
#'   column named after `covariate` containing the stratum labels.
#'
#' @examples
#' set.seed(20240429)
#' dat <- sim_data_se(n = 500)
#' dat_long <- get_ipcw_wgt_se(dat)
#' get_ipcw_km_prob_x_se(dat_long, pre_times = seq(0, 2429, 1))
#'
#' @importFrom survival survfit Surv
#' @importFrom tibble tibble
#' @export
get_ipcw_km_prob_x_se <- function(data, covariate = "x", weight_var = "wgt", 
                               event_var = "delta",
                               pre_times = seq(0, 50, 1)) {
  km_formula <- as.formula(
    paste("Surv(tstart, tstop,", event_var, ") ~", covariate)
    )
  
  ipcw_km_surv_fit <- survfit(km_formula, data = data,
                              weights = data[[weight_var]], 
                              timefix = FALSE)
  
  result <- tibble(
    time   = summary(ipcw_km_surv_fit, times = pre_times)$time,
    surv   = summary(ipcw_km_surv_fit, times = pre_times)$surv,
    strata = summary(ipcw_km_surv_fit, times = pre_times)$strata
  )
  names(result)[3] <- covariate
  result
}




#' Fit an IPCW-weighted Cox proportional hazards model
#'
#' Fits a weighted Cox model with a robust sandwich variance and returns
#' log-hazard ratio estimates together with both model-based and robust
#' standard errors, plus the exponentiated hazard ratio and 95% confidence
#' interval based on the robust standard error.
#'
#' @param data A data frame in long (counting-process) format, as returned by
#'   [get_ipcw_wgt_se()]. Must contain columns `tstart`, `tstop`, `delta`, `id`,
#'   the covariate named by `covariate`, and the weight column named by
#'   `weight`.
#' @param covariate Character string. Name of the predictor covariate column.
#'   Default is `"x"`.
#' @param weight Character string. Name of the weight column. Default is
#'   `"wgt"`.
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
#' set.seed(20240429)
#' dat <- sim_data_se(n = 500)
#' dat_long <- get_ipcw_wgt_se(dat)
#' get_ipcw_cox_fit_se(dat_long, weight = "wgt")
#'
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom dplyr full_join rename
#' @export
get_ipcw_cox_fit_se <- function(data, covariate = "x", weight = "wgt") {
  cox_formula <- as.formula(
    paste("Surv(tstart, tstop, delta) ~", covariate, "+ cluster(id)")
  )
  ipcw_cox_fit <- coxph(cox_formula, data = data,
                        weights = data[[weight]], timefix = FALSE)
  full_join(
    tidy(ipcw_cox_fit)[, 1:4] |>
      rename(log_hr = estimate, log_hr_se = std.error, log_hr_rob_se = robust.se),
    tidy(ipcw_cox_fit, exponentiate = TRUE, conf.int = TRUE)[, c(1, 2, 7, 8)] |>
      rename(hr = estimate, hr_ci_low = conf.low, hr_ci_high = conf.high),
    by = "term"
  )
}





#' Fit a standard (unweighted) Cox proportional hazards model
#'
#' Convenience wrapper around [survival::coxph()] using counting-process
#' (`tstart`, `tstop`) time variables. Returns log-hazard ratio estimates with
#' standard errors and the exponentiated hazard ratio with 95% confidence
#' intervals.
#'
#' @param data A data frame in long (counting-process) format containing
#'   columns `tstart`, `tstop`, `delta`, and the covariate named by
#'   `covariate`.
#' @param covariate Character string. Name of the predictor covariate column.
#'   Default is `"x"`.
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
#' set.seed(20240429)
#' dat <- sim_data_se(n = 500)
#' dat_long <- get_ipcw_wgt_se(dat)
#' get_cox_fit_se(dat_long)
#'
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom dplyr full_join rename
#' @export
get_cox_fit_se <- function(data, covariate = "x") {
  cox_formula <- as.formula(paste("Surv(tstart, tstop, delta) ~", covariate))
  cox_fit <- coxph(cox_formula, data = data, timefix = FALSE)
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
#'   [get_ipcw_cox_fit_se()] or [get_cox_fit_se()].
#' @param B Integer. Number of bootstrap samples used.
#'
#' @return A single numeric value: the bootstrap variance estimate.
#'
#' @examples
#' # Toy example with 10 simulated bootstrap log HRs
#' boot_results <- data.frame(log_hr = rnorm(10, mean = -0.5, sd = 0.2))
#' get_boot_var_se(boot_results, B = 10)
#'
#' @export
get_boot_var_se <- function(data, B) {
  sum((data$log_hr - mean(data$log_hr))^2) / B
}


#' Compute bootstrap percentile confidence interval for the log hazard ratio
#'
#' @param data A data frame with a column `log_hr` containing bootstrap
#'   log hazard ratio estimates, as returned by multiple calls to
#'   [get_ipcw_cox_fit_se()] or [get_cox_fit_se()].
#'
#' @return A named numeric vector of length 2 giving the 2.5th and 97.5th
#'   percentiles of the bootstrap distribution.
#'
#' @examples
#' boot_results <- data.frame(log_hr = rnorm(500, mean = -0.5, sd = 0.2))
#' get_boot_pci_se(boot_results)
#'
#' @importFrom stats quantile
#' @export
get_boot_pci_se <- function(data) {
  quantile(data$log_hr, c(0.025, 0.975))
}


#' Bootstrap IPCW weighted data for single-event survival analysis
#'
#' Draws `B` bootstrap samples from the original wide-format data, fits IPCW
#' weights to each sample via [get_ipcw_wgt_se()], and returns the results as a
#' list of long-format data frames. The list can be passed to
#' [plot_ipcw_km_boot_ci_se()] for plotting or used directly to bootstrap other
#' quantities (e.g. Cox hazard ratios via [get_ipcw_cox_fit_se()]).
#'
#' @param data A wide-format data frame containing the columns specified by
#'   `time_var`, `event_var`, and `cens_cov`.
#' @param B Integer. Number of bootstrap samples. Default is `500`.
#' @param time_var Character string. Name of the observed follow-up time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column
#'   (1 = event, 0 = censored). Default is `"delta"`.
#' @param cens_cov Character string. Name of the covariate column used to model
#'   the censoring distribution. Default is `"W2"`.
#' @param seed Optional integer seed for reproducibility. Default is `NULL`.
#'
#' @return A list of `B` data frames, each in long (counting-process) format
#'   as returned by [get_ipcw_wgt_se()].
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- sim_data_se(n = 200)
#' boot_list <- get_ipcw_boot_se(dat, B = 50)
#' }
#'
#' @export
get_ipcw_boot_se <- function(data, B = 500, time_var = "t", event_var = "delta",
                           cens_cov = "W2", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  lapply(seq_len(B), function(i) {
    boot_dat <- data[sample(nrow(data), nrow(data), replace = TRUE), ]
    get_ipcw_wgt_se(boot_dat, time_var = time_var, event_var = event_var,
                  cens_cov = cens_cov)
  })
}


#' Plot IPCW Kaplan-Meier curves with bootstrap percentile confidence intervals
#'
#' Combines bootstrap survival curves from [get_ipcw_boot_se()] with point
#' estimates from the original IPCW-weighted data to produce a Kaplan-Meier
#' plot with 95% percentile confidence intervals.
#'
#' @param boot_data A list of long-format data frames as returned by
#'   [get_ipcw_boot_se()].
#' @param orig_data A single long-format data frame for the original
#'   (non-bootstrapped) dataset, as returned by [get_ipcw_wgt_se()]. Used for
#'   the point estimates.
#' @param pre_times Numeric vector of times at which to evaluate survival
#'   probabilities. Default is `seq(0, 50, 1)`.
#' @param covariate Character string. Name of the stratification covariate
#'   column. Default is `"x"`.
#' @param weight_var Character string. Name of the IPCW weight column. Default
#'   is `"wgt"`.
#' @param event_var Character string. Name of the event indicator column in
#'   the long-format data. Default is `"delta"`.
#'
#' @return A [ggplot2::ggplot()] object. Add layers or themes to customise.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat      <- sim_data_se(n = 200)
#' dat_long <- get_ipcw_wgt_se(dat)
#' boot_list <- get_ipcw_boot_se(dat, B = 50)
#' plot_ipcw_km_boot_ci_se(boot_list, dat_long, pre_times = seq(0, 500, 10))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_step geom_ribbon labs theme_minimal
#' @importFrom dplyr bind_rows group_by across all_of summarise left_join
#' @importFrom stats quantile
#' @export
plot_ipcw_km_boot_ci_se <- function(boot_data, orig_data, pre_times = seq(0, 50, 1),
                                  covariate = "x", weight_var = "wgt",
                                  event_var = "delta") {
  boot_results <- lapply(boot_data, function(long_dat) {
    get_ipcw_km_prob_x_se(
      long_dat, 
      covariate = covariate,
      weight_var = weight_var, 
      pre_times = pre_times)
  })

  ci_summ <- bind_rows(boot_results) |>
    group_by(across(all_of(c("time", covariate)))) |>
    summarise(
      lpci = quantile(surv, 0.025),
      upci = quantile(surv, 0.975),
      .groups = "drop"
    )

  orig_est <- get_ipcw_km_prob_x_se(orig_data, covariate = covariate,
                                   weight_var = weight_var,
                                   pre_times = pre_times)

  plot_dat <- left_join(orig_est, ci_summ, by = c("time", covariate))

  ggplot(plot_dat, aes(x = time,
                       color = .data[[covariate]],
                       fill  = .data[[covariate]])) +
    geom_ribbon(aes(ymin = lpci, ymax = upci), alpha = 0.2, linetype = "blank") +
    geom_step(aes(y = surv)) +
    labs(x = "Time", y = "Survival Probability") +
    theme_minimal()
}
