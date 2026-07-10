#' Simulate competing risks survival data
#'
#' Generates a dataset with two competing events under a sub-distribution
#' hazard model parameterized by a four-level baseline covariate `z1`.
#' Optionally adds independent or covariate-dependent censoring.
#'
#' @param n Integer. Number of subjects to simulate. Default is 100.
#' @param censoring Character string specifying the censoring mechanism.
#'   One of `"none"`, `"independent"`, or `"baseline"` (censoring rate depends
#'   on `z1`). Default is `"none"`.
#' @param beta1 Numeric. Log sub-distribution hazard ratio for `z1 = 1` vs
#'   `z1 = 0` for event type 1. Default is `log(1.5)`.
#' @param beta2 Numeric. Log sub-distribution hazard ratio for `z1 = 2` vs
#'   `z1 = 0` for event type 1. Default is `log(2.25)`.
#' @param beta3 Numeric. Log sub-distribution hazard ratio for `z1 = 3` vs
#'   `z1 = 0` for event type 1. Default is `log(3.4)`.
#' @param p Numeric in (0, 1). Baseline probability of event type 1 in the
#'   reference group (`z1 = 0`). Default is `0.3`.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{z1}{Four-level factor covariate (levels 0-3).}
#'     \item{delta}{Factor event indicator with levels `"censor"`, `"event_1"`,
#'       `"event_2"`.}
#'     \item{t}{Event or censoring time.}
#'   }
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "baseline")
#' table(dat$delta)
#'
#' @importFrom survival survfit Surv
#' @importFrom stats rexp runif uniroot
#' @export
sim_data_cr <- function(n = 100, censoring = "none",
                        beta1 = log(1.5), beta2 = log(2.25), beta3 = log(3.4),
                        p = 0.3) {
  if (!(censoring %in% c("none", "independent", "baseline"))) {
    stop("censoring must be one of none, independent, or baseline")
  }

  z1 <- as.factor(sample(0:3, n, replace = TRUE))

  U <- runif(n)
  lp <- ifelse(z1 == 0, 0,
               ifelse(z1 == 1, beta1,
                      ifelse(z1 == 2, beta2, beta3)))

  sdf <- function(x, lp) 1 - (1 - p * (1 - exp(-x)))^exp(lp)

  delta <- ifelse(U < sapply(lp, function(lp) sdf(Inf, lp)), 1, 2)
  t <- sapply(seq_len(n), function(i) {
    if (delta[i] == 1) {
      uniroot(function(x) sdf(x, lp[i]) - U[i],
              c(0, 100), extendInt = "upX")$root
    } else {
      rexp(1, exp(-1))
    }
  })

  if (censoring == "independent") {
    cens  <- rexp(n, rate = 1 / 8)
    delta <- ifelse(t < cens, delta, 0)
    t     <- ifelse(t < cens, t, cens)
  }

  if (censoring == "baseline") {
    cens  <- rexp(n, rate = (1 / 4) * exp(-1 * (as.numeric(as.character(z1)) - 1.5)))
    delta <- ifelse(t < cens, delta, 0)
    t     <- ifelse(t < cens, t, cens)
  }

  delta <- factor(delta, 0:2, labels = c("censor", "event_1", "event_2"))
  t[t == 0] <- min(t[t > 0]) / 2

  data.frame(z1 = z1, delta = delta, t = t)
}


#' Convert wide competing risks data to long (counting-process) format
#'
#' Splits the wide dataset at every censoring time, creating one or more rows
#' per subject. The resulting long dataset is suitable for fitting Cox models
#' for the censoring distribution and for use with [add_ipcw_weights_cr()].
#'
#' @param dat A wide-format competing risks data frame. Must contain the columns
#'   specified by `time_var`, `event_var`, and `covariate`.
#' @param time_var Character string. Name of the event/censoring time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column
#'   (factor with levels for censoring and the two event types). Default is
#'   `"delta"`.
#' @param covariate Character string. Name of the covariate column. Default is
#'   `"z1"`.
#' @param cens_level Character string. Factor level in `event_var` representing
#'   censoring. Default is `"censor"`.
#' @param event2_level Character string. Factor level in `event_var`
#'   representing the competing event (event type 2). Default is `"event_2"`.
#'
#' @return A data frame in long format with columns `id`, `delta`, `censor`,
#'   `event2_time`, `tstart`, `tstop`, and the covariate column.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 100, censoring = "baseline")
#' dat_long <- wide_to_long_cr(dat)
#' head(dat_long)
#'
#' @importFrom survival survSplit Surv
#' @export
wide_to_long_cr <- function(dat, time_var = "t", event_var = "delta",
                             covariate = "z1", cens_level = "censor",
                             event2_level = "event_2") {
  dat$id          <- seq_len(nrow(dat))
  dat$censor      <- (dat[[event_var]] == cens_level)
  dat$event2_time <- ifelse(dat[[event_var]] == event2_level, dat[[time_var]], NA)
  dat$tstart      <- 0

  times <- sort(unique(dat[[time_var]][dat$censor]))

  data_long1 <- survSplit(dat, cut = times, end = time_var, start = "tstart",
                          event = event_var)
  data_long1 <- data_long1[order(data_long1$id, data_long1[[time_var]]), ]

  data_long2 <- survSplit(dat, cut = times, end = time_var, start = "tstart",
                          event = "censor")
  data_long2 <- data_long2[order(data_long2$id, data_long2[[time_var]]), ]

  data_long1$censor <- data_long2$censor

  names(data_long1)[names(data_long1) == time_var] <- "tstop"
  if (event_var != "delta")
    names(data_long1)[names(data_long1) == event_var] <- "delta"

  data_long1
}


#' Add IPCW weights to competing risks long-format data
#'
#' Estimates the probability of remaining uncensored and appends it as column
#' `p_notcens` to the dataset. Supports two estimation strategies: a Cox
#' proportional hazards model (`strat = "no"`) or non-parametric (KM) estimates
#' within each level of the covariate (`strat = "yes"`).
#'
#' @param data_long A data frame in long format, as returned by
#'   [wide_to_long_cr()]. Must contain columns `tstart`, `tstop`, `censor`, and
#'   the covariate named by `covariate`.
#' @param covariate Character string. Name of the covariate column used to
#'   model the censoring distribution. Default is `"z1"`.
#' @param strat Character. `"no"` (default) fits a single Cox model for the
#'   censoring distribution using `covariate`. `"yes"` estimates the censoring
#'   distribution non-parametrically within each stratum of `covariate`.
#' @param new_data Optional data frame to which weights are applied. If `NULL`
#'   (default), weights are computed for `data_long` itself.
#' @param by.start Logical. If `TRUE` (default), the weight is
#'   P(not censored by `tstart`). If `FALSE`, the weight is
#'   P(not censored during the interval from `tstart` to `tstop`).
#'
#' @return `new_data` (or `data_long` if `new_data` is `NULL`) with an
#'   additional column `p_notcens`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 100, censoring = "baseline")
#' dat_long <- wide_to_long_cr(dat)
#' dat_long <- add_ipcw_weights_cr(dat_long, strat = "no")
#' summary(dat_long$p_notcens)
#'
#' @importFrom survival coxph survfit Surv
#' @importFrom stats predict
#' @export
add_ipcw_weights_cr <- function(data_long, covariate = "z1", strat = "no",
                              new_data = NULL, by.start = TRUE) {
  if (is.null(new_data)) new_data <- data_long

  temp <- new_data[, c(covariate, "tstart", "tstop", "censor")]
  if (by.start) {
    temp$tstop  <- temp$tstart
    temp$tstart <- 0
  }
  temp$p_notcens <- NA

  if (strat == "yes") {
    cov_levels <- levels(data_long[[covariate]])
    cens_mod <- lapply(cov_levels, function(lv) {
      survfit(Surv(tstart, tstop, censor) ~ 1,
              data = data_long[data_long[[covariate]] == lv, ], timefix = FALSE)
    })
    for (i in seq_along(cov_levels)) {
      lv  <- cov_levels[i]
      idx <- temp[[covariate]] == lv
      temp$p_notcens[idx] <-
        summary(cens_mod[[i]], temp$tstop[idx],  extend = TRUE)$surv /
        summary(cens_mod[[i]], temp$tstart[idx], extend = TRUE)$surv
    }
  }

  if (strat == "no") {
    cens_formula <- as.formula(paste("Surv(tstart, tstop, censor) ~", covariate))
    cens_mod <- coxph(cens_formula, data = data_long, timefix = FALSE)
    temp$p_notcens[temp$tstop == 0] <- 1
    temp$p_notcens[is.na(temp$p_notcens)] <-
      exp(-predict(cens_mod,
                   newdata = temp[is.na(temp$p_notcens), ],
                   type = "expected"))
  }

  new_data$p_notcens <- temp$p_notcens
  new_data
}


#' Naive (unweighted) cumulative incidence estimate
#'
#' Estimates the cumulative incidence of event type 1 using the standard
#' Aalen-Johansen estimator without any IPCW adjustment. Serves as a
#' comparison method when censoring is independent of covariates.
#'
#' @param dat A wide-format competing risks data frame containing the columns
#'   specified by `time_var` and `event_var`.
#' @param esttimes Numeric vector of times at which to return estimates.
#' @param time_var Character string. Name of the event/censoring time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column.
#'   Default is `"delta"`.
#'
#' @return A numeric vector of cumulative incidence estimates at `esttimes`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "independent")
#' cuminc_naive_cr(dat, esttimes = seq(0, 5, 0.5))
#'
#' @importFrom survival survfit Surv
#' @export
cuminc_naive_cr <- function(dat, esttimes, time_var = "t", event_var = "delta") {
  naive_formula <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ 1")
  )
  fit <- survfit(naive_formula, data = dat)
  est <- summary(fit, times = esttimes, extend = TRUE)
  est$pstate[, 2]
}


#' Weighted-average (non-parametric) IPCW cumulative incidence estimate
#'
#' Estimates the marginal cumulative incidence of event type 1 by computing
#' stratum-specific Aalen-Johansen estimates and then combining them with
#' sample-proportion weights. This is the non-parametric IPCW approach.
#'
#' @param dat A wide-format competing risks data frame containing the columns
#'   specified by `time_var`, `event_var`, and `covariate`.
#' @param esttimes Numeric vector of times at which to return estimates.
#'   Defaults to 100 equally spaced points from 0 to 10.
#' @param time_var Character string. Name of the event/censoring time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column.
#'   Default is `"delta"`.
#' @param covariate Character string. Name of the stratification covariate
#'   column. Default is `"z1"`.
#'
#' @return A numeric vector of weighted-average cumulative incidence estimates
#'   at `esttimes`. Values beyond the minimum of the stratum-specific maximum
#'   follow-up times are set to `NA`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "baseline")
#' cuminc_waverage_cr(dat, esttimes = seq(0, 5, 0.5))
#'
#' @importFrom survival strata survfit Surv
#' @export
cuminc_waverage_cr <- function(dat,
                             esttimes = seq(from = 0, to = 10, length.out = 100),
                             time_var = "t", event_var = "delta",
                             covariate = "z1") {
  wavg_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ strata(", covariate, ")"
  ))
  strat_fit  <- survfit(wavg_formula, data = dat)
  est        <- summary(strat_fit, times = esttimes, extend = TRUE)
  cov_levels <- levels(dat[[covariate]])
  n_levels   <- length(cov_levels)
  pstate     <- matrix(est$pstate[, 2], nrow = length(esttimes), ncol = n_levels)
  n          <- nrow(dat)
  weighted   <- rowSums(
    sapply(seq_along(cov_levels), function(i) {
      pstate[, i] * sum(dat[[covariate]] == cov_levels[i])
    })
  ) / n
  tau <- min(sapply(cov_levels, function(lv) {
    max(dat[[time_var]][dat[[covariate]] == lv])
  }))
  weighted[esttimes >= tau] <- NA
  weighted
}


#' Cox model IPCW cumulative incidence estimate
#'
#' Estimates the marginal cumulative incidence of event type 1 using a weighted
#' Aalen-Johansen estimator, where the weights are the IPCW weights stored in
#' the `p_notcens` column (as added by [add_ipcw_weights_cr()]).
#'
#' @param data_long A long-format data frame with IPCW weights, as returned by
#'   [add_ipcw_weights_cr()]. Must contain columns `tstart`, `tstop`, `delta`,
#'   `id`, and `p_notcens`.
#' @param esttimes Numeric vector of times at which to return estimates.
#'   Defaults to 100 equally spaced points from 0 to 10.
#' @param extend Logical. If `FALSE`, estimates beyond the minimum of the
#'   stratum-specific maximum follow-up times are set to `NA`. Default is
#'   `TRUE`.
#' @param covariate Character string. Name of the covariate column, used only
#'   when `extend = FALSE` to compute the truncation time. Default is `"z1"`.
#'
#' @return A numeric vector of IPCW cumulative incidence estimates at
#'   `esttimes`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "baseline")
#' dat_long <- wide_to_long_cr(dat)
#' dat_long <- add_ipcw_weights_cr(dat_long, strat = "no")
#' cuminc_ipcw_cr(dat_long, esttimes = seq(0, 5, 0.5))
#'
#' @importFrom survival survfit Surv
#' @export
cuminc_ipcw_cr <- function(data_long,
                         esttimes = seq(from = 0, to = 10, length.out = 100),
                         extend = TRUE, covariate = "z1") {
  fit  <- survfit(Surv(tstart, tstop, delta) ~ 1,
                  data = data_long, weights = 1 / p_notcens,
                  id = id, timefix = FALSE)
  ipcw <- summary(fit, times = esttimes, extend = TRUE)

  if (!extend) {
    cov_levels <- levels(data_long[[covariate]])
    tau <- min(sapply(cov_levels, function(lv) {
      max(data_long$tstop[data_long[[covariate]] == lv])
    }))
    ipcw$pstate[, 2][esttimes >= tau] <- NA
  }
  ipcw$pstate[, 2]
}


#' Prepare long-format data for Fine-Gray regression
#'
#' After a type-2 event, subjects are artificially re-entered into the risk set
#' (as in the Fine-Gray sub-distribution hazard model). This function appends
#' those additional rows, split at every censoring time.
#'
#' @param data_long A long-format data frame, as returned by [wide_to_long_cr()].
#' @param covariate Character string. Name of the covariate column. Default is
#'   `"z1"`.
#' @param event2_level Character string. Factor level in the `delta` column
#'   representing the competing event. Default is `"event_2"`.
#'
#' @return A data frame in Fine-Gray format, with additional rows for subjects
#'   who experienced the competing event, sorted by `id` and `tstart`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 100, censoring = "baseline")
#' dat_long <- wide_to_long_cr(dat)
#' dat_long_fg <- fg_split_cr(dat_long)
#' nrow(dat_long_fg) > nrow(dat_long)
#'
#' @export
fg_split_cr <- function(data_long, covariate = "z1", event2_level = "event_2") {
  times <- sort(unique(data_long$tstop[data_long$censor == 1]))

  event2_dat <- data_long[data_long$delta == event2_level,
                           c("id", covariate, "event2_time")]
  cens_dat <- data.frame(
    tstart = c(0, times),
    tstop  = c(times, Inf),
    censor = 0,
    delta  = "censor"
  )

  fulldat <- merge(event2_dat, cens_dat)
  fulldat <- fulldat[fulldat$event2_time < fulldat$tstop, ]
  fulldat$tstart[fulldat$event2_time > fulldat$tstart] <-
    fulldat$event2_time[fulldat$event2_time > fulldat$tstart]

  missing_cols <- setdiff(names(data_long), names(fulldat))
  fulldat[missing_cols] <- NA

  fulldat <- rbind(data_long, fulldat[names(data_long)])
  fulldat[order(fulldat$id, fulldat$tstart), ]
}


#' Add Fine-Gray weights to Fine-Gray split data
#'
#' Computes the probability of remaining uncensored after a type-2 event
#' and appends it as column `p_notcens_after_death`.
#'
#' @param data_long_fg A data frame in Fine-Gray format, as returned by
#'   [fg_split_cr()].
#' @param covariate Character string. Name of the covariate column. Default is
#'   `"z1"`.
#' @param strat Character. Passed to [add_ipcw_weights_cr()]. `"no"` (default)
#'   uses a Cox model; `"yes"` uses stratum-specific KM estimates.
#'
#' @return `data_long_fg` with an additional column `p_notcens_after_death`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 100, censoring = "baseline")
#' dat_long    <- wide_to_long_cr(dat)
#' dat_long_fg <- fg_split_cr(dat_long)
#' dat_long_fg <- add_fg_weights_cr(dat_long_fg, strat = "no")
#' summary(dat_long_fg$p_notcens_after_death)
#'
#' @export
add_fg_weights_cr <- function(data_long_fg, covariate = "z1", strat = "no") {
  temp <- data.frame(
    tstart = data_long_fg$event2_time,
    tstop  = data_long_fg$tstart,
    censor = NA
  )
  temp[[covariate]] <- data_long_fg[[covariate]]
  temp$p_notcens_after_death <- NA
  temp$p_notcens_after_death[is.na(temp$tstart)] <- 1
  temp$p_notcens_after_death[!is.na(temp$tstart) & temp$tstart >= temp$tstop] <- 1

  still_na <- is.na(temp$p_notcens_after_death)
  if (any(still_na)) {
    ref_data <- data_long_fg[data_long_fg$tstart < data_long_fg$event2_time |
                               is.na(data_long_fg$event2_time), ]
    temp$p_notcens_after_death[still_na] <-
      add_ipcw_weights_cr(ref_data, covariate = covariate, strat = strat,
                       new_data = temp[still_na, ],
                       by.start = FALSE)$p_notcens
  }

  data_long_fg$p_notcens_after_death <- temp$p_notcens_after_death
  data_long_fg
}


#' Naive Fine-Gray sub-distribution hazard regression
#'
#' Fits a Fine-Gray model using the standard [survival::finegray()] approach,
#' without any IPCW adjustment. Serves as a comparison to [fg_weighted_cr()].
#'
#' @param dat A wide-format competing risks data frame containing the columns
#'   specified by `time_var`, `event_var`, and `covariate`.
#' @param time_var Character string. Name of the event/censoring time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column.
#'   Default is `"delta"`.
#' @param covariate Character string. Name of the covariate column. Default is
#'   `"z1"`.
#'
#' @return A matrix with one row per term and two columns: the log
#'   sub-distribution hazard ratio and its standard error.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "independent")
#' fg_naive_cr(dat)
#'
#' @importFrom survival finegray coxph Surv
#' @export
fg_naive_cr <- function(dat, time_var = "t", event_var = "delta",
                     covariate = "z1") {
  fg_formula  <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ."))
  pdata <- finegray(fg_formula, data = dat, timefix = FALSE)
  m1 <- coxph(
    as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", covariate)),
    weights = fgwt, data = pdata, timefix = FALSE
  )
  summary(m1)$coef[, c(1, 3)]
}


#' FG-weighted Fine-Gray sub-distribution hazard regression
#'
#' Fits a Fine-Gray model on the Fine-Gray split dataset, weighting by the
#' probability of remaining uncensored after the competing event
#' (`p_notcens_after_death`), and uses a robust sandwich variance via
#' `cluster(id)`.
#'
#' @param data_long_fg A data frame in Fine-Gray format with weights, as
#'   returned by [add_fg_weights_cr()].
#' @param covariate Character string. Name of the covariate column. Default is
#'   `"z1"`.
#' @param extend Logical. If `FALSE`, data are truncated at the minimum of the
#'   stratum-specific maximum follow-up times before fitting. Default is `TRUE`.
#' @param event1_level Character string. Factor level in the `delta` column
#'   representing the primary event. Default is `"event_1"`.
#'
#' @return A matrix with one row per term and two columns: the log
#'   sub-distribution hazard ratio and its robust standard error.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_cr(n = 300, censoring = "baseline")
#' dat_long    <- wide_to_long_cr(dat)
#' dat_long_fg <- fg_split_cr(dat_long)
#' dat_long_fg <- add_fg_weights_cr(dat_long_fg, strat = "no")
#' fg_weighted_cr(dat_long_fg)
#'
#' @importFrom survival coxph Surv
#' @export
fg_weighted_cr <- function(data_long_fg, covariate = "z1", extend = TRUE,
                         event1_level = "event_1") {
  if (!extend) {
    cov_levels <- levels(data_long_fg[[covariate]])
    tau <- min(
      vapply(cov_levels, function(lv) {
        idx <- data_long_fg[[covariate]] == lv &
          (is.na(data_long_fg$event2_time) |
             data_long_fg$tstop <= data_long_fg$event2_time)
        max(data_long_fg$tstop[idx])
      }, numeric(1))
    )
    data_long_fg$delta[data_long_fg$tstop > tau] <- "censor"
    data_long_fg$tstop[data_long_fg$tstop > tau] <- tau
    data_long_fg <- data_long_fg[data_long_fg$tstop > data_long_fg$tstart, ]
  }
  fg_formula <- as.formula(sprintf(
    'Surv(tstart, tstop, delta == "%s") ~ %s + cluster(id)',
    event1_level, covariate
  ))
  m1 <- coxph(
    fg_formula,
    weights = p_notcens_after_death,
    data    = data_long_fg[data_long_fg$p_notcens_after_death > 0, ],
    timefix = FALSE
  )
  summary(m1)$coef[, c(1, 4)]
}


#' Compute bootstrap percentile confidence intervals for competing risks estimates
#'
#' Computes the 2.5th and 97.5th percentiles of the bootstrap distribution for
#' each column of a matrix of bootstrap estimates. Intended for competing
#' risks quantities (e.g. cumulative incidence at a set of time points, or
#' Fine-Gray regression terms), where each bootstrap replicate produces more
#' than one estimate, unlike the single scalar log hazard ratio handled by
#' [get_boot_pci_se()].
#'
#' @param boot_mat A numeric matrix of bootstrap estimates, with one row per
#'   bootstrap replicate and one column per estimand (e.g. time point or
#'   regression term), as returned by row-binding the results of multiple
#'   calls to [cuminc_waverage_cr()], [cuminc_ipcw_cr()], or [fg_weighted_cr()].
#'
#' @return A 2-row numeric matrix with rows `"lower"` and `"upper"` giving the
#'   2.5th and 97.5th percentiles of the bootstrap distribution for each
#'   column of `boot_mat`. Columns with any `NA` bootstrap estimate return
#'   `NA` for both bounds.
#'
#' @examples
#' set.seed(1)
#' boot_mat <- matrix(rnorm(500 * 3), nrow = 500, ncol = 3)
#' get_boot_pci_cr(boot_mat)
#'
#' @importFrom stats quantile
#' @export
get_boot_pci_cr <- function(boot_mat) {
  result <- apply(boot_mat, 2, function(x) {
    if (any(is.na(x))) c(NA_real_, NA_real_) else quantile(x, c(0.025, 0.975), names = FALSE)
  })
  rownames(result) <- c("lower", "upper")
  result
}


#' Bootstrap IPCW weighted data for competing risks survival analysis
#'
#' Draws `B` bootstrap samples from the original wide-format competing risks
#' data, converts each sample to long (counting-process) format via
#' [wide_to_long_cr()], and appends IPCW weights via [add_ipcw_weights_cr()].
#' Returns the results as a list of long-format data frames that can be used
#' directly with [cuminc_ipcw_cr()], or further processed with
#' [fg_split_cr()] and [add_fg_weights_cr()] for a Fine-Gray bootstrap.
#'
#' @param data A wide-format competing risks data frame containing the
#'   columns specified by `time_var`, `event_var`, and `covariate`.
#' @param B Integer. Number of bootstrap samples. Default is `500`.
#' @param time_var Character string. Name of the event/censoring time column.
#'   Default is `"t"`.
#' @param event_var Character string. Name of the event indicator column.
#'   Default is `"delta"`.
#' @param covariate Character string. Name of the covariate column. Default
#'   is `"z1"`.
#' @param cens_level Character string. Factor level in `event_var`
#'   representing censoring. Default is `"censor"`.
#' @param event2_level Character string. Factor level in `event_var`
#'   representing the competing event (event type 2). Default is
#'   `"event_2"`.
#' @param strat Character. Passed to [add_ipcw_weights_cr()]. `"no"`
#'   (default) fits a single Cox model for the censoring distribution;
#'   `"yes"` estimates the censoring distribution non-parametrically within
#'   each stratum of `covariate`.
#' @param seed Optional integer seed for reproducibility. Default is `NULL`.
#'
#' @return A list of `B` data frames, each in long (counting-process) format
#'   with an IPCW weight column `p_notcens`, as returned by
#'   [add_ipcw_weights_cr()].
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "baseline")
#' boot_list <- get_ipcw_boot_cr(dat, B = 50)
#' }
#'
#' @export
get_ipcw_boot_cr <- function(data, B = 500, time_var = "t", event_var = "delta",
                              covariate = "z1", cens_level = "censor",
                              event2_level = "event_2", strat = "no",
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  lapply(seq_len(B), function(i) {
    boot_dat <- data[sample(nrow(data), nrow(data), replace = TRUE), ]
    dat_long <- wide_to_long_cr(boot_dat, time_var = time_var,
                                 event_var = event_var, covariate = covariate,
                                 cens_level = cens_level,
                                 event2_level = event2_level)
    add_ipcw_weights_cr(dat_long, covariate = covariate, strat = strat)
  })
}


#' Plot IPCW cumulative incidence with bootstrap percentile confidence intervals
#'
#' Combines bootstrap cumulative incidence curves from [get_ipcw_boot_cr()]
#' with the point estimate from the original IPCW-weighted data to produce a
#' cumulative incidence plot with 95% percentile confidence intervals,
#' mirroring [plot_ipcw_km_boot_ci_se()] for the single-event case.
#'
#' @param boot_data A list of long-format data frames with IPCW weights, as
#'   returned by [get_ipcw_boot_cr()].
#' @param orig_data A single long-format data frame with IPCW weights for the
#'   original (non-bootstrapped) dataset, as returned by
#'   [add_ipcw_weights_cr()]. Used for the point estimate.
#' @param esttimes Numeric vector of times at which to evaluate cumulative
#'   incidence. Defaults to 100 equally spaced points from 0 to 10.
#' @param extend Logical. Passed to [cuminc_ipcw_cr()]. Default is `TRUE`.
#' @param covariate Character string. Name of the covariate column, used only
#'   when `extend = FALSE` to compute the truncation time. Default is `"z1"`.
#'
#' @return A [ggplot2::ggplot()] object. Add layers or themes to customise.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' dat <- sim_data_cr(n = 200, censoring = "baseline")
#' dat_long  <- add_ipcw_weights_cr(wide_to_long_cr(dat), strat = "no")
#' boot_list <- get_ipcw_boot_cr(dat, B = 50)
#' plot_ipcw_cuminc_boot_ci_cr(boot_list, dat_long, esttimes = seq(0, 5, 0.1))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_step geom_ribbon labs theme_minimal
#' @export
plot_ipcw_cuminc_boot_ci_cr <- function(boot_data, orig_data,
                                         esttimes = seq(from = 0, to = 10, length.out = 100),
                                         extend = TRUE, covariate = "z1") {
  boot_mat <- t(vapply(boot_data, function(x) {
    cuminc_ipcw_cr(x, esttimes = esttimes, extend = extend, covariate = covariate)
  }, numeric(length(esttimes))))

  pci <- get_boot_pci_cr(boot_mat)

  orig_est <- cuminc_ipcw_cr(orig_data, esttimes = esttimes, extend = extend,
                              covariate = covariate)

  plot_dat <- data.frame(
    time = esttimes,
    est  = orig_est,
    lpci = pci["lower", ],
    upci = pci["upper", ]
  )

  ggplot(plot_dat, aes(x = time)) +
    geom_ribbon(aes(ymin = lpci, ymax = upci), alpha = 0.2, linetype = "blank") +
    geom_step(aes(y = est)) +
    labs(x = "Time", y = "Cumulative incidence") +
    theme_minimal()
}
