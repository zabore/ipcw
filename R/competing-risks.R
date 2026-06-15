#' Simulate competing risks survival data
#'
#' Generates a dataset with two competing events under a sub-distribution
#' hazard model parameterized by a four-level baseline covariate `z1`.
#' Optionally adds independent or covariate-dependent administrative censoring.
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
#' dat <- sim_data_CR(n = 200, censoring = "baseline")
#' table(dat$delta)
#'
#' @importFrom survival survfit Surv
#' @importFrom stats rexp runif uniroot
#' @export
sim_data_CR <- function(n = 100, censoring = "none",
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
#' for the censoring distribution and for use with [add_ipcw_weights()].
#'
#' @param dat A data frame with columns `t` (event/censoring time), `delta`
#'   (factor with levels `"censor"`, `"event_1"`, `"event_2"`), and `z1`
#'   (covariate). Typically the output of [sim_data_CR()].
#'
#' @return A data frame in long format with columns `id`, `z1`, `delta`,
#'   `censor`, `event2_time`, `tstart`, and `tstop`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 100, censoring = "baseline")
#' dat_long <- wide_to_long_CR(dat)
#' head(dat_long)
#'
#' @importFrom survival survSplit Surv
#' @importFrom dplyr rename
#' @export
wide_to_long_CR <- function(dat) {
  dat$id          <- seq_len(nrow(dat))
  dat$censor      <- (dat$delta == "censor")
  dat$event2_time <- ifelse(dat$delta == "event_2", dat$t, NA)
  dat$tstart      <- 0

  times <- sort(unique(dat$t[dat$censor]))

  data_long1 <- survSplit(dat, cut = times, end = "t", start = "tstart",
                          event = "delta")
  data_long1 <- data_long1[order(data_long1$id, data_long1$t), ]

  data_long2 <- survSplit(dat, cut = times, end = "t", start = "tstart",
                          event = "censor")
  data_long2 <- data_long2[order(data_long2$id, data_long2$t), ]

  data_long1$censor <- data_long2$censor
  rename(data_long1, tstop = t)
}


#' Add IPCW weights to competing risks long-format data
#'
#' Estimates the probability of remaining uncensored and appends it as column
#' `p_notcens` to the dataset. Supports two estimation strategies: a Cox
#' proportional hazards model (`strat = "no"`) or non-parametric (KM) estimates
#' within each level of `z1` (`strat = "yes"`).
#'
#' @param data_long A data frame in long format, as returned by
#'   [wide_to_long_CR()]. Must contain columns `z1`, `tstart`, `tstop`, and
#'   `censor`.
#' @param strat Character. `"no"` (default) fits a single Cox model for the
#'   censoring distribution using `z1` as a covariate. `"yes"` estimates
#'   the censoring distribution non-parametrically within each stratum of `z1`.
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
#' dat <- sim_data_CR(n = 100, censoring = "baseline")
#' dat_long <- wide_to_long_CR(dat)
#' dat_long <- add_ipcw_weights(dat_long, strat = "no")
#' summary(dat_long$p_notcens)
#'
#' @importFrom survival coxph survfit Surv
#' @importFrom stats predict
#' @export
add_ipcw_weights <- function(data_long, strat = "no", new_data = NULL,
                             by.start = TRUE) {
  if (is.null(new_data)) new_data <- data_long

  temp <- new_data[, c("z1", "tstart", "tstop", "censor")]
  if (by.start) {
    temp$tstop  <- temp$tstart
    temp$tstart <- 0
  }
  temp$p_notcens <- NA

  if (strat == "yes") {
    cens_mod <- lapply(levels(data_long$z1), function(x) {
      survfit(Surv(tstart, tstop, censor) ~ 1,
              data = data_long[data_long$z1 == x, ], timefix = FALSE)
    })
    for (i in seq_along(levels(data_long$z1))) {
      lv  <- levels(data_long$z1)[i]
      idx <- temp$z1 == lv
      temp$p_notcens[idx] <-
        summary(cens_mod[[i]], temp$tstop[idx],  extend = TRUE)$surv /
        summary(cens_mod[[i]], temp$tstart[idx], extend = TRUE)$surv
    }
  }

  if (strat == "no") {
    cens_mod <- coxph(Surv(tstart, tstop, censor) ~ z1,
                      data = data_long, timefix = FALSE)
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
#' @param dat A wide-format competing risks data frame with columns `t` and
#'   `delta` (factor with levels `"censor"`, `"event_1"`, `"event_2"`).
#' @param esttimes Numeric vector of times at which to return estimates.
#'
#' @return A numeric vector of cumulative incidence estimates at `esttimes`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 200, censoring = "independent")
#' cuminc_naive(dat, esttimes = seq(0, 5, 0.5))
#'
#' @importFrom survival survfit Surv
#' @export
cuminc_naive <- function(dat, esttimes) {
  fit <- survfit(Surv(t, delta) ~ 1, data = dat)
  est <- summary(fit, times = esttimes, extend = TRUE)
  est$pstate[, 2]
}


#' Weighted-average (non-parametric) IPCW cumulative incidence estimate
#'
#' Estimates the marginal cumulative incidence of event type 1 by computing
#' stratum-specific Aalen-Johansen estimates (stratified by `z1`) and then
#' combining them with sample-proportion weights. This is the non-parametric
#' IPCW approach.
#'
#' @param dat A wide-format competing risks data frame with columns `t`,
#'   `delta`, and `z1`.
#' @param esttimes Numeric vector of times at which to return estimates.
#'   Defaults to 100 equally spaced points from 0 to 10.
#'
#' @return A numeric vector of weighted-average cumulative incidence estimates
#'   at `esttimes`. Values beyond the minimum of the stratum-specific maximum
#'   follow-up times are set to `NA`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 200, censoring = "baseline")
#' cuminc_waverage(dat, esttimes = seq(0, 5, 0.5))
#'
#' @importFrom survival strata survfit Surv
#' @export
cuminc_waverage <- function(dat,
                            esttimes = seq(from = 0, to = 10, length.out = 100)) {
  strat  <- survfit(Surv(t, delta) ~ strata(z1), data = dat)
  est    <- summary(strat, times = esttimes, extend = TRUE)
  pstate <- matrix(est$pstate[, 2], nrow = length(esttimes), ncol = 4)
  n      <- nrow(dat)
  weighted <- (pstate[, 1] * sum(dat$z1 == 0) +
               pstate[, 2] * sum(dat$z1 == 1) +
               pstate[, 3] * sum(dat$z1 == 2) +
               pstate[, 4] * sum(dat$z1 == 3)) / n
  tau <- min(max(dat$t[dat$z1 == 0]),
             max(dat$t[dat$z1 == 1]),
             max(dat$t[dat$z1 == 2]),
             max(dat$t[dat$z1 == 3]))
  weighted[esttimes >= tau] <- NA
  weighted
}


#' Cox model IPCW cumulative incidence estimate
#'
#' Estimates the marginal cumulative incidence of event type 1 using a weighted
#' Aalen-Johansen estimator, where the weights are the IPCW weights stored in
#' the `p_notcens` column (as added by [add_ipcw_weights()]).
#'
#' @param data_long A long-format data frame with IPCW weights, as returned by
#'   [add_ipcw_weights()]. Must contain columns `tstart`, `tstop`, `delta`,
#'   `id`, and `p_notcens`.
#' @param esttimes Numeric vector of times at which to return estimates.
#'   Defaults to 100 equally spaced points from 0 to 10.
#' @param extend Logical. If `FALSE`, estimates beyond the minimum of the
#'   stratum-specific maximum follow-up times are set to `NA`. Default is
#'   `TRUE`.
#'
#' @return A numeric vector of IPCW cumulative incidence estimates at
#'   `esttimes`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 200, censoring = "baseline")
#' dat_long <- wide_to_long_CR(dat)
#' dat_long <- add_ipcw_weights(dat_long, strat = "no")
#' cuminc_ipcw(dat_long, esttimes = seq(0, 5, 0.5))
#'
#' @importFrom survival survfit Surv
#' @export
cuminc_ipcw <- function(data_long,
                        esttimes = seq(from = 0, to = 10, length.out = 100),
                        extend = TRUE) {
  fit  <- survfit(Surv(tstart, tstop, delta) ~ 1,
                  data = data_long, weights = 1 / p_notcens,
                  id = id, timefix = FALSE)
  ipcw <- summary(fit, times = esttimes, extend = TRUE)

  if (!extend) {
    tau <- min(max(data_long$tstop[data_long$z1 == 0]),
               max(data_long$tstop[data_long$z1 == 1]),
               max(data_long$tstop[data_long$z1 == 2]),
               max(data_long$tstop[data_long$z1 == 3]))
    ipcw$pstate[, 2][esttimes >= tau] <- NA
  }
  ipcw$pstate[, 2]
}


#' Prepare long-format data for Fine-Gray weighted regression
#'
#' After a type-2 event (event_2), subjects are artificially re-entered into
#' the risk set (as in the Fine-Gray sub-distribution hazard model). This
#' function appends those additional rows, split at every censoring time.
#'
#' @param data_long A long-format data frame, as returned by [wide_to_long_CR()].
#'
#' @return A data frame in Fine-Gray format, with additional rows for subjects
#'   who experienced event_2, sorted by `id` and `tstart`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 100, censoring = "baseline")
#' dat_long <- wide_to_long_CR(dat)
#' dat_long_fg <- fg_split(dat_long)
#' nrow(dat_long_fg) > nrow(dat_long)
#'
#' @importFrom dplyr select
#' @export
fg_split <- function(data_long) {
  times <- sort(unique(data_long$tstop[data_long$censor == 1]))

  event2_dat <- data_long[data_long$delta == "event_2", ] |>
    select(id, z1, event2_time)
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
  fulldat <- rbind(data_long, fulldat[names(data_long)])
  fulldat[order(fulldat$id, fulldat$tstart), ]
}


#' Add Fine-Gray IPCW weights to Fine-Gray split data
#'
#' Computes the probability of remaining uncensored after a type-2 event
#' and appends it as column `p_notcens_after_death`.
#'
#' @param data_long_fg A data frame in Fine-Gray format, as returned by
#'   [fg_split()].
#' @param strat Character. Passed to [add_ipcw_weights()]. `"no"` (default)
#'   uses a Cox model; `"yes"` uses stratum-specific KM estimates.
#'
#' @return `data_long_fg` with an additional column `p_notcens_after_death`.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 100, censoring = "baseline")
#' dat_long    <- wide_to_long_CR(dat)
#' dat_long_fg <- fg_split(dat_long)
#' dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")
#' summary(dat_long_fg$p_notcens_after_death)
#'
#' @export
add_fg_weights <- function(data_long_fg, strat = "no") {
  temp <- data.frame(
    tstart = data_long_fg$event2_time,
    tstop  = data_long_fg$tstart,
    censor = NA,
    z1     = data_long_fg$z1
  )
  temp$p_notcens_after_death <- NA
  temp$p_notcens_after_death[is.na(temp$tstart)] <- 1
  temp$p_notcens_after_death[!is.na(temp$tstart) & temp$tstart >= temp$tstop] <- 1

  still_na <- is.na(temp$p_notcens_after_death)
  if (any(still_na)) {
    ref_data <- data_long_fg[data_long_fg$tstart < data_long_fg$event2_time |
                               is.na(data_long_fg$event2_time), ]
    temp$p_notcens_after_death[still_na] <-
      add_ipcw_weights(ref_data, strat = strat,
                       new_data = temp[still_na, ],
                       by.start = FALSE)$p_notcens
  }

  data_long_fg$p_notcens_after_death <- temp$p_notcens_after_death
  data_long_fg
}


#' Naive Fine-Gray sub-distribution hazard regression
#'
#' Fits a Fine-Gray model using the standard [survival::finegray()] approach,
#' without any IPCW adjustment. Serves as a comparison to [fg_weighted()].
#'
#' @param dat A wide-format competing risks data frame with columns `t`,
#'   `delta`, and `z1`.
#'
#' @return A matrix with one row per term and two columns: the log
#'   sub-distribution hazard ratio and its standard error.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 200, censoring = "independent")
#' fg_naive(dat)
#'
#' @importFrom survival finegray coxph Surv
#' @export
fg_naive <- function(dat) {
  pdata <- finegray(Surv(t, delta) ~ ., data = dat, timefix = FALSE)
  m1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ z1,
              weights = fgwt, data = pdata, timefix = FALSE)
  summary(m1)$coef[, c(1, 3)]
}


#' IPCW-weighted Fine-Gray sub-distribution hazard regression
#'
#' Fits a Fine-Gray model on the Fine-Gray split dataset, weighting by the
#' inverse of the probability of remaining uncensored after event_2
#' (`p_notcens_after_death`), and uses a robust sandwich variance via
#' `cluster(id)`.
#'
#' @param data_long_fg A data frame in Fine-Gray format with weights, as
#'   returned by [add_fg_weights()].
#' @param extend Logical. If `FALSE`, data are truncated at the minimum of the
#'   stratum-specific maximum follow-up times before fitting. Default is `TRUE`.
#'
#' @return A matrix with one row per term and two columns: the log
#'   sub-distribution hazard ratio and its robust standard error.
#'
#' @examples
#' set.seed(42)
#' dat <- sim_data_CR(n = 200, censoring = "baseline")
#' dat_long    <- wide_to_long_CR(dat)
#' dat_long_fg <- fg_split(dat_long)
#' dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")
#' fg_weighted(dat_long_fg)
#'
#' @importFrom survival coxph Surv
#' @export
fg_weighted <- function(data_long_fg, extend = TRUE) {
  if (!extend) {
    tau <- min(
      vapply(0:3, function(lv) {
        idx <- data_long_fg$z1 == lv &
          (is.na(data_long_fg$event2_time) |
             data_long_fg$tstop <= data_long_fg$event2_time)
        max(data_long_fg$tstop[idx])
      }, numeric(1))
    )
    data_long_fg$delta[data_long_fg$tstop > tau]  <- "censor"
    data_long_fg$tstop[data_long_fg$tstop > tau]  <- tau
    data_long_fg <- data_long_fg[data_long_fg$tstop > data_long_fg$tstart, ]
  }
  m1 <- coxph(
    Surv(tstart, tstop, delta == "event_1") ~ z1 + cluster(id),
    weights = p_notcens_after_death,
    data    = data_long_fg[data_long_fg$p_notcens_after_death > 0, ],
    timefix = FALSE
  )
  summary(m1)$coef[, c(1, 4)]
}
