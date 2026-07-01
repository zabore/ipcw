# Suppress R CMD check notes for column names used in NSE contexts
#' @importFrom stats as.formula
#' @importFrom rlang .data
NULL

utils::globalVariables(c(
  # shared structural columns (used as bare names in dplyr/survival calls)
  "id", "tstart", "tstop", "delta", "censor",
  # weight columns used as bare names in survfit/coxph weights= argument
  "inv_wgt", "p_notcens", "p_notcens_after_death",
  # competing risks structural columns
  "event2_time", "fgwt",
  # broom column names returned as data frames
  "estimate", "std.error", "robust.se", "conf.low", "conf.high",
  # plot_ipcw_km_boot_ci_se / plot_ipcw_cuminc_boot_ci_cr columns used as
  # bare names in dplyr/ggplot2
  "surv", "time", "lpci", "upci", "est"
))
