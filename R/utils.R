# Suppress R CMD check notes for column names used in NSE contexts
utils::globalVariables(c(
  # shared
  "id", "tstart", "tstop", "delta", "censor",
  # single-event
  "wgt", "inv_wgt",
  # competing risks
  "z1", "p_notcens", "p_notcens_after_death",
  "event2_time", "fgwt",
  # broom column names returned as data frames
  "estimate", "std.error", "robust.se", "conf.low", "conf.high"
))
