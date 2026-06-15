dat_cr <- local({
  set.seed(1)
  sim_data_CR(n = 200, censoring = "baseline")
})

test_that("wide_to_long_CR returns long-format data with required columns", {
  dat_long <- wide_to_long_CR(dat_cr)

  expect_s3_class(dat_long, "data.frame")
  expect_true(nrow(dat_long) > nrow(dat_cr))
  expect_true(all(c("id", "tstart", "tstop", "delta", "censor",
                    "event2_time", "z1") %in% names(dat_long)))
})

test_that("wide_to_long_CR tstart < tstop for all rows", {
  dat_long <- wide_to_long_CR(dat_cr)
  expect_true(all(dat_long$tstart < dat_long$tstop))
})

test_that("add_ipcw_weights (Cox) adds p_notcens in valid range", {
  dat_long   <- wide_to_long_CR(dat_cr)
  dat_weights <- add_ipcw_weights(dat_long, strat = "no")

  expect_true("p_notcens" %in% names(dat_weights))
  expect_true(all(!is.na(dat_weights$p_notcens)))
  expect_true(all(dat_weights$p_notcens > 0 & dat_weights$p_notcens <= 1))
})

test_that("add_ipcw_weights (stratified) adds p_notcens in valid range", {
  dat_long    <- wide_to_long_CR(dat_cr)
  dat_weights <- add_ipcw_weights(dat_long, strat = "yes")

  expect_true("p_notcens" %in% names(dat_weights))
  expect_true(all(!is.na(dat_weights$p_notcens)))
  expect_true(all(dat_weights$p_notcens > 0 & dat_weights$p_notcens <= 1))
})

test_that("cuminc_naive returns correct-length numeric in [0, 1]", {
  times  <- seq(0, 5, 0.5)
  result <- cuminc_naive(dat_cr, esttimes = times)

  expect_length(result, length(times))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("cuminc_waverage returns correct-length numeric in [0, 1]", {
  times  <- seq(0, 5, 0.5)
  result <- cuminc_waverage(dat_cr, esttimes = times)

  expect_length(result, length(times))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("cuminc_ipcw returns correct-length numeric in [0, 1]", {
  dat_long <- wide_to_long_CR(dat_cr)
  dat_long <- add_ipcw_weights(dat_long, strat = "no")
  times    <- seq(0, 5, 0.5)
  result   <- cuminc_ipcw(dat_long, esttimes = times)

  expect_length(result, length(times))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("fg_split returns more rows than wide_to_long_CR", {
  dat_long    <- wide_to_long_CR(dat_cr)
  dat_long_fg <- fg_split(dat_long)

  expect_true(nrow(dat_long_fg) >= nrow(dat_long))
  expect_true(all(c("id", "tstart", "tstop", "delta", "z1",
                    "event2_time") %in% names(dat_long_fg)))
})

test_that("add_fg_weights adds p_notcens_after_death in valid range", {
  dat_long    <- wide_to_long_CR(dat_cr)
  dat_long_fg <- fg_split(dat_long)
  dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")

  expect_true("p_notcens_after_death" %in% names(dat_long_fg))
  expect_true(all(!is.na(dat_long_fg$p_notcens_after_death)))
  expect_true(all(dat_long_fg$p_notcens_after_death > 0 &
                    dat_long_fg$p_notcens_after_death <= 1))
})

test_that("fg_naive returns a matrix with correct dimensions", {
  result <- fg_naive(dat_cr)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 3)   # three non-reference levels of z1
})

test_that("fg_weighted returns a matrix with correct dimensions", {
  dat_long    <- wide_to_long_CR(dat_cr)
  dat_long_fg <- fg_split(dat_long)
  dat_long_fg <- add_fg_weights(dat_long_fg, strat = "no")
  result      <- fg_weighted(dat_long_fg)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 3)
})
