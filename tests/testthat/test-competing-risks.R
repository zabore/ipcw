dat_cr <- local({
  set.seed(1)
  sim_data_cr(n = 200, censoring = "baseline")
})

test_that("wide_to_long_cr returns long-format data with required columns", {
  dat_long <- wide_to_long_cr(dat_cr)

  expect_s3_class(dat_long, "data.frame")
  expect_true(nrow(dat_long) > nrow(dat_cr))
  expect_true(all(c("id", "tstart", "tstop", "delta", "censor",
                    "event2_time", "z1") %in% names(dat_long)))
})

test_that("wide_to_long_cr tstart < tstop for all rows", {
  dat_long <- wide_to_long_cr(dat_cr)
  expect_true(all(dat_long$tstart < dat_long$tstop))
})

test_that("add_ipcw_weights_cr (Cox) adds p_notcens in valid range", {
  dat_long   <- wide_to_long_cr(dat_cr)
  dat_weights <- add_ipcw_weights_cr(dat_long, strat = "no")

  expect_true("p_notcens" %in% names(dat_weights))
  expect_true(all(!is.na(dat_weights$p_notcens)))
  expect_true(all(dat_weights$p_notcens > 0 & dat_weights$p_notcens <= 1))
})

test_that("add_ipcw_weights_cr (stratified) adds p_notcens in valid range", {
  dat_long    <- wide_to_long_cr(dat_cr)
  dat_weights <- add_ipcw_weights_cr(dat_long, strat = "yes")

  expect_true("p_notcens" %in% names(dat_weights))
  expect_true(all(!is.na(dat_weights$p_notcens)))
  expect_true(all(dat_weights$p_notcens > 0 & dat_weights$p_notcens <= 1))
})

test_that("cuminc_naive_cr returns correct-length numeric in [0, 1]", {
  times  <- seq(0, 5, 0.5)
  result <- cuminc_naive_cr(dat_cr, esttimes = times)

  expect_length(result, length(times))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("cuminc_waverage_cr returns correct-length numeric in [0, 1]", {
  times  <- seq(0, 5, 0.5)
  result <- cuminc_waverage_cr(dat_cr, esttimes = times)

  expect_length(result, length(times))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("cuminc_ipcw_cr returns correct-length numeric in [0, 1]", {
  dat_long <- wide_to_long_cr(dat_cr)
  dat_long <- add_ipcw_weights_cr(dat_long, strat = "no")
  times    <- seq(0, 5, 0.5)
  result   <- cuminc_ipcw_cr(dat_long, esttimes = times)

  expect_length(result, length(times))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("fg_split_cr returns more rows than wide_to_long_cr", {
  dat_long    <- wide_to_long_cr(dat_cr)
  dat_long_fg <- fg_split_cr(dat_long)

  expect_true(nrow(dat_long_fg) >= nrow(dat_long))
  expect_true(all(c("id", "tstart", "tstop", "delta", "z1",
                    "event2_time") %in% names(dat_long_fg)))
})

test_that("add_fg_weights_cr adds p_notcens_after_death in valid range", {
  dat_long    <- wide_to_long_cr(dat_cr)
  dat_long_fg <- fg_split_cr(dat_long)
  dat_long_fg <- add_fg_weights_cr(dat_long_fg, strat = "no")

  expect_true("p_notcens_after_death" %in% names(dat_long_fg))
  expect_true(all(!is.na(dat_long_fg$p_notcens_after_death)))
  expect_true(all(dat_long_fg$p_notcens_after_death > 0 &
                    dat_long_fg$p_notcens_after_death <= 1))
})

test_that("fg_naive_cr returns a matrix with correct dimensions", {
  result <- fg_naive_cr(dat_cr)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 3)   # three non-reference levels of z1
})

test_that("fg_weighted_cr returns a matrix with correct dimensions", {
  dat_long    <- wide_to_long_cr(dat_cr)
  dat_long_fg <- fg_split_cr(dat_long)
  dat_long_fg <- add_fg_weights_cr(dat_long_fg, strat = "no")
  result      <- fg_weighted_cr(dat_long_fg)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 3)
})

test_that("get_ipcw_boot_cr returns a list of long-format weighted data frames", {
  boot_list <- get_ipcw_boot_cr(dat_cr, B = 3, strat = "no", seed = 1)

  expect_length(boot_list, 3)
  expect_true(all(vapply(boot_list, is.data.frame, logical(1))))
  expect_true(all(vapply(boot_list, function(x) "p_notcens" %in% names(x), logical(1))))
})

test_that("get_boot_pci_cr returns a 2-row matrix with lower < upper", {
  set.seed(1)
  boot_mat <- matrix(rnorm(500 * 3, mean = 0.5, sd = 0.2), nrow = 500, ncol = 3)
  result   <- get_boot_pci_cr(boot_mat)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 3))
  expect_equal(rownames(result), c("lower", "upper"))
  expect_true(all(result["lower", ] < result["upper", ]))
})

test_that("get_boot_pci_cr returns NA for columns with any NA", {
  boot_mat <- matrix(c(rnorm(500), rep(NA, 500)), nrow = 500, ncol = 2)
  result   <- get_boot_pci_cr(boot_mat)

  expect_true(all(is.na(result[, 2])))
  expect_true(all(!is.na(result[, 1])))
})

test_that("plot_ipcw_cuminc_boot_ci_cr returns a ggplot object", {
  dat_long  <- add_ipcw_weights_cr(wide_to_long_cr(dat_cr), strat = "no")
  boot_list <- get_ipcw_boot_cr(dat_cr, B = 3, strat = "no", seed = 1)
  esttimes  <- seq(0, 5, 1)
  result    <- plot_ipcw_cuminc_boot_ci_cr(boot_list, dat_long, esttimes = esttimes)

  expect_s3_class(result, "ggplot")
})
