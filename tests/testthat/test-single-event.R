dat_se <- local({
  set.seed(1)
  sim_data_se(n = 200)
})

test_that("get_ipcw_wgt returns long-format data with weight column", {
  dat_long <- get_ipcw_wgt(dat_se)

  expect_s3_class(dat_long, "data.frame")
  expect_true(nrow(dat_long) > nrow(dat_se))
  expect_true("wgt" %in% names(dat_long))
  expect_true("tstart" %in% names(dat_long))
  expect_true("tstop"  %in% names(dat_long))
  expect_true("id"     %in% names(dat_long))
})

test_that("get_ipcw_wgt produces valid weights", {
  dat_long <- get_ipcw_wgt(dat_se)

  expect_true(all(!is.na(dat_long$wgt)))
  expect_true(all(dat_long$wgt >= 1))
})

test_that("get_ipcw_wgt tstart < tstop for all rows", {
  dat_long <- get_ipcw_wgt(dat_se)
  expect_true(all(dat_long$tstart < dat_long$tstop))
})

test_that("get_ipcw_cox_fit returns expected columns", {
  dat_long <- get_ipcw_wgt(dat_se)
  result   <- get_ipcw_cox_fit(dat_long, weight = "wgt")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("term", "log_hr", "log_hr_se", "log_hr_rob_se",
                    "hr", "hr_ci_low", "hr_ci_high") %in% names(result)))
  expect_true(all(result$hr > 0))
})

test_that("get_cox_fit returns expected columns", {
  dat_long <- get_ipcw_wgt(dat_se)
  result   <- get_cox_fit(dat_long)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("term", "log_hr", "log_hr_se",
                    "hr", "hr_ci_low", "hr_ci_high") %in% names(result)))
  expect_true(all(result$hr > 0))
})

test_that("get_ipcw_km_prob_x returns tibble with correct columns", {
  dat_long <- get_ipcw_wgt(dat_se)
  result   <- get_ipcw_km_prob_x(dat_long, pre_times = seq(0, 100, 10))

  expect_true(all(c("time", "surv", "x") %in% names(result)))
  expect_true(all(result$surv >= 0 & result$surv <= 1))
})

test_that("get_boot_var returns a single non-negative numeric", {
  boot_data <- data.frame(log_hr = rnorm(50, mean = -0.5, sd = 0.2))
  result <- get_boot_var(boot_data, B = 50)

  expect_length(result, 1)
  expect_true(is.numeric(result))
  expect_gte(result, 0)
})

test_that("get_boot_pci returns a two-element vector in correct order", {
  boot_data <- data.frame(log_hr = rnorm(500, mean = -0.5, sd = 0.2))
  result <- get_boot_pci(boot_data)

  expect_length(result, 2)
  expect_lt(result[1], result[2])
})
