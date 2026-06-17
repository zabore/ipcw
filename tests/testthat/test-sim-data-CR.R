test_that("sim_data_cr returns a data frame with correct structure", {
  set.seed(1)
  dat <- sim_data_cr(n = 100)

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100)
  expect_named(dat, c("z1", "delta", "t"))
})

test_that("sim_data_cr delta has correct factor levels", {
  set.seed(1)
  dat <- sim_data_cr(n = 200)

  expect_s3_class(dat$delta, "factor")
  expect_equal(levels(dat$delta), c("censor", "event_1", "event_2"))
})

test_that("sim_data_cr column values are in valid ranges", {
  set.seed(1)
  dat <- sim_data_cr(n = 200)

  expect_true(all(dat$t > 0))
  expect_s3_class(dat$z1, "factor")
  expect_setequal(levels(dat$z1), c("0", "1", "2", "3"))
})

test_that("sim_data_cr censoring argument works", {
  set.seed(1)
  dat_none <- sim_data_cr(n = 500, censoring = "none")
  expect_equal(sum(dat_none$delta == "censor"), 0)

  set.seed(1)
  dat_ind <- sim_data_cr(n = 500, censoring = "independent")
  expect_gt(sum(dat_ind$delta == "censor"), 0)

  set.seed(1)
  dat_bl <- sim_data_cr(n = 500, censoring = "baseline")
  expect_gt(sum(dat_bl$delta == "censor"), 0)
})

test_that("sim_data_cr errors on invalid censoring argument", {
  expect_error(sim_data_cr(censoring = "random"), "censoring must be one of")
})

test_that("sim_data_cr is reproducible with set.seed", {
  set.seed(99)
  dat1 <- sim_data_cr(n = 100, censoring = "baseline")
  set.seed(99)
  dat2 <- sim_data_cr(n = 100, censoring = "baseline")
  expect_equal(dat1, dat2)
})
