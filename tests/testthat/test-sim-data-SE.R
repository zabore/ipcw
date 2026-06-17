test_that("sim_data_se returns a data frame with correct structure", {
  set.seed(1)
  dat <- sim_data_se(n = 100)

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100)
  expect_named(dat, c("S", "t", "delta", "x", "W2"))
})

test_that("sim_data_se column values are in valid ranges", {
  set.seed(1)
  dat <- sim_data_se(n = 200)

  expect_true(all(dat$S > 0))
  expect_true(all(dat$t > 0))
  expect_true(all(dat$t <= dat$S))
  expect_setequal(unique(dat$delta), c(0L, 1L))
  expect_setequal(unique(dat$x),     c(0L, 1L))
  expect_true(all(dat$W2 > 0) && all(dat$W2 < 1))
})

test_that("sim_data_se respects the n argument", {
  set.seed(1)
  expect_equal(nrow(sim_data_se(n = 50)),  50)
  expect_equal(nrow(sim_data_se(n = 500)), 500)
})

test_that("sim_data_se is reproducible with set.seed", {
  set.seed(42)
  dat1 <- sim_data_se(n = 100)
  set.seed(42)
  dat2 <- sim_data_se(n = 100)
  expect_equal(dat1, dat2)
})

test_that("sim_data_se x_prop parameter shifts treatment prevalence", {
  set.seed(1)
  dat_low  <- sim_data_se(n = 2000, x_prop = 0.1)
  dat_high <- sim_data_se(n = 2000, x_prop = 0.9)
  expect_lt(mean(dat_low$x),  0.3)
  expect_gt(mean(dat_high$x), 0.7)
})
