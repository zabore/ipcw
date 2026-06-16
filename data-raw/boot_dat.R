## code to prepare `boot_dat` dataset goes here
set.seed(20240429)
dat <- sim_data_se(n = 500)

boot_dat <- 
  get_ipcw_boot(
    data = dat, 
    B = 500, 
    time_var = "t", 
    event_var = "delta",
    cens_cov = "W2", 
    seed = 20240917)

usethis::use_data(boot_dat, overwrite = TRUE)
