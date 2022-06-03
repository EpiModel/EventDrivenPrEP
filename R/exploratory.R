# Exploring EpiModelHIV modules to incorporate EDP

at <- 100

dat <- initialize_msm(est, param, init, control)
dat <- aging_msm(dat, at)
dat <- departure_msm(dat, at)
dat <- arrival_msm(dat, at)
dat <- hivtest_msm(dat, at)
dat <- hivtx_msm(dat, at)
dat <- hivprogress_msm(dat, at)
dat <- hivvl_msm(dat, at)
dat <- simnet_msm(dat, at)
dat <- acts_msm(dat, at)
dat <- condoms_msm(dat, at)
dat <- position_msm(dat, at)
dat <- prep_msm(dat, at)
dat <- hivtrans_msm(dat, at)
dat <- stitrans_msm(dat, at)
dat <- stirecov_msm(dat, at)
dat <- prevalence_msm(dat, at)
