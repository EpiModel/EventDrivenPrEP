
##
## 04. Small-scale epidemic simulation for testing/debugging
##

## Packages
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

# Load the `NETSIZE` value and the formatted `netsize_string`
# NETSIZE <- 1e4     # to override (before sourcing the file)
source("R/utils-netsize.R")

## Parameters
epistats <- readRDS("data/input/epistats.rds")
netstats <- readRDS(paste0("data/input/netstats-", netsize_string, ".rds"))
est <- readRDS(paste0("data/input/netest-", netsize_string, ".rds"))

param <- param_msm(
  netstats               = netstats,
  epistats               = epistats,
  a.rate                 = 0.00049,
  hiv.test.rate          = c(0.00385, 0.00380, 0.00690),
  tx.init.rate           = c(0.1775, 0.190, 0.2521),
  tx.halt.partial.rate   = c(0.0062, 0.0055, 0.0031),
  tx.reinit.partial.rate = c(0.00255, 0.00255, 0.00255),
  hiv.trans.scale        = c(2.44, 0.424, 0.270),
  riskh.start            = 52,
  prep.start             = 2,
  prep.start.prob        = rep(0.66, 3)
)
init <- init_msm()
control <- control_msm(
  simno = 1,
  nsteps = 100,
  nsims = 1,
  ncores = 1,
  verbose = TRUE
)

sim <- netsim(est, param, init, control)
