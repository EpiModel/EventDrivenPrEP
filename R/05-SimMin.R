
##
## 05. Medium-level epidemic simulation for testing and parameter exploration
##

# Required variables:
#   - NETSIZE
#   - nsims
#   - ncores

## Packages
suppressMessages({
  library("EpiModelHIV")
  library("EpiModelHPC")
})

## Parameters
# Load the shared variables required for this run of the project
source("R/utils-netsize.R")
## To override the valu of NETSIZE in "shared_vars.R"
# NETSIZE <- 102000
# netsize_string <- format(NETSIZE, scientific = FALSE)

epistats <- readRDS("data/input/epistats.rds")
netstats <- readRDS(paste0("data/input/netstats-", netsize_string, ".rds"))
est <- readRDS(paste0("data/input/netest-", netsize_string, ".rds"))

param <- param_msm(
  netstats = netstats,
  epistats = epistats,
  a.rate = 0.00049,
  hiv.test.rate = c(0.00385, 0.00380, 0.00690),
  tx.init.rate = c(0.1775, 0.190, 0.2521),
  tx.halt.partial.rate = c(0.0062, 0.0055, 0.0031),
  tx.reinit.partial.rate = c(0.00255, 0.00255, 0.00255),
  hiv.trans.scale = c(2.44, 0.424, 0.270),
  riskh.start = 52 * 59,
  prep.start = (52 * 60) + 1,
  prep.start.prob = rep(0.66, 3)
)

init <- init_msm()

control <- control_msm(
  nsteps = 52 * 60,
  nsims = ncores,
  ncores = ncores,
  verbose = TRUE,
  verbose.int = 250,
  verbose.FUN = verbose.hpc.net
)

## Simulation
sim <- netsim(est, param, init, control)

## Save-Min
savesim(sim, save.min = TRUE, save.max = FALSE, compress = TRUE, data.dir = "data/output/")
