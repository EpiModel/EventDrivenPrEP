# Settings ---------------------------------------------------------------------
context <- "hpc"
source("R/utils-0_project_settings.R")

# Run the simulations ----------------------------------------------------------
library("EpiModelHIV")

# Necessary files
source("R/utils-default_inputs.R") # generate `path_to_restart`, `param`, `init`

# Controls
source("R/utils-targets.R")
control <- control_msm(
  start               = restart_time,
  nsteps              = restart_time + 10,
  nsims               = 1,
  ncores              = 1,
  initialize.FUN      = reinit_msm,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
  verbose             = TRUE
)

est <- readRDS(path_to_restart)
netsim(est, param, init, control)

