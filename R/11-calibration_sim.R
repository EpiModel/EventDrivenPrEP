##
## 11. Epidemic Model Parameter Calibration, Simulations
##

# Required variables:
#   - `scenario`
#   - `batch_num`
#   - `ncores`
if (interactive()) {
  scenario <- NULL
  batch_num <- 1
  ncores <- 2
}

# Setup ------------------------------------------------------------------------
suppressMessages({
  library("EpiModelHIV")
  library("EpiModelHPC")
})

# Load the `NETSIZE` value and the formatted `netsize_string`
# NETSIZE <- 1e4 # to override (before sourcing the file)
source("R/utils-netsize.R")

output_dir <- "data/output/calib"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source("R/utils-netsim_inputs.R")
source("R/utils-targets.R")

control <- control_msm(
  nsteps = calibration_length,
  nsims = ncores,
  ncores = ncores,
  cumulative.edgelist = TRUE,
  truncate.el.cuml = 0,
  verbose = TRUE,
  tracker.list = calibration_trackers # created in R/utils-targets.R
)

# apply the scenario to the parameters
param <- use_scenario(param, scenario)

# Simulation -------------------------------------------------------------------
print(paste0("Starting simulation for scenario: ", scenario[["id"]]))
print(paste0("Batch number: ", batch_num))
sim <- netsim(est, param, init, control)

## Save the simulation
#  I am not very familiar with the `savesim` workflow so I am doing it this
#  way for now

file_name <- paste0("simcalib__", scenario[["id"]], "__", batch_num, ".rds")

print(paste0("Saving simulation in file: ", file_name))

saveRDS(sim, paste0(output_dir, "/", file_name))

print("Done!")

