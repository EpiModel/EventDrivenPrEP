##
## 00. Shared variables setup
##

est_dir <- "./data/intermediate/estimates/"
diag_dir <- "./data/intermediate/diagnostics/"
calib_dir <- "./data/intermediate/calibration/"
scenarios_dir <- "./data/intermediate/scenarios/"

# Information for the HPC workflows
current_git_branch <- "alg_calib"
mail_user <- "aleguil@emory.edu" # or any other mail provider

# Relevant time steps for the simulation
calibration_end    <- 364 * 45
restart_time       <- calibration_end + 1
prep_start         <- restart_time + 364 * 2
intervention_start <- prep_start + 364 * 5
intervention_end   <- intervention_start + 364 * 10
