##
## 00. Shared variables setup
##

est_dir <- "./data/intermediate/estimates/"
diag_dir <- "./data/intermediate/diagnostics/"
calib_dir <- "./data/intermediate/calibration/"
scenarios_dir <- "./data/intermediate/scenarios/"

# Information for the HPC workflows
# current_git_branch <- "main"
# mail_user <- "clchand@emory.edu" # or any other mail provider

# Relevant time steps for the simulation
calibration_end    <- 364 * 60
restart_time       <- calibration_end + 1
prep_start         <- restart_time + 364 * 5
edp_start          <- prep_start + 364 * 10
intervention_start <- prep_start + 364 * 10
intervention_end   <- intervention_start + 364 * 10
