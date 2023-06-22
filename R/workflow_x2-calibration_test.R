##
## Epidemic Model Parameter Calibration, HPC setup
##

# Libraries --------------------------------------------------------------------
library("slurmworkflow")
library("EpiModelHPC")
library("EpiModelHIV")

# Settings ---------------------------------------------------------------------
source("./R/utils-0_project_settings.R")
context <- "hpc"
max_cores <- 8

source("./R/utils-default_inputs.R") # make `path_to_est`, `param` and `init`
source("./R/utils-hpc_configs.R") # creates `hpc_configs`

# ------------------------------------------------------------------------------

# Workflow creation
wf <- create_workflow(
  wf_name = "model_calibration",
  default_sbatch_opts = hpc_configs$default_sbatch_opts
)

# Update RENV on the HPC
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_renv_restore(
    git_branch = current_git_branch,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = hpc_configs$renv_sbatch_opts
)

# Controls
source("./R/utils-targets.R")
control <- control_msm(
  nsteps              = calibration_end,
  nsims               = 1,
  ncores              = 1,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
 .checkpoint.dir     = "./temp/cp_calib",
 .checkpoint.clear   = TRUE,
 .checkpoint.steps   = 5 * 52 * 7,
  verbose             = FALSE
)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_est, param, init, control,
    scenarios_list = NULL, # scenarios_list,
    output_dir = "./data/intermediate/calibration",
    libraries = "EpiModelHIV",
    n_rep = 64,
    n_cores = max_cores,
    save_pattern = "simple",
    max_array_size = 999,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL,TIME_LIMIT",
    "cpus-per-task" = max_cores,
    "time" = "48:00:00",
    "mem-per-cpu" = "10G"
  )
)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_merge_netsim_scenarios_tibble(
      sim_dir = "data/intermediate/calibration",
      output_dir = "data/intermediate/calibration/merged_tibbles",
      steps_to_keep = 52 * 7 * 60,
      cols = dplyr::everything(),
      n_cores = 8,
      setup_lines = hpc_configs$r_loader
    ),
    sbatch_opts = list(
      "mail-type" = "END",
      "cpus-per-task" = 8,
      "time" = "04:00:00",
      "mem-per-cpu" = "5G"
    )
)


# Send the workflow folder to the <HPC> and run it
#
# $ scp -r ./workflows/model_calibration <HPC>:<project_dir>/workflows/
#
# on the HPC:
# chmod +x workflows/model_calibration/start_workflow.sh
# $ ./workflows/model_calibration/start_workflow.sh

# Once the worfklow is finished download the data from the HPC
#
# $ scp -r <HPC>:<project_dir>/data/intermediate/calibration/assessments.rds ./data/intermediate/calibration/
#
# and analyse them locally using: "./R/12-calibration_eval.R"

# rm -rf workflows/model_calibration temp/cp_calib data/intermediate/calibration
