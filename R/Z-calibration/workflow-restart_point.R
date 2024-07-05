## HPC Workflow: Restart Point
##
## Define a workflow to run a few hundred replications of the default paramater
## and choose the best restart point

# Restart R before running this script (Ctrl_Shift_F10 / Cmd_Shift_0)

# Settings ---------------------------------------------------------------------
library(slurmworkflow)
library(EpiModelHPC)
library(EpiModelHIV)
library(dplyr)

context <- "hpc"
source("./R/hpc_configs.R")
source("R/shared_variables.R", local = TRUE)
source("./R/netsim_settings.R")
source("./R/utils-targets.R")
max_cores <- 8

est <- readRDS(path_to_est)
control <- control_msm(
  .checkpoint.dir = "tmp/ckpt/",
  .checkpoint.step = 364,
  nsteps              = calibration_end,
  .tracker.list       = calibration_trackers,
  verbose             = FALSE
)

wf <- make_em_workflow("restart_point", override = TRUE)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_est, param, init, control,
    scenarios_list = NULL,
    output_dir = calib_dir,
    save_pattern = "all",
    n_rep = 64,
    n_cores = max_cores,
    max_array_size = 500,
    setup_lines = hpc_node_setup
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL",
    "cpus-per-task" = max_cores,
    "time" = "24:00:00",
    "mem-per-cpu" = "5G"
  )
)

# Process calibrations
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_merge_netsim_scenarios_tibble(
      sim_dir = calib_dir,
      output_dir = fs::path(calib_dir, "merged_tibbles"),
      steps_to_keep = Inf,
      cols = dplyr::everything(),
      n_cores = max_cores,
      setup_lines = hpc_node_setup
    ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "02:00:00",
    "mem-per-cpu" = "5G"
  )
)

# wf <- add_workflow_step(
#   wf_summary = wf,
#   step_tmpl = step_tmpl_do_call_script(
#     r_script = "R/Z-calibration/choose_restart.R",
#     args = list(hpc_context = TRUE),
#     setup_lines = hpc_node_setup
#   ),
#   sbatch_opts = list(
#     "cpus-per-task" = max_cores,
#     "time" = "02:00:00",
#     "mem-per-cpu" = "5G"
#   )
# )

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "R/Z-calibration/process_calib_plots.R",
    args = list(hpc_context = TRUE, scenario = "empty_scenario"),
    setup_lines = hpc_node_setup
  ),
  sbatch_opts = list(
    "mail-type" = "END",
    "cpus-per-task" = max_cores,
    "time" = "02:00:00",
    "mem-per-cpu" = "5G"
  )
)
