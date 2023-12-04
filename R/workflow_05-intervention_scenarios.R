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
max_cores <- 30

source("./R/utils-default_inputs.R") # make `path_to_est`, `param` and `init`
source("./R/utils-hpc_configs.R") # creates `hpc_configs`

# ------------------------------------------------------------------------------

# Workflow creation
wf <- create_workflow(
  wf_name = "intervention_scenarios",
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
  start               = restart_time,
  nsteps              = intervention_end,
  nsims               = 1,
  ncores              = 1,
  initialize.FUN      = reinit_msm,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
  verbose             = FALSE
)

scenarios_df <- readr::read_csv("./data/input/scenarios.csv")
#scenarios_df <- tibble(
#  .scenario.id    = c("scenario_1", "scenario_2", "scenario_3", "scenario_4"),
#  .at             = 1,
#  edp.start.scenario = c(1, 2, 3, 4)
#)
scenarios_list <- EpiModel::create_scenario_list(scenarios_df)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_restart, param, init, control,
    scenarios_list = scenarios_list,
    output_dir = "./data/intermediate/scenarios",
    libraries = "EpiModelHIV",
    save_pattern = "simple",
    n_rep = 120,
    n_cores = max_cores,
    max_array_size = 500,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL,TIME_LIMIT",
    "cpus-per-task" = max_cores,
    "time" = "24:00:00",
    "mem" = 0
  )
)

# Process calibrations
#
# produce a data frame with the calibration targets for each scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "./R/41-intervention_scenarios_process.R",
    args = list(
      context = "hpc",
      ncores = 15
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "04:00:00",
    "mem-per-cpu" = "4G",
    "mail-type" = "END"
  )
)

# Contour plot scenarios: Adherence vs. Coverage--------------------------------
# Workflow creation
wf <- create_workflow(
  wf_name = "contour_plots1",
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
  start               = restart_time,
  nsteps              = intervention_end,
  nsims               = 1,
  ncores              = 1,
  initialize.FUN      = reinit_msm,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
  verbose             = FALSE
)

daily.switch <- seq(0.1, 1, .1)
edp.switch <- seq(0.1, 1, 0.1)

s <- expand.grid(daily = daily.switch, edp = edp.switch)

sc_contour_plots <- tibble(
  .scenario.id = paste0("daily_switch_", s$daily, "_edp_switch_", s$edp),
  .at = 1,
  # contour plot specific
  prep.switch.prob = s$daily,
  edp.swith.prob = s$edp,
  # shared
  edp.start.scenario = 1
)

sc_contour_plots_list <- EpiModel::create_scenario_list(sc_contour_plots)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_restart, param, init, control,
    scenarios_list = sc_contour_plots,
    output_dir = "./data/intermediate/scenarios/contourplots2",
    libraries = "EpiModelHIV",
    save_pattern = "simple",
    n_rep = 120,
    n_cores = max_cores,
    max_array_size = 500,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL,TIME_LIMIT",
    "cpus-per-task" = max_cores,
    "time" = "24:00:00",
    "mem" = 0
  )
)

# Process calibrations
#
# produce a data frame with the calibration targets for each scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "./R/41-intervention_scenarios_process.R",
    args = list(
      context = "hpc",
      ncores = 15
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "04:00:00",
    "mem-per-cpu" = "4G",
    "mail-type" = "END"
  )
)

# Contour plot scenarios: PrEP switching --------------------------------------
# Workflow creation
wf <- create_workflow(
  wf_name = "contour_plots2",
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
  start               = restart_time,
  nsteps              = intervention_end,
  nsims               = 1,
  ncores              = 1,
  initialize.FUN      = reinit_msm,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
  verbose             = FALSE
)

daily.switch <- seq(0.1, 1, .1)
edp.switch <- seq(0.1, 1, 0.1)

s <- expand.grid(daily = daily.switch, edp = edp.switch)

sc_contour_plots <- tibble(
  .scenario.id = paste0("daily_switch_", s$daily, "_edp_switch_", s$edp),
  .at = 1,
  # contour plot specific
  prep.switch.prob = s$daily,
  edp.switch.prob = s$edp,
  # shared
  edp.start.scenario = 1
)

sc_contour_plots_list <- EpiModel::create_scenario_list(sc_contour_plots)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_restart, param, init, control,
    scenarios_list = sc_contour_plots_list,
    output_dir = "./data/intermediate/scenarios/contourplots2",
    libraries = "EpiModelHIV",
    save_pattern = "simple",
    n_rep = 120,
    n_cores = max_cores,
    max_array_size = 500,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL,TIME_LIMIT",
    "cpus-per-task" = max_cores,
    "time" = "24:00:00",
    "mem" = 0
  )
)

# Process calibrations
#
# produce a data frame with the calibration targets for each scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "./R/41-intervention_scenarios_process.R",
    args = list(
      context = "hpc",
      ncores = 15
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "04:00:00",
    "mem-per-cpu" = "4G",
    "mail-type" = "END"
  )
)


# rm -rf workflows/intervention_scenarios

# Send the workflow folder to the <HPC> and run it
#
# $ scp -r ./workflows/intervention_scenarios <HPC>:<project_dir>/workflows/
#
# on the HPC:
# chmod +x workflows/intervention_scenarios/start_workflow.sh
# $ ./workflows/intervention_scenarios/start_workflow.sh

# download the data
# scp -r <HPC>:<project_dir>/data/intermediate/scenarios ./data/intermediate/



