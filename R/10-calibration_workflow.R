##
## 10. Epidemic Model Parameter Calibration, HPC setup
##

# Setup ------------------------------------------------------------------------
library("slurmworkflow")
library("EpiModelHPC")

# hpc_configs <- swf_configs_hyak(hpc = "mox", partition = "csde")
hpc_configs <- swf_configs_rsph(
  partition = "preemptable",
  mail_user = "user@emory.edu"
)
max_cores <- 28

# Workflow creation ------------------------------------------------------------
wf <- create_workflow(
  wf_name = "calibration",
  default_sbatch_opts = hpc_configs$default_sbatch_opts
)

# Update RENV on the HPC -------------------------------------------------------
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_renv_restore(
    git_branch = "main",
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = hpc_configs$renv_sbatch_opts
)

# # Run the simulations ----------------------------------------------------------
source("R/utils-netsize.R")
source("R/utils-netsim_inputs.R")
source("R/utils-targets.R")

control <- control_msm(
  nsteps = calibration_length,
  nsims = 1, ncores = 1,
  cumulative.edgelist = TRUE,
  truncate.el.cuml = 0,
  verbose = FALSE,
  tracker.list = calibration_trackers # created in R/utils-targets.R
)

scenarios.df <- read.csv("data/input/calib_scenarios.csv")
scenarios.list <- EpiModel::create_scenario_list(scenarios.df)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    est, param, init, control,
    scenarios_list = scenarios.list,
    output_dir = "data/output/calib",
    libraries = "EpiModelHIV",
    n_rep = 50,
    n_cores = max_cores,
    max_array_size = 999,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "01:00:00",
    "mem" = "0" # special: all mem on node
  )
)

# Process calibrations ---------------------------------------------------------
# produce a data frame with the calibration targets for each scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "R/12-calibration_process.R",
    args = list(
      ncores = 15,
      nsteps = 52
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
# to send the workflows on the HPC
# scp -r workflows klone.hyak.uw.edu:gscratch/BigNets/
# from the BigNets folder on Klone: workflows/calibration/start_workflow.sh

# to get the data back
# scp klone.hyak.uw.edu:gscratch/BigNets/data/output/calib/assessments.rds data/output/calib/
#
# here I only download the processed data. Not all the simulations
