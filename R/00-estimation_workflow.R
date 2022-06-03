##
## 00. Network Model Estimation, HPC setup
##

# Setup ------------------------------------------------------------------------
library("slurmworkflow")
library("EpiModelHPC")
# The size of network to be used is defined in "R/utils-netsize.R"

# hpc_configs <- swf_configs_hyak(hpc = "mox", partition = "csde")
hpc_configs <- swf_configs_rsph(
  partition = "epimodel",
  mail_user = "user@emory.edu"
)
max_cores <- 10

# Workflow creation ------------------------------------------------------------
wf <- create_workflow(
  wf_name = "estimation",
  default_sbatch_opts = hpc_configs$default_sbatch_opts
)

# Update RENV on the HPC -------------------------------------------------------
#
# this template will run:
#   git pull
#   renv::restore()
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_renv_restore(
    git_branch = "main",
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = hpc_configs$renv_sbatch_opts
)

# Estimate the networks --------------------------------------------------------
#
# this template uses the syntax of `base::do.call`
# the arguments in `args` will be made available as variables to the script on
# the HPC
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "R/01-estimation.R",
    args = list(
      ncores = max_cores
   ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "24:00:00",
    "mem" = "0" # special: all mem on node
  )
)

# Generate the diagnostics data ------------------------------------------------
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "R/02-diagnostics.R",
    args = list(
      ncores = max_cores,
      nsims = 50,
      nsteps = 500
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "04:00:00",
    "mem-per-cpu" = "4G",
    "mail-type" = "END" # to get a mail upon completion
  )
)

# to send the workflows on the HPC
# scp -r ~/git/BigNets/workflows/estimation sph:/projects/epimodel/sjenness/BigNets/workflows/estimation
# from the BigNets folder on Mox: workflows/estimation/start_workflow.sh

# to get the data back
# scp -r sph:projects/epimodel/sjenness/BigNets/data/input data/input/
