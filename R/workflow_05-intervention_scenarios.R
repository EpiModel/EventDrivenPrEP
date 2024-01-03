##
## Epidemic Model Intervention Scenarios, HPC setup
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

# Sensitivity analysis: EDP adherence --------------------------------
# Workflow creation
wf <- create_workflow(
  wf_name = "adhr_sens",
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

sc_no_list <- list()
sc_no_list[["no_edp_sc"]] <- tibble(
  .scenario.id = "0_no_edp",
  .at = 1,
  edp.start = intervention_end + 1
)

hi_adhr <- c(seq(-0.74, 0.26, 0.1), 0)

df <- data.frame()

for (i in hi_adhr) {
  prep.adhr.edp.dist <- reallocate_pcp(reall = i)

  df <- rbind(df, prep.adhr.edp.dist)
}

names(df)[1] <- "prep.adhr.edp.dist_1"
names(df)[2] <- "prep.adhr.edp.dist_2"
names(df)[3] <- "prep.adhr.edp.dist_3"
names(df)[4] <- "prep.adhr.edp.dist_4"

df <- df[rep(seq_len(nrow(df)), each = 4), ]
df$edp.start.scenario <- rep(1:4, 12)

sc_interv_list <- list()
sc_interv_list[["adhr_sens"]] <- tibble(
  .scenario.id = paste0("hi_adhr_", df$prep.adhr.edp.dist_4, "_scenario_", df$edp.start.scenario),
  .at = 1,
  prep.start = intervention_start,
  prep.edp.start = intervention_start,
  prep.adhr.edp.dist_1 = df$prep.adhr.edp.dist_1,
  prep.adhr.edp.dist_2 = df$prep.adhr.edp.dist_2,
  prep.adhr.edp.dist_3 = df$prep.adhr.edp.dist_3,
  prep.adhr.edp.dist_4 = df$prep.adhr.edp.dist_4,
  edp.start.scenario = df$edp.start.scenario
)

sc_df_list <- c(
  sc_no_list,
  sc_interv_list
)

scenarios_list <- purrr::reduce(
  sc_df_list,
  \(out, d_sc) c(out, EpiModel::create_scenario_list(d_sc)),
  .init =
)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_restart, param, init, control,
    scenarios_list = scenarios_list,
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

# output one data file per scenario
#wf <- add_workflow_step(
#  wf_summary = wf,
#  step_tmpl = step_tmpl_do_call(
#    what = EpiModelHPC::merge_netsim_scenarios_tibble,
#    args = list(
#      sim_dir = "data/intermediate/scenarios/contourplots1",
#      output_dir = "adhr_sens_output",
#      steps_to_keep = 364 * 10,
#      cols = rlang::quo(dplyr::matches("^doxy")) # rlang::quo required for lazy eval
#    ),
#    setup_lines = hpc_configs$r_loader
#  ),
#  sbatch_opts = list(
#    "mail-type" = "END",
#    "cpus-per-task" = 1,
#    "time" = "02:00:00",
#    "mem" = "15G"
#  )
#)



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


hi_adhr <- c(seq(-0.74, 0.26, 0.1), 0)

df <- data.frame()

for (i in hi_adhr) {
  prep.adhr.edp.dist <- reallocate_pcp(reall = i)

  df <- rbind(df, prep.adhr.edp.dist)
}

names(df)[1] <- "prep.adhr.edp.dist_1"
names(df)[2] <- "prep.adhr.edp.dist_2"
names(df)[3] <- "prep.adhr.edp.dist_3"
names(df)[4] <- "prep.adhr.edp.dist_4"

edp.prep.start.prob <- c(.000827328, 0.000636321, 0.000995523)

prop_change <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2.0)

df.cov <- data.frame(rep(NA, length(prop_change)))

for (j in edp.prep.start.prob) {
  edp.prep.start <- j*prop_change

  df.cov <- cbind(df.cov, edp.prep.start)
}

df.cov <- df.cov[,-1]
names(df.cov)[1] <- "edp.prep.start.prob_1"
names(df.cov)[2] <- "edp.prep.start.prob_2"
names(df.cov)[3] <- "edp.prep.start.prob_3"
df.cov$prop_change <- prop_change

df.merge <- merge(df, df.cov)

contour_plots1 <- tibble(
  .scenario.id = paste0("edp_startprob_", df.merge$prop_change, "_hi_adhr_", df.merge$prep.adhr.edp.dist_4),
  .at = 1,
  prep.start = intervention_start,
  prep.edp.start = intervention_start,
  prep.adhr.edp.dist_1 = df.merge$prep.adhr.edp.dist_1,
  prep.adhr.edp.dist_2 = df.merge$prep.adhr.edp.dist_2,
  prep.adhr.edp.dist_3 = df.merge$prep.adhr.edp.dist_3,
  prep.adhr.edp.dist_4 = df.merge$prep.adhr.edp.dist_4,
  edp.prep.start.prob_1 = df.merge$edp.prep.start.prob_1,
  edp.prep.start.prob_2 = df.merge$edp.prep.start.prob_2,
  edp.prep.start.prob_3 = df.merge$edp.prep.start.prob_3,
  edp.start.scenario = 1
)

sc_contour_plots_list <- EpiModel::create_scenario_list(contour_plots1)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_restart, param, init, control,
    scenarios_list = sc_contour_plots_list,
    output_dir = "./data/intermediate/scenarios/contourplots1",
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

# output one datafile per scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call(
    what = EpiModelHPC::merge_netsim_scenarios_tibble,
    args = list(
      sim_dir = "data/intermediate/scenarios/contourplots1",
      output_dir = "contourplots1_output",
      steps_to_keep = 364 * 10,
      cols = rlang::quo(dplyr::matches("^doxy")) # rlang::quo required for lazy eval
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "END",
    "cpus-per-task" = 1,
    "time" = "02:00:00",
    "mem" = "15G"
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
  edp.start.scenario = 1,
  prep.start = intervention_start,
  prep.edp.start = intervention_start
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



### Workflow for just renv and merge steps

# Workflow creation
wf <- create_workflow(
  wf_name = "cp1_tibbles",
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

# output one datafile per scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call(
    what = EpiModelHPC::merge_netsim_scenarios_tibble,
    args = list(
      sim_dir = "data/intermediate/scenarios/contourplots1",
      output_dir = "contourplots1_output",
      steps_to_keep = 364 * 10,
      cols = rlang::quo(dplyr::matches("^doxy")) # rlang::quo required for lazy eval
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "END",
    "cpus-per-task" = 1,
    "time" = "02:00:00",
    "mem" = "15G"
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



