##
## Epidemic Model Intervention Scenarios, HPC setup
##

# Libraries --------------------------------------------------------------------
library("slurmworkflow")
library("EpiModelHPC")
library("EpiModelHIV")
library("dplyr")

# Settings ---------------------------------------------------------------------
source("./R/utils-0_project_settings.R")
context <- "hpc"
max_cores <- 30

source("./R/utils-default_inputs.R") # make `path_to_est`, `param` and `init`
source("./R/utils-hpc_configs.R") # creates `hpc_configs`

#################################### Main analysis (Table 3) ######################################

# Workflow creation
wf <- create_workflow(
  wf_name = "table3",
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

# Intervention scenarios

sc_no_list <- list()
sc_no_list[["no_edp_sc"]] <- tibble(
  .scenario.id = "0_no_edp",
  .at = 1,
  edp.start = intervention_end + 1
)

edp.prep.start.prob <- c(9.78E-05, 0.0001, 0.000791814)

prop_change <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)

df <- data.frame(rep(NA, length(prop_change)))

for (j in edp.prep.start.prob) {
  edp.prep.start <- j*prop_change

  df <- cbind(df, edp.prep.start)
}

names(df) <- make.unique(names(df))

df <- df |>
  select(-1) |>
  rename(edp.prep.start.prob_1 = 1,
         edp.prep.start.prob_2 = 2,
         edp.prep.start.prob_3 = 3) |>
  slice(rep(1:n(), 2)) |>
  mutate(
    prop_change = rep(prop_change, 2),
    elig_scenario = rep(c(1,4), each = 10, length.out = length(prop_change))
  ) |>
  filter(!(elig_scenario == 4 & prop_change < 1))

df.elig <- data.frame(
  elig_scenario = 2:3,
  prop_change = rep(1,2),
  edp.prep.start.prob_1 = rep(edp.prep.start.prob[1], 2),
  edp.prep.start.prob_2 = rep(edp.prep.start.prob[2], 2),
  edp.prep.start.prob_3 = rep(edp.prep.start.prob[3], 2)
)

df <- rbind(df, df.elig)

df <- df |>
  arrange(elig_scenario, prop_change)

sc_interv_list <- list()
sc_interv_list[["table3"]] <- tibble(
  .scenario.id = paste0("scenario_", df$elig_scenario, "_prop_change_", df$prop_change),
  .at = 1,
  prep.start = intervention_start,
  prep.edp.start = intervention_start,
  edp.prep.start.prob_1 = df$edp.prep.start.prob_1,
  edp.prep.start.prob_2 = df$edp.prep.start.prob_2,
  edp.prep.start.prob_3 = df$edp.prep.start.prob_3,
  edp.start.scenario = df$elig_scenario
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
    output_dir = "./data/intermediate/scenarios/table3",
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
      sim_dir = "data/intermediate/scenarios/table3",
      output_dir = "data/output/table3",
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

######################## Sensitivity analysis: EDP adherence (Table 4) ##############################
# Workflow creation

wf <- create_workflow(
  wf_name = "table4",
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

df <- df |>
  rename(
    prep.adhr.edp.dist_1 = 1,
    prep.adhr.edp.dist_2 = 2,
    prep.adhr.edp.dist_3 = 3,
    prep.adhr.edp.dist_4 = 4
  ) |>
  slice(rep(1:n(), each = 4)) |>
  mutate(edp.start.scenario = rep(1:4, 12)) |>
  arrange(edp.start.scenario, prep.adhr.edp.dist_4)

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
    output_dir = "./data/intermediate/scenarios/table4",
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
      sim_dir = "data/intermediate/scenarios/table4",
      output_dir = "data/output/table4",
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

# Contour plot scenarios: Adherence vs. Coverage (Figure 2) --------------------------------
# Workflow creation
wf <- create_workflow(
  wf_name = "figure2",
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

edp.prep.start.prob <- c(9.78E-05, 0.0001, 0.000791814)

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
    output_dir = "./data/intermediate/scenarios/figure2",
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
      sim_dir = "data/intermediate/scenarios/figure2",
      output_dir = "data/output/figure2",
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

################################### Contour plot scenarios: PrEP switching ####################################
# Workflow creation
wf <- create_workflow(
  wf_name = "figure3",
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
    output_dir = "./data/intermediate/scenarios/figure3",
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
      sim_dir = "data/intermediate/scenarios/figure3",
      output_dir = "data/output/figure3",
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

############################# Contour plot scenarios: EDP vs. Daily PrEP ################################
# Workflow creation
wf <- create_workflow(
  wf_name = "contour_plots3",
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

prep.start.prob <- c(9.78E-05, 0.0001, 0.000791814)

prop_change <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2.0)
do_prop_change <- seq(1,2, 0.1)

df.edp.cov <- data.frame(rep(NA, length(prop_change)))
df.do.cov <- data.frame(rep(NA, length(do_prop_change)))

for (j in prep.start.prob) {
  edp.prep.start <- j*prop_change
  do.prep.start <- j*do_prop_change

  df.edp.cov <- cbind(df.edp.cov, edp.prep.start)
  df.do.cov <- cbind(df.do.cov, edp.prep.start)
}

df.edp.cov <- df.edp.cov[,-1]
names(df.edp.cov)[1] <- "edp.prep.start.prob_1"
names(df.edp.cov)[2] <- "edp.prep.start.prob_2"
names(df.edp.cov)[3] <- "edp.prep.start.prob_3"
df.edp.cov$edp_prop_change <- prop_change

df.do.cov <- df.do.cov[,-1]
names(df.do.cov)[1] <- "do.prep.start.prob_1"
names(df.do.cov)[2] <- "do.prep.start.prob_2"
names(df.do.cov)[3] <- "do.prep.start.prob_3"
df.do.cov$do_prop_change <- do_prop_change

df.merge <- merge(df.edp.cov, df.do.cov)

contour_plots3 <- tibble(
  .scenario.id = paste0("edp_startprob_", df.merge$edp_prop_change, "_do_startprob_", df.merge$do_prop_change),
  .at = 1,
  prep.start = intervention_start,
  prep.edp.start = intervention_start,
  edp.prep.start.prob_1 = df.merge$edp.prep.start.prob_1,
  edp.prep.start.prob_2 = df.merge$edp.prep.start.prob_2,
  edp.prep.start.prob_3 = df.merge$edp.prep.start.prob_3,
  prep.start.prob_1 = df.merge$do.prep.start.prob_1,
  prep.start.prob_2 = df.merge$do.prep.start.prob_2,
  prep.start.prob_3 = df.merge$do.prep.start.prob_3,
  edp.start.scenario = 1
)

sc_contour_plots_list <- EpiModel::create_scenario_list(contour_plots3)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_restart, param, init, control,
    scenarios_list = sc_contour_plots_list,
    output_dir = "./data/intermediate/scenarios/contourplots3",
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
      sim_dir = "data/intermediate/scenarios/contourplots3",
      output_dir = "data/output/contourplots3",
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

## Test workflow ----------------------------------------------------------------------------------------

# Workflow creation
wf <- create_workflow(
  wf_name = "test",
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

# Intervention scenarios

sc_no_list <- list()
sc_no_list[["no_edp_sc"]] <- tibble(
  .scenario.id = "0_no_edp",
  .at = 1,
  edp.start = intervention_end + 1
)

sc_interv_list <- list()
sc_interv_list[["test"]] <- tibble(
  .scenario.id = paste0("scenario_", df_test$edp.start.scenario, "_prop_change_", df_test$prep.adhr.edp.dist_4),
  .at = 1,
  prep.start = intervention_start,
  prep.edp.start = intervention_start,
  edp.prep.start.prob_1 = ifelse(is.na(df_test$edp.prep.start.prob_1), param$edp.prep.start.prob[1],
                                 df_test$edp.prep.start.prob_1),
  edp.prep.start.prob_2 = ifelse(is.na(df_test$edp.prep.start.prob_2), param$edp.prep.start.prob[2],
                                 df_test$edp.prep.start.prob_2),
  edp.prep.start.prob_3 = ifelse(is.na(df_test$edp.prep.start.prob_3), param$edp.prep.start.prob[3],
                                 df_test$edp.prep.start.prob_3),
  edp.start.scenario = 1,
  prep.adhr.edp.dist_1 = ifelse(is.na(df_test$prep.adhr.edp.dist_1), param$prep.adhr.edp.dist[1],
                                 df_test$prep.adhr.edp.dist_1),
  prep.adhr.edp.dist_2 = ifelse(is.na(df_test$prep.adhr.edp.dist_2), param$prep.adhr.edp.dist[2],
                                 df_test$prep.adhr.edp.dist_2),
  prep.adhr.edp.dist_3 = ifelse(is.na(df_test$prep.adhr.edp.dist_3), param$prep.adhr.edp.dist[3],
                                 df_test$prep.adhr.edp.dist_3),
  prep.adhr.edp.dist_4 = ifelse(is.na(df_test$prep.adhr.edp.dist_4), param$prep.adhr.edp.dist[4],
                                df_test$prep.adhr.edp.dist_4),
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

# output one datafile per scenario
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call(
    what = EpiModelHPC::merge_netsim_scenarios_tibble,
    args = list(
      sim_dir = "data/intermediate/scenarios/test",
      output_dir = "output/test",
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



