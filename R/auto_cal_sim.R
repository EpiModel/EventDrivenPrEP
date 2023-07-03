model_fun <- function(proposal) {
  # Libraries ------------------------------------------------------------------
  library("EpiModelHIV")
  library("dplyr")

  # Settings -------------------------------------------------------------------
  source("./R/utils-0_project_settings.R")
  context <- "hpc"

  # Inputs ---------------------------------------------------------------------
  source("./R/utils-default_inputs.R", local = TRUE)
  est <- readRDS(path_to_est)

  source("R/utils-targets.R")
  control <- control_msm(
    nsteps              = calibration_end,
    nsims               = 1,
    ncores              = 1,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = 0,
    .tracker.list       = calibration_trackers,
    # .checkpoint.dir     = "./temp/cp_calib",
    # .checkpoint.clear   = TRUE,
    # .checkpoint.steps   = 15 * 364,
    verbose             = FALSE
  )

  # Proposal to scenario -------------------------------------------------------
  scenario_df <- proposal

  scenario_df[[".scenario.id"]] <- scenario_df[[".proposal_index"]]
  scenario_df[[".at"]] <- 1
  scenario_df[[".proposal_index"]] <- NULL
  scenario_df[[".wave"]] <- NULL
  scenario_df[[".iteration"]] <- NULL
  scenario <- EpiModel::create_scenario_list(scenario_df)[[1]]

  param_sc <- EpiModel::use_scenario(param, scenario)

  # Simulation and processing --------------------------------------------------
  sim <- netsim(est, param_sc, init, control)

  as_tibble(sim) %>%
    mutate_calibration_targets() %>%
    filter(time >= max(time) - 364) %>%
    select(c(sim, any_of(names(targets)))) %>%
    group_by(sim) %>%
    summarise(across( everything(), ~ mean(.x, na.rm = TRUE)))
}

calibration2_fun <- function(proposal) {
  # Libraries ------------------------------------------------------------------
  library("EpiModelHIV")
  library("dplyr")

  # Settings -------------------------------------------------------------------
  source("./R/utils-0_project_settings.R")
  context <- "hpc"

  # Inputs ---------------------------------------------------------------------
  source("./R/utils-default_inputs.R", local = TRUE)
  est <- readRDS(path_to_restart)

  source("R/utils-targets.R")
  control <- control_msm(
    start               = restart_time,
    nsteps              = intervention_start,
    nsims               = 1,
    ncores              = 1,
    initialize.FUN      = reinit_msm,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = 0,
    .tracker.list       = calibration_trackers,
    verbose             = FALSE
  )

  # Proposal to scenario -------------------------------------------------------
  scenario <- EpiModelHPC::swfcalib_proposal_to_scenario(proposal)
  param_sc <- EpiModel::use_scenario(param, scenario)

  # Simulation and processing --------------------------------------------------
  sim <- netsim(est, param_sc, init, control)

  as_tibble(sim) %>%
    mutate_calibration_targets() %>%
    filter(time >= max(time) - 364) %>%
    select(c(sim, any_of(names(targets)))) %>%
    group_by(sim) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
}
