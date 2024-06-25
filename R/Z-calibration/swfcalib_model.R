## swfcalib Model Function
##
## Define helper functions to create the `model` function for an swfcalib
## process.
##
## This script should not be run directly. But `sourced` from the
## swfcalib_config scripts.

make_model_fn <- function(calib_steps) {
  force(calib_steps)
  function(proposal) {
    # Libraries ------------------------------------------------------------------
    library("EpiModelHIV")
    library("dplyr")

    # Inputs ---------------------------------------------------------------------
    source("./R/utils-0_project_settings.R")
    context <- "hpc"
    source("./R/utils-default_inputs.R")
    source("./R/utils-targets.R")

    est <- readRDS(path_to_est)
    control <- control_msm(
      nsteps              = calibration_end,
      .tracker.list       = calibration_trackers,
      verbose             = FALSE
    )

    # init$init_attr <- readRDS("./d_init_attr.rds")

    # Proposal to scenario -------------------------------------------------------
    scenario <- EpiModelHPC::swfcalib_proposal_to_scenario(proposal)
    param_sc <- EpiModel::use_scenario(param, scenario)

    param_sc$rgc.prob <- plogis(qlogis(param_sc$ugc.prob) + log(1.25))
    param_sc$rct.prob <- plogis(qlogis(param_sc$uct.prob) + log(1.25))

    # Simulation and processing --------------------------------------------------
    sim <- netsim(est, param_sc, init, control)

    as_tibble(sim) |>
      mutate_calibration_targets() |>
      filter(time >= max(time) - calib_steps) |>
      select(c(sim, any_of(names(targets)), num)) |>
      group_by(sim) |>
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  }
}
