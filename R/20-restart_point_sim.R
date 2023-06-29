##
## 20. Epidemic Model Restart Point, Local simulation runs
##

# Settings ---------------------------------------------------------------------
context <- "local"
source("R/utils-0_project_settings.R")

# Run the simulations ----------------------------------------------------------
library("EpiModelHIV")

# Necessary files
source("R/utils-default_inputs.R") # generate `path_to_est`, `param` and `init`

# Controls
source("R/utils-targets.R")
control <- control_msm(
  nsteps              = 364*5, #calibration_end,
  nsims               = 1,
  ncores              = 1,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
  verbose             = FALSE
)

param$tx.init.rate <- rep(0.1, 3)
# No scenarios are used here

EpiModelHPC::netsim_scenarios(
  path_to_est, param, init, control,
  scenarios_list = NULL,
  n_rep = 3,
  n_cores = 3,
  output_dir = "data/intermediate/calibration",
  libraries = "networkLite",
  save_pattern = "restart" # more data is required to allow restarting
)

# Check the files produced
list.files("data/intermediate/calibration")

d <- readRDS("./data/intermediate/calibration/sim__empty_scenario__1.rds") |>
  as_tibble()
tail(d$i__B)
tail(d$i_dx__B)
tail(d$i_tx__B)

# tx.init.rate_1,0.2982484,numeric
# tx.init.rate_2,0.3677387,numeric
# tx.init.rate_3,0.357965,numeric
#
# tx.init.rate_1,0.042600143,numeric,
# tx.init.rate_2,0.052528857,numeric,
# tx.init.rate_3,0.051131543,numeric,

# P2p <- function(P, n) 1 - (1 - P)^(1/n)
# p2P <- function(p, n) 1 - (1 - p)^n
#
simp <- function(p, n) any(runif(n) < p)
replicate(1e4, simp(0.1, 28)) |> mean()
