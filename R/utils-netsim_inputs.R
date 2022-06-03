##
## utils. Default Inputs for EpiModel::nestim
##

# Required variables: (provided by the calling script)
# - NETSIZE
# - netsize_string

# Setup ------------------------------------------------------------------------
suppressMessages({
  library("EpiModelHIV")
  library("EpiModelHPC")
})

epistats <- readRDS("data/input/epistats.rds")
netstats <- readRDS(paste0("data/input/netstats-", netsize_string, ".rds"))
est <- readRDS(paste0("data/input/netest-", netsize_string, ".rds"))

# Relevant times
calibration_length <- 52 * 10
prep_start         <- calibration_length + 52 * 5 + 1
interv_start       <- prep_start + 52 * 5
nsteps             <- interv_start + 52 * 10 - 1

# Parameters -------------------------------------------------------------------
df_params <- readr::read_csv("data/input/params.csv", )
param <- param.net(
  data.frame.params = df_params,
  netstats = netstats,
  epistats = epistats,
  riskh.start = prep_start - 53,
  prep.start = prep_start
)

# Initial conditions -----------------------------------------------------------
init <- init_msm()

