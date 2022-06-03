##
## 02. Network Model Diagnostics
##

# Required variables:
#   - ncores
#   - nsims
#   - nsteps
if (interactive()) {
  ncores <- 2
  nsims <- 10
  nsteps <- 500
}

# Setup ------------------------------------------------------------------------
suppressMessages({
  library("EpiModelHIV")
})

# Load the `NETSIZE` value and the formatted `netsize_string`
# NETSIZE <- 1e4 # to override (before sourcing the file)
source("R/utils-netsize.R")

fn <- paste0("data/input/netest-", netsize_string, ".rds")
est <- readRDS(fn)

# Main -------------------------------------------------------------------------

fit_main <- est[["fit_main"]]

model_main_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = TRUE) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = TRUE) +
  nodefactor("deg.casl", levels = TRUE) +
  degrange(from = 3) +
  concurrent +
  nodematch("role.class", diff = TRUE) +
  degree(0:3)

dx_main <- netdx(
  fit_main,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  nwstats.formula = model_main_dx,
  skip.dissolution = TRUE,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
  set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 2e5)
)

dx_main_static <- EpiModel::netdx(
  fit_main,
  dynamic = FALSE,
  nsims = 10000,
  nwstats.formula = model_main_dx,
  skip.dissolution = TRUE,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

dx <- list(dx_main = dx_main, dx_main_static = dx_main_static)
fn <- paste0("data/input/netdx-main-", netsize_string, ".rds")
saveRDS(dx, file = fn)
rm(dx, dx_main, dx_main_static)

# Casual -----------------------------------------------------------------------

fit_casl <- est[["fit_casl"]]

model_casl_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = TRUE) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = TRUE) +
  nodefactor("deg.main", levels = TRUE) +
  degrange(from = 4) +
  concurrent +
  nodematch("role.class", diff = TRUE) +
  degree(0:4)

dx_casl <- netdx(
  fit_casl,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  nwstats.formula = model_casl_dx,
  skip.dissolution = TRUE,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
  set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 2e5)
)

dx_casl_static <- netdx(
  fit_casl,
  dynamic = FALSE,
  nsims = 10000,
  nwstats.formula = model_casl_dx,
  skip.dissolution = TRUE,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

dx <- list(dx_casl = dx_casl, dx_casl_static = dx_casl_static)
fn <- paste0("data/input/netdx-casl-", netsize_string, ".rds")
saveRDS(dx, file = fn)
rm(dx, dx_casl, dx_casl_static)

# One-Off ----------------------------------------------------------------------

fit_inst <- est[["fit_inst"]]

model_inst_dx <- ~edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = TRUE) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = TRUE) +
  nodefactor("risk.grp", levels = TRUE) +
  nodefactor("deg.tot", levels = TRUE) +
  nodematch("role.class", diff = TRUE) +
  degree(0:4)

dx_inst <- netdx(
  fit_inst,
  nsims = 10000,
  dynamic = FALSE,
  nwstats.formula = model_inst_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

dx <- list(dx_inst = dx_inst)
fn <- paste0("data/input/netdx-inst-", netsize_string, ".rds")
saveRDS(dx, file = fn)
