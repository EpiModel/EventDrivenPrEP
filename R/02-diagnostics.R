##
## 02. Network Model Diagnostics
##

## Packages ##
rm(list = ls())
library("methods")
suppressMessages(library("EpiModelHIV"))

est <- readRDS("data/input/netest-10k.rds")

#ncores <- 30
ncores <- 1
nsims <- 10
nsteps <- 1000

# Main --------------------------------------------------------------------

fit_main <- est[[1]]

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
  fit_main, nsims = nsims, ncores = ncores, nsteps = nsteps,
  nwstats.formula = model_main_dx, skip.dissolution = TRUE,
  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5),
  set.control.stergm = control.simulate.network(MCMC.burnin.min = 2e5))

dx_main_static <- EpiModel::netdx(
  fit_main, dynamic = FALSE, nsims = 10,
  nwstats.formula = model_main_dx, skip.dissolution = TRUE,
  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))


# Casual ------------------------------------------------------------------

fit_casl <- est[[2]]

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
  fit_casl, nsims = nsims, ncores = ncores, nsteps = nsteps,
  nwstats.formula = model_casl_dx, skip.dissolution = TRUE,
  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5),
  set.control.stergm = control.simulate.network(MCMC.burnin.min = 2e5))

dx_casl_static <- netdx(
  fit_casl, dynamic = FALSE, nsims = 10,
  nwstats.formula = model_casl_dx, skip.dissolution = TRUE,
  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))


# One-Off -----------------------------------------------------------------

fit_inst <- est[[3]]

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
  fit_inst, nsims = 10, dynamic = FALSE,
  nwstats.formula = model_inst_dx,
  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5)) #other options can improve fit

dx <- list(dx_main = dx_main, dx_main_static = dx_main_static,
           dx_casl = dx_casl, dx_casl_static = dx_casl_static,
           dx_inst = dx_inst)

saveRDS(dx, file = "data/input/netdx.rds")


# Interactive Dx Analysis -------------------------------------------------

if (interactive()) {

  netstats <- readRDS("data/input/netstats-10k.rds")
  dx <- readRDS("data/input/netdx.rds")

  # Main
  print(dx$dx_main, digits = 2)
  plot(dx$dx_main)

  netstats$main

  print(dx$dx_main_static, digits = 2)
  plot(dx$dx_main_static)

  # Casual
  print(dx$dx_casl, digits = 2)
  plot(dx$dx_casl)

  netstats$casl

  print(dx$dx_casl_static, digits = 2)
  plot(dx$dx_casl_static)

  # Inst
  print(dx$dx_inst, digits = 2)
  plot(dx$dx_inst)

}

mcmc.diagnostics(fit_main$fit)
mcmc.diagnostics(fit_casl$fit)
mcmc.diagnostics(fit_inst$fit)

