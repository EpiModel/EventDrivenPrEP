##
## 01. Network Model Estimation
##

# Required variables:
#   - ncores
if (interactive()) {
  ncores <- 2
}

## Packages ##
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
suppressMessages(library("ARTnet"))

# Load the `NETSIZE` value and the formatted `netsize_string`
# NETSIZE <- 1e4 # to override (before sourcing the file)
source("R/utils-netsize.R")

# 0. Initialize Network --------------------------------------------------------


epistats <- build_epistats(
  geog.lvl = "city",
  geog.cat = "Atlanta",
  init.hiv.prev = c(0.33, 0.137, 0.084),
  race = TRUE,
  time.unit = 7
)
saveRDS(epistats, file = "data/input/epistats.rds")

netparams <- build_netparams(
  epistats = epistats,
  smooth.main.dur = TRUE
)

netstats <- build_netstats(
  epistats,
  netparams,
  expect.mort = 0.000478213,
  network.size = NETSIZE
)
saveRDS(netstats, file = paste0("data/input/netstats-", netsize_string, ".rds"))

num <- netstats$demog$num
nw <- EpiModel::network_initialize(num, directed = FALSE)

attr.names <- names(netstats$attr)
attr.values <- netstats$attr
nw <- EpiModel::set_vertex_attribute(nw, attr.names, attr.values)
nw_main <- nw_casl <- nw_inst <- nw

# 1. Main Model ----------------------------------------------------------------

# Formula
model_main <- ~ edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  concurrent +
  degrange(from = 3) +
  nodematch("role.class", diff = TRUE, levels = 1:2)

# Target Stats
netstats_main <- c(
  edges = netstats$main$edges,
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  nodefactor_age.grp = netstats$main$nodefactor_age.grp[-1],
  nodematch_race = netstats$main$nodematch_race_diffF,
  nodefactor_race = netstats$main$nodefactor_race[-1],
  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
  concurrent = netstats$main$concurrent,
  degrange = 0,
  nodematch_role.class = c(0, 0)
)
netstats_main <- unname(netstats_main)

# Fit model
fit_main <- netest(
  nw_main,
  formation = model_main,
  target.stats = netstats_main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm = control.ergm(
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 10000,
    MCMC.interval = 5000,
    parallel = ncores
  ),
  verbose = FALSE
)
fit_main <- trim_netest(fit_main)

# 2. Casual Model ---------------------------------------------------------

# Formula
model_casl <- ~ edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = c(-1, -5)) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.main", levels = -3) +
  concurrent +
  degrange(from = 4) +
  nodematch("role.class", diff = TRUE, levels = 1:2)

# Target Stats
netstats_casl <- c(
  edges = netstats$casl$edges,
  nodematch_age.grp = netstats$casl$nodematch_age.grp,
  nodefactor_age.grp = netstats$casl$nodefactor_age.grp[-c(1, 5)],
  nodematch_race = netstats$casl$nodematch_race_diffF,
  nodefactor_race = netstats$casl$nodefactor_race[-1],
  nodefactor_deg.main = netstats$casl$nodefactor_deg.main[-3],
  concurrent = netstats$casl$concurrent,
  degrange = 0,
  nodematch_role.class = c(0, 0)
)
netstats_casl <- unname(netstats_casl)

# Fit model
fit_casl <- netest(
  nw_casl,
  formation = model_casl,
  target.stats = netstats_casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm = control.ergm(
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 10000,
    MCMC.interval = 5000,
    parallel = ncores
  ),
  verbose = FALSE
)
fit_casl <- trim_netest(fit_casl)

# 3. One-Off Model -------------------------------------------------------------

# Formula
model_inst <- ~ edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("risk.grp", levels = -5) +
  nodefactor("deg.tot", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = 1:2)

# Target Stats
netstats_inst <- c(
  edges = netstats$inst$edges,
  nodematch_age.grp = sum(netstats$inst$nodematch_age.grp),
  nodefactor_age.grp = netstats$inst$nodefactor_age.grp[-1],
  nodematch_race = netstats$inst$nodematch_race_diffF,
  nodefactor_race = netstats$inst$nodefactor_race[-1],
  nodefactor_risk.grp = netstats$inst$nodefactor_risk.grp[-5],
  nodefactor_deg.tot = netstats$inst$nodefactor_deg.tot[-1],
  nodematch_role.class = c(0, 0)
)
netstats_inst <- unname(netstats_inst)

# Fit model
fit_inst <- netest(
  nw_inst,
  formation = model_inst,
  target.stats = netstats_inst,
  coef.diss = dissolution_coefs(~ offset(edges), 1),
  set.control.ergm = control.ergm(
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 10000,
    MCMC.interval = 5000,
    parallel = ncores
  ),
  verbose = FALSE
)
fit_inst <- trim_netest(fit_inst)

# 4. Save Data -----------------------------------------------------------------

out <- list(fit_main = fit_main, fit_casl = fit_casl, fit_inst = fit_inst)
saveRDS(out, file = paste0("data/input/netest-", netsize_string, ".rds"))


# Extra Control Args ------------------------------------------------------

# fit_main <- netest(nw_main,
#                    formation = model_main,
#                    target.stats = netstats_main,
#                    coef.diss = netstats$main$diss.byage,
#                    set.control.ergm = control.ergm(init.method = "MPLE",
#                                            MCMLE.effectiveSize = NULL,
#                                            MCMC.burnin = 1e6,
#                                            MCMC.interval = 1e5,
#                                            MCMC.samplesize = 10000,
#                                            init.MPLE.samplesize = 2e8,
#                                            MPLE.constraints.ignore = TRUE,
#                                            parallel = 10,
#                                            SAN.nsteps = 2e8),
#                    verbose = TRUE)

