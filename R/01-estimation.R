##
## 01. Network Model Estimation
##

## Packages ##
rm(list = ls())
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("ARTnetData"))
library("ARTnet")
library("tidyverse")

ncores <- parallel::detectCores() - 1

# 0. Initialize Network ---------------------------------------------------

#install older version of ARTnet
#library(remotes)
#install_github("EpiModel/ARTnet", ref = 'v2.5.0')

#install local version of ARTnet for any changes
#pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\ARTnet")

epistats <- build_epistats(
  geog.lvl = "city",
  geog.cat = "Atlanta",
  init.hiv.prev = c(0.33, 0.137, 0.084),
  race = TRUE,
  #time.unit = 7, #change to 1 for daily
  browser = TRUE
)
saveRDS(epistats, file = "data/input/epistats.rds")

netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
netstats <- build_netstats(
  epistats,
  netparams,
  expect.mort = 0.000478213, #divide by 7 for daily
  network.size = 10000
)
saveRDS(netstats, file = "data/input/netstats-10k.rds")

num <- netstats$demog$num
nw <- EpiModel::network_initialize(num, directed = FALSE)


attr.names <- names(netstats$attr)
attr.values <- netstats$attr
nw <- EpiModel::set_vertex_attribute(nw, attr.names, attr.values)
nw_main <- nw_casl <- nw_inst <- nw
#replicates networks x3 for partnership types: main, casual, one-time


# 1. Main Model -----------------------------------------------------------

# Formula
model_main <- ~edges +
  nodematch("age.grp", diff = TRUE) + #heterogeneity
  nodefactor("age.grp", levels = -1) + #mean degree
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  concurrent +
  degrange(from = 3) + #max degree
  nodematch("role.class", diff = TRUE, levels = 1:2)
# role classes related to sexual positioning for MSM
# insertive man won't mix with insertive man
# versatile 3rd cat can mix with either insertive or receptive

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
cbind(netstats_main)
netstats_main <- unname(netstats_main)

# Fit model
fit_main <- netest(nw_main, #empty network
                   formation = model_main, #formula for frmation model
                   target.stats = netstats_main, #data for formation model
                   coef.diss = netstats$main$diss.byage, #frmula and data for the dissolution model
                   #allows avg duration in weeks for younger pairing to be shorter than
                   #avg duration for older pairings
                   set.control.ergm = control.ergm(MCMLE.maxit = 500, #controls MCMC estimation process
                                                   SAN.maxit = 4,
                                                   SAN.nsteps.times = 10,
                                                   MCMC.samplesize = 10000,
                                                   MCMC.interval = 5000,
                                                   parallel = ncores),
                   verbose = FALSE)

summary(fit_main)

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



# 2. Casual Model ---------------------------------------------------------

# Formula
model_casl <- ~edges +
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
cbind(netstats_casl)
netstats_casl <- unname(netstats_casl)

# Fit model
fit_casl <- netest(nw_casl,
                   formation = model_casl,
                   target.stats = netstats_casl,
                   coef.diss = netstats$casl$diss.byage,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 4,
                                                   SAN.nsteps.times = 10,
                                                   MCMC.samplesize = 10000,
                                                   MCMC.interval = 5000,
                                                   parallel = ncores),
                   verbose = FALSE)


# 3. One-Off Model --------------------------------------------------------

# Formula
model_inst <- ~edges +
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
cbind(netstats_inst)
netstats_inst <- unname(netstats_inst)

# Fit model
fit_inst <- netest(nw_inst,
                   formation = model_inst,
                   target.stats = netstats_inst,
                   coef.diss = dissolution_coefs(~offset(edges), 1),
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 4,
                                                   SAN.nsteps.times = 10,
                                                   MCMC.samplesize = 10000, # can be increased to improve fit
                                                   MCMC.interval = 5000,
                                                   parallel = ncores),
                   verbose = FALSE)


# 4. Save Data ------------------------------------------------------------

fit_main$fit$newnetworks <- NULL
fit_casl$fit$newnetworks <- NULL
fit_inst$fit$newnetworks <- NULL

out <- list(fit_main = fit_main, fit_casl = fit_casl, fit_inst = fit_inst)
saveRDS(out, file = "data/input/netest-10k.rds")

#refit model to daily time steps instead of weekly time steps
#will need to go in to ARTnet
#dissolution model is hardcoded in weekly time steps
#unclear if we need to manually

#main and casl partnerships don't need to be changed
#dissolution models fr main and casl need to be adjusted for daily time steps
#models need to be refit accordingly

#for inst models - no formation/dissolution
#all interpreted as incidence of partnerships
#need to be adjusted for by size of time step

#look at ARTnet package dcumentation that builds model and data out
#e.g., epistats_build (see what needs to be updated)


