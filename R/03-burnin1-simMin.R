
##
## 03. Epidemic Model Burnin, Stage 1, Parameter Calibration
##

#library(remotes)
#remotes::install_github("EpiModel/EpiModelHIV-p@EDP-mod-timeupdate")
#renv::snapshot()

#pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModelHIV-p")
#pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModel")

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

pull_env_vars()

network.size <- 10000
time.unit <- 7

fn.fx <- function(file, tu = 7, ns = 10000)
  paste0("data/input/", file, "-time", tu, "-size", ns, ".rds")

## Parameters
epistats <- readRDS(fn.fx("epistats"))
netstats <- readRDS(fn.fx("netstats"))
est <- readRDS(fn.fx("netest"))

param <- param_msm(netstats = netstats,
                   epistats = epistats,
                   a.rate = 0.00049,
                   hiv.test.rate = c(0.00385, 0.00380, 0.00690),
                   tx.init.prob = c(0.1775, 0.190, 0.2521),
                   tx.halt.partial.prob = c(0.0062, 0.0055, 0.0031),
                   tx.reinit.partial.prob = c(0.00255, 0.00255, 0.00255),
                   trans.scale = c(2.44, 0.424, 0.270),
                   riskh.start = 52 * 59,
                   prep.start = (52 * 60) + 1,
                   prep.start.prob = 0.66)

init <- init_msm()

control <- control_msm(simno = fsimno,
                       nsteps = 52 * 60,
                       nsims = ncores,
                       ncores = ncores,
                       cumulative.edgelist = TRUE,
                       truncate.el.cuml = 0,
                       verbose = TRUE)
## Simulation
sim <- netsim(est, param, init, control)

## Save-Min
savesim(sim, save.min = TRUE, save.max = FALSE, compress = TRUE, data.dir = "data/output/")
