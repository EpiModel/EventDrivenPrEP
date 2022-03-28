
##
## 03. Epidemic Model Burnin, Stage 1, Parameter Calibration
##

library(remotes)
remotes::install_github("EpiModel/EpiModelHIV-p@EDP-mod-timeupdate")
remotes::install_github("EpiModel/ARTnet")
renv::snapshot()

#pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModelHIV-p")
#pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModel")

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

pull_env_vars()

network.size <- 10000
time.unit <- 7

fn.fx <- function(file, tu = time.unit, ns = network.size)
  paste0("data/input/", file, "-time", tu, "-size", ns, ".rds")

## Parameters
epistats <- readRDS(fn.fx("epistats"))
netstats <- readRDS(fn.fx("netstats"))
est <- readRDS(fn.fx("netest"))

param <- param_msm(netstats = netstats,
                   epistats = epistats,

                   # Clinical
                   hiv.test.rate = c((0.00385/7)*time.unit,
                                     (0.00380/7)*time.unit,
                                     (0.00690/7)*time.unit),
                   test.window.int = 21 / time.unit,
                   tx.init.prob = c((0.1775/7)*time.unit,
                                    (0.190/7)*time.unit,
                                    (0.2521/7)*time.unit),
                   tx.halt.partial.prob = c((0.0062/7)*time.unit,
                                            (0.0055/7)*time.unit,
                                            (0.0031/7)*time.unit),
                   tx.reinit.partial.prob = c((0.00255/7)*time.unit,
                                              (0.00255/7)*time.unit,
                                              (0.00255/7)*time.unit),

                   # HIV natural history
                   max.time.off.tx.full.int = (364/time.unit) * 15,
                   max.time.on.tx.partial.int = (364/time.unit) * 10,
                   max.time.off.tx.partial.int = (364/time.unit) * 10,
                   vl.acute.rise.int = (6.4/7)*time.unit,
                   vl.acute.fall.int = (6.4/7)*time.unit,
                   vl.aids.onset.int = (520/7)*time.unit,
                   vl.aids.int = (104/7)*time.unit,
                   vl.tx.down.slope = (0.25/7)*time.unit,
                   vl.tx.aids.down.slope = (0.25/7)*time.unit,
                   vl.tx.up.slope = (0.25/7)*time.unit,

                   # Demographic
                   a.rate = 0.00049,

                   # HIV transmission prob
                   hiv.trans.scale = c(2.44, 0.424, 0.270),

                   # STI epi
                   rgc.ntx.int = (16.8/7)*time.unit,
                   ugc.ntx.int = (16.8/7)*time.unit,
                   gc.tx.int = (1.4/7)*time.unit,
                   rct.ntx.int = (32/7)*time.unit,
                   uct.ntx.int = (32/7)*time.unit,
                   ct.tx.int = (1.4/7)*time.unit,

                   # PrEP
                   riskh.start = (364/time.unit)*59,
                   prep.start = ((364/time.unit)*60) + 1,
                   prep.start.prob = 0.66,
                   prep.tst.int = 90 / time.unit,
                   prep.risk.int = 182 / time.unit,
                   prep.sti.screen.int = 182 / time.unit,
                   prep.risk.reassess.int = 364/time.unit,

                   # Partner notification
                   part.ident.main.window.int = (12/7)*time.unit,
                   part.ident.casl.window.int = (12/7)*time.unit,
                   part.ident.ooff.window.int = (12/7)*time.unit,
                   part.prep.start.prob = c((0.5/7)*time.unit,
                                            (0.5/7)*time.unit,
                                            (0.5/7)*time.unit),
                   part.tx.init.prob = c((0.6/7)*time.unit,
                                         (0.6/7)*time.unit,
                                         (0.8/7)*time.unit),
                   part.tx.halt.prob = c((0.00102/7)*time.unit,
                                         (0.00102/7)*time.unit,
                                         (0.00071/7)*time.unit),
                   part.tx.reinit.prob = c((0.5/7)*time.unit,
                                           (0.5/7)*time.unit,
                                           (0.5/7)*time.unit)
                   )

init <- init_msm()

control <- control_msm(simno = fsimno,
                       nsteps = (364/time.unit) * 60,
                       nsims = ncores,
                       ncores = ncores,
                       verbose = TRUE)



# Weekly Test (before changes)
library(remotes)
remotes::install_github("EpiModel/EpiModelHIV-p")
remotes::install_github("EpiModel/ARTnet", ref = 'v2.5.0')

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

pull_env_vars()

network.size <- 10000
time.unit <- 7

fn.fx <- function(file, tu = time.unit, ns = network.size)
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
                       verbose = TRUE,
                       verbose.int = 250,
                       verbose.FUN = verbose.hpc.net)

# Testing time-dependent modules

at <- 2
dat <- initialize_msm(est, param, init, control)

dat <- aging_msm(dat, at)
dat <- arrival_msm(dat, at)
dat <- departure_msm(dat, at)
dat <- hivtest_msm(dat, at)
dat <- hivtx_msm(dat, at)
dat <- hivprogress_msm(dat, at)
dat <- hivvl_msm(dat, at)
dat <- simnet_msm(dat, at)
dat <- acts_msm(dat, at)
dat <- condoms_msm(dat, at)
dat <- position_msm(dat, at)
dat <- prep_msm(dat, at)
dat <- hivtrans_msm(dat, at)
dat <- stitrans_msm(dat, at)
dat <- stirecov_msm(dat, at)
dat <- stitx_msm(dat, at)
dat <- prevalence_msm(dat, at)

#check rates
debug(acts_msm)
acts_msm(dat, at)
undebug(acts_msm)

debug(prep_msm)
prep_msm(dat, at)
undebug(prep_msm)

debug(riskhist_msm)
riskhist_msm(dat, at)
undebug(riskhist_msm)

# Check ir100 by race and by STI
debug(prevalence_msm)
prevalence_msm(dat, at)
undebug(riskhist_msm)



## Simulation
sim <- netsim(est, param, init, control)

## Inspect simulation
sim_df <- as.data.frame(sim)
summary(sim_df)
sim_df

## Save-Min
savesim(sim, save.min = TRUE, save.max = FALSE, compress = TRUE, data.dir = "data/output/")
