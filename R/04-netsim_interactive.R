
##
## 04. Small-scale epidemic simulation for testing/debugging
##

## Packages
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

# Load the `NETSIZE` value and the formatted `netsize_string`
# NETSIZE <- 1e4     # to override (before sourcing the file)
source("R/utils-netsize.R")

## Parameters
epistats <- readRDS("data/input/epistats.rds")
netstats <- readRDS(paste0("data/input/netstats-", netsize_string, ".rds"))
est <- readRDS(paste0("data/input/netest-", netsize_string, ".rds"))

time.unit <- 1

param <- param_msm(netstats = netstats,
                   epistats = epistats,

                   # Clinical
                   hiv.test.rate = c((0.00385/7)*time.unit,
                                     (0.00380/7)*time.unit,
                                     (0.00690/7)*time.unit),
                   test.window.int = 21 / time.unit,
                   tx.init.rate = c((0.1775/7)*time.unit,
                                  (0.190/7)*time.unit,
                                  (0.2521/7)*time.unit),
                   tx.halt.partial.rate = c((0.0062/7)*time.unit,
                                            (0.0055/7)*time.unit,
                                            (0.0031/7)*time.unit),
                   tx.reinit.partial.rate = c((0.00255/7)*time.unit,
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
                                           (0.5/7)*time.unit),

                   riskh.start            = 52*7,
                   prep.start             = 2*7,
                   prep.start.prob        = rep(0.66, 3)
)

init <- init_msm()

control <- control_msm(
  simno = 1,
  nsteps = 100*7,
  nsims = 1,
  ncores = 1,
  verbose = TRUE
)

sim <- netsim(est, param, init, control)
