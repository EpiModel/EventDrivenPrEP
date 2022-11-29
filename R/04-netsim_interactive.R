
##
## 04. Small-scale epidemic simulation for testing/debugging
##

## Packages
pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModelHIV-p")
suppressMessages(library("EpiModelHIV"))
library(dplyr)

## Parameters
epistats <- readRDS("data/intermediate/estimates/epistats.rds")
netstats <- readRDS("data/intermediate/estimates/netstats.rds")
est      <- readRDS("data/intermediate/estimates/netest.rds")

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
                   prep.tst.int = 90 / time.unit,
                   prep.risk.int = 182 / time.unit,
                   prep.sti.screen.int = 182 / time.unit,
                   prep.risk.reassess.int = 364/time.unit,
                   prep.discont.rate = rep(1 - (2 ^ (-1 / (224.4237))), 3), # divide 224.4237 by 7 for weekly time steps

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

                   riskh.start            = 1,
                   prep.start             = 26*7, # needs to start at least 6 months after riskh.start
                   prep.start.prob        = rep(0.66, 3),

                   # New EDP parameters
                   prep.edp.start     = 400, # the timeline used for LAI PrEP
                   prep.daily.prob    = 0.5, # the probability of starting daily vs. EDP
                   prep.adhr.dist.edp = c(0.11, 0.06, 0.09, 0.74), # 4 adherence classes for EDP
                   prep.adhr.rr.edp   = c(1, 0.24, 0.14, 0.03) # relative risk based on adherence class

                   # set LNT parameter to false

)

init <- init_msm()

control <- control_msm(
  simno = 1,
  nsteps = 364*2,
  nsims = 1,
  ncores = 7,
  verbose = TRUE,
  raw.output = FALSE # will output raw data including raw attribute vectors up until that time step
)

debug(hivtrans_msm)
undebug(hivtrans_msm)

options(error = recover)

sim <- netsim(est, param, init, control)

# Explore sim object

## Explore the number of people starting daily oral PrEP vs. EDP
sim$epi$edp.starters
sim$epi$daily.starters

x <- 1:728

prepDailyStart <- ifelse(is.na(sim$epi$daily.starters$sim1), 0, sim$epi$daily.starters$sim1)

cum.prepDailyStart <- cumsum(prepDailyStart)
cum.prepDailyStart

prepEDPStart <- ifelse(is.na(sim$epi$edp.starters$sim1), 0, sim$epi$edp.starters$sim1)
prepEDPStart

cum.prepEDPStart <- cumsum(prepEDPStart)
cum.prepEDPStart

plot(x, y = cum.prepDailyStart, type = "l", col = "red", xlab = "Day", ylab = "Cumulative Number of PrEP Starters",
     main = paste("Probability of starting Daily PrEP =", param$prep.daily.prob),
     sub = paste("Total MSM starting PrEP by day 600 =", cum.prepDailyStart[600] + cum.prepEDPStart[600]))
lines(x, y = cum.prepEDPStart, type = "l", col = "blue")
legend("topleft", legend = c("Daily PrEP", "Event-Driven PrEP"), col = c("red", "blue"), lty = 1)

## Explore the change in prep.daily.prob over time
plot(x, y = sim$epi$prep.daily.prob$sim1, xlab = "Day", ylab = "Probability of Daily Oral PrEP vs. EDP")

## Explore the distribution of EDP PrEP classes among those starting PrEP
### set raw.output = TRUE in control settings

a <- table(sim[[1]]$attr$prepClass.edp)
a/sum(a)

b <- table(sim[[1]]$attr$prepClass)
b/sum(b)

sim[[1]]$epi$edp.class.1
sim[[1]]$epi$edp.class.2
sim[[1]]$epi$edp.class.3
sim[[1]]$epi$edp.class.4

sim[[1]]$epi$incid.edp.1
sim[[1]]$epi$incid.edp.2
sim[[1]]$epi$incid.edp.3
sim[[1]]$epi$incid.edp.4

## Explore EDP adherence class distribution over time
plot(x, y = sim[[1]]$epi$edp.class.1, type = "l", col = "red",
     xlab = "Day", ylab = "Number of EDP Users",
     xlim = c(400, 728),
     ylim = c(0, max(unlist(sim[[1]]$epi$edp.class.4), na.rm = T)),
     main = "EDP Users by Adherence Class")
lines(x, sim[[1]]$epi$edp.class.2, type = "l", col = "blue")
lines(x, sim[[1]]$epi$edp.class.3, type = "l", col = "green")
lines(x, sim[[1]]$epi$edp.class.4, type = "l", col = "black")
legend("topleft", legend = c("None", "Bad", "Good", "Excellent"),
       col = c("red", "blue", "green", "black"), lty = 1)

# Explore HIV incidence by EDP use
sum(sim$epi$incid$sim1, na.rm = T)
sum(sim$epi$incid.edp$sim1, na.rm = T)
sum(sim$epi$incid.edp.1$sim1, na.rm = T)
sum(sim$epi$incid.edp.2$sim1, na.rm = T)
sum(sim$epi$incid.edp.3$sim1, na.rm = T)
sum(sim$epi$incid.edp.4$sim1, na.rm = T)


## Explore history of prepClass.edp attribute for EDP users

attr_history <- get_attr_history(sim)
attr_history

unique(attr_history$prepClass.edp$uids)
length(unique(attr_history$prepClass.edp$uids))

attr_history_merged <- left_join(attr_history$prepClass.edp, attr_history$sex.edp, by = c("time", "uids")) %>%
  left_join(., attr_history$prepTimeLastPill, by = c("time", "uids")) %>%
  select(time, uids, values.x, values.y, values) %>%
  rename("prepClass" = "values.x",
         "sex" = "values.y",
         "prepTimeLastPill" = "values")

id_7031 <- filter(attr_history_merged, uids == 7031)
id_37 <- filter(attr_history_merged, uids == 37)
