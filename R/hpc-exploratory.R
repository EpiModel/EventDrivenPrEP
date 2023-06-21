## Explore HPC results
library(EpiModelHPC)

df_raw <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023/assessments_raw.rds")
df <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023/assessments.rds")
sim1_1 <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023/sim__EDPStart1__1.rds")
sim2_1 <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023/sim__EDPStart2__1.rds")
sim3_1 <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023/sim__EDPStart3__1.rds")
sim4_1 <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023/sim__EDPStart4__1.rds")


## Merge scenarios together
merge_netsim_scenarios("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/data/intermediate/calibration/7JUNE2023",
                       "test_merged")

sim1 <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/test_merged/merged__EDPStart2.rds")
sim2 <- readRDS("C:/Users/clchand/OneDrive - Emory University/EpiModel-repos/EventDrivenPrEP/test_merged/merged__EDPStart4.rds")

plot(sim1_1, y = "ir100")
plot(sim2_1, y = "ir100")
plot(sim4_1, y = "ir100")

sim <- sim4_1

plot(sim, y = c("prep.switch", "edp.switch"),
     sim.lines = TRUE,
     mean.line = FALSE,
     qnts = FALSE)

plot(sim, y = c("daily.starters", "edp.starters"),
     xlab = "Day", ylab = "Number of PrEP Starters",
     main = "PrEP Starters by Regimen",
     xlim = c(0, 4000),
     lty = 1,
     sim.lines = TRUE,
     sim.col = c("#B2182B", "#2166AC"),
     mean.line = FALSE,
     qnts = FALSE,
     sim.alpha = 0.1)
legend(2500, 250, legend = c("Daily PrEP",
                             "Event-Driven PrEP"),
       lty = 1,
       lwd = 2,
       col = c("#B2182B", "#2166AC"))

mean(sapply(sim$epi$edp.starters, sum, na.rm = T))
mean(sapply(sim$epi$daily.starters, sum, na.rm = T))

plot(sim, y = c("edp.class.1", "edp.class.2", "edp.class.3", "edp.class.4"),
     xlab = "Day", ylab = "Number of EDP Users",
     main = "EDP Users by Adherence Class",
     xlim = c(0, 4000),
     lty = 1,
     sim.lines = TRUE,
     sim.col = c("#B2182B", "#339966", "orange", "#2166AC"),
     mean.line = FALSE,
     qnts = FALSE)
     #mean.col = c("#B2182B", "#339966", "orange", "#2166AC"),
     #qnts.col = c("#B2182B", "#339966", "orange", "#2166AC"))
legend(1, 0.35, legend = c("None",
                               "Poor",
                               "Good",
                               "Excellent"),
       lty = 1,
       lwd = 2,
       col = c("#B2182B", "#339966", "orange", "#2166AC"))

mean(sapply(sim$epi$edp.class.1, mean, na.rm = TRUE))
mean(sapply(sim$epi$edp.class.2, mean, na.rm = TRUE))
mean(sapply(sim$epi$edp.class.3, mean, na.rm = TRUE))
mean(sapply(sim$epi$edp.class.4, mean, na.rm = TRUE))

# Explore HIV incidence by EDP use
incid.total <- sapply(sim$epi$incid, sum, na.rm = TRUE)
mean(incid.total)

incid.daily <- sapply(sim$epi$incid.daily, sum, na.rm = TRUE)
mean(incid.daily)

incid.edp <- sapply(sim$epi$incid.edp, sum, na.rm = TRUE)
mean(incid.edp)

incid.edp1 <- sapply(sim$epi$incid.edp.1, sum, na.rm = TRUE)
mean(incid.edp1)

incid.edp2 <- sapply(sim$epi$incid.edp.2, sum, na.rm = TRUE)
mean(incid.edp2)

incid.edp3 <- sapply(sim$epi$incid.edp.3, sum, na.rm = TRUE)
mean(incid.edp3)

incid.edp4 <- sapply(sim$epi$incid.edp.4, sum, na.rm = TRUE)
mean(incid.edp4)

d <- data.frame(incid.total, incid.daily, incid.edp, incid.edp1, incid.edp2, incid.edp3, incid.edp4)
d$total.edp <- rowSums(d[4:7])
d

write.csv(d, "incidence.csv")

# Explore number of pills per day per person
pillspp <- sapply(sim$epi$pillspp, mean, na.rm = TRUE)
mean(pillspp)
mean(pillspp)*7

pills <- sapply(sim$epi$pills, sum, na.rm = TRUE)
mean(pills)

hi <- mean(sapply(sim$epi$prepCurr.hi, sum, na.rm = TRUE)/7*5.5)
med <- mean(sapply(sim$epi$prepCurr.med, sum, na.rm = TRUE)/7*2.5)
low <- mean(sapply(sim$epi$prepCurr.lo, sum, na.rm = TRUE)/7)

daily.pills <- floor(hi+med+low)
daily.pills


# Plotting multiple simulations with ggplot -----------------------------------
