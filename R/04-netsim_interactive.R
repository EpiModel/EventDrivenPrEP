
##
## 04. Small-scale epidemic simulation for testing/debugging
##

## Packages
pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModelHIV-p")
#suppressMessages(library("EpiModelHIV"))
library(dplyr)
library(ggplot2)

## Parameters
epistats <- readRDS("data/intermediate/estimates/epistats.rds")
netstats <- readRDS("data/intermediate/estimates/netstats.rds")
est      <- readRDS("data/intermediate/estimates/netest.rds")

time.unit <- 1

param <- param_msm(netstats = netstats,
                   epistats = epistats,

                   # Clinical
                   hiv.test.rate = c((0.01325/7)*time.unit,
                                     (0.0125/7)*time.unit,
                                     (0.0124/7)*time.unit),
                   test.window.int = 21 / time.unit,
                   tx.init.rate = c((0.092/7)*time.unit,
                                  (0.092/7)*time.unit,
                                  (0.127/7)*time.unit),
                   tx.halt.partial.rate = c((0.0102/7)*time.unit,
                                            (0.0102/7)*time.unit,
                                            (0.0071/7)*time.unit),
                   tx.reinit.partial.rate = c((0.00066/7)*time.unit,
                                              (0.00066/7)*time.unit,
                                              (0.00291/7)*time.unit),

                   # HIV natural history
                   max.time.off.tx.full.int = (364/time.unit) * 15,
                   max.time.on.tx.partial.int = (364/time.unit) * 10,
                   max.time.off.tx.partial.int = (364/time.unit) * 10,
                   vl.acute.rise.int = (3/7)*time.unit,
                   vl.acute.fall.int = (3/7)*time.unit,
                   vl.aids.onset.int = (520/7)*time.unit,
                   vl.aids.int = (104/7)*time.unit,
                   vl.tx.down.rate = (0.25/7)*time.unit,
                   vl.tx.aids.down.rate = (0.25/7)*time.unit,
                   vl.tx.up.rate = (0.25/7)*time.unit,

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
                   prep.discont.int = c(33.42*7, 57.48*7, 57.39*7),
                   edp.start.scenario = 2,

                   # Partner notification
                   part.ident.main.window.int = (12/7)*time.unit,
                   part.ident.casl.window.int = (12/7)*time.unit,
                   part.ident.ooff.window.int = (12/7)*time.unit,
                   part.prep.start.prob = c((0.5/7)*time.unit,
                                            (0.5/7)*time.unit,
                                            (0.5/7)*time.unit),
                   part.tx.init.rate = c((0.6/7)*time.unit,
                                         (0.6/7)*time.unit,
                                         (0.8/7)*time.unit),
                   part.tx.reinit.rate = c((0.5/7)*time.unit,
                                           (0.5/7)*time.unit,
                                           (0.5/7)*time.unit),

                   # PrEP start
                   riskh.start = 1,
                   prep.start = 182

)

init <- init_msm()

control <- control_msm(
  simno = 1,
  nsteps = 364*2,
  nsims = 1,
  ncores = 1,
  verbose = TRUE,
  raw.output = FALSE # will output raw data including raw attribute vectors up until that time step
)

#debug(prep_msm)
#undebug(prep_msm)

#options(error = recover)

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
     main = paste("EDP Start Scenario =", param$edp.start.scenario),
     sub = paste("Total MSM starting PrEP by day 600 =", cum.prepDailyStart[600] + cum.prepEDPStart[600]))
lines(x, y = cum.prepEDPStart, type = "l", col = "blue")
legend("topleft", legend = c("Daily PrEP", "Event-Driven PrEP"), col = c("red", "blue"), lty = 1)

## Explore the distribution of EDP PrEP classes among those starting PrEP
### set raw.output = TRUE in control settings

b <- table(sim$attr$sim1$prepClass)
prop.table(b)

sim$epi$edp.class.1
sim$epi$edp.class.2
sim$epi$edp.class.3
sim$epi$edp.class.4

sum(sim$epi$edp.class.1$sim1, na.rm = TRUE)
sum(sim$epi$edp.class.2$sim1, na.rm = TRUE)
sum(sim$epi$edp.class.3$sim1, na.rm = TRUE)
sum(sim$epi$edp.class.4$sim1, na.rm = TRUE)

## Explore EDP adherence class distribution over time
plot(x, y = as.vector(unlist(sim$epi$edp.class.1)), type = "l", col = "red",
     xlab = "Day", ylab = "Number of EDP Users",
     xlim = c(400, 728),
     ylim = c(0, max(as.vector(unlist(sim$epi$edp.class.4)), na.rm = T)),
     main = "EDP Users by Adherence Class")
lines(x, as.vector(unlist(sim$epi$edp.class.2)), type = "l", col = "blue")
lines(x, as.vector(unlist(sim$epi$edp.class.3)), type = "l", col = "green")
lines(x, as.vector(unlist(sim$epi$edp.class.4)), type = "l", col = "black")
legend("topleft", legend = c("None", "Bad", "Good", "Excellent"),
       col = c("red", "blue", "green", "black"), lty = 1)

## Stacked bar chart for adherence

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

unique(attr_history$edp.prepClass$uids)
length(unique(attr_history$edp.prepClass$uids))

prop.table(table(attr_history$edp.prepClass$values))

attr_history_merged <- left_join(attr_history$edp.prepClass, attr_history$sex.edp, by = c("time", "uids")) %>%
  left_join(., attr_history$lastPrepCombo, by = c("time", "uids")) %>%
  left_join(., attr_history$timeLastSex, by = c("time", "uids")) %>%
  select(time, uids, values.x, values.y, values.x.x, values.y.y) %>%
  rename("prepClass" = "values.x",
         "sex" = "values.y",
         "lastPrepCombo" = "values.x.x",
         "timeLastSex" = "values.y.y") %>%
  mutate(timeSinceLastSex = time- timeLastSex)

id_a <- filter(attr_history_merged, uids == 4493)
id_b <- filter(attr_history_merged, uids == 107)
id_c <- filter(attr_history_merged, uids == 400)
id_d <- filter(attr_history_merged, uids == 3169)
id_e <- filter(attr_history_merged, uids == 7013)

#write.csv(id_d, "id_d.csv")

# Plotting PrEP adherence class over time

ggplot(data = attr_history_merged, aes(x = time, y = factor(uids), group = prepClass)) +
  geom_point(aes(color = factor(prepClass)))

attr_history_trunc_time <- attr_history_merged %>%
  filter(time >= 600) %>%
  filter(3000<= uids & uids <= 5000)

ggplot(data = attr_history_trunc_time, aes(x = time, y = factor(uids), group = prepClass)) +
  geom_point(aes(color = factor(prepClass))) +
  theme_minimal()

attr_history_trunc_uids <- attr_history_merged %>%
  filter(600 <= time) %>%
  filter(2500 <= uids & uids <= 2600)

ggplot(data = attr_history_trunc_uids, aes(x = time, y = factor(uids))) +
  geom_point(aes(color = factor(prepClass))) +
  labs(color="EDP Adherence\nClass") +
  xlab("Day") +
  ylab("EDP User ID") +
  ggtitle("Adherence Class Trajectories among EDP Users") +
  scale_color_discrete(labels=c("None", "Poor", "Good", "Excellent", "NA"))


# Stacked bar chart

prepClass <- attr_history_merged %>%
  select(time, prepClass) %>%
  group_by(time, prepClass) %>%
  summarise(value = n())

time <- unique(prepClass$time)
length(time)

ggplot(data = prepClass,
       aes(x = time, y = value, fill = factor(prepClass))) +
         geom_bar(position="fill", stat="identity")

ggplot(data = prepClass,
       aes(x = time, y = value, fill = factor(prepClass))) +
  geom_bar(position="stack", stat="identity") +
  labs(fill="EDP Adherence\nClass") +
  xlab("Day") +
  ylab("EDP Users") +
  ggtitle("Distribution of Adherence Classes among EDP Users") +
  scale_fill_discrete(labels=c("None", "Poor", "Good", "Excellent", "NA"))

prepClass_4 <- attr_history_merged %>%
  select(time, prepClass) %>%
  group_by(time, prepClass) %>%
  summarise(value = n()) %>%
  filter(prepClass == 4)

prepClass_3 <- attr_history_merged %>%
  select(time, prepClass) %>%
  group_by(time, prepClass) %>%
  summarise(value = n()) %>%
  filter(prepClass == 3)

a <- a %>%
  mutate(time = 1:728) %>%
  left_join(prepClass_4, by = "time")

a <- a %>%
  select(sim1, time, value) %>%
  filter(time >= 400)

## Explore number of people switching from daily to EDP and EDP to daily

sim$epi$prep.switch
sim$epi$edp.switch

y1 <- as.vector(sim$epi$prep.switch$sim1)
y2 <- as.vector(sim$epi$edp.switch$sim1)

plot(x, y1, type = "l",
     xlab = "Day", ylab = "PrEP Users",
     xlim = c(400, 728),
     main = "Number of PrEP Users Switching Regimens")
lines(x, y2, type = "l", col = "blue")
legend("topright", legend = c("Daily to EDP", "EDP to Daily"),
       col = c("black", "blue"), lty = 1)


