
##
## 04. Small-scale epidemic simulation for testing/debugging
##

## Libraries ------------------------------------------------------------------
#pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModelHIV-p")
library("EpiModelHIV")
library(dplyr)
library(ggplot2)

# Settings ---------------------------------------------------------------------
source("R/utils-0_project_settings.R")

# Necessary files
epistats <- readRDS("data/intermediate/estimates/epistats-local.rds")
netstats <- readRDS("data/intermediate/estimates/netstats-local.rds")
est      <- readRDS("data/intermediate/estimates/netest-local.rds")

param <- param.net(
  data.frame.params   = read.csv("data/input/params.csv"),
  netstats            = netstats,
  epistats            = epistats,
  prep.start          = 182,
  riskh.start         = 1
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
     main = paste("EDP Start Scenario =", sim$param$edp.start.scenario),
     sub = paste("Total MSM starting PrEP by day 728 =", cum.prepDailyStart[728] + cum.prepEDPStart[728]))
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

# Pill count tracker
summary(sim$epi$pills, na.rm = T)
sum(sim$epi$pills, na.rm = T)
plot(x, unlist(sim$epi$pills))

summary(sim$epi$pillspp)
mean(sim$epi$pillspp$sim1, na.rm = T)*7
plot(x, unlist(sim$epi$pillspp))

summary(sim$epi$pills7d)
plot(x, unlist(sim$epi$pills7d))

RcppRoll::roll_meanr(unlist(sim$epi$pills7d), 7)

# HIV incidence for daily PrEP
sum(sim$epi$incid$sim1, na.rm = T)
sum(sim$epi$incid.daily$sim1, na.rm = T)

# Pill counts for daily PrEP users
hi <- sum(sim$epi$prepCurr.hi$sim1, na.rm = T)/7*5.5
med <- sum(sim$epi$prepCurr.med$sim1, na.rm = T)/7*2.5
low <- sum(sim$epi$prepCurr.low$sim1, na.rm = T)/7

daily.pills <- floor(hi+med+low)
daily.pills

## Explore history of prepClass.edp attribute for EDP users

attr_history <- get_attr_history(sim)
attr_history

unique(attr_history$edp.prepClass$uids)
length(unique(attr_history$edp.prepClass$uids))

prop.table(table(attr_history$edp.prepClass$values))

attr_history_merged <- left_join(attr_history$edp.prepClass, attr_history$sex.edp, by = c("time", "uids")) %>%
  left_join(., attr_history$lastPrepCombo, by = c("time", "uids")) %>%
  left_join(., attr_history$timeLastSex, by = c("time", "uids")) %>%
  left_join(., attr_history$pillCount, by = c("time", "uids")) %>%
  left_join(., attr_history$pillCount7d, by = c("time", "uids")) %>%
  select(time, uids, values.x, values.y, values.x.x, values.y.y, values.x.x.x, values.y.y.y) %>%
  rename("prepClass" = "values.x",
         "sex" = "values.y",
         "lastPrepCombo" = "values.x.x",
         "timeLastSex" = "values.y.y",
         "pillCount" = "values.x.x.x",
         "pillCount7d" = "values.y.y.y") %>%
  mutate(timeSinceLastSex = time - timeLastSex)

id_a <- filter(attr_history_merged, uids == 4560)
id_b <- filter(attr_history_merged, uids == 6113)
id_c <- filter(attr_history_merged, uids == 7203)
id_d <- filter(attr_history_merged, uids == 314)
id_e <- filter(attr_history_merged, uids == 8530)

write.csv(id_a, "sample_id.csv")

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


