# Analysis for EDP manuscript #################################################

# Settings --------------------------------------------------------------------
source("./R/utils-0_project_settings.R")
library(here)
library(EpiModelHPC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gt)
library(writexl)

# Process Data ----------------------------------------------------------------
## Merge files as tibbles
merge_netsim_scenarios_tibble(
  "data/intermediate/scenarios",
  "test_tibble", # saves dataframes in folder called "test_tibble"
  3640 # only includes the last 3640 time steps / 10 years
)

## Read in files
read_scenarios <- function(file_name) {
  out <- readRDS(here("test_tibble", paste0("df__", file_name, ".rds")))
  return(out)
}

df_baseline <- read_scenarios("baseline")
df_1 <- read_scenarios("scenario_1")
df_2 <- read_scenarios("scenario_2")
df_3 <- read_scenarios("scenario_3")
df_4 <- read_scenarios("scenario_4")
df_m10 <- read_scenarios("scenario_minus10")
df_m20 <- read_scenarios("scenario_minus20")
df_m30 <- read_scenarios("scenario_minus30")
df_m40 <- read_scenarios("scenario_minus40")
df_m50 <- read_scenarios("scenario_minus50")
df_m60 <- read_scenarios("scenario_minus60")
df_m70 <- read_scenarios("scenario_minus70")
df_m80 <- read_scenarios("scenario_minus80")
df_m90 <- read_scenarios("scenario_minus90")
df_p10 <- read_scenarios("scenario_plus10")
df_p20 <- read_scenarios("scenario_plus20")
df_p30 <- read_scenarios("scenario_plus30")
df_p40 <- read_scenarios("scenario_plus40")
df_p50 <- read_scenarios("scenario_plus50")
df_p60 <- read_scenarios("scenario_plus60")
df_p70 <- read_scenarios("scenario_plus70")
df_p80 <- read_scenarios("scenario_plus80")
df_p90 <- read_scenarios("scenario_plus90")
df_p100 <- read_scenarios("scenario_plus100")
outcomes <- readRDS(here("data", "intermediate", "scenarios", "outcomes.rds"))

# Functions -------------------------------------------------------------------
last_year <- intervention_end - 364

outcome_fun <- function(simulation) {
  simulation <- simulation %>%
    mutate(edp.prepCurr = coalesce(edp.prepCurr, 0),
           do.proportion = prepCurr/prepElig,
           edp.proportion = edp.prepCurr/edp.prepElig,
           total.prep = prepCurr + edp.prepCurr,
           do.covered.sex = sex.do.hi/sex.do,
           do.covered.discordant = sex.do.disc.hi/sex.do.disc,
           edp.covered.sex = sex.edp.4/sex.edp,
           edp.covered.discordant = sex.edp.disc.4/sex.edp.disc,
           do.prep = prepCurr/(prepCurr + edp.prepCurr),
           edp.prep = edp.prepCurr/(prepCurr + edp.prepCurr))

  df <- simulation %>%
    filter(time %in% last_year:intervention_end) %>%
    group_by(batch_number, sim) %>%
    summarise(prevalence = mean(i.prev),
              incidence = mean(ir100),
              do.prop = mean(do.proportion),
              do.covered = mean(do.covered.sex),
              do.covered.disc = mean(do.covered.discordant),
              prepEligible = mean(prepElig),
              edp.prepEligible = mean(edp.prepElig),
              edp.covered.avg = mean(edp.covered.sex),
              edp.covered.disc = mean(edp.covered.discordant),
              edp.prop = mean(edp.proportion),
              prepCurr.avg = mean(prepCurr)
              ) %>%
    ungroup() %>%

              # HIV prevalence
    summarise(prev = round(median(prevalence)*100, 2),
              prev.ll = round(quantile(prevalence, 0.025, names = FALSE)*100, 2),
              prev.ul = round(quantile(prevalence, 0.975, names = FALSE)*100, 2),

              # HIV incidence
              ir100.med = round(median(incidence), 2),
              ir100.ll = round(quantile(incidence, 0.025, names = FALSE), 2),
              ir100.ul = round(quantile(incidence, 0.975, names = FALSE), 2),

              # Proportion on daily oral PrEP
              do.med = round(median(do.prop)*100, 2),
              do.ll = round(quantile(do.prop*100, 0.025, names = FALSE), 2),
              do.ul = round(quantile(do.prop*100, 0.975, names = FALSE), 2),

              # Proportion of sex acts covered by DO PrEP
              do.sex.med = round(median(do.covered)*100, 2),
              do.sex.ll = round(quantile(do.covered, 0.025, names = FALSE)*100, 2),
              do.sex.ul = round(quantile(do.covered, 0.975, names = FALSE)*100, 2),

              # Proportion of discordant sex acts covered by DO PrEP
              do.disc.med = round(median(do.covered.disc)*100, 2),
              do.disc.ll = round(quantile(do.covered.disc, 0.025, names = FALSE)*100, 2),
              do.disc.ul = round(quantile(do.covered.disc, 0.975, names = FALSE)*100, 2),

              # Number eligible for DO PrEP
              do.elig = round(median(prepEligible)),
              do.elig.ll = round(quantile(prepEligible, 0.025, names = FALSE)),
              do.elig.ul = round(quantile(prepEligible, 0.975, names = FALSE)),

              # Proportion of those eligible on EDP
              edp.med = round(median(edp.prop)*100, 2),
              edp.ll = round(quantile(edp.prop, 0.025, na.rm = TRUE, names = FALSE)*100, 2),
              edp.ul = round(quantile(edp.prop, 0.975, na.rm = TRUE, names = FALSE)*100, 2),

              # Number eligible for EDP
              edp.elig = round(median(edp.prepEligible)),
              edp.elig.ll = round(quantile(edp.prepEligible, 0.025, na.rm = TRUE, names = FALSE)),
              edp.elig.ul = round(quantile(edp.prepEligible, 0.975, na.rm = TRUE, names = FALSE)),

              # Proportion of sex acts covered by EDP
              edp.sex.med = round(median(edp.covered.avg)*100, 2),
              edp.sex.ll = round(quantile(edp.covered.avg, 0.025, na.rm = TRUE, names = FALSE)*100, 2),
              edp.sex.ul = round(quantile(edp.covered.avg, 0.975, na.rm = TRUE, names = FALSE)*100, 2),

              # Proportion of discordant sex acts covered by EDP
              edp.disc.med = round(median(edp.covered.disc), 2),
              edp.disc.ll = round(quantile(edp.covered.disc, 0.025, na.rm = TRUE, names = FALSE), 2),
              edp.disc.ul = round(quantile(edp.covered.disc, 0.975, na.rm = TRUE, names = FALSE), 2))

  # Total PrEP users by end of simulation
  df1 <- simulation %>%
    filter(time == intervention_end) %>%

    summarise(# Proportion of PrEP users on DO PrEP
              do.prep.med = round(median(do.prep)*100, 2),
              do.prep.ll = round(quantile(do.prep, 0.025, names = FALSE)*100, 2),
              do.prep.ul = round(quantile(do.prep, 0.975, names = FALSE)*100, 2),

              # Proportion of PrEP users on EDP
              edp.prep.med = round(median(edp.prep)*100, 2),
              edp.prep.ll = round(quantile(edp.prep, 0.025, names = FALSE)*100, 2),
              edp.prep.ul = round(quantile(edp.prep, 0.975, names = FALSE)*100, 2),

              # Total PrEP users
              total.prep.med = round(median(total.prep)),
              total.prep.ll = round(quantile(total.prep, 0.025, names = FALSE)),
              total.prep.ul = round(quantile(total.prep, 0.975, names = FALSE)))

  # Percent infections averted
  base <- df_baseline %>%
    group_by(batch_number, sim) %>%
    summarise(incidence = sum(incid, na.rm = TRUE))

  incid.base <- as.vector(base$incidence)

  sim.stats <- simulation %>%
    group_by(batch_number, sim) %>%
    summarise(incidence = sum(incid, na.rm = TRUE),

              edp.pills = sum(pills, na.rm = TRUE),
              prepCurr.hi = sum(prepCurr.hi, na.rm = T),
              prepCurr.med = sum(prepCurr.med, na.rm= T),
              prepCurr.low = sum(prepCurr.low, na.rm = T))

  incid.comp <- as.vector(sim.stats$incidence)

  vec.nia <- incid.base - incid.comp

  vec.pia <- vec.nia/incid.base
  out.pia <- round(data.frame(pia.median = median(vec.pia)*100,
                              pia.ll = quantile(vec.pia, 0.025, names = FALSE)*100,
                              pia.ul = quantile(vec.pia, 0.975, names = FALSE)*100), 2)

  # Number of pills needed to prevent 1 infection
  edp.pills <- mean(sim.stats$edp.pills)

  hi <- mean(sim.stats$prepCurr.hi)/7*5.5
  med <- mean(sim.stats$prepCurr.med)/7*2.5
  low <- mean(sim.stats$prepCurr.low)/7

  daily.pills <- floor(hi+med+low)
  daily.pills

  total.pills <- edp.pills + daily.pills
  total.pills

  vec.nnt <- total.pills/vec.nia
  nnt.out <- round(data.frame(nnt.median = median(vec.nnt),
                              nnt.ll = quantile(vec.nnt, 0.025, names = FALSE),
                              nnt.ul = quantile(vec.nnt, 0.975, names = FALSE)))

  return(list(df1 = df1,
              df = df,
              pia = out.pia,
              pills = data.frame(total.pills, edp.pills, daily.pills),
              nnt = nnt.out))
}

days_fun <- function(simulation) {
  simulation %>%
    mutate(day = rep(0:3640, 120))
}

plotdf_fun <- function(simulation, counterfactual) {
  simulation %>%
    group_by(day) %>%
    summarise(ir100.mean = mean(ir100),
              ir100.ll = quantile(ir100, 0.025),
              ir100.ul = quantile(ir100, 0.975)) %>%
    mutate(scenario = rep(counterfactual, 3641))
}

# Analysis---------------------------------------------------------------------

report_fun <- function(outcome_list) {
  df <- cbind(outcome_list$df1,
              outcome_list$df,
              outcome_list$pia,
              outcome_list$pills,
              outcome_list$nnt)

  return(df)
}

# List of scenario dataframes
df_list <- list(df_m90 = df_m90, df_m80 = df_m80, df_m70 = df_m70,
                df_m60 = df_m60, df_m50 = df_m50, df_m40 = df_m40,
                df_m30 = df_m30, df_m20 = df_m20,
                df_m10 = df_m10, df_1 = df_1, df_p10 = df_p10, df_p20 = df_p20,
                df_p30 = df_p30, df_p40 = df_p40, df_p50 = df_p50,
                df_p60 = df_p60, df_p70 = df_p70, df_p80 = df_p80,
                df_p90 = df_p90, df_p100 = df_p100,
                df_2 = df_2, df_3 = df_3, df_4 = df_4)

# Create df for Table 1
baseline <- outcome_fun(df_baseline)
df <- report_fun(baseline)

for (i in df_list) {
  outcomes_list <- outcome_fun(i)
  report_output <- report_fun(outcomes_list)

  df <- rbind(df, report_output)
}

df$scenarios <- c("Baseline", names(df_list))
df <- df %>% relocate(scenarios)

# Save df as excel file
write_xlsx(df, "raw_output.xlsx")

# Publication-ready table
df |>
  gt(rowname_col = "scenarios") |>
  tab_header(
    title = "Table 1",
  )
# Plots -----------------------------------------------------------------------

# Create dataframe for plot
df.list <- list(df_baseline, df_1, df_2, df_3, df_4, df_5, df_6, df_7,
                df_8, df_9, df_11, df_12, df_13, df_14, df_15)

df_baseline <- days_fun(df_baseline)
df_1 <- days_fun(df_1)
df_2 <- days_fun(df_2)
df_3 <- days_fun(df_3)
df_4 <- days_fun(df_4)
df_5 <- days_fun(df_5)
df_6 <- days_fun(df_6)
df_7 <- days_fun(df_7)
df_8 <- days_fun(df_8)
df_9 <- days_fun(df_9)
df_11 <- days_fun(df_11)
df_12 <- days_fun(df_12)
df_13 <- days_fun(df_13)
df_14 <- days_fun(df_14)
df_15 <- days_fun(df_15)

test <- df_1 %>%
  filter(day == 1)

df_baselinep <- plotdf_fun(df_baseline, "Baseline")
df_1p <- plotdf_fun(df_1, "Scenario 1")
df_2p <- plotdf_fun(df_2, "Scenario 2")
df_3p <- plotdf_fun(df_3, "Scenario 3")
df_4p <- plotdf_fun(df_4, "Scenario 4")
df_5p <- plotdf_fun(df_5, "-50%")
df_6p <- plotdf_fun(df_6, "-40%")
df_7p <- plotdf_fun(df_7, "-30%")
df_8p <- plotdf_fun(df_8, "-20%")
df_9p <- plotdf_fun(df_9, "-10%")
df_11p <- plotdf_fun(df_11, "+10%")
df_12p <- plotdf_fun(df_12, "+20%")
df_13p <- plotdf_fun(df_13, "+30%")
df_14p <- plotdf_fun(df_14, "+40%")
df_15p <- plotdf_fun(df_15, "+50%")

df_plot1 <- bind_rows(df_baselinep, df_5p, df_6p, df_7p,
                     df_8p, df_9p, df_11p, df_12p, df_13p, df_14p, df_15p)

df_plot2 <- bind_rows(df_baselinep, df_1p, df_2p, df_3p, df_4p)

## ggplot

ggplot(df_plot1, aes(x = day, y = ir100.mean, color = scenario)) +
  geom_line() +
  theme_minimal()

ggplot(df_plot2, aes(x = day, y = ir100.mean, color = scenario)) +
  geom_line() +
  theme_minimal()



