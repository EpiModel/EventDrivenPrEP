# Analysis for EDP manuscript #################################################

# Settings --------------------------------------------------------------------
source("./R/utils-0_project_settings.R")
source("./R/utils_process_functions.R")
source("./R/utils-format.R")
library(EpiModelHPC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gt)
library(writexl)
library(zoo)
library(metR)
library(viridis)
library(stringr)
library(patchwork)
library(scales)

# Table 3 ---------------------------------------------------------------------

# Process Data
## Merge files as tibbles
#merge_netsim_scenarios_tibble(
#  "data/intermediate/table3",
#  "data/output/table3", # saves dataframes in folder called "test_tibble"
#  3640 # only includes the last 3640 time steps / 10 years
#)

## Read in files
sc_dir <- "data/output/table3"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)

## Process each scenario
d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

# Create df and calculate PIA and NNT
d_sc_raw <- bind_rows(d_ls)
d_sc_raw <- pia_nnt_calc(d_sc_raw, no_scenarios = 19)

table3 <- format_table(d_sc_raw, var_labels, format_patterns)
readr::write_csv(table3, "table3.csv")

# Table 3 Sensitivity Analysis ---------------------------------------------------------------------

# Process Data
## Merge files as tibbles
#merge_netsim_scenarios_tibble(
#  "data/intermediate/table3_sens",
#  "data/output/table3_sens", # saves dataframes in folder called "table3_sens"
#  3640 # only includes the last 3640 time steps / 10 years
#)

## Read in files
sc_dir <- "data/output/table3_sens"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)
sc_infos_tbl_baseline <- EpiModelHPC::get_scenarios_tibble_infos("data/output/table3")
sc_infos_tbl <- rbind(sc_infos_tbl_baseline[1, ], sc_infos_tbl[-1,])

## Process each scenario
d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

# Create df and calculate PIA and NNT
d_sc_raw <- bind_rows(d_ls)
d_sc_raw <- pia_nnt_calc(d_sc_raw, no_scenarios = 6)

table3_sens <- format_table(d_sc_raw, var_labels, format_patterns)
readr::write_csv(table3_sens, "table3_sens.csv")

## Explore HIV incidence by PrEP regimen
d_sc_raw2 <- d_sc_raw |>
  mutate(do.incid.user.ratio = incid.daily_ly/prepCurr_ly,
         edp.incid.user.ratio = incid.edp_ly/edp.prepCurr_ly)

table3_sens <- format_table(d_sc_raw2, var_labels, format_patterns)


# Table 4 ---------------------------------------------------------------------

# Process Data
## Merge files as tibbles
merge_netsim_scenarios_tibble(
  "data/intermediate/table4",
  "data/output/table4_output",
  3640 # only includes the last 3640 time steps / 10 years
)

## Read in files
sc_dir <- "data/output/table4"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)
sc_infos_tbl_baseline <- EpiModelHPC::get_scenarios_tibble_infos("data/output/table3")
sc_infos_tbl <- rbind(sc_infos_tbl_baseline[1, ], sc_infos_tbl[-1,])

## Read in files
sc_dir <- "data/output/table4"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)

## Process each scenario
d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

## Create df and calculate PIA and NNT
d_sc_raw <- bind_rows(d_ls)
d_sc_raw <- pia_nnt_calc(d_sc_raw, 49)

table4 <- format_table(d_sc_raw, var_labels, format_patterns)
table4 <- table4[order(as.integer(str_sub(table4$scenario_name, -1, -1))), ]
table4 <- table4 |>
  slice(49, 1:48)

readr::write_csv(table4, "table4.csv")

# Contour plot 1 --------------------------------------------------------------

## Process Data
### Merge files as tibbles
#merge_netsim_scenarios_tibble(
#  "data/intermediate/contourplots1",
#  "cp1_tibble", # saves dataframes in folder called "adhr_sens_tibble"
#  3640 # only includes the last 3640 time steps / 10 years
#)

sc_dir <- "data/output/figure2"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)
sc_infos_tbl_baseline <- EpiModelHPC::get_scenarios_tibble_infos("data/output/table3")
sc_infos_tbl <- rbind(sc_infos_tbl_baseline[1, ], sc_infos_tbl)

sc_info <- sc_infos_tbl[1, ]

d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

## Create df and calculate PIA and NNT
d_sc_raw <- bind_rows(d_ls)

new_rows <- data.frame(
  scenario_name = rep("edp_startprob_1.4_hi_adhr_0.6", 30)
)
d_sc_raw <- bind_rows(d_sc_raw, new_rows)
d_sc_raw <- sort(d_sc_raw)

d_sc_raw <- pia_nnt_calc(d_sc_raw, 133)

d_sc_raw <- d_sc_raw |>
  filter(!is.na(batch_number))

readr::write_csv(d_sc_raw, "figure2_raw.csv")

## Create table output
figure2 <- format_table(d_sc_raw, var_labels, format_patterns)
readr::write_csv(figure2, "figure2.csv")

## Use table output for editing
d_sc_raw <- read.csv("figure2.csv")

## Dfs for loess model
out.pia <- d_sc_raw |>
  filter(scenario_name != "0_no_edp") |>
  mutate(edp.prep.start.prob = as.numeric(str_match(scenario_name, "edp_startprob_(.*?)_hi_adhr")[,2]),
         hi.adhr.prob = as.numeric(str_replace(scenario_name, '(.*?)hi_adhr_(.*?)', ''))) |>
  select(!scenario_name) |>
  select(c(edp.prep.start.prob, hi.adhr.prob, pia))

summary(out.pia)

out.nnt <- d_sc_raw |>
  filter(scenario_name != "0_no_edp") |>
  mutate(edp.prep.start.prob = as.numeric(str_match(scenario_name, "edp_startprob_(.*?)_hi_adhr")[,2]),
         hi.adhr.prob = as.numeric(str_replace(scenario_name, '(.*?)hi_adhr_(.*?)', '')),
         nnt = ifelse(is.infinite(nnt), 0, nnt)) |>
  select(!scenario_name) |>
  select(c(edp.prep.start.prob, hi.adhr.prob, nnt))

test <- out.nnt |>
  group_by(edp.prep.start.prob, hi.adhr.prob) |>
  summarise(average_nnt = mean(nnt, na.rm = TRUE))

summary(out.nnt)

## Loess model
loess_pia <- loess(pia ~ edp.prep.start.prob * hi.adhr.prob, data = out.pia, span = 0.25)
fit_pia <- expand.grid(list(edp.prep.start.prob = seq(min(out.pia$edp.prep.start.prob), max(out.pia$edp.prep.start.prob), length.out = 132),
                            hi.adhr.prob = seq(min(out.pia$hi.adhr.prob), max(out.pia$hi.adhr.prob), length.out = 132)))
fit_pia$pia <- as.numeric(predict(loess_pia, newdata = fit_pia))

loess_nnt <- loess(nnt ~ edp.prep.start.prob * hi.adhr.prob, data = out.nnt, span = 0.25)
fit_nnt <- expand.grid(list(edp.prep.start.prob = seq(min(out.nnt$edp.prep.start.prob), max(out.nnt$edp.prep.start.prob), length.out = 132),
                            hi.adhr.prob = seq(min(out.nnt$hi.adhr.prob), max(out.nnt$hi.adhr.prob), length.out = 132)))
fit_nnt$nnt <- as.numeric(predict(loess_nnt, newdata = fit_nnt))

## Plots
pia <- ggplot(fit_pia, aes(edp.prep.start.prob, hi.adhr.prob)) +
  geom_raster(aes(fill = pia), interpolate = TRUE) +
  geom_contour(aes(z = pia), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Probability of Excellent EDP Adherence", x = "EDP Initiation Probability", fill = "PIA") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1,
                     labels = scales::label_percent())

pia

fit_nnt <- fit_nnt |>
  filter(edp.prep.start.prob >= 1)

nnt <- ggplot(fit_nnt, aes(edp.prep.start.prob, hi.adhr.prob)) +
  geom_raster(aes(fill = nnt), interpolate = TRUE) +
  geom_contour(aes(z = nnt), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Probability of Excellent EDP Adherence", x = "EDP Initiation Probability", fill = "NPNT") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1,
                     labels = label_number())

nnt

nested <- (pia|nnt) +
  plot_annotation(tag_levels = 'A')
nested

png('figure2.png', width=3600, height=2048, res=300)
nested
dev.off()


# Contour plot 2 --------------------------------------------------------------

# Process Data
## Merge files as tibbles
#merge_netsim_scenarios_tibble(
#  "data/intermediate/contourplots2",
#  "cp2_tibble", # saves dataframes in folder called "adhr_sens_tibble"
#  3640 # only includes the last 3640 time steps / 10 years
#)

sc_dir <- "data/output/figure3"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)
sc_infos_tbl_baseline <- EpiModelHPC::get_scenarios_tibble_infos("data/output/table3")
sc_infos_tbl <- rbind(sc_infos_tbl_baseline[1, ], sc_infos_tbl)

d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

# Create df and calculate PIA and NNT
d_sc_raw <- bind_rows(d_ls)

d_sc_raw <- pia_nnt_calc(d_sc_raw, 101)
readr::write_csv(d_sc_raw, "figure3_raw.csv")

out.pia <- d_sc_raw |>
  filter(scenario_name != "0_no_edp") |>
  mutate(daily.switch = as.numeric(str_match(scenario_name, "daily_switch_(.*?)_edp_switch")[,2]),
         edp.switch = as.numeric(str_replace(scenario_name, '(.*?)edp_switch_(.*?)', ''))) |>
  select(!scenario_name) |>
  select(c(daily.switch, edp.switch, pia))

summary(out.pia)

out.nnt <- d_sc_raw |>
  filter(scenario_name != "0_no_edp") |>
  mutate(daily.switch = as.numeric(str_match(scenario_name, "daily_switch_(.*?)_edp_switch")[,2]),
         edp.switch = as.numeric(str_replace(scenario_name, '(.*?)edp_switch_(.*?)', ''))) |>
  select(!scenario_name) |>
  select(c(daily.switch, edp.switch, nnt))

summary(out.nnt)

# Loess model
loess_pia <- loess(pia ~ daily.switch * edp.switch, data = out.pia, span = 0.25)
fit_pia <- expand.grid(list(daily.switch = seq(min(out.pia$daily.switch), max(out.pia$daily.switch), length.out = 100),
                            edp.switch = seq(min(out.pia$edp.switch), max(out.pia$edp.switch), length.out = 100)))
fit_pia$pia <- as.numeric(predict(loess_pia, newdata = fit_pia))

loess_nnt <- loess(nnt ~ daily.switch * edp.switch, data = out.nnt, span = 0.25)
fit_nnt <- expand.grid(list(daily.switch = seq(min(out.nnt$daily.switch), max(out.nnt$daily.switch), length.out = 100),
                            edp.switch = seq(min(out.nnt$edp.switch), max(out.nnt$edp.switch), length.out = 100)))
fit_nnt$nnt <- as.numeric(predict(loess_nnt, newdata = fit_nnt))

pia <- ggplot(fit_pia, aes(daily.switch, edp.switch)) +
  geom_raster(aes(fill = pia), interpolate = TRUE) +
  geom_contour(aes(z = pia), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Event-Driven PrEP to Daily PrEP Switch (Probability)",
       x = "Daily PrEP to Event-Driven PrEP Switch (Probability)",
       fill = "PIA") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1,
                     labels = scales::label_percent())

nnt <- ggplot(fit_nnt, aes(daily.switch, edp.switch)) +
  geom_raster(aes(fill = nnt), interpolate = TRUE) +
  geom_contour(aes(z = nnt), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Event-Driven PrEP to Daily PrEP Switch (Probability)",
       x = "Daily PrEP to Event-Driven PrEP Switch (Probability)",
       fill = "NPNT") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1)

nested <- (pia|nnt) +
  plot_annotation(tag_levels = 'A')
nested

png('figure3.png', width=3600, height=2048, res=300)
nested
dev.off()

# Test --------------------------------------------------------------

# Process Data
## Merge files as tibbles
merge_netsim_scenarios_tibble(
  "data/intermediate/test6",
  "data/output/test_tibble6", # saves dataframes in folder called "test_tibble"
  3640 # only includes the last 3640 time steps / 10 years
)

## Read in files
sc_dir <- "data/output/test_tibble2"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)

## Process each scenario
d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

# Create df and calculate PIA and NNT
d_sc_raw <- bind_rows(d_ls)
d_sc_raw <- pia_nnt_calc(d_sc_raw, no_scenarios = 4)

table_test2 <- format_table(d_sc_raw, var_labels, format_patterns)
readr::write_csv(test5, "test.csv")
