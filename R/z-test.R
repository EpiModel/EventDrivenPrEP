# Scratchpad for interactive testing before integration in a script

library(EpiModelHIV)
library(dplyr)
library(tidyr)
source("./R/utils-targets.R")

d <- readRDS("./data/intermediate/calibration/merged_tibbles/df__empty_scenario.rds")

glimpse(d)

d_tar <- mutate_calibration_targets(d) |>
  filter(time >= max(time) - (52 * 7)) %>%
  select(sim, any_of(names(targets))) %>%
  group_by(sim) %>%
  summarise(across(
    everything(),
    ~ mean(.x, na.rm = TRUE)
  )) %>%
  ungroup() |>
  select(- c(sim)) %>%
  summarise(across(
    everything(),
    list(
      q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
      q2 = ~ quantile(.x, 0.50, na.rm = TRUE),
      q3 = ~ quantile(.x, 0.75, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"
  ))

as.list(d_tar)

ckpt <- readRDS("./1.rds")
