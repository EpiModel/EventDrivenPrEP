## Contour plots
library(dplyr)
library(here)
library(ggplot2)
library(zoo)
library(metR)
library(viridis)
library(stringr)
library(patchwork)

sc_dir <- "cp2_tibble"
sc_infos_tbl <- EpiModelHPC::get_scenarios_tibble_infos(sc_dir)
sc_infos_tbl_baseline <- EpiModelHPC::get_scenarios_tibble_infos("test_tibble")
sc_infos_tbl <- rbind(sc_infos_tbl, sc_infos_tbl_baseline[1, ])

sc_info <- sc_infos_tbl[1, ]

# Create processing function
process_one_scenario_tibble <- function(sc_info) {
  # loading the file
  d_sc <- readRDS(sc_info$file_path)
  d_sc <- mutate(d_sc, scenario_name = sc_info$scenario_name) |>
    select(scenario_name, batch_number, sim, time, everything())

  # global mutate
  d_sc <- d_sc |>
    mutate(
      prev = i.num / (i.num + s.num),
      incid.sti = incid.gc + incid.ct
    )

  # last year summaries
  d_sc_ly <- d_sc |>
    filter(time > max(time) - 364) |>
    group_by(scenario_name, batch_number, sim) |>
    summarise(
      across(
        c(prev, disease.mr100),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}_ly"
      ),
      .groups = "drop" # ungroup the tibble after the summary
    )

  # cummulative summaries
  d_sc_cml <- d_sc |>
    filter(time > max(time) - 10 * 364) |>
    group_by(scenario_name, batch_number, sim) |>
    summarise(
      across(
        starts_with("incid"),
        ~ sum(.x, na.rm = TRUE),
        .names = "{.col}_cml"
      ),
      across(
        starts_with("prepCurr."),
        ~ sum(.x, na.rm = TRUE),
        .names = "{.col}_cml"
      ),
      edp.pills = sum(pills, na.rm = TRUE),
      .groups = "drop" # ungroup the tibble after the summary
    )

  # joining
  d_cmb <- left_join(
    d_sc_ly, d_sc_cml,
    by = c("scenario_name", "batch_number", "sim")
  )

  return(d_cmb)
}

d_ls <- lapply(
  seq_len(nrow(sc_infos_tbl)),
  \(i) process_one_scenario_tibble(sc_infos_tbl[i, ])
)

d_sc_raw <- bind_rows(d_ls)
readr::write_csv(d_sc_raw, "sc_raw.csv")

d_sc_raw <- read.csv("sc_raw.csv")

# Calculate PIA
incid.base.tbl <- d_sc_raw |>
  filter(scenario_name == "baseline")
incid.base <- as.vector(incid.base.tbl$incid_cml)

d_sc_raw$incid.base <- rep(incid.base, 101)
d_sc_raw$nia <- d_sc_raw$incid.base - d_sc_raw$incid_cml
d_sc_raw$pia <- d_sc_raw$nia / d_sc_raw$incid.base

out.pia <- d_sc_raw |>
  filter(scenario_name != "baseline") |>
  #group_by(scenario_name) |>
  #summarise(pia.median = median(pia)) |>
  mutate(daily.switch = as.numeric(str_match(scenario_name, "daily_switch_(.*?)_edp_switch")[,2]),
         edp.switch = as.numeric(str_replace(scenario_name, '(.*?)edp_switch_(.*?)', ''))) |>
  select(!scenario_name) |>
  select(c(daily.switch, edp.switch, pia))

summary(out.pia)

# Calculate NNT
d_sc_raw <- d_sc_raw |>
  mutate(hi.pills = prepCurr.hi_cml/7*5.5,
         med.pills = prepCurr.med_cml/7*2.5,
         low.pills = prepCurr.low_cml/7) |>
  mutate(daily.pills = hi.pills + med.pills + low.pills) |>
  mutate(total.pills = edp.pills + daily.pills) |>
  mutate(nnt = total.pills/nia)

out.nnt <- d_sc_raw |>
  filter(scenario_name != "baseline") |>
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
  labs(y = "Event-Driven PrEP to Daily PrEP Switch (Probability)", x = "Daily PrEP to Event-Driven PrEP Switch (Probability)", fill = "PIA") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1,
                     labels = scales::label_percent())

nnt <- ggplot(fit_nnt, aes(daily.switch, edp.switch)) +
  geom_raster(aes(fill = nnt), interpolate = TRUE) +
  geom_contour(aes(z = nnt), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Event-Driven PrEP to Daily PrEP Switch (Probability)", x = "Daily PrEP to Event-Driven PrEP Switch (Probability)", fill = "NNT") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1)

nested <- (pia|nnt) +
  plot_annotation(title = "120 simulations per scenario",
                  tag_levels = 'A')
nested

# sample first two batches only---------------------------------------------------------------------------
d_sc_raw_sample <- d_sc_raw |>
  filter(batch_number %in% 1:2)

d_sc_raw <- d_sc_raw_sample

# Calculate PIA
incid.base.tbl <- d_sc_raw |>
  filter(scenario_name == "baseline")
incid.base <- as.vector(incid.base.tbl$incid_cml)

d_sc_raw$incid.base <- rep(incid.base, 101)
d_sc_raw$nia <- d_sc_raw$incid.base - d_sc_raw$incid_cml
d_sc_raw$pia <- d_sc_raw$nia / d_sc_raw$incid.base

out.pia <- d_sc_raw |>
  filter(scenario_name != "baseline") |>
  #group_by(scenario_name) |>
  #summarise(pia.median = median(pia)) |>
  mutate(daily.switch = as.numeric(str_match(scenario_name, "daily_switch_(.*?)_edp_switch")[,2]),
         edp.switch = as.numeric(str_replace(scenario_name, '(.*?)edp_switch_(.*?)', ''))) |>
  select(!scenario_name) |>
  select(c(daily.switch, edp.switch, pia))

summary(out.pia)

# Calculate NNT
d_sc_raw <- d_sc_raw |>
  mutate(hi.pills = prepCurr.hi_cml/7*5.5,
         med.pills = prepCurr.med_cml/7*2.5,
         low.pills = prepCurr.low_cml/7) |>
  mutate(daily.pills = hi.pills + med.pills + low.pills) |>
  mutate(total.pills = edp.pills + daily.pills) |>
  mutate(nnt = total.pills/nia)

out.nnt <- d_sc_raw |>
  filter(scenario_name != "baseline") |>
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
  labs(y = "Event-Driven PrEP to Daily PrEP Switch (Probability)", x = "Daily PrEP to Event-Driven PrEP Switch (Probability)", fill = "PIA") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1,
                     labels = scales::label_percent())

nnt <- ggplot(fit_nnt, aes(daily.switch, edp.switch)) +
  geom_raster(aes(fill = nnt), interpolate = TRUE) +
  geom_contour(aes(z = nnt), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Event-Driven PrEP to Daily PrEP Switch (Probability)", x = "Daily PrEP to Event-Driven PrEP Switch (Probability)", fill = "NNT") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1)

nested_sample <- (pia|nnt) +
  plot_annotation(title = "Sample 60 simulations",
                  tag_levels = 'A')
nested_sample

