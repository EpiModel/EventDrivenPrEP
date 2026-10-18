# test processing

d_baseline <- readRDS(here("test_tibble", "df__baseline.rds"))
d_sc <- readRDS(sc_info$file_path)
d_sc <- mutate(d_sc, scenario_name = sc_info$scenario_name) |>
  select(scenario_name, batch_number, sim, time, everything())

glimpse(d_sc)

d_sc <- d_sc |>
  mutate(
    prev = i.num / (i.num + s.num),
    incid.sti = incid.gc + incid.ct
  )

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

glimpse(d_sc_ly)

d_sc_cml <- d_sc |>
  filter(time > max(time) - 10 * 364) |>
  group_by(scenario_name, batch_number, sim) |>
  summarise(
    across(
      starts_with("incid."),
      ~ sum(.x, na.rm = TRUE),
      .names = "{.col}_cml"
    ),
    .groups = "drop" # ungroup the tibble after the summary
  )

glimpse(d_sc_cml)

d_cmb <- left_join(
  d_sc_ly, d_sc_cml,
  by = c("scenario_name", "batch_number", "sim")
)

glimpse(d_cmb)
