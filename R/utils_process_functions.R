# utils_process_functions.R

# Functions -------------------------------------------------------------------
# Create processing function
process_one_scenario_tibble <- function(sc_info) {
  # loading the file
  d_sc <- readRDS(sc_info$file_path)
  d_sc <- mutate(d_sc, scenario_name = sc_info  $scenario_name) |>
    select(scenario_name, batch_number, sim, time, everything())

  # global mutate
  d_sc <- d_sc |>
    mutate(
      prev = i.num / (i.num + s.num),
      incid.sti = incid.gc + incid.ct,
      edp.prepCurr = coalesce(edp.prepCurr, 0),
      do.proportion = prepCurr/prepElig,
      edp.proportion = edp.prepCurr/edp.prepElig,
      total.users = prepCurr + edp.prepCurr,
      do.covered.sex = sex.do.hi/sex.do,
      do.covered.discordant = sex.do.disc.hi/sex.do.disc,
      edp.covered.sex = sex.edp.4/sex.edp,
      edp.covered.discordant = sex.edp.disc.4/sex.edp.disc,
      do.prep = prepCurr/(prepCurr + edp.prepCurr),
      edp.prep = edp.prepCurr/(prepCurr + edp.prepCurr)
    )

  d_sc_baseline <- d_sc |>
    filter(grepl("baseline", scenario_name)) |>
    mutate(total.proportion = prepCurr/prepElig)

  d_sc_2_3 <- d_sc |>
    filter(grepl("scenario_2|scenario_3", scenario_name)) |>
    mutate(total.proportion = (prepCurr + edp.prepCurr)/(prepElig + edp.prepElig))

  d_sc_1_4 <- d_sc |>
    filter(!grepl("scenario_2|scenario_3|baseline", scenario_name)) |>
    mutate(total.proportion = (prepCurr + edp.prepCurr)/prepElig)

  d_sc <- rbind(d_sc_baseline, d_sc_2_3, d_sc_1_4)

  # last year summaries
  d_sc_ly <- d_sc |>
    filter(time > max(time) - 364) |>
    group_by(scenario_name, batch_number, sim) |>
    summarise(
      across(
        c(prev, disease.mr100, do.proportion, do.covered.sex,
          do.covered.discordant, prepElig, edp.prepElig,
          edp.covered.sex, edp.covered.discordant, edp.proportion,
          edp.prepCurr, prepCurr, do.prep, edp.prep, total.users, total.proportion),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}_ly"
      ),
      across(
        starts_with("ir100"),
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

  # PrEP coverage on last day of simulation
  d_sc_end <- d_sc |>
    filter(time == max(time)) |>
    select(do.prep, edp.prep, total.users,
           do.proportion, edp.proportion, prepElig, edp.prepElig) |>
    rename_with(~ paste("end", .x, sep = "_"))

  # joining
  d_cmb <- left_join(
    d_sc_ly, d_sc_cml, d_sc_end,
    by = c("scenario_name", "batch_number", "sim")
  )

  return(d_cmb)
}

# Function for calculating PIA and NNT
pia_nnt_calc <- function(d_sc_raw, no_scenarios) {
  # Calculate PIA
  incid.base.tbl <- d_sc_raw |>
    filter(scenario_name == "baseline")
  incid.base <- as.vector(incid.base.tbl$incid_cml)

  d_sc_raw$incid.base <- rep(incid.base, no_scenarios)
  d_sc_raw$nia <- d_sc_raw$incid.base - d_sc_raw$incid_cml
  d_sc_raw$pia <- d_sc_raw$nia / d_sc_raw$incid.base

  # Calculate NNT & Percent pills difference from baseline
  d_sc_raw <- d_sc_raw |>
    mutate(hi.pills = prepCurr.hi_cml/7*5.5,
           med.pills = prepCurr.med_cml/7*2.5,
           low.pills = prepCurr.low_cml/7) |>
    mutate(daily.pills = hi.pills + med.pills + low.pills) |>
    mutate(total.pills = edp.pills + daily.pills) |>
    mutate(nnt = total.pills/nia)

  # Percent pills difference
  total.pills.tbl <- d_sc_raw |>
    filter(scenario_name == "baseline")
  total.pills.base <- as.vector(total.pills.tbl$total.pills)

  d_sc_raw$total.pills.base <- rep(total.pills.base, no_scenarios)
  d_sc_raw <- d_sc_raw |>
    mutate(pills.perc.diff = (total.pills - total.pills.base)/total.pills.base)

  return(d_sc_raw)
}
