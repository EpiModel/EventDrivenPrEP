##
## 12. Epidemic Model Parameter Calibration, Processing
##

# Required variables:
#   - `ncores`
#   - `nsteps`
if (interactive()) {
  ncores <- 2
  nsteps <- 52
}

# Setup ------------------------------------------------------------------------
suppressMessages({
  library("EpiModelHIV")
  library("future.apply")
})

future::plan(future::multicore, workers = ncores)
calib_dir <- "data/output/calib"

# Process each file in parallel ------------------------------------------------
calib_files <- list.files(
  calib_dir,
  pattern = "^sim__.*rds$",
  full.names = TRUE
)

source("R/utils-targets.R")
assessments <- future.apply::future_lapply(
  calib_files,
  process_one_calibration, # in R/utils-targets.R
  nsteps = nsteps
)

# Merge all and combine --------------------------------------------------------
assessments <- bind_rows(assessments)

assessments <- assessments %>%
  select(- c(sim, batch)) %>%
  group_by(scenario_name) %>%
  summarise(across(
    everything(),
    list(
      q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
      q2 = ~ quantile(.x, 0.50, na.rm = TRUE),
      q3 = ~ quantile(.x, 0.75, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"
  ))
  # this last bloc calculate the q1, median and q3 for all of the variables

# Save the result --------------------------------------------------------------
# small rds to be downloaded and assessed localy
saveRDS(assessments, paste0(calib_dir, "/assessments.rds"))
