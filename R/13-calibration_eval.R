##
## 13. Epidemic Model Parameter Calibration, Local evaluation
##
#

# Setup ------------------------------------------------------------------------
suppressMessages({
  library("EpiModel")
  library("dplyr")
})

d <- readRDS("data/output/calib/assessments.rds")

glimpse(d)

