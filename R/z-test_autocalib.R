library(dplyr)
library(tidyr)

results <- readRDS("./results.rds")

results |>
  group_by(.iteration) |>
  summarize(across(starts_with("hiv.trans.scale_3"), list(
        min = ~min(.x),
        max = ~max(.x)
        )))
