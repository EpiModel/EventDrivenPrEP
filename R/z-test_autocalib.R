library(dplyr)
library(tidyr)

results <- readRDS("./results.rds")

results |>
  group_by(.iteration) |>
  summarize(across(starts_with("prep.start.prob_2"), list(
        min = ~min(.x),
        max = ~max(.x)
        ))) |> print(n = 100)

# targets_val = c(0.33, 0.127, 0.09),

results |>
  filter(
    abs(cc.prep.B - 0.33) < 0.02,
    abs(cc.prep.H - 0.127) < 0.02,
    abs(cc.prep.W - 0.09) < 0.02
  ) |>
  select(starts_with("hiv.trans.scale")) |>
  summary()

results |>
  select(starts_with("aids.off.tx")) |> summary()

results |>
  group_by(.iteration) |>
  filter(
    abs(i.prev.dx.B - 0.33) < 0.02,
    abs(i.prev.dx.H - 0.127) < 0.02,
    abs(i.prev.dx.W - 0.09) < 0.02
  ) |>
  tally() |> print(n = 100)
