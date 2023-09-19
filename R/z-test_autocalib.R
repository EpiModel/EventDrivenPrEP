library(dplyr)
library(tidyr)
library(ggplot2)

results <- readRDS("./results.rds")

results |>
  group_by(.iteration) |>
  summarize(across(starts_with("prep.start.prob_2"), list(
        min = ~min(.x),
        max = ~max(.x)
        ))) |> print(n = 100)

results |>
  ggplot(aes(x = prep.start.prob_2, y = cc.prep.H)) +
  geom_smooth() +
  geom_point() +
  geom_hline(yintercept = 0.229)

# targets_val = c(0.33, 0.127, 0.09),

results |>
  filter(
    abs(cc.prep.B - 0.199) < 0.02,
    abs(cc.prep.H - 0.229) < 0.02,
    abs(cc.prep.W - 0.321) < 0.02
  )

d1 <- results |>
  filter(
    abs(cc.prep.B - 0.199) < 0.01
  ) |>
  pull(prep.start.prob_1)
length(d1)
summary(d1)

d2 <- results |>
  filter(
    abs(cc.prep.H - 0.229) < 0.01
  ) |>
  pull(prep.start.prob_2)
length(d2)
summary(d2)

d3 <- results |>
  filter(
    abs(cc.prep.H - 0.331) < 0.01
  ) |>
  pull(prep.start.prob_3)
length(d3)
summary(d3)


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
