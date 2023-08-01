library(dplyr)
library(tidyr)

results <- readRDS("./results.rds")
retain_prop <- 0.3
thresholds <- rep(0.02, 3)
n_enough <- 50

job <- list()
job$params <- c("hiv.trans.scale_1", "hiv.trans.scale_2", "hiv.trans.scale_3")
job$targets <- c("i.prev.dx.B", "i.prev.dx.H", "i.prev.dx.W")
job$targets_val <- c(0.33, 0.127, 0.084)

new_ranges <- list()
for (i in seq_along(job$targets)) {
  params <- results[[ job$params[i] ]]
  values <- results[[ job$targets[i] ]]
  target <- job$targets_val[i]

  d <- dplyr::tibble(
    params = params,
    score = abs(values - target)
  )
  d <- dplyr::arrange(d, score)
  d <- head(d, ceiling(nrow(d) * retain_prop))
  new_ranges[[i]] <- range(d$params)
}

p_ok <- results[, c(job$params, job$targets)]
for (j in seq_along(job$targets)) {
  values <- p_ok[[ job$targets[j] ]]
  target <- job$targets_val[j]
  thresh <- thresholds[j]

  p_ok <- p_ok[abs(values - target) < thresh, ]
}

if (nrow(p_ok) > n_enough) {
  res <- p_ok[, job$params]
  # get the n_tuple where all values are the closest to the median
  best <- dplyr::summarise(res, dplyr::across(
      dplyr::everything(),
      ~ abs(.x - median(.x)))
  )
  best <- which.min(rowSums(best))
  res_best <- res[best, ]
}

res_best

res_fmt <- res %>%
  select(c(starts_with("i.prev.dx"), starts_with("hiv.trans.scale"))) %>%
  mutate(
    i.prev.dx.B = abs(i.prev.dx.B - 0.33),
    i.prev.dx.H = abs(i.prev.dx.H - 0.127),
    i.prev.dx.W = abs(i.prev.dx.W - 0.084),
    score = i.prev.dx.B + i.prev.dx.H + i.prev.dx.W
  ) %>%
  arrange(score)

p_ok <- res_fmt %>%
  filter(
    i.prev.dx.B < 0.02,
    i.prev.dx.H < 0.02,
    i.prev.dx.W < 0.02
  )

results |>
  select(.iteration, starts_with("ir100.ct")) |>
  group_by(.iteration) |>
  summarise(across(everything(), list(
    min = min,
    mean = mean,
    max = max
  )))
