### utils-format.R
library(dplyr)
library(tidyr)

format_table <- function(d, var_labels, format_patterns) {
  formatters <- make_formatters(var_labels, format_patterns)

  d_out <- d |>
    sum_quants(0.025, 0.5, 0.975) |>
    pivot_longer(-scenario_name) |>
    separate(name, into = c("name", "quantile"), sep = "_/_") |>
    pivot_wider(names_from = quantile, values_from = value) |>
    filter(name %in% names(var_labels)) |>
    mutate(
      clean_val = purrr::pmap_chr(
        list(name, l, m, h),
        ~ common_format(formatters, ..1, ..2, ..3, ..4))
    ) |>
    select(-c(l, m, h)) |>
    mutate(
      name = var_labels[name]
    ) |>
    pivot_wider(names_from = name, values_from = clean_val) |>
    arrange(scenario_name)

  reorder_cols(d_out, var_labels)
}

make_formatters <- function(var_labels, format_patterns) {
  fmts <- vector(mode = "list", length = length(var_labels))
  for (nms in names(var_labels)) {
    for (fp in format_patterns) {
      if (any(stringr::str_detect(nms, fp$patterns))) {
        fmts[[nms]] <- fp$fun
        break()
      }
    }
  }
  fmts
}


sum_quants <- function(d, ql = 0.025, qm = 0.5, qh = 0.975) {
  d |>
    ungroup() |>
    select(-c(batch_number, sim)) |>
    group_by(scenario_name) |>
    summarise(across(
      everything(),
      list(
        l = ~ quantile(.x, ql, na.rm = TRUE),
        m = ~ quantile(.x, qm, na.rm = TRUE),
        h = ~ quantile(.x, qh, na.rm = TRUE)
      ),
      .names = "{.col}_/_{.fn}"
    ),
    .groups = "drop"
    )
}


reorder_cols <- function(d, var_labels) {
  missing_cols <- setdiff(names(d), var_labels)
  cols_order <- c(missing_cols, intersect(var_labels, names(d)))
  d[, cols_order]
}

common_format <- function(formatters, name, ql, qm, qh) {
  if (is.na(qm)) {
    "-"
  } else {
    paste0(
      formatters[[name]](qm), " (", formatters[[name]](ql),
      ", ", formatters[[name]](qh), ")"
    )
  }
}
