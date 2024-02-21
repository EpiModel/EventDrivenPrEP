### utils-format.R
library(dplyr)
library(tidyr)

var_labels <- c(
  "total.proportion_ly"       = "Overall PrEP Coverage",
  "do.proportion_ly"          = "DO-PrEP Coverage",
  "edp.proportion_ly"         = "EDP Coverage",
  "do.prep_ly"                = "% DO-PrEP (proportion of all PrEP users using DO-PrEP",
  "edp.prep_ly"               = "% EDP (proportion of all PrEP users using DO-PrEP",
  "total.users_ly"            = "Total PrEP Users",
  "prev_ly"                   = "HIV Prevalence",
  "ir100_ly"                  = "HIV Cumulative Incidence",
  "pia"                       = "Percent Infections Averted (PIA)",
  "nnt"                       = "Number Needed to Treat (NNT)",
  "pills.perc.diff"           = "Percent Difference in Total Pills Taken",
  "do.covered.sex_ly"         = "% of Acts Covered Under DO-PrEP",
  "edp.covered.sex_ly"        = "% of Acts Covered Under EDP",
  "do.covered.discordant_ly"  = "% of Discordant Acts Covered Under DO-PrEP",
  "edp.covered.discordant_ly" = "% of Discordant Acts Covered Under EDP"
)

format_patterns <- list(
  small_num = list(
    patterns = c("^ir100"),
    fun = scales::label_number(0.01)
  ),
  perc = list(
    patterns = c("^prev", "pia", "^do.covered", "^edp.covered", ".proportion_ly", ".prep_ly", "^pills."),
    fun = scales::label_percent(0.1)
  ),
  # formatter with a catch all pattern. Must be last.
  default = list(
    patterns = ".*",
    fun = scales::label_number(1)
  )
)

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
