### utils-format.R
library(dplyr)
library(tidyr)

var_labels <- c(
  "total.proportion_ly"       = "Overall PrEP Coverage",
  "total.prop.neg_ly"         = "Overall PrEP Coverage (among HIV negative)",
  "do.proportion_ly"          = "DO-PrEP Coverage",
  "do.prop.neg_ly"            = "DO-PrEP Coverage (among HIV negative)",
  "edp.proportion_ly"         = "EDP Coverage",
  "edp.prop.neg_ly"           = "EDP Coverage (among HIV negative)",
  "do.prep_ly"                = "% DO-PrEP (proportion of all PrEP users using DO-PrEP)",
  "edp.prep_ly"               = "% EDP (proportion of all PrEP users using DO-PrEP)",
  "total.users_ly"            = "Total PrEP Users",
  "users.perc.diff"           = "Percent Difference in Total PrEP Users",
  "prev_ly"                   = "HIV Prevalence",
  "ir100_ly"                  = "HIV Cumulative Incidence",
  "pia"                       = "Percent Infections Averted (PIA)",
  "nnt"                       = "Number of Pills Needed to Prevent One Infections (NPNT)",
  "pills.perc.diff"           = "Percent Difference in Total Pills Taken",
  "pills_pp"                  = "Number of Pills Needed to Prevent Per PrEP User",
  "do.covered.sex_ly"         = "% of Acts Covered Under DO-PrEP",
  "edp.covered.sex_ly"        = "% of Acts Covered Under EDP",
  "do.covered.discordant_ly"  = "% of Discordant Acts Covered Under DO-PrEP",
  "edp.covered.discordant_ly" = "% of Discordant Acts Covered Under EDP",
  "incid.daily_cml"           = "Incident infections among Daily PrEP Users",
  "incid.edp_cml"             = "Incident infections among Event-Driven PrEP users",
  "incid.daily_ly"            = "Incident Infections in the Last Year (DO PrEP)",
  "incid.edp_ly"              = "Incident Infections in the Last Year (EDP)",
  "do.incid.user.ratio"       = "Ratio of incident HIV among DO PrEP users to DO PrEP users (last year)",
  "edp.incid.user.ratio"      = "Ratio of incident HIV among EDP users to EDP users (last year)"
)

format_patterns <- list(
  small_num = list(
    patterns = c("^ir100"),
    fun = scales::label_number(0.01)
  ),
  perc = list(
    patterns = c("^prev", "pia", "^do.covered", "^edp.covered", ".proportion_ly",
                 ".prep_ly", ".neg_ly", ".perc.diff"),
    fun = scales::label_percent(0.1)
  ),
  perc2 = list(
    patterns = ".ratio",
    fun = scales::label_number(0.01)
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
