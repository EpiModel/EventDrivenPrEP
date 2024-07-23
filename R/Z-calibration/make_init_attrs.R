# Scratchpad for interactive testing before integration in a script

source("./R/shared_variables.R")
library(dplyr)

sample_join <- function(x, y, by) {
  x <- dplyr::mutate(x, .orig.pos = seq_len(dplyr::n()))
  d_combs <- lapply(x[by], unique) |> expand.grid() |> dplyr::as_tibble()

  d_final_list <- vector(mode = "list", length = nrow(d_combs))
  for (i in seq_len(nrow(d_combs))) {
    dt_start <- semi_join(x, d_combs[i, ])
    dt_target <- semi_join(y, d_combs[i, ])
    missing_cols <- setdiff(names(dt_start), names(dt_target))

    dt_target <- sample_n(dt_target, nrow(dt_start), replace = TRUE)
    d_final_list[[i]] <- dplyr::bind_cols(dt_target, dt_start[missing_cols])
  }

  dplyr::bind_rows(d_final_list) |>
    dplyr::arrange(.orig.pos) |>
    dplyr::select(-.orig.pos)
}

# get attributes from the network esitmate
est <- readRDS("./data/intermediate/estimates/netest-hpc.rds")
nw <- est[[1]][["newnetwork"]]
otha <- network::list.vertex.attributes(nw)
otha <- setdiff(
  otha,
  c("na", "vertex.names", "active", "testatus.active", "tergm_pid")
)

d_attrs <- lapply(
  otha, network::get.vertex.attribute,
  x = nw, na.omit = FALSE, null.na = TRUE, unlist = TRUE
)
names(d_attrs) <- otha
d_attrs <- as_tibble(d_attrs)

# prepare the "target" population
sim <- readRDS("./data/intermediate/estimates/restart-hpc.rds")
d_tar <- as_tibble(sim$attr[[1]])

d_tar <- select(d_tar, -c(starts_with("deg."))) |>
  mutate(across(
    c(ends_with(".time"), ends_with("infTime"), starts_with("prep.ind."),
      last.neg.test),
    function(x) x - calibration_end
  ))

# format `diag.status` to be able to join the dfs
d_attrs$dg2 <- d_attrs$diag.status
d_tar$dg2 <- ifelse(is.na(d_tar$diag.status), 0, d_tar$diag.status)

# sample join and save
d_final <- sample_join(d_attrs, d_tar, by = c("race", "dg2", "role.class"))
d_final$dg2 <- NULL
saveRDS(d_final, "d_init_attr.rds")
