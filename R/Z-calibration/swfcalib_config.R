## swfcalib Configuration 1 (pre-prep)
##
## Set up the configuration for the first calibration. This takes place before
## the restart point and before PrEP is started.
##
## This script should not be run directly. But `sourced` from the swfcalib
## workflow

source("R/Z-calibration/swfcalib_model.R", local = TRUE)
model_fn <- make_model_fn(restart = FALSE, calib_steps = year_steps)
targets <- EpiModelHIV::get_calibration_targets()

params_df <- readr::read_csv("data/input/params.csv") |>
  dplyr::select(-type, -source) |>
  dplyr::mutate(value = as.numeric(value)) |>
  tidyr::pivot_wider(names_from = param)

n_sims <- 512

calib_object <- list(
  # state = list() # managed internally
  config = list(
    simulator = model_fn,
    root_directory = "data/calib",
    max_iteration = 100,
    n_sims = n_sims,
    default_proposal = dplyr::select(
      params_df,
      prep.start.prob_1, prep.start.prob_2, prep.start.prob_3,
      a.rate,
      aids.off.tx.mort.rate,
      hiv.test.rate_1, hiv.test.rate_2, hiv.test.rate_3,
      tx.init.rate_1, tx.init.rate_2, tx.init.rate_3,
      ugc.prob,
      uct.prob,
      tx.halt.partial.rate_1, tx.halt.partial.rate_2, tx.halt.partial.rate_3,
      hiv.trans.scale_1, hiv.trans.scale_2, hiv.trans.scale_3,
      mddcoef.hiv_1, mddcoef.hiv_2,
      mdd.diag.gen.prob_1, mdd.diag.gen.prob_2,
      mdd.prob.base
    )
  ),
  waves = list(
    # wave1 = list(
    #   job0 = list(
    #     targets = "disease.mr100",
    #     targets_val = targets["disease.mr100"],
    #     params = c("aids.off.tx.mort.rate"), # target: 0.00385
    #     initial_proposals = dplyr::tibble(
    #       aids.off.tx.mort.rate = sample(seq(0.0001, 0.001, length.out = n_sims)),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job1 = list(
    #     targets = "cc.prep.B",
    #     targets_val = targets["cc.prep.B"],
    #     params = c("prep.start.prob_1"),
    #     initial_proposals = dplyr::tibble(
    #       prep.start.prob_1 = seq(0.001, 0.01, length.out = n_sims),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job2 = list(
    #     targets = "cc.prep.H",
    #     targets_val = targets["cc.prep.H"],
    #     params = c("prep.start.prob_2"),
    #     initial_proposals = dplyr::tibble(
    #       prep.start.prob_2 = seq(0.001, 0.01, length.out = n_sims),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job3 = list(
    #     targets = "cc.prep.W",
    #     targets_val = targets["cc.prep.W"],
    #     params = c("prep.start.prob_3"),
    #     initial_proposals = dplyr::tibble(
    #       prep.start.prob_3 = seq(0.005, 0.02, length.out = n_sims),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job4 = list(
    #     targets = "mdd.prphiv0",
    #     targets_val = targets["mdd.prphiv0"],
    #     params = c("mdd.prob.base"),
    #     initial_proposals = dplyr::tibble(
    #       mdd.prob.base = seq(0.1, 0.3, length.out = n_sims),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job5 = list(
    #     targets = c("mdd.diag.prphiv0"),
    #     targets_val = targets[c("mdd.diag.prphiv0")],
    #     params = c("mdd.diag.gen.prob_1"),
    #     initial_proposals = dplyr::tibble(
    #       mdd.diag.gen.prob_1 = seq(0.08, 0.3, length.out = n_sims)
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job5 = list(
    #     targets = c("mdd.diag.prphiv1"),
    #     targets_val = targets[c("mdd.diag.prphiv1")],
    #     params = c("mdd.diag.gen.prob_2"),
    #     initial_proposals = dplyr::tibble(
    #       mdd.diag.gen.prob_2 = seq(0.08, 0.3, length.out = n_sims),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   )
    # ),

    wave2 = list(
      job1 = list(
        targets = "cc.dx.B",
        targets_val = targets["cc.dx.B"],
        params = c("hiv.test.rate_1"), # target: 0.00385
        initial_proposals = dplyr::tibble(
          hiv.test.rate_1 = sample(seq(0.002, 0.0035, length.out = n_sims)),
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),
      job2 = list(
        targets = "cc.dx.H",
        targets_val = targets["cc.dx.H"],
        params = c("hiv.test.rate_2"), # target: 0.0038
        initial_proposals = dplyr::tibble(
          hiv.test.rate_2 = sample(seq(0.002, 0.0035, length.out = n_sims)),
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
        ),
      job3 = list(
        targets = "cc.dx.W",
        targets_val = targets["cc.dx.W"],
        params = c("hiv.test.rate_3"), # target: 0.0069
        initial_proposals = dplyr::tibble(
          hiv.test.rate_3 = sample(seq(0.002, 0.0035, length.out = n_sims)),
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      )
    # ,
    #   job6 = list(
    #     targets = paste0("cc.linked1m.", c("B", "H", "W")),
    #     targets_val = targets[paste0("cc.linked1m.", c("B", "H", "W"))],
    #     params = paste0("tx.init.rate_", 1:3),
    #     initial_proposals = dplyr::tibble(
    #       tx.init.rate_1 = sample(seq(0.1, 0.5, length.out = n_sims)),
    #       tx.init.rate_2 = sample(tx.init.rate_1),
    #       tx.init.rate_3 = sample(tx.init.rate_1)
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 3)
    #   ),
    #   job7 = list(
    #     targets = "mdd.prphiv1",
    #     targets_val = targets["mdd.prphiv1"],
    #     params = c("mddcoef.hiv_2"),
    #     initial_proposals = dplyr::tibble(
    #       mddcoef.hiv_2 = seq(1.1, 1.5, length.out = n_sims)
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   )
    ),
    wave3 = list(
      job1 = list(
        targets = paste0("cc.vsupp.", c("B", "H", "W")),
        targets_val = targets[paste0("cc.vsupp.", c("B", "H", "W"))],
        params = paste0("tx.halt.partial.rate_", 1:3),
        initial_proposals = dplyr::tibble(
          tx.halt.partial.rate_1 = sample(seq(0.002, 0.006, length.out = n_sims)),
          tx.halt.partial.rate_2 = sample(tx.halt.partial.rate_1),
          tx.halt.partial.rate_3 = sample(tx.halt.partial.rate_1)
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 3)
      )
    ),

    wave4 = list(
      job4 = list(
        targets = "ir100.gc",
        targets_val = targets["ir100.gc"],
        params = c("ugc.prob"), # target:
        initial_proposals = dplyr::tibble(
          ugc.prob = sample(seq(0.23, 0.26,length.out = n_sims)),
          ),
        make_next_proposals = swfcalib::make_proposer_se_range(n_sims, retain_prop = 0.3),
        get_result = swfcalib::determ_end_thresh(
          thresholds = 1,
          n_enough = 100
        )
      ),
      job5 = list(
        targets = "ir100.ct",
        targets_val = targets["ir100.ct"],
        params = c("uct.prob"), # target:
        initial_proposals = dplyr::tibble(
          uct.prob = sample(seq(0.16, 0.2, length.out = n_sims)),
        ),
        make_next_proposals = swfcalib::make_proposer_se_range(n_sims, retain_prop = 0.3),
        get_result = swfcalib::determ_end_thresh(
          thresholds = 1,
          n_enough = 100
        )
      )
    ),
    wave5 = list(
      job1 = list(
        targets = paste0("i.prev.dx.", c("B", "H", "W")),
        targets_val = targets[paste0("i.prev.dx.", c("B", "H", "W"))],
        params = paste0("hiv.trans.scale_", 1:3),
        initial_proposals = dplyr::tibble(
          hiv.trans.scale_1 = sample(seq(2, 4, length.out = n_sims)),
          hiv.trans.scale_2 = sample(seq(0.2, 0.6, length.out = n_sims)),
          hiv.trans.scale_3 = sample(seq(0.2, 0.6, length.out = n_sims))
        ),
        make_next_proposals = swfcalib::make_proposer_se_range(n_sims, retain_prop = 0.3),
        get_result = swfcalib::determ_end_thresh(
          thresholds = c(0.02, 0.02, 0.01),
          n_enough = 100
        )
      )
    )
    # ,
    # wave6 = list(
    #   job0 = list(
    #     targets = "disease.mr100",
    #     targets_val = targets["disease.mr100"],
    #     params = c("aids.off.tx.mort.rate"), # target: 0.00385
    #     initial_proposals = dplyr::tibble(
    #       aids.off.tx.mort.rate = sample(seq(0.0001, 0.001, length.out = n_sims)),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
    #   )
    # ),
    # wave7 = list(
    #   job0 = list(
    #     targets = "num",
    #     targets_val = targets["num"],
    #     params = c("a.rate"),
    #     initial_proposals = dplyr::tibble(
    #       aids.off.tx.mort.rate = sample(seq(0.0004, 0.0005, length.out = n_sims)),
    #     ),
    #     make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
    #     get_result = swfcalib::determ_poly_end(1, poly_n = 5)
    #   )
    # )
  )
)
