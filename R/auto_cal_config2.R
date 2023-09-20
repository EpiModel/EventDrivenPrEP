source("R/auto_cal_fns.R")
source("R/auto_cal_sim.R")

calib_object <- list(
  waves = list(
    # wave1 = list(
    #   job1 = list(
    #     targets = "cc.prep.B",
    #     targets_val = 0.199,
    #     params = c("prep.start.prob_1"),
    #     initial_proposals = dplyr::tibble(
    #       prep.start.prob_1 = seq(1e-8, 1e-2, length.out = n_sims),
    #     ),
    #     make_next_proposals = make_noisy_proposer(n_sims, floor(n_sims / 3)),
    #     get_result = determ_noisy_end(0.01, 100)
    #     # make_next_proposals = make_shrink_proposer(n_sims),
    #     # get_result = determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job2 = list(
    #     targets = "cc.prep.H",
    #     targets_val = 0.229,
    #     params = c("prep.start.prob_2"),
    #     initial_proposals = dplyr::tibble(
    #       prep.start.prob_2 = seq(1e-8, 1e-2, length.out = n_sims),
    #     ),
    #     make_next_proposals = make_noisy_proposer(n_sims, floor(n_sims / 3)),
    #     get_result = determ_noisy_end(0.01, 100)
    #     # make_next_proposals = make_shrink_proposer(n_sims),
    #     # get_result = determ_poly_end(0.001, poly_n = 5)
    #   ),
    #   job3 = list(
    #     targets = "cc.prep.W",
    #     targets_val = 0.321,
    #     params = c("prep.start.prob_3"),
    #     initial_proposals = dplyr::tibble(
    #       prep.start.prob_3 = seq(1e-8, 1e-2, length.out = n_sims),
    #     ),
    #     make_next_proposals = make_noisy_proposer(n_sims, floor(n_sims / 3)),
    #     get_result = determ_noisy_end(0.01, 100)
    #     # make_next_proposals = make_shrink_proposer(n_sims),
    #     # get_result = determ_poly_end(0.001, poly_n = 5)
    #   )
    # ),
    wave2 = list(
      job1 = list(
        targets = "ir100.gc",
        targets_val = 12.81,
        params = "ugc.prob", # target:
        initial_proposals = dplyr::tibble(
          ugc.prob = seq(0.15, 0.45, length.out = n_sims)
        ),
        make_next_proposals = make_shrink_proposer_rm0(n_sims, shrink = 3 / 2),
        get_result = determ_poly_end_rm0(0.05, poly_n = 5)
      ),
      job1 = list(
        targets = "ir100.ct",
        targets_val = 14.59,
        params = "uct.prob", # target:
        initial_proposals = dplyr::tibble(
          uct.prob = seq(0.15, 0.45, length.out = n_sims)
        ),
        make_next_proposals = make_shrink_proposer_rm0(n_sims, shrink = 3 / 2),
        get_result = determ_poly_end_rm0(0.05, poly_n = 5)
      )
    )
  ),
  config = list(
    simulator = calibration2_fun,
    default_proposal = dplyr::tibble(
      prep.start.prob_1 = 9.782055e-05,
      prep.start.prob_2 = 1e-4,
      prep.start.prob_3 = 0.0007918135,
      # remove after
      ugc.prob = 0.2822787,
      uct.prob = 0.2703707
    ),
    root_directory = "data/calib",
    max_iteration = 100,
    n_sims = n_sims
  )
  # state = list() # managed internally
)
