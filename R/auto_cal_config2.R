source("R/auto_cal_fns.R")
source("R/auto_cal_sim.R")

calib_object <- list(
  waves = list(
    wave1 = list(
      job1 = list(
        targets = "cc.prep.B",
        targets_val = 0.199,
        params = c("prep.start.prob_1"),
        initial_proposals = dplyr::tibble(
          prep.start.prob_1 = seq(0.001, 0.01, length.out = n_sims),
        ),
        make_next_proposals = make_shrink_proposer(n_sims),
        get_result = determ_poly_end(0.001, poly_n = 5)
      ),
      job2 = list(
        targets = "cc.prep.H",
        targets_val = 0.229,
        params = c("prep.start.prob_2"),
        initial_proposals = dplyr::tibble(
          prep.start.prob_2 = seq(0.001, 0.01, length.out = n_sims),
        ),
        make_next_proposals = make_shrink_proposer(n_sims),
        get_result = determ_poly_end(0.001, poly_n = 5)
      ),
      job3 = list(
        targets = "cc.prep.W",
        targets_val = 0.321,
        params = c("prep.start.prob_3"),
        initial_proposals = dplyr::tibble(
          prep.start.prob_3 = seq(0.001, 0.01, length.out = n_sims),
        ),
        make_next_proposals = make_shrink_proposer(n_sims),
        get_result = determ_poly_end(0.001, poly_n = 5)
      )
    ),
    wave2 = list(
      job1 = list(
        targets = "ir100.gc",
        targets_val = 12.81,
        params = c("ugc.prob", "rgc.prob"), # target:
        initial_proposals = dplyr::tibble(
          ugc.prob = seq(0.15, 0.45, length.out = n_sims),
          rgc.prob = plogis(qlogis(ugc.prob) + log(1.25))
        ),
        make_next_proposals = make_sti_range_proposer(n_sims),
        get_result = determ_trans_end(
          retain_prop = 0.3,
          thresholds = 1,
          n_enough = 300
        )
      ),
      job2 = list(
        targets = "ir100.ct",
        targets_val = 14.59,
        params = c("uct.prob", "rct.prob"), # target:
        initial_proposals = dplyr::tibble(
          uct.prob = seq(0.15, 0.45, length.out = n_sims),
          rct.prob = plogis(qlogis(uct.prob) + log(1.25))
        ),
        make_next_proposals = make_sti_range_proposer(n_sims),
        get_result = determ_trans_end(
          retain_prop = 0.3,
          thresholds = 1,
          n_enough = 300
        )
      )
    )
  ),
  config = list(
    simulator = calibration2_fun,
    default_proposal = dplyr::tibble(
      prep.start.prob_1 = 0.006,
      prep.start.prob_2 = 0.006,
      prep.start.prob_3 = 0.006,
      # remove after
      ugc.prob = 0.1902003,
      rgc.prob = plogis(qlogis(ugc.prob) + log(1.25)),
      uct.prob = 0.1714429,
      rct.prob = plogis(qlogis(uct.prob) + log(1.25)),
    ),
    root_directory = "data/calib",
    max_iteration = 100,
    n_sims = n_sims
  )
  # state = list() # managed internally
)

