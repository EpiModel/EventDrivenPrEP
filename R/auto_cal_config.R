source("R/auto_cal_fns.R")
source("R/auto_cal_sim.R")

calib_object <- list(
  waves = list(
    wave1 = list(
      job0 = list(
        targets = "disease.mr100",
        targets_val = 0.273,
        params = c("aids.off.tx.mort.rate"), # target: 0.00385
        initial_proposals = dplyr::tibble(
          aids.off.tx.mort.rate = seq(0.0002, 0.0007, length.out = n_sims),
        ),
        make_next_proposals = make_shrink_proposer(n_sims),
        get_result = determ_poly_end(0.001, poly_n = 5)
      ),
      # job1 = list(
      #   targets = "cc.dx.B",
      #   targets_val = 0.847,
      #   params = c("hiv.test.rate_1"), # target: 0.00385
      #   initial_proposals = dplyr::tibble(
      #     hiv.test.rate_1 = seq(0.0004, 0.001, length.out = n_sims),
      #     ),
      #   make_next_proposals = make_shrink_proposer(n_sims),
      #   get_result = determ_poly_end(0.001, poly_n = 5)
      # ),
      # job2 = list(
      #   targets = "cc.dx.H",
      #   targets_val = 0.818,
      #   params = c("hiv.test.rate_2"), # target: 0.0038
      #   initial_proposals = dplyr::tibble(
      #     hiv.test.rate_2 = seq(0.0004, 0.001, length.out = n_sims),
      #     ),
      #   make_next_proposals = make_shrink_proposer(n_sims),
      #   get_result = determ_poly_end(0.001, poly_n = 5)
      #   ),
      # job3 = list(
      #   targets = "cc.dx.W",
      #   targets_val = 0.862,
      #   params = c("hiv.test.rate_3"), # target: 0.0069
      #   initial_proposals = dplyr::tibble(
      #     hiv.test.rate_3 = seq(0.0005, 0.0012, length.out = n_sims),
      #   ),
      #   make_next_proposals = make_shrink_proposer(n_sims),
      #   get_result = determ_poly_end(0.001, poly_n = 5)
      #   ),
      job4 = list(
        targets = "ir100.gc",
        targets_val = 12.81,
        params = c("ugc.prob", "rgc.prob"), # target:
        initial_proposals = dplyr::tibble(
          ugc.prob = seq(0.2, 0.3, length.out = n_sims),
          rgc.prob = plogis(qlogis(ugc.prob) + log(1.25))
        ),
        make_next_proposals = make_sti_range_proposer(n_sims),
        get_result = determ_trans_end(
          retain_prop = 0.3,
          thresholds = 1.5,
          n_enough = 100
        )
      ),
      job5 = list(
        targets = "ir100.ct",
        targets_val = 14.59,
        params = c("uct.prob", "rct.prob"), # target:
        initial_proposals = dplyr::tibble(
          uct.prob = seq(0.15, 0.25, length.out = n_sims),
          rct.prob = plogis(qlogis(uct.prob) + log(1.25))
        ),
        make_next_proposals = make_sti_range_proposer(n_sims),
        get_result = determ_trans_end(
          retain_prop = 0.3,
          thresholds = 1.5,
          n_enough = 100
        )
      ),
    #   job6 = list(
    #     targets = paste0("cc.linked1m.", c("B", "H", "W")),
    #     targets_val = c(0.829, 0.898, 0.890),
    #     params = paste0("tx.init.rate_", 1:3),
    #     initial_proposals = dplyr::tibble(
    #       tx.init.rate_1 = sample(seq(0.01, 0.2, length.out = n_sims)),
    #       tx.init.rate_2 = sample(tx.init.rate_1),
    #       tx.init.rate_3 = sample(tx.init.rate_1),
    #     ),
    #     make_next_proposals = make_ind_shrink_proposer(n_sims),
    #     get_result = determ_ind_poly_end(0.001, poly_n = 3)
    #   )
    # ),
    wave2 = list(
      job1 = list(
        targets = paste0("cc.vsupp.", c("B", "H", "W")),
        targets_val = c(0.602, 0.620, 0.710),
        params = paste0("tx.halt.partial.rate_", 1:3),
        initial_proposals = dplyr::tibble(
          tx.halt.partial.rate_1 = sample(seq(0.0001, 0.001, length.out = n_sims)),
          tx.halt.partial.rate_2 = sample(tx.halt.partial.rate_1),
          tx.halt.partial.rate_3 = sample(tx.halt.partial.rate_1)
        ),
        make_next_proposals = make_ind_shrink_proposer(n_sims),
        get_result = determ_ind_poly_end(0.001, poly_n = 3)
      )
    ) ,
    wave3 = list(
      job1 = list(
        targets = paste0("i.prev.dx.", c("B", "H", "W")),
        targets_val = c(0.33, 0.127, 0.09),
        params = paste0("hiv.trans.scale_", 1:3),
        initial_proposals = {
          job <- list()
          job$params <- paste0("hiv.trans.scale_", 1:3)
          p_ranges <- list(
            c(1, 5),
            c(0.1, 1),
            c(0.1, 1)
          )
          proposals <- lhs::randomLHS(n_sims, length(job$params))
          for (i in seq_along(job$params)) {
            spread <- p_ranges[[i]][2] - p_ranges[[i]][1]
            rmin <- p_ranges[[i]][1]
            proposals[, i] <- proposals[, i] * spread + rmin
          }
          colnames(proposals) <- job$params
          dplyr::as_tibble(proposals)
        },
        make_next_proposals = make_range_proposer_lhs(n_sims),
        get_result = determ_trans_end_rmse(
          retain_prop = 0.4,
          thresholds = rep(0.02, 3),
          n_enough = 100
        )
        # get_result = determ_lin_poly_end(c(0.005, 0.01, 0.01), poly_n = 1)
      )
    )

  ),
  config = list(
    simulator = model_fun,
    default_proposal = dplyr::tibble(
      hiv.test.rate_1 = 0.00072050,
      hiv.test.rate_2 = 0.00061664,
      hiv.test.rate_3 = 0.00085655,
      tx.init.rate_1 = 0.0600,
      tx.init.rate_2 = 0.0747,
      tx.init.rate_3 = 0.0726,
      ugc.prob = 0.01,
      rgc.prob = plogis(qlogis(ugc.prob) + log(1.25)),
      uct.prob = 0.01,
      rct.prob = plogis(qlogis(uct.prob) + log(1.25)),
      tx.halt.partial.rate_1 = 0.0007740,
      tx.halt.partial.rate_2 = 0.0007167,
      tx.halt.partial.rate_3 = 0.0004727,
      hiv.trans.scale_1 = 2.5,
      hiv.trans.scale_2 = 0.5,
      hiv.trans.scale_3 = 0.3,
      aids.off.tx.mort.rate = 0.0006
    ),
    root_directory = "data/calib",
    max_iteration = 100,
    n_sims = n_sims
  )
  # state = list() # managed internally
)
