# Utilities --------------------------------------------------------------------
# makes the trackers race specific
# indiv == TRUE: 3 trackers: tracker___B, tracker___H, tracker___W
# full == TRUE: add one trackers: tracker___ALL
epi_tracker_by_race <- function(ls_funs, races = 1:3,
                                races_names = c("B", "H", "W"),
                                indiv = TRUE, full = TRUE) {

 ls_races <- if (indiv) as.list(races) else list()
 races_names <- if (indiv) races_names else c()

  if (full) {
    ls_races <- c(ls_races, list(races))
    races_names <- c(races_names, "ALL")
  }

  epi_tracker <- lapply(
    ls_races,
    function(race) lapply(ls_funs, do.call, args = list(r_ind = race))
  )

  epi_tracker <- unlist(epi_tracker)
  names(epi_tracker) <- paste0(
    names(epi_tracker), "___",
    unlist(lapply(races_names, rep, times = length(ls_funs)))
  )

  epi_tracker
}

# Trackers ---------------------------------------------------------------------
epi_n <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "active")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & active == 1, na.rm = TRUE)
    })
  }
}

# HIV Trackers
epi_s <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 0, na.rm = TRUE)
    })
  }
}

# eligible to prep
epi_s_prep_elig <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status", "prepElig")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 0 & prepElig == 1, na.rm = TRUE)
    })
  }
}

# on prep
epi_s_prep <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status", "prepStat")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 0 & prepStat == 1, na.rm = TRUE)
    })
  }
}

epi_i <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 1, na.rm = TRUE)
    })
  }
}

epi_i_dx <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status", "diag.status")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 1 & diag.status == 1, na.rm = TRUE)
    })
  }
}

epi_i_tx <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status", "tx.status")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 1 & tx.status == 1, na.rm = TRUE)
    })
  }
}

epi_i_sup <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status", "vl.last.supp")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind & status == 1 & vl.last.supp == at, na.rm = TRUE)
    })
  }
}

epi_i_sup_dur <- function(r_ind) {
  function(dat, at) {
    needed_attributes <- c("race", "status", "vl.last.usupp")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% r_ind &
          status == 1 &
          at - vl.last.usupp >= 52,
          na.rm = TRUE)
    })
  }
}

# linked in less than `weeks` step
epi_linked_time <- function(weeks) {
  function(r_ind) {
    function(dat, at) {
      needed_attributes <- c("race", "tx.init.time", "diag.time")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% r_ind &
          tx.init.time - diag.time <= weeks,
          na.rm = TRUE
        )
      })
    }
  }
}

# STI trackers
epi_gc_i <- function(hiv_status) {
  function(r_ind) {
    function(dat, at) {
      needed_attributes <- c("race", "rGC", "uGC", "status")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% r_ind &
          status %in% hiv_status &
          (rGC == 1 | uGC == 1),
          na.rm = TRUE
        )
      })
    }
  }
}

epi_ct_i <- function(hiv_status) {
  function(r_ind) {
    function(dat, at) {
      needed_attributes <- c("race", "rCT", "uCT", "status")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% r_ind &
          status %in% hiv_status &
          (rCT == 1 | uCT == 1),
          na.rm = TRUE
        )
      })
    }
  }
}

epi_gc_s <- function(hiv_status) {
  function(r_ind) {
    function(dat, at) {
      needed_attributes <- c("race", "rGC", "uGC", "status")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% r_ind &
          status %in% hiv_status &
          (rGC == 0 & uGC == 0),
          na.rm = TRUE
        )
      })
    }
  }
}

epi_ct_s <- function(hiv_status) {
  function(r_ind) {
    function(dat, at) {
      needed_attributes <- c("race", "rCT", "uCT", "status")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% r_ind &
          status %in% hiv_status &
          (rCT == 0 & uCT == 0),
          na.rm = TRUE
        )
      })
    }
  }
}
