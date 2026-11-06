#' ----------------------------------------------------------------------------#
#' Title: Helper Function File                                                 #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 25-03-2024                                                    #
#' -----                                                                       #
#' Last Modified: 13-02-2025                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#



### Functions ------------------------------------------------------------------
simulate <- function(
    gen_model_list, method_list, conditions, repetitions,
    out_dir = getwd()) {
  #' Main simulation function
  #' gen_model_list is a list of objects of class gen_model
  #' method_list is a list of objects inheriting from the class method with
  #' their own fit and calculate_performance measures methods
  #' conditions is a list of simulation conditions
  #' repetitions is a integer giving the number of repetitions per condition

  require(data.table)

  ## Setup-simulation and initialize objects

  # Create sim object by initiating objects from the gen_model list
  sim_grid <- as.data.table(
    expand.grid(c(conditions, list(gen_model = gen_model_list)))
  )
  sim_grid$conditions <- apply(
    subset(sim_grid, select = -gen_model),
    1, function(x) as.list(x)
  )

  # Incorporate simulation conditions into gen_model
  sim_grid[, gen_model := mapply(add_conditions, gen_model, conditions)]
  ncond <- nrow(sim_grid)

  # Replicate the rows of sim_grid repetitions time
  sim_grid <- sim_grid[rep(seq_len(nrow(sim_grid)), each = repetitions), ]

  cat(
    "In total,", nrow(sim_grid), "data sets will be simulated and",
    length(method_list) * nrow(sim_grid), "analyses will be performed."
  )

  ## Simulate the data
  cat("\n\nSimulating Data: ")
  sim_count <- counter()
  sim_grid[, dat := future_lapply(
    gen_model,
    function(x) {
      cat(sim_count(), " ")
      get_tsm_data(sim_tsm(x))
    },
    future.seed = TRUE
  )]

  ## Add analysis methods to the simulation and fill in information
  if (length(method_list) == 1) {
    sim_grid[, method := list(list(method_list))]
  } else {
    sim_grid[, method := list(method_list)]
  }

  sim_grid[, method := mapply(
    function(method, conditions, gen_model) {
      lapply(method, fill_in_method,
        conditions = conditions,
        gen_model = gen_model
      )
    }, method, conditions, gen_model,
    SIMPLIFY = FALSE
  )]

  ## Fit analysis methods and record performance measures
  cat("\n\nModel Fitting:    ")
  fit_count <- counter()
  # Fit model by condition
  for (i in seq_len(ncond)) {
    # Indicator for all repetitions of a condition
    ind <- seq(1, repetitions) + repetitions * (i - 1)
    # Model fitting loop
    sim_grid$method[ind] <- mapply(
      function(method, data) {
        cat(fit_count(), " ")
        lapply(method, fit, data = data)
      },
      method = sim_grid$method[ind], data = sim_grid$dat[ind],
      SIMPLIFY = FALSE
    )
    # Save intermittend results in temp file
    save(sim_grid, file = paste0(out_dir, "/temp.Rdata"))
  }

  c("\n\n")

  # Return simulation grid
  return(sim_grid)
}

restart_simulation <- function(out_dir) {
  # Load partial simulation
  load(paste0(out_dir, "/temp.Rdata"))

  # Determine where the simulation stopped
  not_converged <- sapply(sim_grid$method, function(method) {
    all(!sapply(method, function(x) slot(x, "converged")))
  })

  for (i in seq_along(not_converged)) {
    if (not_converged[i] && all(not_converged[i:length(not_converged)])) {
      start <- i
      break
    }
  }

  # Derive repetitions
  repetitions <- nrow(sim_grid) / length(unique(sim_grid$gen_model))

  # Restart simulation
  cat("\n\nModel Fitting:    ")
  fit_count <- counter()

  for (i in seq_len(start - 1)) {
    cat(fit_count(), " ")
  }

  for (i in ceiling(start / repetitions):length(unique(sim_grid$gen_model))) {
    # Indicator for all repetitions of a condition
    ind <- seq(1, repetitions) + repetitions * (i - 1)

    # Model fitting loop
    sim_grid$method[ind] <- mapply(
      function(method, data) {
        cat(fit_count(), " ")
        lapply(method, fit, data = data)
      },
      method = sim_grid$method[ind], data = sim_grid$dat[ind],
      SIMPLIFY = FALSE
    )
    # Save intermittend results in temp file
    save(sim_grid, file = paste0(out_dir, "/temp.Rdata"))
  }

  c("\n\n")

  # Return simulation grid
  return(sim_grid)
}

add_conditions <- function(model, conditions) {
  #' Convenience function to add the conditions list to the respective slots
  #' in the gen_model objects. Conditions will be added in order and according
  #' to name, first to slots with the same name, then to parameters with the
  #' same name and all remaining conditions will be appended to the paramters
  #' slot.

  # Extract slot names
  slot_names <- slotNames(model)
  pars_names <- names(slot(model, "pars"))

  # Assign conditions that are in slot names to slot
  ind_slotn <- which(names(conditions) %in% slot_names)
  for (ind in ind_slotn) {
    slot(model, names(conditions)[ind]) <- conditions[[ind]]
  }

  # Assign conditions that are in pars names to pars
  ind_parsn <- which(names(conditions) %in% pars_names)
  for (ind in ind_parsn) {
    slot(model, "pars")[[names(conditions)[ind]]] <- conditions[[ind]]
  }

  # Append all other conditions to the pars vector
  slot(model, "pars") <- c(
    slot(model, "pars"),
    conditions[-c(ind_parsn, ind_slotn)]
  )

  return(model)
}

fill_in_method <- function(method, conditions, gen_model) {
  #' Convenience functio to copy over information and fill in the meta slots
  #' of the method objects
  slot(method, "gen_model") <- slot(gen_model, "model_name")
  slot(method, "conditions") <- conditions
  return(method)
}

sim_tsm <- function(model) {
  #' Function to simulate from a dynamic time-series model that is specified
  #' in a gen_model object.

  # Initialize data frame to hold simulation results
  dat <- data.frame(time = seq(
    from = 0, to = slot(model, "time"),
    by = slot(model, "stepsize")
  ))

  # Set initial state value to specified starting value
  dat$state[[1]] <- lapply(slot(model, "start"), function(x) {
    eval(x[[3]], slot(model, "pars"))
  })
  names(dat$state[[1]]) <- lapply(slot(model, "start"), function(x) x[[2]])

  if (slot(model, "model_type") == "SSM") {
    # Apply state function iteratively to get state values starting form the
    # second time point
    for (t in seq(2, nrow(dat))) {
      dat$state[[t]] <- mapply(
        # Apply function to combinations of state and dynamic error expressions
        function(x, y) {
          # Evaluate dynamic error function in parameter environment
          error <- eval(
            str2lang(
              paste(deparse(y[[3]]),
                rnorm(1, 0, 1),
                sep = " * "
              )
            ),
            list2env(c(dat$state[[t - 1]], slot(model, "pars")))
          )
          # Evaluate state function in parameter environment and add error
          state_plus_one <- eval(
            str2lang(paste(deparse(x[[3]]), error,
              sep = " + "
            )),
            list2env(c(dat$state[[t - 1]], slot(model, "pars")))
          )
          return(slot(model, "boundary")(state_plus_one))
        },
        x = slot(model, "state_eq"), y = slot(model, "dynamic_error"),
        SIMPLIFY = FALSE
      )
      # Add names
      names(dat$state[[t]]) <- lapply(
        slot(model, "state_eq"),
        function(x) x[[2]]
      )
    }
  } else if (slot(model, "model_type") == "DE") {
    # Run Euler-Maruyama method starting from second time point
    # Instantiate new state list that holds the oversampled state
    state <- list(dat$state[[1]])
    for (t in seq(2, dat$time[nrow(dat)] / slot(model, "delta") + 1)) {
      state[[t]] <- mapply(
        # Apply function to combinations of state and dynamic error expressions
        function(x, y) {
          # Evaluate the Wiener (diffusion) component in the parameter env.
          wiener <- eval(
            str2lang(
              paste(deparse(y[[3]]),
                rnorm(
                  1, 0,
                  sqrt(slot(model, "delta"))
                ),
                sep = " * "
              )
            ),
            list2env(c(state[[t - 1]], slot(model, "pars")))
          )
          # Evaluate the differential (drift) component in the parameter env
          # and add previous state value und Wiener component
          state_plus_one <- eval(
            str2lang(paste0(
              deparse(x[[2]]), " + ",
              slot(model, "delta"), "*(",
              deparse(x[[3]]), ")", " + ", wiener
            )),
            list2env(c(state[[t - 1]], slot(model, "pars")))
          )
          return(slot(model, "boundary")(state_plus_one))
        },
        x = slot(model, "state_eq"), y = slot(model, "dynamic_error"),
        SIMPLIFY = FALSE
      )
      # Add names
      names(state[[t]]) <- lapply(
        slot(model, "state_eq"),
        function(x) x[[2]]
      )
    }
    # Copy subsample at the desired time points over to data frame
    # The rounding is unelegant but required to make the matching work with
    # floating values.
    dat$state <- state[which(round(seq(
      from = dat$time[1], to = dat$time[nrow(dat)],
      by = slot(model, "delta")
    ), digits = 10) %in% round(seq(
      from = dat$time[1], to = dat$time[nrow(dat)],
      by = slot(model, "stepsize")
    ), digits = 10))]
  }

  # Generate observation values from state values by adding observation errors
  dat$observation <- lapply(
    dat$state,
    function(state) {
      val <- lapply(
        slot(model, "meas_eq"),
        function(x) {
          eval(
            x[[3]],
            list2env(c(state, slot(model, "pars")))
          )
        }
      )
      names(val) <- lapply(slot(model, "meas_eq"), function(x) x[[2]])
      return(val)
    }
  )

  return(dat)
}

get_tsm_data <- function(data) {
  #' Convenience function to extract and reformat data generated by sim_tsm.
  #' This function takes the individual lists created by sim_tsm and reformats
  #' them into a data frame.
  # Validation
  require(data.table)

  # Select data featrues
  res <- list()

  res$time <- data$time
  res$observation <- do.call(
    rbind.data.frame,
    lapply(data$observation, rbind)
  )

  res$state <- do.call(rbind.data.frame, lapply(data$state, rbind))

  nam <- unlist(sapply(res, names)) # Extract names
  nam <- c("time", nam) # Add time name
  res <- do.call(cbind, res) # Reformat data
  names(res) <- nam # Fix names

  res <- as.data.table(sapply(res, unlist)) # Get rid of lists within columns

  return(res)
}

calc_mse <- function(method, data) {
  #' Convenience function to calculate mse between method estimate and data
  return(sqrt(mean((slot(method, "estimate") - data$y)^2)))
}

ci_test <- function(method, data) {
  #' Convenience function ot test if the actual state lies in the confidence
  #' interval

  ci_incl <- mapply(function(y, ub, lb) y < ub & y > lb,
    y = data$y, ub = slot(method, "ci")$ub, lb = slot(method, "ci")$lb
  )

  ci_coverage <- mean(ci_incl)
  return(ci_coverage)
}

extract_results <- function(sim) {
  #' Convenience function to extract the performance measures (i.e. mse, gcv,
  #' and ci coverage) as well as all meta info from the two lists created during
  #' the simulation.


  dat <- subset(sim, select = c(-gen_model, -conditions, -dat, -method))[
    rep(seq_len(nrow(sim)), each = length(sim$method[[1]]))
  ]

  dat[, method := unlist(sim$method)]

  results <- lapply(dat$method, function(x) {
    c(list(
      method = slot(x, "method_name"),
      model = slot(x, "gen_model"),
      mse = slot(x, "mse"),
      gcv = slot(x, "gcv"),
      ci_coverage = slot(x, "ci_coverage")
    ), slot(x, "conditions"))
  })

  results <- do.call(rbind.data.frame, results)

  return(results)
}

gp_post_pred <- function(
    fit, time, obs, state, f_name, draws = 200,
    alpha = 1, lwd = 1, rect = FALSE,
    quant = c(.60, .80, .95), ...) {
  #' Function for the posterior predictive plotting of a gaussian process

  f_draws <- as.data.frame(fit$draws(f_name, format = "draws_df"))
  f_draws <- f_draws[, seq(1, ncol(f_draws) - 3)]

  if (!rect) {
    plot(x = time, y = obs, ...)

    lapply(sample(seq(1, nrow(f_draws)), draws),
      function(x, f_draws, time, lwd, alpha) {
        lines(
          x = time, y = f_draws[x, seq_len(length(time))],
          col = rgb(red = 1, green = 0, blue = 0, alpha = alpha),
          lwd = lwd
        )
      },
      f_draws = f_draws, time = time, alpha = alpha, lwd = lwd
    )

    lines(x = time, y = state)
  } else if (rect) {
    plot_poly <- function(x, quants, time) {
      quant_names <- c((1 - x) / 2, 1 - ((1 - x) / 2))
      quants <- quants[, names(quants) %in% quant_names]
      polygon(
        x = c(time, rev(time)),
        y = c(quants[time, 1], quants[rev(time), 2]),
        col = rgb(red = 1, green = x / 2, blue = 0)
      )
    }

    quants <- lapply(f_draws, function(x, quant) {
      quant <- c((1 - quant) / 2, 1 - ((1 - quant) / 2))
      quant_val <- quantile(x, probs = quant)
      names(quant_val) <- quant
      return(quant_val)
    }, quant = quant)

    quants <- as.data.frame(do.call(rbind, quants))

    plot(x = time, y = obs, type = "n")

    for (x in quant[order(quant, decreasing = TRUE)]) {
      plot_poly(x, quants, time)
    }

    points(x = time, y = obs)

    lines(x = time, y = state)
  }
}

gp_posterior_3d <- function(
    fit, f_name, draws = 200, grid_size = 225,
    state = NULL, obs = NULL, smooth = NULL) {
  #' Function to generate 3d posterior predictive plots of a 2d gaussian
  #' process distribution using kernel smoothers.

  require(plotly)
  require(MASS)
  require(tidyr)

  f_draws <- as.data.frame(fit$draws("f", format = "draws_df"))
  f_draws <- f_draws[
    sample(seq_len(nrow(f_draws)), draws),
    seq(1, ncol(f_draws) - 3)
  ]
  time <- seq_len(ncol(f_draws))
  colnames(f_draws) <- time
  f_draws <- gather(f_draws, key = "time", value = "state")
  f_draws$time <- as.numeric(f_draws$time)

  kd <- kde2d(f_draws$time, f_draws[["state"]], n = 225)

  fig <- plot_ly(x = kd$x, y = kd$y, z = t(kd$z)) %>%
    add_surface()

  if (!is.null(state)) {
    z <- rep(max(kd$z) * 1.1, length(state))
    fig <- add_trace(fig,
      y = ~state,
      x = ~time,
      z = ~z,
      type = "scatter3d", mode = "lines",
      line = list(size = 5, color = "red")
    )
  }

  if (!is.null(smooth)) {
    z <- rep(max(kd$z) * 1.1, length(smooth))
    fig <- add_trace(fig,
      y = ~smooth,
      x = ~time,
      z = ~z,
      type = "scatter3d", mode = "lines",
      line = list(size = 5, color = "#00eeff")
    )
  }

  if (!is.null(obs)) {
    z <- rep(max(kd$z) * 1.1, length(state))
    fig <- add_markers(fig,
      y = ~obs,
      x = ~time,
      z = ~z,
      type = "scatter3d",
      marker = list(size = 3, color = "orange")
    )
  }

  fig <- layout(fig,
    showlegend = FALSE,
    title = "Posterior Gaussian Process Distribution",
    scene = list(
      xaxis = list(title = "Time"),
      yaxis = list(title = "State"),
      zaxis = list(title = "Probability Density"),
      camera = list(eye = list(x = 1.25, y = 1.25, z = 1.25))
    )
  )

  fig
}

make_exemplar_plot <- function(
    gen_model, main, xlab, ylab, cex, expr = NULL, ...) {
  #' Convenience function to generate exemplar plots for each of the methods
  dat_list <- sim_tsm(gen_model)
  dat <- get_tsm_data(dat_list)

  plot(dat$time, dat$y_obs,
    main = main, xlab = xlab, ylab = ylab,
    cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, ...
  )
  lines(dat$time, dat$y, lwd = 2.5)
  eval(expr)
}

plot_results <- function(res, outcome, style = "all", ..., legend = TRUE) {
  #' Convenience function to create discriptive plots of the performance
  #' measures.

  require(data.table)
  require(ggplot2)
  require(ggdist)

  grp_var <- sapply(list(...), function(x) match.arg(x, names(res)))
  if (length(grp_var) > 0) {
    factors <- res[, c("model", grp_var), with = FALSE]
    res[, group := apply(factors, 1, function(x) {
      paste0(names(x), " ", x, collapse = "-")
    })]
  } else {
    factors <- res[, model]
    res[, group := as.character(model)]
  }

  if (style == "all") {
    gg <- ggplot2::ggplot(res, aes(
      x = group, y = !!sym(outcome), color = model
    )) +
      ggdist::stat_halfeye(
        aes(fill = model),
        adjust = .5,
        width = .6,
        .width = 0,
        justification = -.3,
        point_colour = NA
      ) +
      geom_boxplot(
        width = .25,
        outlier.shape = NA
      ) +
      geom_point(
        size = 1.5,
        alpha = .2,
        position = position_jitter(
          seed = 1, width = .1
        )
      ) +
      facet_grid(
        cols = vars(method),
        labeller = labeller(method = label_wrap_gen(20))
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 20)
      ) +
      ylab(outcome)
  } else if (style == "mean") {
    res$outcome <- res[, c(outcome), with = FALSE]
    if (length(grp_var) > 0) {
      res_sum <- res[, .(
        mean = mean(outcome, na.rm = TRUE),
        se = sd(outcome, na.rm = TRUE) / sqrt(.N),
        group = unique(group)
      ), by = c("model", "method", grp_var)]
    } else {
      res_sum <- res[, .(
        mean = mean(outcome, na.rm = TRUE),
        se = sd(outcome, na.rm = TRUE) / sqrt(.N),
        group = unique(group)
      ), by = c("model", "method")]
    }
    gg <- ggplot2::ggplot(res_sum, aes(
      x = group, y = mean, color = model
    )) +
      geom_point(
        size = 1.5
      ) +
      geom_errorbar(
        aes(ymin = mean - se, ymax = mean + se)
      ) +
      facet_grid(
        cols = vars(method),
        labeller = labeller(method = label_wrap_gen(20))
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 20)
      ) +
      ylab(outcome)
  } else if (style == "missing") {
    res$outcome <- res[, c(outcome), with = FALSE]
    res_mis <- res[, .(
      missing = sum(is.na(outcome) / length(outcome)),
      group = unique(group)
    ), by = c("model", "method", grp_var)]

    gg <- ggplot2::ggplot(res_mis, aes(
      x = group, y = missing, fill = model
    )) +
      geom_col() +
      facet_grid(
        cols = vars(method),
        labeller = labeller(method = label_wrap_gen(20))
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 20)
      ) +
      ylab(outcome)
  }
  if (length(grp_var) > 0) {
    g <- suppressWarnings(ggplot_build(gg))
    breaks <- g$layout$panel_params[[1]]$x$breaks
    breaks_new <- sapply(breaks, function(x) {
      x <- strsplit(x, "-")[[1]]
      x <- gsub("[^0-9.-]", "", x)
      return(paste(x[seq(2, length(x))], collapse = " "))
    })
    breaks_new <- unlist(lapply(
      strwrap(breaks_new,
        width = 1,
        simplify = FALSE
      ),
      paste0,
      collapse = "\n"
    ))
    names(breaks_new) <- NULL
    gg <- gg + scale_x_discrete(labels = breaks_new)
    suppressWarnings(print(gg))
  }
  if (!legend) {
    gg <- gg + theme(legend.position = "none")
    suppressWarnings(print(gg))
  }
  return(gg)
}

inspect <- function(data, method_name, sim) {
  methods <- sim$method[[data]]
  method <- methods[
    sapply(methods, function(x) x@method_name) == method_name
  ][[1]]
  print(method)
  plot(method, sim = sim)
}

model_search <- function(fit_fun, predictors, data) {
  #' Function to conduct and exhaustative model search based on AIC and BIC
  #' model weights.

  # Create list of all substest of model terms.
  term_list <- unlist(lapply(seq_along(predictors), function(i) {
    combn(predictors, i, FUN = function(z) {
      lapply(seq_along(z), function(m) {
        combn(z, m, paste, collapse = ":")
      })
    }, simplify = FALSE)
  }), recursive = FALSE)

  # Collapse model terms to rhs of model formula
  res <- data.frame(
    predictors = unlist(lapply(
      term_list,
      function(x) paste(unlist(x), collapse = " + ")
    ))
  )

  # Apply fit_fun to rhs'
  res$model <- lapply(res$predictors, fit_fun, data = data)

  # Calculate AIC, BIC and weights
  res$AIC <- sapply(res$model, AIC)
  res$BIC <- sapply(res$model, BIC)
  res$w_AIC <- model_weight(res$AIC)
  res$w_BIC <- model_weight(res$BIC)

  # Output
  cat(
    "The model preferred by the AIC is model: ",
    which.max(res$w_AIC), ". With a weight of: ", max(res$w_AIC),
    "\n"
  )

  cat(
    "The model preferred by the BIC is model: ",
    which.max(res$w_BIC), ". With a weight of: ", max(res$w_BIC),
    "\n"
  )

  return(res)
}

model_weight <- function(ic) {
  #' Convenience function for information criterion weights

  delta_ic <- ic - min(ic)
  l_m_data <- exp(-delta_ic / 2)
  w_ic <- l_m_data / sum(l_m_data)
  return(w_ic)
}

quiet <- function(x) {
  #' Convenience function to silence functions during the simulation to keep
  #' the console clean.

  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

counter <- function() {
  i <- 0
  function() {
    i <<- i + 1
    return(i)
  }
}

wrapper <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}
