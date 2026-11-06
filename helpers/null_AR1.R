#' ----------------------------------------------------------------------------#
#' Title: AR(1) Null-Condition Data Generator                                  #
#' Purpose: Simulate series where AR(1) is the TRUE model, to estimate the     #
#'          empirical size of the tests (and the ARMA(1,1) effect of adding   #
#'          measurement error) on the same n_obs grid as the main simulation.  #
#' ----------------------------------------------------------------------------#

### Dependencies ---------------------------------------------------------------
library(data.table)

### Functions ------------------------------------------------------------------

sim_ar1 <- function(n_obs, phi, dyn_er = 1, meas_er = 0, burnin = 100) {
  #' Simulate one series from a latent AR(1) plus white measurement noise.
  #'
  #'   y_t     = phi * y_{t-1} + dyn_er  * e_t,   e_t ~ N(0, 1)   (latent)
  #'   y_obs_t = y_t            + meas_er * u_t,  u_t ~ N(0, 1)   (observed)
  #'
  #' Notes:
  #'  - For y, AR(1) is exactly correct -> tests on y estimate SIZE.
  #'  - For y_obs with meas_er > 0, the true model is ARMA(1,1), so an AR(1)
  #'    fit is (mildly, linearly) misspecified -> tests on y_obs estimate
  #'    power against measurement-error-induced misspecification.
  #'  - Initialized from the stationary distribution plus a burn-in, so no
  #'    start-up transient.
  #'
  #' @param n_obs   Series length
  #' @param phi     AR(1) coefficient, |phi| < 1
  #' @param dyn_er  Innovation SD of the latent process (default 1)
  #' @param meas_er Measurement error SD (default 0 -> y_obs identical to y)
  #' @param burnin  Discarded initial observations (default 100)
  #' @return data.table with columns time, y, y_obs
  
  stopifnot(abs(phi) < 1, dyn_er >= 0, n_obs > 1)
  
  n_tot <- n_obs + burnin
  innov <- rnorm(n_tot, 0, dyn_er)
  y <- numeric(n_tot)
  # Draw start from the stationary distribution (SD = dyn_er / sqrt(1 - phi^2))
  y[1] <- rnorm(1, 0, if (dyn_er > 0) dyn_er / sqrt(1 - phi^2) else 0)
  for (t in 2:n_tot) y[t] <- phi * y[t - 1] + innov[t]
  y <- y[(burnin + 1):n_tot]
  
  data.table(
    time  = seq_len(n_obs),
    y     = y,
    y_obs = y + rnorm(n_obs, 0, meas_er)
  )
}

make_ar1_grid <- function(
    n_obs   = c(22, 43, 64, 85, 127, 169, 190, 253),
    phi     = c(0.3, 0.6),
    dyn_er  = 1,
    meas_er = c(0, 1),
    repetitions = 100,
    seed    = NULL
) {
  #' Build a simulation grid in the same shape as the main sim object:
  #' one row per replication, with a `dat` list column holding the series.
  
  if (!is.null(seed)) set.seed(seed)
  
  grid <- as.data.table(expand.grid(
    n_obs = n_obs, phi = phi, dyn_er = dyn_er, meas_er = meas_er,
    rep = seq_len(repetitions)
  ))
  grid[, sim_id := .I]
  grid[, dat := lapply(.I, function(i) {
    sim_ar1(n_obs[i], phi[i], dyn_er[i], meas_er[i])
  })]
  grid[, gen_model := "ar1_null"]
  grid[]
}