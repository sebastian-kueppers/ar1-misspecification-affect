##########################################
### null_type = "rb_refit": AR(1) residual bootstrap WITH re-estimation

###
### Procedure (one call = one series, all requested tests share the same B draws):
###  1. Fit AR(1) by ML (same arima() call as the helpers) or by CLS
###     (estimator = "CLS": lm of y_t on y_{t-1}, complete pairs only, as in
###     _analysis.R); observed statistics on the standardized residuals.
###  2. Innovations e = centered residuals (first residual dropped, NAs dropped).
###  3. For b = 1..B: draw e* (with replacement, length n + burn), rebuild the series
###     by x*_t = phi_hat x*_{t-1} + e*_t, discard burn-in, re-insert NAs at the
###     OBSERVED missing positions, REFIT AR(1) (same estimator), standardize
###     residuals, compute stats.
###  4. p = (1 + #{null >= obs}) / (B_valid + 1)   (same formula as the RP helpers;
###     for S-rho this replaces the package's k/B p-value).
###
### The statistics are NOT defined here. They come from the per-test helpers,
### which must be sourced first:
###   bds_stat()   helpers/bds.R      srho_stat()  helpers/srho.R
###   gspec_stat() helpers/gspec.R    adcf_stat()  helpers/adcf.R
###   hsic_stat()  helpers/hsic.R
### Each takes a standardized residual vector (NAs allowed) as first argument.
### Only the tests named in `tests` are computed. The statistics draw no random
### numbers, so for a given seed the bootstrap series are identical whichever
### tests are requested.
##########################################

# Registry of statistic functions. Referenced by name (not via get()) so that
# `future` detects and exports them to parallel workers.
.rb_stat_fns <- function() {
  list(bds = bds_stat, srho = srho_stat, gspec = gspec_stat,
       adcf = adcf_stat, hsic = hsic_stat)
}

# Default statistic settings; reproduce the statistics used in the validated
# refit simulation (BDS_power/power_sim_null_type.R). Override per test via
# `stat_args`, e.g. stat_args = list(adcf = list(max_lag = 3)).
.rb_stat_defaults <- function() {
  list(bds   = list(m = 2, eps = 1),
       srho  = list(max_lag = 1),
       gspec = list(na_handling = "strip"),
       adcf  = list(max_lag = 1, na_handling = "pairwise"),
       hsic  = list(max_lag = 1, na_handling = "pairwise"))
}

.rb_std <- function(r) (r - mean(r, na.rm = TRUE)) / sd(r, na.rm = TRUE)

# AR(1) fit -> list(phi, res) with res aligned to x, or NULL on failure.
# "ML": exact Gaussian ML via Kalman (arima); NA gaps are bridged by phi^k.
# "CLS": OLS of x_t on x_{t-1}; pairs containing an NA are dropped, so the
#        residual after every gap (incl. the night) is NA.
.rb_fit <- function(x, estimator) {
  if (estimator == "ML") {
    fit <- tryCatch(arima(x, order = c(1, 0, 0), method = "ML"), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    return(list(phi = unname(coef(fit)["ar1"]), res = as.numeric(residuals(fit))))
  }
  x_lag <- c(NA, x[-length(x)])
  fit <- tryCatch(lm(x ~ x_lag, na.action = na.exclude), error = function(e) NULL)
  if (is.null(fit) || !is.finite(coef(fit)["x_lag"])) return(NULL)
  list(phi = unname(coef(fit)["x_lag"]), res = as.numeric(residuals(fit)))
}

# Returns a named list (one element per requested test)
test_rb_refit <- function(x, tests = c("bds", "srho", "gspec", "adcf", "hsic"),
                          B = 99, burn = 50, seed = NULL, stat_args = list(),
                          estimator = c("ML", "CLS")) {
  estimator <- match.arg(estimator)
  all_fns <- .rb_stat_fns()
  unknown <- setdiff(tests, names(all_fns))
  if (length(unknown)) stop("Unknown test(s): ", paste(unknown, collapse = ", "))

  if (!is.null(seed)) set.seed(seed)
  x      <- as.numeric(x)
  n      <- length(x)
  na_pos <- is.na(x)
  fns    <- all_fns[tests]
  args   <- lapply(setNames(tests, tests), function(t)
    utils::modifyList(.rb_stat_defaults()[[t]],
                      if (is.null(stat_args[[t]])) list() else stat_args[[t]]))
  k      <- length(fns)
  na_res <- setNames(rep(list(list(statistic = NA_real_, p.value = NA_real_)), k), tests)
  calc   <- function(z) vapply(tests, function(t)
    tryCatch(as.numeric(do.call(fns[[t]], c(list(z), args[[t]])))[1],
             error = function(e) NA_real_),
    numeric(1))

  # 1. fit + observed statistics
  fit <- .rb_fit(x, estimator)
  if (is.null(fit)) return(na_res)
  phi <- fit$phi
  res <- fit$res
  # CLS is unconstrained: |phi| >= 1 would make the bootstrap series explode
  if (abs(phi) >= 1) return(na_res)
  # same guard as the test_* helpers (also skips B refits for degenerate series)
  if (sum(!is.na(res)) < 10 || sd(res, na.rm = TRUE) == 0) return(na_res)
  obs <- calc(.rb_std(res))

  # 2. innovations
  # (ML: drop the first residual, which is the stationary-start residual;
  #  CLS: position 1 is NA anyway)
  e <- res[-1]; e <- e[!is.na(e)]; e <- e - mean(e)

  # 3. bootstrap null with refit
  null <- vapply(seq_len(B), function(b) {
    eb <- sample(e, n + burn, replace = TRUE)
    xb <- as.numeric(stats::filter(eb, phi, method = "recursive"))[(burn + 1):(n + burn)]
    xb[na_pos] <- NA_real_
    fb <- .rb_fit(xb, estimator)
    if (is.null(fb)) return(rep(NA_real_, k))
    calc(.rb_std(fb$res))
  }, numeric(k))
  null <- matrix(null, nrow = k)

  # 4. p-values
  out <- lapply(seq_len(k), function(i) {
    nk <- null[i, ]; nk <- nk[is.finite(nk)]
    p  <- if (!is.finite(obs[i]) || length(nk) < B / 2) NA_real_
          else (1 + sum(nk >= obs[i])) / (length(nk) + 1)
    list(statistic = unname(obs[i]), p.value = p)
  })
  setNames(out, tests)
}
