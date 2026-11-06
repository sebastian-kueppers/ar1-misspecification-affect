### Dependencies ---------------------------------------------------------------
library(dHSIC)    # HSIC statistic

### HSIC statistic -------------------------------------------------------------
hsic_stat <- function(
    z,
    max_lag     = 3,
    na_handling = c("strip", "pairwise"),
    min_pairs   = 8
) {
  #' HSIC portmanteau on a standardized residual vector: Gaussian-kernel HSIC
  #' between z_t and z_{t-j}, summed over lags 1..max_lag.
  #'
  #' Statistic only (no null). Used by test_hsic (permutation null) and by
  #' test_rb_refit (refit bootstrap null), so both nulls share one definition.
  #'
  #' @param z           Standardized residuals; may contain NAs.
  #' @param max_lag     Highest lag entering the portmanteau.
  #' @param na_handling "strip": drop all NAs before forming lagged pairs.
  #'                    "pairwise": keep NAs; per lag use only complete pairs,
  #'                    so no pair straddles a former gap.
  #' @param min_pairs   Pairwise only: skip a lag with fewer complete pairs.
  #'
  #' @return Numeric scalar, or NA_real_ if no lag is usable.
  na_handling <- match.arg(na_handling)
  z <- as.numeric(z)

  # Gaussian kernel, bandwidth 1 (residuals standardized => unit scale;
  # permutation-invariant per marginal).
  .hsic_pair <- function(a, b)
    dHSIC::dhsic(list(a, b), kernel = "gaussian")$dHSIC

  if (na_handling == "strip") {
    z  <- z[!is.na(z)]
    nx <- length(z)
    if (nx <= max_lag + 2) return(NA_real_)
    return(sum(vapply(seq_len(max_lag), function(j)
      .hsic_pair(z[1:(nx - j)], z[(1 + j):nx]), numeric(1))))
  }

  nx      <- length(z)
  na_mask <- is.na(z)
  lag_masks <- lapply(seq_len(max_lag), function(j) {
    if (nx - j < 1) return(logical(0))
    !na_mask[1:(nx - j)] & !na_mask[(1 + j):nx]
  })
  valid_lags <- which(vapply(lag_masks,
                             function(m) sum(m) >= min_pairs, logical(1)))
  if (length(valid_lags) == 0) return(NA_real_)
  sum(vapply(valid_lags, function(j) {
    m <- lag_masks[[j]]
    .hsic_pair(z[1:(nx - j)][m], z[(1 + j):nx][m])
  }, numeric(1)))
}

### HSIC test ------------------------------------------------------------------
test_hsic <- function(
    ts,
    input       = c("raw", "residual"),
    na_handling = c("strip", "pairwise"),
    max_lag     = 3,
    B           = 499,
    min_pairs   = 8,
    seed        = NULL
) {
  #' Permutation-based HSIC test for serial dependence.
  #'
  #' As test_adcf, but the statistic is hsic_stat(), calibrated by a
  #' residual-permutation null.
  #'
  #' @param input       "raw": fit AR(1) internally (Kalman handles internal
  #'                    NAs), use standardized residuals. "residual": treat
  #'                    `ts` as already-residual (standardized here).
  #' @param na_handling See hsic_stat(). Under "pairwise" the permutation holds
  #'                    NA positions fixed (identical complete-pair set in
  #'                    observed and null).
  #' @param max_lag     Highest lag entering the portmanteau (default 3).
  #' @param B           Permutation replications (default 499).
  #' @param min_pairs   Pairwise only: skip a lag with fewer complete pairs.
  #' @param seed        Optional RNG seed.
  #'
  #' @return list(statistic, p.value); upper-tail p-value. NA on failure.
  input       <- match.arg(input)
  na_handling <- match.arg(na_handling)
  fail <- list(statistic = NA_real_, p.value = NA_real_)
  if (!is.null(seed)) set.seed(seed)

  ## --- residual preparation (inlined; keeps NAs iff pairwise) ---
  if (input == "raw") {
    fit <- tryCatch(arima(ts, order = c(1, 0, 0), method = "ML"),
                    error = function(e) NULL)
    if (is.null(fit)) return(fail)
    res <- as.numeric(residuals(fit))
  } else {
    res <- as.numeric(ts)
  }
  if (sum(!is.na(res)) < 10 || sd(res, na.rm = TRUE) == 0) return(fail)
  ar_res_z <- (res - mean(res, na.rm = TRUE)) / sd(res, na.rm = TRUE)
  if (na_handling == "strip") ar_res_z <- ar_res_z[!is.na(ar_res_z)]

  .stat <- function(x) hsic_stat(x, max_lag = max_lag,
                                 na_handling = na_handling, min_pairs = min_pairs)
  obs <- tryCatch(.stat(ar_res_z), error = function(e) NA_real_)
  if (is.na(obs)) return(fail)

  if (na_handling == "strip") {
    boot <- replicate(B, tryCatch(.stat(sample(ar_res_z)),
                                  error = function(e) NA_real_))
  } else {
    obs_idx <- which(!is.na(ar_res_z))               # NA positions held fixed
    boot <- replicate(B, {
      sh <- ar_res_z; sh[obs_idx] <- sample(ar_res_z[obs_idx])
      tryCatch(.stat(sh), error = function(e) NA_real_)
    })
  }

  boot <- boot[!is.na(boot)]
  if (length(boot) < B / 2) return(list(statistic = obs, p.value = NA_real_))
  list(statistic = obs,
       p.value   = (1 + sum(boot >= obs)) / (length(boot) + 1))
}
