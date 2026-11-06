### Dependencies ---------------------------------------------------------------
# tseries::bds.test

### BDS statistic --------------------------------------------------------------
bds_stat <- function(z, m = 2, eps = 1) {
  #' BDS statistic (tseries::bds.test) on a standardized residual vector.
  #'
  #' Statistic only (no null). Used by test_bds (permutation null) and by
  #' test_rb_refit (refit bootstrap null), so both nulls share one definition.
  #'
  #' @param z   Standardized residuals; NAs are removed (bds.test requires
  #'            NA-free input; former gaps become false adjacencies).
  #' @param m   Embedding dimension (default 2).
  #' @param eps Epsilon, in SDs of the standardized residuals (default 1).
  #' @return    Numeric scalar.
  z <- as.numeric(z)
  z <- z[!is.na(z)]
  tseries::bds.test(z, m = m, eps = eps)$statistic[[1]]
}

### BDS test -------------------------------------------------------------------

test_bds <- function(
    ts,
    input = c("raw", "residual"),
    m     = 2,
    eps   = 1,
    B     = 499,
    seed  = NULL
) {
  #' Permutation-based BDS test for serial dependence.
  #'
  #' Fits AR(1) (input = "raw") or takes residuals directly (input =
  #' "residual"), standardizes, and tests the residuals for departures from
  #' i.i.d. using tseries::bds.test with a permutation null.
  #'
  #' @param input "raw": fit AR(1) internally (Kalman handles internal NAs),
  #'              use standardized residuals. "residual": treat `ts` as an
  #'              already-residual series (standardized here).
  #' @param m     Embedding dimension (default 2).
  #' @param eps   Epsilon bandwidth, in SDs of the standardized residuals
  #'              (default 1).
  #' @param B     Permutation replications (default 499).
  #' @param seed  Optional RNG seed.
  #'
  #' @return list(statistic, p.value); upper-tail p-value. NA on failure.
  #'
  #' NA handling: NAs are always removed from the residuals before the test,
  #' because tseries::bds.test errors on any NA and embeds the m-histories
  #' internally (no seam for a gap-aware / pairwise variant). Former overnight
  #' gaps therefore become false adjacencies in the m-histories; this is a
  #' known, second-order limitation to note in the write-up rather than fix in
  #' code. Use test_adcf / test_hsic / test_gspec for a gap-faithful test.
  input <- match.arg(input)
  fail  <- list(statistic = NA_real_, p.value = NA_real_)
  if (!is.null(seed)) set.seed(seed)

  ## --- residual preparation (inlined) ---
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
  ar_res_z <- ar_res_z[!is.na(ar_res_z)]        # bds.test requires NA-free input

  # Observed BDS statistic on actual residuals
  obs <- tryCatch(bds_stat(ar_res_z, m = m, eps = eps),
                  error = function(e) NA_real_)
  if (is.na(obs)) return(fail)

  # Permutation null: shuffle residuals B times and recompute statistic
  boot_stats <- replicate(B, {
    tryCatch(bds_stat(sample(ar_res_z, replace = FALSE), m = m, eps = eps),
             error = function(e) NA_real_)
  })
  boot_stats <- boot_stats[!is.na(boot_stats)]
  if (length(boot_stats) < B / 2) return(list(statistic = obs, p.value = NA_real_))

  # Upper-tail: residual dependence inflates the BDS statistic
  p_val <- (1 + sum(boot_stats >= obs)) / (length(boot_stats) + 1)
  list(statistic = obs, p.value = p_val)
}
