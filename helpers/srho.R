#' ----------------------------------------------------------------------------#
#' Title: Srho (Entropy) Test for serial dependence                            #
#' Input: raw series (AR(1) fit) or residual series -> Srho.test.ts on the     #
#'        resulting residuals -> list(statistic, p.value)                      #
#' ----------------------------------------------------------------------------#

### Dependencies ---------------------------------------------------------------
library(tseriesEntropy)

### Srho statistic -------------------------------------------------------------
srho_stat <- function(z, max_lag = 1, bw = "reference", method = "integral") {
  #' Srho entropy measure (tseriesEntropy::Srho.ts) on a standardized residual
  #' vector, summed over lags 1..max_lag (with max_lag = 1: Srho at lag 1).
  #'
  #' Statistic only (no null), used by test_rb_refit. test_srho keeps the
  #' package's own permutation test (Srho.test.ts), which computes the same
  #' per-lag Srho with the same bw / method settings.
  #'
  #' @param z       Standardized residuals; NAs are removed (Srho.ts requires
  #'                NA-free input; former gaps become false adjacencies).
  #' @param max_lag Highest lag (default 1).
  #' @param bw      Bandwidth selector passed to Srho.ts (default "reference").
  #' @param method  Integration method passed to Srho.ts (default "integral").
  #' @return        Numeric scalar.
  z <- as.numeric(z)
  z <- z[!is.na(z)]
  s <- tseriesEntropy::Srho.ts(z, lag.max = max_lag, bw = bw,
                               method = method, plot = FALSE)
  sum(as.numeric(s@.Data)[seq_len(max_lag)])
}

### Srho test ------------------------------------------------------------------

test_srho <- function(
    ts,
    input   = c("raw", "residual"),
    max_lag = 3,
    B       = 499,
    combine = c("bonferroni", "lag1", "min"),
    seed    = NULL
) {
  #' Entropy-based (Srho) test for serial dependence in AR(1) residuals.
  #'
  #' Applies Srho.test.ts (permutation null) to a residual series obtained
  #' either from an internal AR(1) fit (input = "raw") or passed directly
  #' (input = "residual").
  #'
  #' @param input   "raw": fit AR(1) internally (Kalman handles internal NAs),
  #'                use standardized residuals. "residual": treat `ts` as an
  #'                already-residual series (standardized here).
  #' @param max_lag Maximum lag tested (default 3).
  #' @param B       Permutation replications.
  #' @param combine Reduce lag-wise results to one:
  #'                "bonferroni" (default): p = min(p) * max_lag, capped at 1
  #'                "lag1"                : lag 1 only
  #'                "min"                 : raw min(p) (anti-conservative)
  #'                The returned statistic is Srho at the SAME lag the p-value
  #'                is taken from (the minimizing lag for bonferroni/min, lag 1
  #'                for lag1), so statistic and p.value stay aligned.
  #' @param seed    Optional RNG seed.
  #'
  #' @return list(statistic, p.value); NA on failure.
  #'
  #' NA handling: NAs are always removed before the test. Srho.test.ts owns its
  #' internal embedding and null, so there is no gap-aware / pairwise variant;
  #' former overnight gaps become false adjacencies (second-order limitation).
  input   <- match.arg(input)
  combine <- match.arg(combine)
  fail <- list(statistic = NA_real_, p.value = NA_real_)
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
  ar_res_z <- ar_res_z[!is.na(ar_res_z)]        # Srho.test.ts requires NA-free input

  out <- tryCatch(
    Srho.test.ts(
      ar_res_z,
      lag.max = max_lag,
      B       = B,
      plot    = FALSE,
      bw      = "reference",
      method  = "integral"
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(fail)

  p_lags   <- as.numeric(out@p.value)   # bootstrap p-value per lag
  srho_lags <- as.numeric(out@.Data)    # Srho statistic per lag
  if (!length(p_lags) || anyNA(p_lags)) return(fail)

  # pick the lag the p-value comes from, and report Srho at that same lag
  idx <- switch(combine,
                bonferroni = which.min(p_lags),
                min        = which.min(p_lags),
                lag1       = 1L
  )
  p_val <- switch(combine,
                  bonferroni = min(1, min(p_lags) * length(p_lags)),
                  min        = min(p_lags),
                  lag1       = p_lags[1]
  )
  stat <- if (idx <= length(srho_lags)) srho_lags[idx] else NA_real_

  list(statistic = stat, p.value = p_val)
}
