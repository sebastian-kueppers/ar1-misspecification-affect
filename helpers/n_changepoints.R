### Dependencies ---------------------------------------------------------------
library(ecp)   # e.divisive()

### Number of change points ----------------------------------------------------
n_changepoints <- function(ts_clean, alpha = 0.05, R = 400, min.size = 30) {
  #' Count regime shifts via E-divisive change-point estimation.
  #'
  #' Runs ecp::e.divisive on the series and counts the estimated change points
  #' whose permutation p-value is below alpha, as a regime-shift metric.
  #'
  #' @param ts_clean Numeric vector, NA-FREE. e.divisive does not accept NAs;
  #'                 pass the series with missing values already removed.
  #' @param alpha    Significance level for accepting a change point (default
  #'                 0.05). Used both as e.divisive's sig.lvl and as the cutoff
  #'                 on the returned p-values.
  #' @param R        Permutations per candidate split (default 400).
  #' @param min.size minimum segment size
  #'
  #' @return Integer count of significant change points, or NA_integer_ on
  #'         failure. e.divisive always reports the series endpoints in its
  #'         estimates; this counts only interior splits with p < alpha, so a
  #'         series with no regime shift returns 0, not a spurious count.
  if (anyNA(ts_clean)) ts_clean <- ts_clean[!is.na(ts_clean)]
  if (length(ts_clean) < 2) return(NA_integer_)

  out <- tryCatch(
    ecp::e.divisive(matrix(ts_clean, ncol = 1),
                    R = R,
                    sig.lvl = alpha,
                    min.size = min.size),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_integer_)

  # p.values has one entry per split e.divisive attempted; the last is the
  # non-significant stop and is typically NA. Count strictly p < alpha.
  pv <- out$p.values
  as.integer(sum(pv < alpha, na.rm = TRUE))
}
