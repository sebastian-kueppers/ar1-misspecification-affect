### Dependencies ---------------------------------------------------------------
library(mgcv)   # gam(), s()

### GAM lag-1 EDF --------------------------------------------------------------
edf_gam <- function(ts_data, k_max = NULL, k_fallback = 10) {
  #' Effective degrees of freedom of a lag-1 smooth as a nonlinearity metric.
  #'
  #' Fits ts_t ~ s(ts_{t-1}) with a thin-plate spline and returns the smooth's
  #' EDF. EDF near 1 => (near-)linear lag-1 relationship; EDF > 1 => curvature,
  #' i.e. nonlinear dependence of the present on the immediate past.
  #'
  #' The lag is built on the RAW series (overnight NAs kept), then any pair
  #' with an NA in either member is dropped, so no lagged pair straddles a
  #' former overnight gap. This matches the pairwise-complete construction used
  #' by the dependence tests.
  #'
  #' @param ts_data    Numeric vector: the series WITH overnight NAs in place.
  #' @param k_max      Basis dimension to try first. Default n_gam - 1, i.e. the
  #'                   richest basis the complete-pair count allows.
  #' @param k_fallback Basis dimension to fall back to if the k_max fit errors
  #'                   (default 10).
  #'
  #' @return Numeric scalar EDF of the lag-1 smooth, or NA_real_ on failure.
  ts_df <- data.frame(
    ts     = ts_data,
    ts_lag = c(NA, ts_data[-length(ts_data)])
  )
  ts_df <- ts_df[!is.na(ts_df$ts) & !is.na(ts_df$ts_lag), , drop = FALSE]
  n_gam <- nrow(ts_df)
  if (n_gam < 4) return(NA_real_)                 # too few pairs to smooth

  if (is.null(k_max)) k_max <- n_gam - 1
  k_max <- min(k_max, n_gam - 1)                  # k cannot exceed n_gam - 1

  .fit_edf <- function(k) {
    fit <- mgcv::gam(ts ~ s(ts_lag, bs = "tp", k = k),
                     method = "ML", data = ts_df)
    unname(sum(summary(fit)$edf))                 # sum: robust if >1 smooth term
  }

  out <- tryCatch(.fit_edf(k_max), error = function(e) NULL)
  if (is.null(out)) {
    kf  <- min(k_fallback, n_gam - 1)
    out <- tryCatch(.fit_edf(kf), error = function(e) NA_real_)
  }
  out
}
