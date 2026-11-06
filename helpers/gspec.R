#' ----------------------------------------------------------------------------#
#' Title: Generalized Spectral Test (Hong & Lee, 2003)                         #
#' Input: raw series (AR(1) fit) or residual series -> M(p) test               #
#' Output: list(statistic, p.value); NA on failure                             #
#'                                                                             #
#' M(p) statistic of Hong & Lee (2003, Econometric Theory 19, 1065-1121),     #
#' Eq. (2.9). Because T is small in ESM data, inference is a permutation null  #
#' rather than the asymptotic N(0,1). The permutation calibration is what      #
#' licenses the NA-aware ("pairwise") statistic below: exact finite-sample     #
#' centering is not required, only that the same statistic is applied to       #
#' observed and permuted data with NA positions held fixed.                    #
#' ----------------------------------------------------------------------------#

### Dependencies: base R only --------------------------------------------------

.gauss_hermite <- function(K) {
  #' Golub-Welsch Gauss-Hermite nodes/weights (weight exp(-x^2)); base R.
  i <- seq_len(K - 1)
  b <- sqrt(i / 2)
  J <- matrix(0, K, K)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b
  ev <- eigen(J, symmetric = TRUE)
  list(nodes = ev$values, weights = sqrt(pi) * ev$vectors[1, ]^2)
}

.daniell <- function(z) ifelse(z == 0, 1, sin(pi * z) / (pi * z))

gspec_stat <- function(e, p = NULL, K = 15, min_pairs = 2,
                       na_handling = c("pairwise", "strip")) {
  #' Hong-Lee M(p) statistic on a (standardized) residual vector e.
  #' NA-aware: rows with NA e are excluded from lag-0 quantities, and each
  #' lag j uses only complete pairs (both endpoints non-NA). With no NAs this
  #' reduces exactly to the original statistic.
  #'
  #' Statistic only (no null). Used by test_gspec (permutation null) and by
  #' test_rb_refit (refit bootstrap null), so both nulls share one definition.
  #'
  #' @param e           Standardized residuals; may contain NAs.
  #' @param p           Bandwidth / effective lag order. NULL (default):
  #'                    max(2, round(1.5 * n^0.2)), n = number of non-NA e.
  #' @param K           Gauss-Hermite quadrature nodes (default 15).
  #' @param min_pairs   Skip a lag with fewer complete pairs than this.
  #' @param na_handling "pairwise" (default): NA-aware statistic as above.
  #'                    "strip": drop NAs first (former gaps become false
  #'                    adjacencies).
  #' @return            Numeric scalar M(p), or NA_real_ on degeneracy.
  na_handling <- match.arg(na_handling)
  e <- as.numeric(e)
  if (na_handling == "strip") e <- e[!is.na(e)]
  if (is.null(p)) p <- max(2, round(1.5 * sum(!is.na(e))^0.2))

  n_all <- length(e)
  na    <- is.na(e)
  idx   <- which(!na)
  n0    <- length(idx)
  if (n0 < 3) return(NA_real_)

  gh <- .gauss_hermite(K)
  u  <- sqrt(2) * gh$nodes
  w  <- gh$weights / sqrt(pi)
  W2 <- outer(w, w)

  E <- exp(1i * outer(e, u))          # rows with NA e are NA

  # --- Lag-0 quantities over non-NA rows ---
  E0       <- E[idx, , drop = FALSE]
  phi_marg <- colMeans(E0)
  Phi0     <- (t(E0) %*% E0) / n0
  s0       <- Phi0 - outer(phi_marg, phi_marg)
  C0 <- (sum(w * (1 - Mod(phi_marg)^2)))^2
  D0 <- (sum(W2 * Mod(s0)^2))^2
  if (!is.finite(C0) || !is.finite(D0) || D0 <= 0) return(NA_real_)

  # --- Sum over lags, complete pairs only ---
  num <- 0; S2 <- 0; S4 <- 0
  for (j in seq_len(n_all - 1)) {
    mask <- !na[(j + 1):n_all] & !na[1:(n_all - j)]
    nj   <- sum(mask)
    if (nj < min_pairs) next
    kj <- .daniell(j / p); k2 <- kj^2
    E1 <- E[(j + 1):n_all, , drop = FALSE][mask, , drop = FALSE]
    E2 <- E[1:(n_all - j), , drop = FALSE][mask, , drop = FALSE]
    Phij <- (t(E1) %*% E2) / nj
    sj   <- Phij - outer(colMeans(E1), colMeans(E2))
    Ij   <- sum(W2 * Mod(sj)^2)
    num  <- num + k2 * nj * Ij
    S2   <- S2 + k2
    if (j <= n_all - 2) S4 <- S4 + kj^4
  }
  if (S4 <= 0) return(NA_real_)

  (num - C0 * S2) / sqrt(2 * D0 * S4)
}

### Test -----------------------------------------------------------------------
test_gspec <- function(
    ts,
    input       = c("raw", "residual"),
    na_handling = c("strip", "pairwise"),
    p           = NULL,
    K           = 15,
    B           = 499,
    min_pairs   = 2,
    asymptotic  = FALSE,
    seed        = NULL
) {
  #' Generalized spectral (Hong-Lee) test for serial dependence.
  #'
  #' @param input       "raw": fit AR(1) internally (Kalman handles internal
  #'                    NAs), use standardized residuals. "residual": treat
  #'                    `ts` as already-residual (standardized here).
  #' @param na_handling "strip": drop all NAs before the statistic.
  #'                    "pairwise": keep NAs; the M(p) statistic uses only
  #'                    complete lagged pairs, and the permutation holds NA
  #'                    positions fixed (identical complete-pair structure in
  #'                    observed and null).
  #' @param p           Bandwidth. Default max(2, round(1.5 * n^0.2)), where n
  #'                    is the number of non-NA residuals. Vary (2, 3, 5) as a
  #'                    sensitivity check.
  #' @param K           Gauss-Hermite nodes (default 15).
  #' @param B           Permutation replications (default 499).
  #' @param min_pairs   Pairwise only: skip a lag with fewer complete pairs.
  #' @param asymptotic  If TRUE, skip permutation and return 1 - pnorm(M).
  #'                    N(0,1) approx unverified at ESM lengths; check size
  #'                    first. Ignored implication under pairwise: the masked
  #'                    statistic is not exactly N(0,1) even asymptotically.
  #' @param seed        Optional RNG seed.
  #'
  #' @return list(statistic, p.value); upper-tail. NA on failure. (asymptotic
  #'         = TRUE returns a bare numeric p-value, matching prior behaviour.)
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
  e <- (res - mean(res, na.rm = TRUE)) / sd(res, na.rm = TRUE)
  if (na_handling == "strip") e <- e[!is.na(e)]

  n <- sum(!is.na(e))                          # effective length for bandwidth
  if (is.null(p)) p <- max(2, round(1.5 * n^0.2))

  M_obs <- gspec_stat(e, p = p, K = K, min_pairs = min_pairs)
  if (!is.finite(M_obs)) return(fail)

  if (asymptotic) return(1 - pnorm(M_obs))

  na_mask <- is.na(e)
  obs_idx <- which(!na_mask)                    # NA positions held fixed
  M_boot <- replicate(B, {
    sh <- e; sh[obs_idx] <- sample(e[obs_idx])
    gspec_stat(sh, p = p, K = K, min_pairs = min_pairs)
  })
  M_boot <- M_boot[is.finite(M_boot)]
  if (length(M_boot) < B / 2) return(list(statistic = M_obs, p.value = NA_real_))

  list(statistic = M_obs,
       p.value   = (1 + sum(M_boot >= M_obs)) / (length(M_boot) + 1))
}
