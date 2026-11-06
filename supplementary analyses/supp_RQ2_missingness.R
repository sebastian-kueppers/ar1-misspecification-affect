# ==============================================================================
# Diagnose NAs in RQ2 helper outputs (edf_gam, optimal_theta_from_SMap,
# n_changepoints)
#
# Run AFTER: packages loaded, df.list / ts_inventory built, rq2_results loaded
# (i.e. after the load("results/rq2_results_...Rdata") line).
#
# Strategy:
#   0. Tabulate NA patterns from rq2_results alone (no series needed).
#   1. Re-run diagnostic mirrors of the helpers ONLY on the NA rows, capturing
#      the stage at which they fail plus the underlying error/warning message.
#   2. Validate the mirrors on a random sample of non-NA rows (they must
#      reproduce the stored values), so the diagnosis can be trusted.
#
# The mirrors re-implement the helper logic stage by stage (the originals
# swallow errors inside tryCatch). If you edit the helpers, edit these too.
# ==============================================================================

library(tidyverse)
library(rlang)


load("results/rq2_results_2026-10-13_112850.Rdata")

# ------------------------------------------------------------------------------
# 0. NA overview from rq2_results alone ----------------------------------------
# ------------------------------------------------------------------------------

rq2_na <- rq2_results %>%
  mutate(na_edf   = is.na(edf_lag),
         na_theta = is.na(theta_opt),
         na_cp    = is.na(changepoint_n),
         n_obs_bin = cut(n_obs, breaks = c(24, 50, 75, 100, 125, 150, Inf),))

# overall NA counts per helper
rq2_na %>% summarise(across(c(na_edf, na_theta, na_cp),
                            list(n = sum, pct = ~ mean(.x) * 100)))

# na_edf_n na_edf_pct na_theta_n na_theta_pct na_cp_n na_cp_pct
# 1     1736   9.809572       1203     6.797762       0         0

# NA rate by series length (is NA missing-at-random w.r.t. n_obs?)
rq2_na %>%
  group_by(n_obs_bin) %>%
  summarise(n = n(),
            pct_na_edf   = mean(na_edf)   * 100,
            pct_na_theta = mean(na_theta) * 100,
            pct_na_cp    = mean(na_cp)    * 100)

# A tibble: 6 × 5
# n_obs_bin     n pct_na_edf pct_na_theta pct_na_cp
# <fct>     <int>      <dbl>        <dbl>     <dbl>
# 1 (24,50]    2973      16.6        20.0           0
# 2 (50,75]   11083       9.86        4.75          0
# 3 (75,100]   1824       3.73        3.23          0
# 4 (100,125]   564       4.61        0.887         0
# 5 (125,150]   289       7.27        2.42          0
# 6 (150,Inf]   964       3.53        1.04          0


# ------------------------------------------------------------------------------
# 1. Utilities -----------------------------------------------------------------
# ------------------------------------------------------------------------------

# Evaluate expr, capturing error message and warnings instead of swallowing them
capture <- function(expr) {
  warns <- character()
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e)
      structure(list(msg = conditionMessage(e)), class = ".diag_error")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  failed <- inherits(val, ".diag_error")
  list(value = if (failed) NULL else val,
       error = if (failed) val$msg else NA_character_,
       warn  = if (length(warns)) paste(unique(warns), collapse = " | ")
       else NA_character_)
}

# Same series extraction as in the RQ2 loop (centered, overnight NAs kept)
ds_names <- sapply(df.list, function(df) unique(df$dataset))

get_series <- function(dataset, user_ID, item) {
  k <- which(ds_names == dataset)
  d <- df.list[[k]] %>% filter(.data$ID == .env$user_ID)
  d %>% pull(paste0(item, "_c"))
}

# Descriptive features likely to explain failures
series_features <- function(ts_data) {
  x <- ts_data[!is.na(ts_data)]
  tibble(
    f_n_complete = length(x),
    f_n_unique   = length(unique(x)),
    f_sd         = if (length(x) > 1) sd(x) else NA_real_,
    f_prop_mode  = if (length(x)) max(table(x)) / length(x) else NA_real_
  )
}

# ------------------------------------------------------------------------------
# 2. Diagnostic mirrors --------------------------------------------------------
# ------------------------------------------------------------------------------

## 2a. edf_gam -----------------------------------------------------------------
diag_edf <- function(ts_data, k_fallback = 10) {
  ts_df <- data.frame(ts = ts_data, ts_lag = c(NA, ts_data[-length(ts_data)]))
  ts_df <- ts_df[!is.na(ts_df$ts) & !is.na(ts_df$ts_lag), , drop = FALSE]
  n_gam <- nrow(ts_df)
  n_uniq_lag <- length(unique(ts_df$ts_lag))

  out <- tibble(edf_status = NA_character_, edf_value = NA_real_,
                edf_k_used = NA_integer_, edf_n_pairs = n_gam,
                edf_n_unique_lag = n_uniq_lag,
                edf_msg_kmax = NA_character_, edf_msg_kfb = NA_character_,
                edf_warn = NA_character_)

  if (n_gam < 4) { out$edf_status <- "too_few_pairs"; return(out) }

  fit_edf <- function(k) {
    fit <- mgcv::gam(ts ~ s(ts_lag, bs = "tp", k = k),
                     method = "ML", data = ts_df)
    unname(sum(summary(fit)$edf))
  }

  k_max <- n_gam - 1
  r1 <- capture(fit_edf(k_max))
  out$edf_msg_kmax <- r1$error
  out$edf_warn     <- r1$warn
  if (is.null(r1$error) || is.na(r1$error)) {
    out$edf_status <- "ok_kmax"; out$edf_value <- r1$value
    out$edf_k_used <- as.integer(k_max)
    return(out)
  }

  kf <- min(k_fallback, n_gam - 1)
  r2 <- capture(fit_edf(kf))
  out$edf_msg_kfb <- r2$error
  if (is.na(r2$error)) {
    out$edf_status <- "ok_fallback"; out$edf_value <- r2$value
    out$edf_k_used <- as.integer(kf)
  } else {
    out$edf_status <- "fail_both"
  }
  out
}

## 2b. optimal_theta_from_SMap -------------------------------------------------
diag_theta <- function(ts_data, colname = "ts", drop_na = FALSE) {
  x <- as.numeric(ts_data)
  if (drop_na) x <- x[!is.na(x)]
  N <- length(x)
  mid <- floor(N * 0.7)

  out <- tibble(theta_status = NA_character_, theta_value = NA_real_,
                theta_E = NA_real_, theta_N = N,
                theta_n_na_x = sum(is.na(x)),
                theta_lib_nonNA  = if (mid >= 1) sum(!is.na(x[1:mid])) else 0L,
                theta_pred_nonNA = if (mid < N) sum(!is.na(x[(mid + 1):N])) else 0L,
                theta_msg = NA_character_,
                theta_would_crash_original = FALSE,
                theta_n_err = NA_integer_, theta_n_fewfinite = NA_integer_,
                theta_n_cor_na = NA_integer_,
                theta_rho_max = NA_real_, theta_rho0 = NA_real_,
                theta_gain = NA_real_, theta_rho_range = NA_real_,
                theta_n_near_max = NA_integer_)

  if (N < 4 || mid < 2) { out$theta_status <- "too_short"; return(out) }

  df <- data.frame(x); names(df) <- colname
  lib <- c(1, mid); pred <- c(mid + 1, N)

  # stage 1: embedding dimension
  r_e <- capture(EmbedDimension(dataFrame = df, lib = lib, pred = pred,
                                maxE = 15, Tp = 1,
                                columns = colname, target = colname,
                                noTime = TRUE, showPlot = FALSE))
  if (!is.na(r_e$error)) {
    out$theta_status <- "embed_error"; out$theta_msg <- r_e$error
    return(out)
  }
  emb <- r_e$value
  if (!any(is.finite(emb$rho))) {
    out$theta_status <- "embed_all_rho_na"
    # original: which.max() -> integer(0) -> if (is.na(e)) errors
    out$theta_would_crash_original <- TRUE
    out$theta_msg <- r_e$warn
    return(out)
  }
  e <- emb$E[which.max(emb$rho)]
  out$theta_E <- e

  # stage 2: S-map over theta
  thetas <- seq(0, 8, by = 0.5)
  per_theta <- lapply(thetas, function(th) {
    r <- capture(SMap(dataFrame = df, lib = lib, pred = pred,
                      columns = colname, target = colname,
                      E = e, Tp = 1, theta = th, noTime = TRUE))
    if (!is.na(r$error)) return(list(status = "error", rho = NA_real_, msg = r$error))
    obs <- r$value$predictions$Observations
    prd <- r$value$predictions$Predictions
    ok  <- is.finite(obs) & is.finite(prd)
    if (sum(ok) < 3) return(list(status = "few_finite", rho = NA_real_,
                                 msg = paste0("n_finite=", sum(ok))))
    rho <- suppressWarnings(cor(obs[ok], prd[ok]))
    if (is.na(rho)) return(list(status = "cor_na", rho = NA_real_,
                                msg = "zero variance in obs or pred"))
    list(status = "ok", rho = rho, msg = NA_character_)
  })
  st  <- vapply(per_theta, `[[`, "", "status")
  rho <- vapply(per_theta, `[[`, 0,  "rho")
  out$theta_n_err       <- sum(st == "error")
  out$theta_n_fewfinite <- sum(st == "few_finite")
  out$theta_n_cor_na    <- sum(st == "cor_na")

  if (all(is.na(rho))) {
    out$theta_status <- paste0("smap_all_failed:", names(sort(table(st), decreasing = TRUE))[1])
    msgs <- na.omit(vapply(per_theta, `[[`, "", "msg"))
    out$theta_msg <- if (length(msgs)) msgs[1] else NA_character_
    return(out)
  }
  rho_fin <- rho[is.finite(rho)]
  out$theta_rho_max    <- max(rho_fin)
  out$theta_rho0       <- rho[1]                    # theta = 0 (global linear)
  out$theta_gain       <- max(rho_fin) - rho[1]     # NA if rho[1] not finite
  out$theta_rho_range  <- diff(range(rho_fin))      # ~1e-17 => flat profile
  out$theta_n_near_max <- sum(rho_fin >= max(rho_fin) - 1e-6)
  out$theta_value  <- thetas[which.max(rho)]
  out$theta_status <- "ok"
  out
}

# ------------------------------------------------------------------------------
# 3. Validate mirrors on non-NA rows -------------------------------------------
# ------------------------------------------------------------------------------

set.seed(1)
val_rows <- rq2_results %>%
  filter(!is.na(edf_lag), !is.na(theta_opt)) %>%
  slice_sample(n = 30)

val <- pmap_dfr(val_rows, function(ID, dataset, user_ID, item, edf_lag,
                                   theta_opt, ...) {
  ts_data <- get_series(dataset, user_ID, item)
  bind_cols(tibble(ID = ID, stored_edf = edf_lag, stored_theta = theta_opt),
            diag_edf(ts_data), diag_theta(ts_data))
})

# should both be TRUE; if not, the mirror has drifted from the helper
all(abs(val$edf_value - val$stored_edf) < 1e-6)
all(val$theta_value == val$stored_theta)

# ------------------------------------------------------------------------------
# 4. Run diagnostics on NA rows only -------------------------------------------
# ------------------------------------------------------------------------------

run_diag <- function(rows, diag_fun, use_clean = FALSE) {
  pmap_dfr(rows, function(ID, dataset, user_ID, item, ...) {
    ts_data <- get_series(dataset, user_ID, item)
    x <- if (use_clean) ts_data[!is.na(ts_data)] else ts_data
    bind_cols(tibble(ID = ID), series_features(ts_data), diag_fun(x))
  })
}

# n_obs / n_na come from rq2_results; keep them for the summaries
rows_edf   <- rq2_na %>% filter(na_edf)
rows_theta <- rq2_na %>% filter(na_theta)

diag_edf_res   <- run_diag(rows_edf,   diag_edf) %>%
  left_join(dplyr::select(rq2_results, ID, n_obs, n_na), by = "ID")
diag_theta_res <- run_diag(rows_theta, diag_theta) %>%
  left_join(dplyr::select(rq2_results, ID, n_obs, n_na), by = "ID")

# ------------------------------------------------------------------------------
# 5. Summaries -----------------------------------------------------------------
# ------------------------------------------------------------------------------

# Reason codes per helper
count(diag_edf_res,   edf_status)

count(diag_theta_res, theta_status)
# # A tibble: 2 × 2
# theta_status                   n
# <chr>                      <int>
# 1 smap_all_failed:cor_na       754
# 2 smap_all_failed:few_finite   449


# Underlying messages
diag_edf_res   %>% count(edf_status, edf_msg_kmax, edf_msg_kfb, sort = TRUE) %>% print(n = 20, width = Inf)
# # A tibble: 1 × 4
# edf_status edf_msg_kmax
# <chr>      <chr>
#   1 fail_both  Ein Term hat weniger einzigartige Kombinationen von Kovariaten als maximal angegebene Freiheitsgrade
# edf_msg_kfb                                                                                              n
# <chr>                                                                                                <int>
#   1 Ein Term hat weniger einzigartige Kombinationen von Kovariaten als maximal angegebene Freiheitsgrade  1736

diag_theta_res %>% count(theta_status, theta_msg, sort = TRUE)                 %>% print(n = 20, width = Inf)
# # A tibble: 5 × 3
# theta_status               theta_msg                        n
# <chr>                      <chr>                        <int>
#   1 smap_all_failed:cor_na     zero variance in obs or pred   754
# 2 smap_all_failed:few_finite n_finite=1                     163
# 3 smap_all_failed:few_finite n_finite=2                     144
# 4 smap_all_failed:few_finite n_finite=0                     133
# 5 smap_all_failed:few_finite zero variance in obs or pred     9


# ------------------------------------------------------------------------------
# 6. Number of thin splices used in GAM ----------------------------------------
# ------------------------------------------------------------------------------

diag_edf_all <- run_diag(rq2_na, diag_edf) %>%
  left_join(dplyr::select(rq2_na, ID, n_obs_bin, affect_type), by = "ID")


# share of all series by which fit was used
diag_edf_all %>%
  count(edf_status) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# A tibble: 3 × 3
# edf_status      n   pct
# <chr>       <int> <dbl>
#   1 fail_both    1736   9.8
# 2 ok_fallback 15826  89.4
# 3 ok_kmax       135   0.8

# the same by length bin
diag_edf_all %>%
  count(n_obs_bin, edf_status) %>%
  group_by(n_obs_bin) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

# A tibble: 13 × 4
# n_obs_bin edf_status      n   pct
# <fct>     <chr>       <int> <dbl>
#   1 (24,50]   fail_both     494  16.6
# 2 (24,50]   ok_fallback  2344  78.8
# 3 (24,50]   ok_kmax       135   4.5
# 4 (50,75]   fail_both    1093   9.9
# 5 (50,75]   ok_fallback  9990  90.1
# 6 (75,100]  fail_both      68   3.7
# 7 (75,100]  ok_fallback  1756  96.3
# 8 (100,125] fail_both      26   4.6
# 9 (100,125] ok_fallback   538  95.4
# 10 (125,150] fail_both      21   7.3
# 11 (125,150] ok_fallback   268  92.7
# 12 (150,Inf] fail_both      34   3.5
# 13 (150,Inf] ok_fallback   930  96.5

# basis dimension actually used (the fallback is min(10, pairs - 1), so check it was 10)
table(diag_edf_all$edf_k_used, useNA = "ifany")

# 8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    28    33  <NA>
#   4    12 15829    22    23     3     9     6     3     4    14     8     3     2     4     7     2     3     1     1     1  1736
