# ==============================================================================
# Supplementary analyses RQ1: misspecification of AR(2) models
#
# Same procedure as RQ1 in _analysis.R, with an AR(2) instead of an AR(1)
# null model:
#   y_t = c + phi1 * y_{t-1} + phi2 * y_{t-2} + e_t
# fitted by CLS (OLS on complete lag sets; no lag spans the night or a missed
# beep). Misspecification is tested on the standardized AR(2) residuals
# against a refit bootstrap null that simulates from the fitted AR(2)
# (test_rb_refit_model(), section 5).
#
# Contents:
#   1. Data (as in _analysis.R)
#   2. Time-series inventory (as in _analysis.R)
#   3. Parameters
#   4. Null model (the only section that differs between the supp_RQ1_* scripts)
#   5. Run tests
#   6. Inspect results, incl. comparison with the AR(1) results of the main text
# ==============================================================================

library(tidyverse)
library(tseries)
library(future)
library(furrr)
library(progressr)
library(rlang)



# ==============================================================================
# 1 Data
# ==============================================================================

# load preprocessed data
load("data/processed/data_processed.Rdata")

# load input sheet for item info
inputsESM <- readxl::read_xlsx("inputs/input_file_ESM.xlsx")

# load helpers (as in _analysis.R)
source("helpers/adcf.R")
source("helpers/bds.R")
source("helpers/gspec.R")
source("helpers/hsic.R")
source("helpers/srho.R")
source("helpers/rb_refit.R")   # statistic registry and defaults, used by test_rb_refit_model() below

# AR(1) results of the main text, for comparison in section 6
rq1_results_main_file <- "results/rq1_results_2026-10-22_125950.Rdata"


# ==============================================================================
# 2 Time-series inventory (copied from _analysis.R)
# ==============================================================================

ts_inventory <- lapply(1:length(df.list), function(k){
  study <- unique(df.list[[k]]$dataset)
  inpstudy <- inputsESM %>%
    filter(dataset_presentation == study) %>%
    filter(!is.na(construct))
  # for each item in this dataset
  lapply(1:nrow(inpstudy), function(j){
    item        <- inpstudy$variable[j]
    construct   <- inpstudy$construct[j]
    affect_type <- inpstudy$affect_type[j]
    item_c      <- paste0(item, "_c")
    # for each person
    df.list[[k]] %>%
      group_by(user_ID = ID, UUID) %>%
      summarise(
        n_obs   = sum(!is.na(.data[[item_c]])),
        n_beeps = n(),
        .groups = "drop"
      ) %>%
      mutate(
        dataset     = study,
        item        = item,
        construct   = construct,
        affect_type = affect_type
      )
  }) %>% bind_rows()
}) %>% bind_rows()

ts_inventory <- ts_inventory %>%
  dplyr::select(dataset, user_ID, UUID, item, construct, affect_type, n_obs, n_beeps) %>%
  arrange(dataset, user_ID, item) %>%
  mutate(ID = row_number()) %>% # unique identifier, identical to _analysis.R
  dplyr::select(ID, dataset, user_ID, UUID, item, construct, affect_type, n_obs, n_beeps)

nrow(ts_inventory)  # person x item series, 17697


# ==============================================================================
# 3 Parameters
# ==============================================================================

# Which tests should be applied?
RQ1.tests <- c(
  "adcf",
  "gspec",
  "hsic",
  "srho",
  "bds"
)

# Should NAs in residuals be deleted pairwise? (FALSE = list-wise deletion)
# BDS and Srho test do not offer pairwise deletion
PAIRWISE_NA_DELETION <- TRUE
na_handling <- if (PAIRWISE_NA_DELETION) "pairwise" else "strip"

# Statistic settings (as in _analysis.R; lag-1 dependence of the residuals)
RQ1.stat_args <- list(
  adcf  = list(max_lag = 1, na_handling = na_handling),
  hsic  = list(max_lag = 1, na_handling = na_handling),
  gspec = list(na_handling = na_handling)
)

# Number of bootstrap replicates (B + 1 a multiple of 20, see _analysis.R)
# B <- 499
B <- 199

# Run on a random subsample of the series? (share in (0, 1]; 1 = all).
# row_id is assigned before subsampling, so each series keeps its bootstrap
# seed (seed = row_id), and row_id is the same as in _analysis.R.
RQ1.subsample      <- 1
RQ1.subsample_seed <- 2026

# label used in the results file name
model_label <- "AR2"

stopifnot(all(RQ1.tests %in% c("adcf", "gspec", "hsic", "srho", "bds")))


# ==============================================================================
# 4 Null model: AR(2), CLS
# ==============================================================================

# Series-specific inputs for fit/simulation; none for the AR(2).
# d = the person's rows of df.list (same rows and order as the series).
null_aux <- function(d) NULL

# CLS fit: OLS of y_t on y_{t-1} and y_{t-2}. Triplets containing an NA are
# dropped, so the first two observed beeps after a gap (incl. the night) are
# not outcomes and their residuals are NA.
null_fit <- function(x, aux) {
  n  <- length(x)
  l1 <- c(NA, x[-n])
  l2 <- c(NA, NA, x[-c(n - 1, n)])
  fit <- lm(x ~ l1 + l2, na.action = na.exclude)
  phi <- unname(coef(fit)[c("l1", "l2")])
  if (any(!is.finite(phi))) return(NULL)
  # stationary iff all roots of 1 - phi1 z - phi2 z^2 lie outside the unit circle
  list(par = phi, res = as.numeric(residuals(fit)),
       stationary = all(Mod(polyroot(c(1, -phi))) > 1))
}

# Simulation: x_t = phi1 x_{t-1} + phi2 x_{t-2} + e_t (the intercept is
# irrelevant because the refit estimates one)
null_sim <- function(par, e, n, burn, aux) {
  eb <- sample(e, n + burn, replace = TRUE)
  as.numeric(stats::filter(eb, par, method = "recursive"))[(burn + 1):(n + burn)]
}

# Descriptives of the observed fit, stored with the results
null_info <- function(x, aux) {
  n  <- length(x)
  l1 <- c(NA, x[-n])
  l2 <- c(NA, NA, x[-c(n - 1, n)])
  fit <- tryCatch(lm(x ~ l1 + l2, na.action = na.exclude), error = function(e) NULL)
  if (is.null(fit)) return(data.frame(phi1 = NA_real_, phi2 = NA_real_, phi2_p = NA_real_,
                                      sigma2 = NA_real_, n_res = NA_integer_))
  ct <- coef(summary(fit))
  data.frame(phi1   = unname(coef(fit)["l1"]),
             phi2   = unname(coef(fit)["l2"]),
             phi2_p = if ("l2" %in% rownames(ct)) ct["l2", "Pr(>|t|)"] else NA_real_,
             sigma2 = sigma(fit)^2,
             n_res  = sum(!is.na(residuals(fit))))
}


# ==============================================================================
# 5 Run tests
# ==============================================================================

# Refit bootstrap for the null model of section 4 (same procedure as
# test_rb_refit() in helpers/rb_refit.R, with the null model passed in):
#   fit_fn(x, aux) -> list(par, res, stationary) or NULL   (null_fit)
#   sim_fn(par, e, n, burn, aux) -> series of length n     (null_sim)
#  1. Fit null model; observed statistics on the standardized residuals.
#  2. Innovations e = centered non-NA residuals.
#  3. For b = 1..B: simulate from the fitted model, re-insert NAs at the
#     OBSERVED missing positions, REFIT, standardize residuals, compute stats.
#  4. p = (1 + #{null >= obs}) / (B_valid + 1)
# Uses .rb_stat_fns(), .rb_stat_defaults() and .rb_std() from helpers/rb_refit.R.
test_rb_refit_model <- function(x, fit_fn, sim_fn, aux = NULL,
                                tests = c("bds", "srho", "gspec", "adcf", "hsic"),
                                B = 99, burn = 50, seed = NULL, stat_args = list()) {
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
  fit <- tryCatch(fit_fn(x, aux), error = function(e) NULL)
  if (is.null(fit)) return(na_res)
  # CLS is unconstrained: a non-stationary fit would make the bootstrap explode
  if (!isTRUE(fit$stationary)) return(na_res)
  res <- fit$res
  # same guard as rb_refit.R
  if (sum(!is.na(res)) < 10 || sd(res, na.rm = TRUE) == 0) return(na_res)
  obs <- calc(.rb_std(res))

  # 2. innovations
  e <- res[!is.na(res)]; e <- e - mean(e)

  # 3. bootstrap null with refit
  null <- vapply(seq_len(B), function(b) {
    xb <- sim_fn(fit$par, e, n, burn, aux)
    xb[na_pos] <- NA_real_
    fb <- tryCatch(fit_fn(xb, aux), error = function(e) NULL)
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

# normalize any test return to c(stat, p) (as in _analysis.R)
.sp <- function(res) {
  if (is.null(res)) return(c(NA_real_, NA_real_))
  if (is.list(res)) return(c(res$statistic %||% NA_real_, res$p.value %||% NA_real_))
  c(NA_real_, as.numeric(res))
}

supp_inventory <- ts_inventory %>%
  filter(n_obs >= 25) %>%
  mutate(row_id = row_number())
if (RQ1.subsample < 1) {
  set.seed(RQ1.subsample_seed)
  supp_inventory <- supp_inventory %>% slice_sample(prop = RQ1.subsample)
}
nrow(supp_inventory)

# set up parallel backend
n_cores <- parallelly::availableCores(omit = 1)
plan(multisession, workers = n_cores)

set.seed(42)
start <- Sys.time()
with_progress({
  p <- progressor(steps = nrow(supp_inventory))
  supp_results <- supp_inventory %>%
    furrr::future_pmap(
      .options = furrr_options(seed = TRUE),
      function(dataset, ID, user_ID, UUID, item, construct, affect_type, n_obs, n_beeps, row_id){
        p()
        k  <- which(sapply(df.list, function(df) unique(df$dataset)) == dataset)
        d  <- df.list[[k]] %>% filter(.data$ID == .env$user_ID)
        ts <- pull(d, paste0(item, "_c"))
        stopifnot(length(ts) == n_beeps)

        aux <- null_aux(d)

        # every test gets a slot; NA if not requested or fit failed
        R <- list(adcf  = c(NA_real_, NA_real_), gspec = c(NA_real_, NA_real_),
                  hsic  = c(NA_real_, NA_real_), srho  = c(NA_real_, NA_real_),
                  bds   = c(NA_real_, NA_real_))

        rb <- tryCatch(test_rb_refit_model(ts, fit_fn = null_fit, sim_fn = null_sim,
                                           aux = aux, tests = RQ1.tests, B = B,
                                           seed = row_id, stat_args = RQ1.stat_args),
                       error = function(e) NULL)
        if (!is.null(rb)) for (t in names(rb)) R[[t]] <- .sp(rb[[t]])

        data.frame(
          ID, dataset, UUID, item, construct, affect_type,
          n_obs, n_na = sum(is.na(ts)), row_id,
          null_info(ts, aux),
          adcf_stat  = R$adcf[1],  adcf_p  = R$adcf[2],
          gspec_stat = R$gspec[1], gspec_p = R$gspec[2],
          hsic_stat  = R$hsic[1],  hsic_p  = R$hsic[2],
          srho_stat  = R$srho[1],  srho_p  = R$srho[2],
          bds_stat   = R$bds[1],   bds_p   = R$bds[2],
          stringsAsFactors = FALSE
        )
      }
    ) %>% bind_rows()
})
Sys.time() - start

# reset to sequential
plan(sequential)

# save results
ts_stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
sub_tag  <- if (RQ1.subsample < 1) paste0("_sub", round(100 * RQ1.subsample)) else ""
save(supp_results, file = paste0("results/supp_rq1_", model_label, "_results_",
                                 ts_stamp, sub_tag, ".Rdata"))


# ==============================================================================
# 6 Inspect results
# ==============================================================================

load("results/supp_rq1_AR2_results_2026-10-31_165240.Rdata")

tests_p <- paste0(RQ1.tests, "_p")

# 1. Misspecification rate per test ============================================
supp_results %>%
  dplyr::select(all_of(tests_p)) %>%
  pivot_longer(everything(), names_to = "test", values_to = "p") %>%
  filter(!is.na(p)) %>%
  group_by(test) %>%
  summarise(
    n_total     = n(),
    n_misspec   = sum(p <= 0.05),
    pct_misspec = mean(p <= 0.05) * 100,
    .groups = "drop"
  )

# A tibble: 5 × 4
# test    n_total n_misspec pct_misspec
# <chr>     <int>     <int>       <dbl>
# 1 adcf_p    16732      2963        17.7
# 2 bds_p     17064      3011        17.6
# 3 gspec_p   17077      3088        18.1
# 4 hsic_p    16732      2973        17.8
# 5 srho_p    17077      2970        17.4

# series without a result (fit failed, non-stationary fit, < 10 residuals)
supp_results %>%
  summarise(across(all_of(tests_p), ~ sum(is.na(.x))))

#   adcf_p gspec_p hsic_p srho_p bds_p
# 1    965     620    965    620   633

# 2. Comparison with the AR(1) results (same series, paired) ===================
main <- local({ load(rq1_results_main_file); rq1_results })

cmp <- supp_results %>%
  dplyr::select(ID, n_obs, all_of(tests_p)) %>%
  inner_join(dplyr::select(main, ID, all_of(tests_p)), by = "ID",
             suffix = c("_supp", "_ar1"))

for (t in tests_p) {
  a <- cmp[[paste0(t, "_ar1")]]  <= 0.05
  s <- cmp[[paste0(t, "_supp")]] <= 0.05
  ok <- !is.na(a) & !is.na(s)
  cat("\n==", t, "== (n =", sum(ok), ")\n")
  tab <- table(AR1  = factor(a[ok], levels = c(FALSE, TRUE)),
               supp = factor(s[ok], levels = c(FALSE, TRUE)))
  print(tab)
  cat(sprintf("AR(1): %.1f%%   %s: %.1f%%\n",
              100 * mean(a[ok]), model_label, 100 * mean(s[ok])))
  print(mcnemar.test(tab))
}

# rejection rates by series length
cmp %>%
  mutate(n_obs_bin = cut(n_obs,
                         breaks = c(25, 50, 75, 100, 125, 150, Inf),
                         labels = c("25-50","51-75","76-100","101-125","126-150",">150"),
                         include.lowest = TRUE)) %>%
  group_by(n_obs_bin) %>%
  summarise(n = n(),
            across(ends_with(c("_ar1", "_supp")), ~ mean(.x <= 0.05, na.rm = TRUE) * 100),
            .groups = "drop")

# 3. Model-specific descriptives ===============================================
# share of series with a significant lag-2 coefficient, and distribution of phi2
summary(supp_results$phi2)
mean(supp_results$phi2_p <= 0.05, na.rm = TRUE)
