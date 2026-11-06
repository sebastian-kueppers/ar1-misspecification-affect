# ==============================================================================
# Supplementary analyses RQ1: series without a refit-bootstrap result
#
# Counts, per test, why a series analysed in RQ1 (n_obs >= 25) has no p-value.
# Reasons follow the order of the checks in test_rb_refit() (helpers/rb_refit.R,
# estimator = "CLS"); each series is assigned the FIRST reason that applies:
#   2a. AR(1) fit failed (lm error or non-finite phi)
#   2b. non-stationary fit, |phi| >= 1
#   3a. fewer than 10 non-missing residuals
#   3b. residual SD = 0
#   4.  observed statistic not computable (*_stat NA), e.g. ADCF with fewer
#       than min_pairs = 8 complete lag-1 residual pairs
#   5.  statistic computed, but p NA (fewer than B/2 valid bootstrap replicates,
#       or an error inside test_rb_refit)
# Reasons 2-3 are recomputed by refitting the CLS AR(1) (no bootstrap);
# reasons 4-5 are read from the saved rq1_results.
#
# Contents:
#   1. Data (as in _analysis.R) and saved RQ1 results
#   2. Time-series inventory (as in _analysis.R)
#   3. Parameters
#   4. Refit and classify
#   5. Counts
# ==============================================================================

library(tidyverse)


# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------

# CLS AR(1) as in _analysis.R / rb_refit.R (.rb_fit, estimator = "CLS"); returns
# the quantities the guards in test_rb_refit() check
refit_checks <- function(ts) {
  ts_lag <- c(NA, ts[-length(ts)])
  fit <- tryCatch(lm(ts ~ ts_lag, na.action = na.exclude), error = function(e) NULL)
  phi <- if (is.null(fit)) NA_real_ else unname(coef(fit)["ts_lag"])
  if (is.null(fit) || !is.finite(phi))
    return(data.frame(fit_ok = FALSE, phi_refit = NA_real_, n_res = NA_integer_,
                      sd_res = NA_real_, n_pairs_lag1 = NA_integer_))
  res <- as.numeric(residuals(fit))
  ok  <- !is.na(res)
  data.frame(fit_ok       = TRUE,
             phi_refit    = phi,
             n_res        = sum(ok),
             sd_res       = sd(res, na.rm = TRUE),
             # complete lag-1 residual pairs (ADCF/HSIC with pairwise deletion)
             n_pairs_lag1 = sum(ok[-1] & ok[-length(ok)]))
}

# first applicable reason, in the order of test_rb_refit()
classify <- function(fit_ok, phi, n_res, sd_res, stat, p) {
  case_when(
    !fit_ok                   ~ "2a fit failed",
    abs(phi) >= 1             ~ "2b |phi| >= 1",
    n_res < 10                ~ "3a < 10 residuals",
    sd_res == 0               ~ "3b residual SD = 0",
    is.na(stat)               ~ "4 statistic NA",
    is.na(p)                  ~ "5 p NA (bootstrap)",
    TRUE                      ~ "valid"
  )
}


# ==============================================================================
# 1 Data and saved RQ1 results
# ==============================================================================

# load preprocessed data
load("data/processed/data_processed.Rdata")

# load input sheet for item info
inputsESM <- readxl::read_xlsx("inputs/input_file_ESM.xlsx")

# RQ1 results reported in the main text (same file as in _analysis.R)
load("results/rq1_results_2026-10-22_125950.Rdata")


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

# Which tests should be checked? (those run in the loaded rq1_results)
RQ1.tests <- c(
  "adcf",
  # "gspec",
  # "hsic",
  # "srho",
  "bds"
)

stopifnot(all(paste0(RQ1.tests, "_p") %in% names(rq1_results)))


# ==============================================================================
# 4 Refit and classify
# ==============================================================================

# Step 1 exclusion: fewer than 25 non-missing observations
table(analysed = ts_inventory$n_obs >= 25)

# analysed
# TRUE
# 17697
# --> was also part of the inclusion criteria for the main _analysis

rq1_inventory <- ts_inventory %>% filter(n_obs >= 25)

# all analysed series should be in rq1_results (full run, no subsample)
stopifnot(nrow(rq1_results) == nrow(rq1_inventory),
          setequal(rq1_results$ID, rq1_inventory$ID))

# person data, split once (faster than filtering df.list per series)
person_data <- lapply(df.list, function(d) split(d, d$ID))
names(person_data) <- sapply(df.list, function(d) unique(d$dataset))

checks <- rq1_inventory %>%
  pmap(function(dataset, ID, user_ID, item, n_beeps, ...) {
    ts <- person_data[[dataset]][[as.character(user_ID)]][[paste0(item, "_c")]]
    stopifnot(length(ts) == n_beeps)
    data.frame(ID = ID, refit_checks(ts))
  }) %>%
  bind_rows()

excl <- rq1_results %>%
  dplyr::select(ID, dataset, n_obs, phi, ends_with("_stat"), ends_with("_p")) %>%
  left_join(checks, by = "ID")

# the refit must reproduce the saved phi (same CLS fit as in _analysis.R)
stopifnot(isTRUE(all.equal(excl$phi, excl$phi_refit)))

for (t in RQ1.tests) {
  excl[[paste0(t, "_reason")]] <- classify(excl$fit_ok, excl$phi_refit, excl$n_res,
                                           excl$sd_res, excl[[paste0(t, "_stat")]],
                                           excl[[paste0(t, "_p")]])
}

# sanity: series flagged by reasons 2-3 must have no p-value
for (t in RQ1.tests) {
  r <- excl[[paste0(t, "_reason")]]
  stopifnot(all(is.na(excl[[paste0(t, "_p")]][substr(r, 1, 1) %in% c("2", "3")])))
}


# ==============================================================================
# 5 Counts
# ==============================================================================

# per test: number of series per reason
excl_counts <- excl %>%
  dplyr::select(ID, all_of(paste0(RQ1.tests, "_reason"))) %>%
  pivot_longer(-ID, names_to = "test", values_to = "reason") %>%
  mutate(test = sub("_reason$", "", test)) %>%
  count(test, reason) %>%
  pivot_wider(names_from = test, values_from = n, values_fill = 0) %>%
  arrange(reason)
excl_counts

# A tibble: 6 × 3
# reason              adcf   bds
# <chr>              <int> <int>
# 1 2a fit failed          6     6
# 2 2b |phi| >= 1         41    41
# 3 3a < 10 residuals      7     7
# 4 4 statistic NA       203     1
# 5 5 p NA (bootstrap)     0     4
# 6 valid              17440 17638


# excluded series per test (should equal nrow(rq1_results) - n valid)
excl_counts %>%
  filter(reason != "valid") %>%
  summarise(across(all_of(RQ1.tests), sum))

# A tibble: 1 × 2
# adcf   bds
# <int> <int>
#   1   257    59

# reason 4: ADCF statistic NA vs. number of complete lag-1 residual pairs
# (adcf_stat is NA with pairwise deletion if < min_pairs = 8 pairs)
if ("adcf" %in% RQ1.tests) {
  excl %>%
    filter(adcf_reason == "4 statistic NA") %>%
    count(pairs_lt_8 = n_pairs_lag1 < 8)
}

# pairs_lt_8   n
# 1       TRUE 203


# ==============================================================================
# 6 Inspect series excluded at the refit stage (reasons 2a, 2b, 3a)
# ==============================================================================

# raw series by ID
get_ts <- function(id) {
  r <- rq1_inventory[rq1_inventory$ID == id, ]
  person_data[[r$dataset]][[as.character(r$user_ID)]][[paste0(r$item, "_c")]]
}

# properties of the CLS sample (complete (y_t, y_{t-1}) pairs)
diag_series <- function(ts) {
  ts_lag <- c(NA, ts[-length(ts)])
  ok  <- !is.na(ts) & !is.na(ts_lag)
  fit <- tryCatch(lm(ts ~ ts_lag, na.action = na.exclude),
                  error = function(e) conditionMessage(e))
  err <- is.character(fit)
  data.frame(
    n_rows        = length(ts),
    n_obs         = sum(!is.na(ts)),
    n_pairs       = sum(ok),                          # CLS sample size
    n_unique      = length(unique(na.omit(ts))),      # distinct response values
    sd_y          = sd(ts, na.rm = TRUE),
    sd_y_pairs    = if (sum(ok) > 1) sd(ts[ok])     else NA_real_,
    sd_lag_pairs  = if (sum(ok) > 1) sd(ts_lag[ok]) else NA_real_,
    lm_error      = if (err) fit else NA_character_,
    phi           = if (err) NA_real_ else unname(coef(fit)["ts_lag"]),
    phi_se        = if (err || is.na(coef(fit)["ts_lag"])) NA_real_
    else coef(summary(fit))["ts_lag", "Std. Error"]
  )
}

reason_col <- paste0(RQ1.tests[1], "_reason")  # reasons 2-3 are identical across tests

refit_excl <- excl %>%
  filter(substr(.data[[reason_col]], 1, 1) %in% c("2", "3")) %>%
  dplyr::select(ID, reason = all_of(reason_col)) %>%
  left_join(dplyr::select(rq1_inventory, ID, dataset, user_ID, item), by = "ID")

refit_diag <- refit_excl %>%
  bind_cols(bind_rows(lapply(refit_excl$ID, function(id) diag_series(get_ts(id)))))

# overview of all excluded series
refit_diag %>%
  group_by(reason) %>%
  summarise(n = n(),
            n_obs_mdn   = median(n_obs),
            n_pairs_mdn = median(n_pairs),
            n_pairs_max = max(n_pairs),
            n_unique_mdn = median(n_unique),
            .groups = "drop")

# by dataset and item (clustering points to design/item causes)
count(refit_diag, reason, dataset, item)

# --- 2a: the 6 series where the AR(1) could not be estimated -----------------
fail_2a <- refit_diag %>% filter(reason == "2a fit failed")
print(as_tibble(fail_2a), width = Inf)
# expected: lm_error NA and sd_lag_pairs == 0 (constant predictor, slope NA)

# raw values
for (id in fail_2a$ID) {
  cat("\n== ID", id, "==\n")
  print(get_ts(id))
}

# plot (lines break at NAs, i.e. at missed beeps and nights)
bind_rows(lapply(fail_2a$ID, function(id) {
  y <- get_ts(id)
  tibble(ID = id, t = seq_along(y), y = y)
})) %>%
  ggplot(aes(t, y)) +
  geom_line() +
  geom_point(size = 0.8) +
  facet_wrap(~ ID, scales = "free_x") +
  labs(x = "Row (incl. overnight NA rows)", y = "Centered item score")

# --- 2b: non-stationary fits --------------------------------------------------
refit_diag %>%
  filter(reason == "2b |phi| >= 1") %>%
  arrange(desc(abs(phi))) %>%
  dplyr::select(ID, dataset, item, n_obs, n_pairs, n_unique, phi, phi_se) %>%
  as_tibble() %>%
  print(n = Inf)


# --- 3a: too few non-NA lag pairs ---------------------------------------------
refit_diag %>% filter(reason == "3a < 10 residuals") %>% count(dataset, user_ID)
# dataset user_ID n
# 1 Emotion Regulation Effort (2020)      73 7
