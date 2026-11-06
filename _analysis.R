# ==============================================================================
# Analysis
# ==============================================================================

# load packages
library(tidyverse)
library(tseries)
library(rEDM)
library(mgcv)
library(ecp)
library(future)
library(furrr)
library(tseries)
library(progressr)
library(data.table)
library(glmmTMB)
library(rlang)
library(ordinal)
library(ggtext)
library(ComplexUpset)
library(ggplot2)
library(dplyr)


# load preprocessed data
load("data/processed/data_processed.Rdata")
load("results/excluded_subjects_data.Rdata")

# load input sheet for item info
inputsESM <- readxl::read_xlsx("inputs/input_file_ESM.xlsx")

# load helpers
# RQ1
source("helpers/adcf.R")
source("helpers/bds.R")
source("helpers/gspec.R")
source("helpers/hsic.R")
source("helpers/srho.R")
source("helpers/rb_refit.R")   # refit bootstrap null; uses the *_stat() functions above

# RQ2
source("helpers/edf_gam.R")
source("helpers/optimal_theta_from_SMap.R")
source("helpers/n_changepoints.R")

# ==============================================================================
## DESCRIPTIVES

## 1. Participants -----------------------------------------------------
# ID restarts at 1 in every dataset, so key on UUID
ids <- rbindlist(lapply(df.list, \(d) d[, c("dataset", "ID", "UUID")]))

uniqueN(ids$UUID)                      # participants, 2108
uniqueN(ids, by = c("dataset", "ID"))  # sanity check


## 2. Observations -----------------------------------------------------
sum(sapply(df.list, nrow))   # all rows, incl. inserted overnight NA rows

# completed beeps: rows with at least one non-missing affect item
sum(sapply(df.list, \(d) {
  cols <- grep("_ES_c$", names(d), value = TRUE)
  sum(rowSums(!is.na(d[, cols, drop = FALSE])) > 0)
}))

# individual datapoints across all person × item series
sum(sapply(df.list, \(d) {
  cols <- grep("_ES_c$", names(d), value = TRUE)
  sum(!is.na(d[, cols, drop = FALSE]))
}))

aff_cols <- function(d) {
  items <- inputsESM %>%
    filter(dataset_presentation == unique(d$dataset), !is.na(construct)) %>%
    pull(variable)
  intersect(paste0(items, "_c"), names(d))
}

n_rows       <- sum(sapply(df.list, nrow))   # includes inserted overnight NA rows
n_beeps      <- sum(sapply(df.list, \(d) sum(rowSums(!is.na(d[, aff_cols(d), drop = FALSE])) > 0)))
n_datapoints <- sum(sapply(df.list, \(d) sum(!is.na(d[, aff_cols(d), drop = FALSE]))))

counts <- c(rows = n_rows, completed_beeps = n_beeps, datapoints = n_datapoints)
counts


## 3. Average Number of Observations per Day ----------------------------

# 0. Which day variable do the datasets have?
Reduce(intersect, lapply(df.list, names))
day_var <- "dayno"   # day number within person (see _data_preprocessing.R);
# inserted overnight NA rows have dayno = NA and are dropped below

# 1. Row level: one row per beep
beep_level <- lapply(df.list, function(d) {
  stopifnot(day_var %in% names(d))
  cols <- aff_cols(d)
  stopifnot(length(cols) > 0)
  answered <- !is.na(d[, cols, drop = FALSE])
  tibble(
    dataset   = d$dataset,
    UUID      = d$UUID,
    day       = d[[day_var]],
    beep_done = rowSums(answered) > 0,   # completed beep
    items_ans = rowSums(answered),       # items answered at this beep
    n_items   = length(cols)
  )
}) %>% bind_rows()

# 2. Person-day level ----------------------------------------------------------
person_day <- beep_level %>%
  group_by(dataset, UUID, day) %>%
  summarise(beeps     = sum(beep_done),
            items_ans = sum(items_ans),
            n_items   = dplyr::first(n_items),
            .groups = "drop") %>%
  mutate(obs_per_item = items_ans / n_items)

# Keep only days with at least one completed beep (inserted overnight NA rows
# and fully missed days would otherwise count as days with 0 observations).
# To include them, skip this filter and compare.
person_day_obs <- person_day %>% filter(beeps > 0)

# 3. Summaries -----------------------------------------------------------------
by_dataset <- person_day_obs %>%
  group_by(dataset) %>%
  summarise(persons          = n_distinct(UUID),
            person_days      = n(),
            beeps_per_day    = mean(beeps),
            sd_beeps         = sd(beeps),
            obs_per_item_day = mean(obs_per_item),
            median_beeps     = median(beeps),
            .groups = "drop") %>%
  arrange(desc(beeps_per_day))
print(by_dataset, n = Inf, width = Inf)

# overall <- tibble(
#   # pooled over all person-days (large datasets weigh more)
#   beeps_per_day_pooled       = mean(person_day_obs$beeps),
#   obs_per_item_day_pooled    = mean(person_day_obs$obs_per_item),
#   # mean of dataset means (each dataset weighs the same)
#   beeps_per_day_by_dataset   = mean(by_dataset$beeps_per_day),
#   obs_per_item_day_by_dataset = mean(by_dataset$obs_per_item_day),
#   person_days = nrow(person_day_obs)
# )
# print(overall, width = Inf)

round(median(by_dataset$beeps_per_day)) # 8


# ==============================================================================
# 1 Build time series and person inventories
# ==============================================================================

## Time series inventory
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
  mutate(ID = row_number()) %>% # unique identifier
  dplyr::select(ID, dataset, user_ID, UUID, item, construct, affect_type, n_obs, n_beeps)

sum(ts_inventory$n_obs) == n_datapoints                       # is TRUE
nrow(ts_inventory)                                            # person x item series, 17697
sum(ts_inventory$n_obs[ts_inventory$n_obs >= 25])             # datapoints in analyzed series, 1241129

## Person inventory
person_inventory <- lapply(seq_along(df.list), function(k){
  d <- df.list[[k]]
  grab <- function(var) if (var %in% names(d)) d[[var]] else NA_real_
  d %>%
    mutate(.DEP = grab("DEP"), .NEUR = grab("NEUR"), .SWL = grab("SWL")) %>%
    group_by(user_ID = ID, UUID) %>%
    summarise(
      dataset = unique(dataset),
      DEP  = if (all(is.na(.DEP)))  NA_real_ else first(na.omit(.DEP)),
      NEUR = if (all(is.na(.NEUR))) NA_real_ else first(na.omit(.NEUR)),
      SWL  = if (all(is.na(.SWL)))  NA_real_ else first(na.omit(.SWL)),
      .groups = "drop"
    )
}) %>%
  bind_rows() %>%
  # within-dataset standardization (over persons)
  group_by(dataset) %>%
  mutate(
    DEP_z  = as.numeric(scale(DEP)),
    NEUR_z = as.numeric(scale(NEUR)),
    SWL_z  = as.numeric(scale(SWL))
  ) %>%
  ungroup()

# Descriptives for time series lengths
summary(ts_inventory$n_obs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 26.00   55.00   62.00   70.13   71.00  204.00

table(cut(ts_inventory$n_obs, breaks = c(24, 50, 75, 100, 125, 150, Inf)))
# (24,50]   (50,75]  (75,100] (100,125] (125,150] (150,Inf]
# 2973     11083      1824       564       289       964

# ==============================================================================
# 2 RQ1: AR(1) misspecification
# ==============================================================================

# --------------------------------------
# PARAMETERS ---------------------------
# --------------------------------------

# Which tests should be applied for RQ1?
RQ1.tests <- c(
  "adcf",
  "gspec",
  "hsic",
  "srho",
  "bds"
)

# Should NAs in residuals (e.g., due to overnight lags) be deleted pairwise?
# (FALSE = list-wise deletion)
# BDS and Srho test do not offer pairwise deletion
PAIRWISE_NA_DELETION <- TRUE
na_handling <- if (PAIRWISE_NA_DELETION) "pairwise" else "strip"

# Statistic settings for the refit bootstrap (helpers/rb_refit.R). Settings not
# listed here keep the rb_refit defaults used in the size simulation.
RQ1.stat_args <- list(
  adcf  = list(max_lag = 1, na_handling = na_handling),
  hsic  = list(max_lag = 1, na_handling = na_handling),
  gspec = list(na_handling = na_handling)
)

# Number of bootstrap replicates (refit bootstrap; B + 1 should be a multiple
# of 20 so that rejecting at p <= .05 has size exactly .05)
B <- 499
# B <- 199

stopifnot(all(RQ1.tests %in% c("adcf", "gspec", "hsic", "srho", "bds")))


# --------------------------------------
# HELPERS ------------------------------
# --------------------------------------

# normalize any test return to c(stat, p)
.sp <- function(res) {
  if (is.null(res)) return(c(NA_real_, NA_real_))
  if (is.list(res)) return(c(res$statistic %||% NA_real_, res$p.value %||% NA_real_))
  c(NA_real_, as.numeric(res))
}


# --------------------------------------
# PERFORM CLASSIFICATION ---------------
# --------------------------------------

# Run RQ1 on a random subsample of the series? (share in (0, 1]; 1 = all).
# row_id is assigned BEFORE subsampling, so each series keeps its bootstrap
# seed (seed = row_id) whether or not it is run as part of a subsample.
RQ1.subsample      <- 1
RQ1.subsample_seed <- 2026

rq1_inventory <- ts_inventory %>%
  filter(n_obs >= 25) %>%
  mutate(row_id = row_number())
if (RQ1.subsample < 1) {
  set.seed(RQ1.subsample_seed)
  rq1_inventory <- rq1_inventory %>% slice_sample(prop = RQ1.subsample)
}
nrow(rq1_inventory)

# set up parallel backend
n_cores <- parallelly::availableCores(omit = 1)
plan(multisession, workers = n_cores)

set.seed(42)
start <- Sys.time()
with_progress({
  p <- progressor(steps = nrow(rq1_inventory))
  rq1_results <- rq1_inventory %>%
    furrr::future_pmap(
      .options = furrr_options(seed = TRUE),
      function(dataset, ID, user_ID, UUID, item, construct, affect_type, n_obs, n_beeps, row_id){
        p()
        k  <- which(sapply(df.list, function(df) unique(df$dataset)) == dataset)
        ts <- df.list[[k]] %>%
          filter(.data$ID == .env$user_ID) %>%
          # arrange(beep) %>% pull(paste0(item, "_c")) # arrange puts NAs to the end of the TS
          pull(paste0(item, "_c"))

        stopifnot(length(ts) == n_beeps)

        # ar_fit <- tryCatch(arima(ts, order = c(1,0,0), method = "ML"),
        #                    error = function(e) NULL)
        # phi    <- if (is.null(ar_fit)) NA_real_ else unname(coef(ar_fit)["ar1"])
        # sigma2 <- if (is.null(ar_fit)) NA_real_ else ar_fit$sigma2

        # CLS: OLS of y_t on y_{t-1}. Pairs containing an NA (overnight NA rows,
        # missed beeps) are dropped, so no lag spans the night and the first
        # observed beep after a gap is not an outcome.
        ts_lag <- c(NA, ts[-length(ts)])
        ar_fit <- tryCatch(lm(ts ~ ts_lag, na.action = na.exclude),
                           error = function(e) NULL)
        phi    <- if (is.null(ar_fit)) NA_real_ else unname(coef(ar_fit)["ts_lag"])
        sigma2 <- if (is.null(ar_fit)) NA_real_ else sigma(ar_fit)^2

        # every test gets a slot; NA if not requested or fit failed
        R <- list(adcf  = c(NA_real_, NA_real_), gspec = c(NA_real_, NA_real_),
                  hsic  = c(NA_real_, NA_real_), srho  = c(NA_real_, NA_real_),
                  bds   = c(NA_real_, NA_real_))

        # Refit bootstrap null (NEW): takes the RAW series and refits the AR(1)
        # by CLS on the observed series and on every bootstrap series; all
        # requested tests share the same B draws
        if (!is.null(ar_fit)) {
          rb <- tryCatch(test_rb_refit(ts, tests = RQ1.tests, B = B, seed = row_id,
                                       stat_args = RQ1.stat_args, estimator = "CLS"),
                         error = function(e) NULL)
          if (!is.null(rb)) for (t in names(rb)) R[[t]] <- .sp(rb[[t]])
        }

        # Residual permutation null (OLD, no refit) on the CLS residuals.
        # na.exclude keeps res aligned to ts: NA wherever y_t or y_{t-1} is NA,
        # incl. the first observed beep after the night. Under pairwise
        # deletion the permutation holds these NA positions fixed.
        # if (!is.null(ar_fit)) {
        #   res <- as.numeric(residuals(ar_fit))
        #   rp <- list(
        #     adcf  = function() test_adcf (res, input = "residual", max_lag = 1,
        #                                   na_handling = na_handling, B = B, seed = row_id),
        #     gspec = function() test_gspec(res, input = "residual",
        #                                   na_handling = na_handling, B = B, seed = row_id),
        #     hsic  = function() test_hsic (res, input = "residual", max_lag = 1,
        #                                   na_handling = na_handling, B = B, seed = row_id),
        #     srho  = function() test_srho (res, input = "residual", max_lag = 1,
        #                                   B = B, seed = row_id),
        #     bds   = function() test_bds  (res, input = "residual", m = 2, eps = 1,
        #                                   B = B, seed = row_id)
        #   )
        #   for (t in RQ1.tests) R[[t]] <- .sp(tryCatch(rp[[t]](), error = function(e) NULL))
        # }

        data.frame(
          ID, dataset, UUID, item, construct, affect_type,
          n_obs, n_na = sum(is.na(ts)), row_id, phi, sigma2,
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
save(rq1_results, file = paste0("results/rq1_results_", ts_stamp, sub_tag, ".Rdata"))


# --------------------------------------
# INSPECT RESULTS ----------------------
# --------------------------------------

load("results/rq1_results_2026-10-22_125950.Rdata")

# 1. Overall misspecification rate ============================================

# ADCF
rq1_results %>%
  filter(!is.na(adcf_p)) %>%
  summarise(
    n_total       = n(),
    n_misspec     = sum(adcf_p <= 0.05),
    pct_misspec   = mean(adcf_p <= 0.05) * 100
  )

# n_total n_misspec pct_misspec
# 1   17440      4295    24.62729

# BDS
rq1_results %>%
  filter(!is.na(bds_p)) %>%
  summarise(
    n_total       = n(),
    n_misspec     = sum(bds_p <= 0.05),
    pct_misspec   = mean(bds_p <= 0.05) * 100
  )

# n_total n_misspec pct_misspec
# 1   17638      3970    22.50822

# confidence intervals
cluster_ci <- function(d, pcol, B = 2000, seed = 2026) {
  pers <- d %>%
    filter(!is.na(.data[[pcol]])) %>%
    group_by(dataset, UUID) %>%
    summarise(k = sum(.data[[pcol]] <= .05), n = n(), .groups = "drop")
  set.seed(seed)
  boot <- replicate(B, {
    s <- pers %>% group_by(dataset) %>%
      slice_sample(prop = 1, replace = TRUE) %>% ungroup()
    sum(s$k) / sum(s$n)
  })
  round(100 * c(est = sum(pers$k) / sum(pers$n),
                quantile(boot, c(.025, .975))), 1)
}

cluster_ci(rq1_results, "adcf_p")
# est  2.5% 97.5%
#   24.6  23.8  25.5
cluster_ci(rq1_results, "bds_p")
# est  2.5% 97.5%
#   22.5  21.7  23.3

# PA / NA:
cluster_ci(filter(rq1_results, affect_type == "NEG"), "adcf_p")
# est  2.5% 97.5%
#   27.7  26.7  28.8

cluster_ci(filter(rq1_results, affect_type == "POS"), "adcf_p")
# est  2.5% 97.5%
#   20.4  19.4  21.4

cluster_ci(filter(rq1_results, affect_type == "NEG"), "bds_p")
# est  2.5% 97.5%
#   25.2  24.2  26.2

cluster_ci(filter(rq1_results, affect_type == "POS"), "bds_p")
# est  2.5% 97.5%
#   18.9  17.9  19.8

# misspecification per dataset
tab_len <- rq1_results %>%
  group_by(dataset) %>%
  summarise(
    n_persons    = n_distinct(UUID),
    n_series     = n(),
    T_median     = median(n_obs),
    T_min        = min(n_obs),
    T_max        = max(n_obs),
    pct_T_gt100  = round(100 * mean(n_obs > 100), 1),
    adcf_n       = sum(!is.na(adcf_p)),
    adcf_pct     = round(100 * mean(adcf_p <= .05, na.rm = TRUE), 1),
    .groups = "drop") %>%
  arrange(desc(T_median))
print(tab_len, n = Inf, width = Inf)

# dataset                                                n_persons n_series T_median T_min T_max pct_T_gt100 adcf_n adcf_pct
# <chr>                                                      <int>    <int>    <dbl> <int> <int>       <dbl>  <int>    <dbl>
#   1 FEEL Study 1                                                 179     1074    173      67   204        97.8   1073     47.5
# 2 Leuven emotions in daily life 2008                            69      759    121      27   175        82.6    756     34.4
# 3 Trier social stress study                                     72      504     93      55   115        27.8    504     38.5
# 4 Emotional events in daily life                               101      505     90      45    98         0      505     22
# 5 Exam results study                                           101     1010     83      58    90         0     1004     55.4
# 6 Leuven emotions in daily life 2012                           100     1000     71      43    77         0      999     21.7
# 7 ACU emotions in daily life                                    69      690     66      43    82         0      689     19
# 8 Leuven couples study 2014                                     98      882     65      53    72         0      879     23.3
# 9 Leuven clinical study                                        111      999     63      35    70         0      999     19.2
# 10 Leuven 3-wave longitudinal study                             197     1773     62      33    76         0     1772     22.5
# 11 Leuven emotions in daily life 2011                            94      658     61      39    73         0      658     18.4
# 12 Objectification in women's daily lives Study 1                78      624     61      29    75         0      622     19.6
# 13 Objectification in women's daily lives Study 2                58      522     60      27    75         0      512     23.2
# 14 Objectification in women's daily lives Study 3                79      711     59      29    74         0      699     23.6
# 15 Leuven couples study 2016                                    174     1044     57.5    43    70         0     1038     27.6
# 16 Goals, Emotions, Moods, and Strategies (GEMS) 2024           228     2508     57      28    70         0     2431     14.6
# 17 Assessing, Understanding, and Regulating Affect (AURA)       167     1503     50      26    63         0     1404     14.7
# 18 Emotion Regulation Effort (2020)                             133      931     46      26    56         0      896     16.2


# all
rq1_results %>%
  dplyr::select(adcf_p, srho_p, hsic_p, gspec_p, bds_p) %>%
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
# 1 adcf_p    17440      4295        24.6
# 2 bds_p     17638      3970        22.5
# 3 gspec_p   17643      4227        24.0
# 4 hsic_p    17440      3699        21.2
# 5 srho_p    17643      4046        22.9

tests <- c("adcf_p", "srho_p", "hsic_p", "gspec_p", "bds_p")

rq1_results %>%
  mutate(n_sig = rowSums(across(all_of(tests), ~ .x <= 0.05), na.rm = TRUE)) %>%
  count(n_sig) %>%
  mutate(pct = n / sum(n) * 100) %>%
  arrange(desc(n_sig))

# n_sig     n       pct
# 1     5  1539  8.696389
# 2     4  1193  6.741256
# 3     3  1021  5.769339
# 4     2  1247  7.046392
# 5     1  2213 12.504944
# 6     0 10484 59.241679

# Visualization in UpSet plot

flags <- rq1_results %>%
  filter(if_all(all_of(tests), ~ !is.na(.x))) %>%
  mutate(across(all_of(tests), ~ .x <= 0.05)) %>%
  dplyr::select(ADCF = adcf_p, S_rho = srho_p, HSIC = hsic_p,
         Gen_Spec = gspec_p, BDS = bds_p)

flags_plot <- flags[rowSums(flags) > 0, ]

n_all  <- nrow(flags)
n_none <- sum(rowSums(flags) == 0)

col_full <- c(
  ADCF     = "#0072B2",
  BDS      = "#D55E00",
  S_rho    = "#009E73",
  HSIC     = "#CC79A7",
  Gen_Spec = "#E69F00"
)

p <- upset(
  flags_plot,
  intersect = names(col_full),
  mode = "exclusive_intersection",
  min_size = 100,
  sort_intersections_by = "cardinality",
  sort_intersections = "descending",
  width_ratio = 0.2,
  name = "Test combination",
  base_annotations = list(
    "Share of all series (%)" = intersection_size(
      mode = "exclusive_intersection",
      bar_number_threshold = 1,
      text = list(size = 3.2),                    # was 2.5
      text_colors = c(on_background = "black", on_bar = "black"),
      text_mapping = aes(
        label = paste0(round(!!get_size_mode("exclusive_intersection") / n_all * 100, 1), "%")
      )
    ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  ),
  queries = lapply(names(col_full), function(s)
    upset_query(set = s, fill = col_full[[s]])),
  set_sizes = upset_set_size(position = "right") +
    scale_y_continuous(
      limits = c(0, n_all),
      breaks = n_all * c(0, 0.25, 0.5, 0.75, 1),
      labels = function(x) paste0(round(x / n_all * 100), "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    ylab("Flagged (of all series)")
)

p

ggsave("figures/upset_misspecification.png", p, width = 9, height = 6, units = "in", dpi = 300, bg = "white")


rq1_long <- rq1_results %>%
  filter(!is.na(affect_type)) %>%
  dplyr::select(ID, n_obs, affect_type, ADCF = adcf_p, BDS = bds_p) %>%
  pivot_longer(c(ADCF, BDS), names_to = "Test", values_to = "p") %>%
  filter(!is.na(p)) %>%                            # testable series per test
  mutate(misspec = as.numeric(p <= .05),
         affect  = factor(recode(affect_type, POS = "Positive affect",
                                 NEG = "Negative affect"),
                          levels = c("Positive affect", "Negative affect")))

rq1_long %>% count(Test, affect)                   # n per facet for the note

p_len <- ggplot(rq1_long, aes(x = n_obs, y = misspec,
                              colour = affect, fill = affect)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              method.args = list(family = binomial, method = "REML"),
              alpha = 0.2, linewidth = 1) +
  geom_rug(data = rq1_long %>% filter(misspec == 0), sides = "b",
           alpha = 0.02, length = unit(0.02, "npc"), show.legend = FALSE) +
  geom_rug(data = rq1_long %>% filter(misspec == 1), sides = "t",
           alpha = 0.02, length = unit(0.02, "npc"), show.legend = FALSE) +
  geom_hline(yintercept = .05, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ Test) +
  scale_colour_manual(values = col_full, name = NULL) +
  scale_fill_manual(values = col_full, name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
  labs(x = "Time series length (non-missing observations)",
       y = "AR(1)-misspecified series") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"))

p_len
ggsave("figures/misspecification_by_length_smooth_ADCF_BDS.png", p_len,
       width = 11, height = 6, dpi = 500, bg = "white")



# 3. Misspecification rate by affect type =====================================
rq1_results %>%
  filter(!is.na(adcf_p <= 0.05)) %>%
  group_by(affect_type) %>%
  summarise(
    n_total     = n(),
    n_misspec   = sum(adcf_p <= 0.05),
    pct_misspec = mean(adcf_p <= 0.05) * 100,
    .groups     = "drop"
  )

# A tibble: 2 × 4
# affect_type n_total n_misspec pct_misspec
# <chr>         <int>     <int>       <dbl>
# 1 NEG           10061      2787        27.7
# 2 POS            7379      1508        20.4

# per dataset
bind_rows(
  rate_tab(rq1_results, dataset, affect_type),
  rate_tab(rq1_results, affect_type) %>% mutate(dataset = "ALL")
) %>% print(n = 50)

# 3. Misspecification rate by construct and dataset ===========================
rate_tab <- function(d, ...) d %>%
  group_by(...) %>%
  summarise(n_series = n(),
            ADCF = round(mean(adcf_p <= .05, na.rm = TRUE) * 100, 1),
            BDS  = round(mean(bds_p  <= .05, na.rm = TRUE) * 100, 1),
            .groups = "drop")

tab_item    <- rate_tab(rq1_results, dataset, item, affect_type)   # item table
tab_dataset <- rate_tab(rq1_results, dataset)                      # dataset table
tab_constr  <- rate_tab(rq1_results, construct, affect_type)


# 4. Misspecification rate by n_obs ===========================================
# bin n_obs into groups and check misspecification rate per bin


## -- Plot per test, time series length, and affect type

# colorblind-safe
col_full <- c("Positive affect" = "#009E73",   # green
              "Negative affect" = "#7B3FA0")    # purple

d2 <- rq1_results %>%
  filter(!is.na(adcf_p), !is.na(bds_p), !is.na(affect_type)) %>%
  mutate(n_obs_bin = cut(n_obs,
                         breaks = c(25, 50, 75, 100, 125, 150, Inf),
                         labels = c("25-50","51-75","76-100","101-125","126-150",">150"))) %>%
  group_by(n_obs_bin, affect_type) %>%
  summarise(ADCF = mean(adcf_p <= .05) * 100,
            BDS  = mean(bds_p  <= .05) * 100,
            n    = n(), .groups = "drop") %>%
  pivot_longer(c(ADCF, BDS), names_to = "Test", values_to = "pct") %>%
  mutate(affect = recode(affect_type, POS = "Positive affect", NEG = "Negative affect"),
         affect = factor(affect, levels = c("Positive affect", "Negative affect")),
         n_obs_bin = factor(n_obs_bin)) %>%
  droplevels()

lab2 <- d2 %>% group_by(n_obs_bin, affect, Test) %>%
  summarise(n = first(n), .groups = "drop")

ggplot(d2, aes(x = n_obs_bin, group = affect)) +
  geom_col(aes(y = 100, fill = affect), alpha = 0.35,
           position = position_dodge(0.8), width = 0.7, show.legend = FALSE) +
  geom_col(aes(y = pct, fill = affect),
           position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = col_full, name = NULL) +
  geom_text(aes(y = pct, label = sprintf("%.1f", pct), group = affect),
            position = position_dodge(0.8), vjust = -0.4, size = 2.6) +
  geom_text(data = lab2,
            aes(x = n_obs_bin, y = 104, label = paste0("italic(n)==", n),
                group = affect),
            position = position_dodge(0.8),
            angle = 45, hjust = 0, vjust = 0,
            inherit.aes = FALSE, size = 2.6, colour = "grey30", parse = TRUE) +
  facet_wrap(~ Test) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "grey50") +
  labs(x = "Time series length (observations)",
       y = "AR(1)-misspecified series (%)") +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("figures/misspecification_by_length_and_affect_ADCF_BDS.png",
       width = 11, height = 7, dpi = 500)


# -- Plot per-participant proportion

pp <- rq1_results %>%
  filter(!is.na(adcf_p), !is.na(bds_p)) %>%
  group_by(dataset, UUID) %>%
  summarise(n_series = n(),
            prop_adcf = mean(adcf_p <= .05),
            prop_bds  = mean(bds_p  <= .05),
            .groups = "drop") %>%
  pivot_longer(c(prop_adcf, prop_bds),
               names_to = "Test", values_to = "prop") %>%
  mutate(Test = recode(Test, prop_adcf = "ADCF", prop_bds = "BDS"))

pp_binned <- pp %>%
  mutate(bin = cut(prop, breaks = seq(0, 1, 0.1), include.lowest = TRUE,
                   labels = c("0–10%","10-20%", "20–30%","30-40%", "40–50%", "50-60%", "60-70%", "70–80%","80–90%", "90-100%"))) %>%
  count(bin, Test)

col_full <- c(ADCF = "#0072B2",
              BDS = "#D55E00")

ggplot(pp_binned, aes(bin, n, fill = Test)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = col_full, name = NULL) +
  labs(x = "Proportion of a participant's series classified as AR(1) misspecified",
       y = "Number of participants") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 10)))

ggsave("figures/misspecification_by_participant_ADCF_BDS.png",
       width = 11, height = 7, dpi = 500)



# also check correlation between n_obs and bds_result
cor.test(rq1_results$n_obs, as.numeric(rq1_results$adcf_p <= 0.05), method = "spearman")
# Spearman's rank correlation rho
#
# data:  rq1_results$n_obs and as.numeric(rq1_results$adcf_p <= 0.05)
# S = 6.8549e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.2246178

cor.test(rq1_results$n_obs, as.numeric(rq1_results$bds_p <= 0.05), method = "spearman")
# data:  rq1_results$n_obs and as.numeric(rq1_results$bds_p <= 0.05)
# S = 6.7646e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho
# 0.2603147


# ==============================================================================
# 3 RQ2: Complexity markers
# ==============================================================================

n_cores <- parallel::detectCores() - 1
plan(multisession, workers = n_cores)

set.seed(42)
start <- Sys.time()
with_progress({
  p <- progressor(steps = nrow(ts_inventory %>% filter(n_obs >= 25)))
  rq2_results <- ts_inventory %>%
    filter(n_obs >= 25) %>%
    mutate(row_id = row_number()) %>%
    furrr::future_pmap(
      .options = furrr_options(seed = TRUE),
      function(ID, dataset, user_ID, UUID, item, construct, affect_type,
               n_obs, n_beeps, row_id){
        p()
        # locate dataset, extract the centered series (overnight NAs kept)
        k <- which(sapply(df.list, function(df) unique(df$dataset)) == dataset)
        ts_data <- df.list[[k]] %>%
          filter(.data$ID == .env$user_ID) %>%
          # arrange(beep) %>%
          pull(paste0(item, "_c"))

        stopifnot(length(ts_data) == n_beeps)

        ts_clean <- ts_data[!is.na(ts_data)]     # for e.divisive only

        edf_lag <- edf_gam(ts_data)

        theta_opt <- optimal_theta_from_SMap(ts_data, drop_na = FALSE)

        # e.divisive cannot take NAs; must use stripped series.
        changepoint_n <- n_changepoints(ts_clean,
                                        R = 400,
                                        min.size = 8) # median num of obs per day

        data.frame(
          ID            = ID,
          dataset       = dataset,
          user_ID       = user_ID,
          UUID          = UUID,
          item          = item,
          construct     = construct,
          affect_type   = affect_type,
          n_obs         = n_obs,
          n_na          = sum(is.na(ts_data)),
          edf_lag       = edf_lag,
          theta_opt     = theta_opt,
          changepoint_n = changepoint_n,
          stringsAsFactors = FALSE
        )
      }
    ) %>% bind_rows()
})
end <- Sys.time()
print(end - start)

plan(sequential)

# save results
ts_stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
save(rq2_results, file = paste0("results/rq2_results_", ts_stamp, ".Rdata"))

# --------------------------------------
# INSPECT RESULTS ----------------------
# --------------------------------------

# --------------------------------------
# DESCRIPTVES --------------------------

rq2_results %>%
  summarise(
    n_series   = n(),
    n_persons  = n_distinct(dataset, UUID),
    n_datasets = n_distinct(dataset),

    # --- Missingness (% of ALL series) ---
    edf_na_pct   = mean(is.na(edf_lag)) * 100,
    theta_na_pct = mean(is.na(theta_opt)) * 100,
    cp_na_pct    = mean(is.na(changepoint_n)) * 100,

    # --- EDF ---
    edf_n   = sum(!is.na(edf_lag)),
    edf_na  = sum(is.na(edf_lag)),
    edf_mdn = median(edf_lag, na.rm = TRUE),
    edf_gt1_mean = mean(edf_lag[edf_lag > 1.001 & !is.na(edf_lag)], na.rm = TRUE),
    edf_gt1_est = mean(edf_lag > 1.001, na.rm = TRUE) * 100,
    edf_gt1_n_pa  = sum(edf_lag > 1.001 & affect_type == "POS", na.rm = TRUE),
    edf_gt1_m_pa  = mean(edf_lag[which(edf_lag > 1.001 & affect_type == "POS")]),
    edf_gt1_sd_pa = sd(edf_lag[which(edf_lag > 1.001 & affect_type == "POS")]),
    edf_gt1_n_na  = sum(edf_lag > 1.001 & affect_type == "NEG", na.rm = TRUE),
    edf_gt1_m_na  = mean(edf_lag[which(edf_lag > 1.001 & affect_type == "NEG")]),
    edf_gt1_sd_na = sd(edf_lag[which(edf_lag > 1.001 & affect_type == "NEG")]),


    # --- Theta ---
    theta_n   = sum(!is.na(theta_opt)),
    theta_na  = sum(is.na(theta_opt)),
    theta_m   = mean(theta_opt, na.rm = TRUE),
    theta_sd  = sd(theta_opt, na.rm = TRUE),
    theta_mdn = median(theta_opt, na.rm = TRUE),
    theta_0     = mean(theta_opt == 0 & !is.na(theta_opt)) * 100,
    theta_0_est = mean(theta_opt == 0, na.rm = TRUE) * 100,
    theta_8     = mean(theta_opt == 8 & !is.na(theta_opt)) * 100,
    theta_8_est = mean(theta_opt == 8, na.rm = TRUE) * 100,

    # --- Change points ---
    cp_m   = mean(changepoint_n, na.rm = TRUE),
    cp_sd  = sd(changepoint_n, na.rm = TRUE),
    cp_ge1     = mean(changepoint_n >= 1 & !is.na(changepoint_n)) * 100,
    cp_ge1_est = mean(changepoint_n >= 1, na.rm = TRUE) * 100
  ) %>%
  tidyr::pivot_longer(everything()) %>%
  print(n = Inf)

#   n_series n_persons n_datasets edf_n edf_na  edf_mdn  edf_gt1 theta_n theta_na  theta_m theta_sd theta_mdn  theta_0  theta_8      cp_m    cp_sd   cp_ge1
# 1    17697      2108         18 15961   1736 1.000015 18.14423   16494     1203 3.004668 3.242397       1.5 32.07833 22.72948 0.9800531 1.173615 55.33706

# -- Helper functions --

col_full <- c("Positive affect" = "#009E73",
              "Negative affect" = "#7B3FA0")

theme_len <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

plot_prop <- function(data, flag, ylab, ref_line = NULL) {
  d <- data %>%
    filter(!is.na({{ flag }})) %>%
    group_by(n_obs_bin, affect) %>%
    summarise(pct = mean({{ flag }}) * 100, n = n(), .groups = "drop")

  p <- ggplot(d, aes(x = n_obs_bin, group = affect)) +
    geom_col(aes(y = 100, fill = affect), alpha = 0.35,
             position = position_dodge(0.8), width = 0.7, show.legend = FALSE) +
    geom_col(aes(y = pct, fill = affect), position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(y = pct, label = sprintf("%.1f", pct)),
              position = position_dodge(0.8), vjust = -0.4, size = 2.6) +
    geom_text(aes(y = 104, label = paste0("italic(n)==", n)),
              position = position_dodge(0.8), angle = 45, hjust = 0, vjust = 0,
              size = 2.6, colour = "grey30", parse = TRUE) +
    scale_fill_manual(values = col_full, name = NULL) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25)) +
    labs(x = "Time series length (non-missing observations)", y = ylab) +
    theme_len
  if (!is.null(ref_line))
    p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", colour = "grey50")
  p
}

plot_count <- function(data, value, ylab) {
  d <- data %>%
    filter(!is.na({{ value }})) %>%
    group_by(n_obs_bin, affect) %>%
    summarise(m = mean({{ value }}), n = n(), .groups = "drop")

  top <- max(d$m) * 1.08                    # panel ends just above the highest bar
  br  <- pretty(c(0, max(d$m)))
  br  <- br[br <= max(d$m)]

  ggplot(d, aes(x = n_obs_bin, y = m, fill = affect, group = affect)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", m)),
              position = position_dodge(0.8), vjust = -0.4, size = 2.6) +
    # n labels sit above the panel, outside the grid
    geom_text(aes(y = top * 1.02, label = paste0("italic(n)==", n)),
              position = position_dodge(0.8), angle = 45, hjust = 0, vjust = 0,
              size = 2.6, colour = "grey30", parse = TRUE) +
    scale_fill_manual(values = col_full, name = NULL) +
    scale_y_continuous(breaks = br, expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, top), clip = "off") +
    labs(x = "Time series length (non-missing observations)", y = ylab) +
    theme_len +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(t = 40, r = 25, b = 5, l = 5))
}


# -- Preprocess Data --

# load RQ1 results for ADCF outcome
# load()

# load RQ2 results
load("results/rq2_results_2026-10-13_112850.Rdata")

rq2_d <- rq2_results %>%
  filter(!is.na(affect_type)) %>%
  mutate(
    n_obs_bin = cut(n_obs,
                    breaks = c(24, 50, 75, 100, 125, 150, Inf),
                    labels = c("25-50", "51-75", "76-100", "101-125", "126-150", ">150")),
    affect = factor(recode(affect_type, POS = "Positive affect", NEG = "Negative affect"),
                    levels = c("Positive affect", "Negative affect")),
    theta_nonlin = theta_opt > 0,
    edf_nonlin   = edf_lag > 1.001
  )

# EDF of nonlinear series (overall, NA, PA)
rq2_results %>% filter(edf_lag > 1.001) %>%
  summarise(m = mean(edf_lag), .by = NULL)
# m
# 1 2.921784

rq2_results %>% filter(edf_lag > 1.001) %>%
  group_by(affect_type) %>% summarise(n = n(), m = mean(edf_lag), sd = sd(edf_lag))
# A tibble: 2 × 4
# affect_type     n     m    sd
# <chr>       <int> <dbl> <dbl>
#   1 NEG          1737  3.05  1.70
# 2 POS          1159  2.72  1.23

# Change points per length bin and affect (for 0.57/0.54 and 2.30/1.95)
rq2_d %>% group_by(n_obs_bin, affect) %>%
  summarise(m = mean(changepoint_n, na.rm = TRUE), .groups = "drop")

# # A tibble: 12 × 3
# n_obs_bin affect              m
# <fct>     <fct>           <dbl>
#   1 25-50     Positive affect 0.567
# 2 25-50     Negative affect 0.545
# 3 51-75     Positive affect 0.943
# 4 51-75     Negative affect 0.779
# 5 76-100    Positive affect 1.64
# 6 76-100    Negative affect 1.49
# 7 101-125   Positive affect 2.09
# 8 101-125   Negative affect 1.73
# 9 126-150   Positive affect 1.38
# 10 126-150   Negative affect 1.20
# 11 >150      Positive affect 2.30
# 12 >150      Negative affect 1.95

theme_len <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))



# -- Plotting --

# EDF
p_edf <- plot_prop(rq2_d, edf_nonlin, "Series with nonlinear lag-1 relation (%)\n(EDF > 1.001)")
p_edf

ggsave("figures/EDF_by_length_and_affect.png", p_edf,
       width = 7, height = 7, dpi = 500)

# Change points
p_cp <- plot_count(rq2_d, changepoint_n, "Mean number of change points")
p_cp

ggsave("figures/changepoints_by_length_and_affect.png", p_cp,
       width = 7, height = 7, dpi = 500)

# optimal theta
p_th_mean <- plot_count(rq2_d, theta_opt, expression("Mean optimal " * theta))
p_th_mean

ggsave("figures/opt_theta_by_length_and_affect.png", p_th_mean,
       width = 7, height = 7, dpi = 500)



# --------------------------------------
# MULTILEVEL MODELS --------------------


# Merge rq1_results and rq2_results
d2 <- rq1_results %>%
  transmute(ID, adcf = adcf_stat, n_obs) %>%
  left_join(transmute(rq2_results, ID, dataset, user_ID,
                      edf_lag, theta_opt, changepoint_n),
            by = "ID") %>%
  mutate(logT   = log(n_obs),
         person = interaction(dataset, user_ID, drop = TRUE))


# -- 1. Nonlinearity: ADCF ~ EDF_GAM

sum(!is.na(rq2_results$edf_lag)) # 15961
sum(rq2_results$edf_lag > 1.001, na.rm = TRUE) # 2896
sum(rq2_results$edf_lag <= 1.001, na.rm = TRUE) # 13065

# -- 1. Nonlinearity: ADCF ~ EDF_GAM (two-part)

d2_edf <- d2 %>% filter(!is.na(edf_lag), !is.na(adcf))

# 1a) Is the lag-1 relation nonlinear (EDF > 1.001)?
m_edf_bin <- glmmTMB(I(edf_lag > 1.001) ~ scale(adcf) + (1 | dataset/person),
                     family = binomial, data = d2_edf)

edf_bin_ci <- confint(m_edf_bin, parm = "theta_", method = "Wald")
edf_bin_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)
# 2.5 %    97.5 %  Estimate
# Std.Dev.(Intercept)|person:dataset 0.3150293 0.5123371 0.4017477
# Std.Dev.(Intercept)|dataset        0.1671717 0.6520988 0.3301703

summary(m_edf_bin)
# Family: binomial  ( logit )
# Formula:          I(edf_lag > 1.001) ~ scale(adcf) + (1 | dataset/person)
# Data: d2_edf
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 13065.1   13095.8   -6528.5   13057.1     15810
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.4017   0.6338
# dataset        (Intercept) 0.3302   0.5746
# Number of obs: 15814, groups:  person:dataset, 2077; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.76548    0.13945  -12.66   <2e-16 ***
#   scale(adcf)  0.78248    0.02533   30.89   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 1b) Among nonlinear series: how nonlinear (EDF - 1)?
d2_edf_nl <- d2_edf %>% filter(edf_lag > 1.001)

m_edf_nl <- glmmTMB(I(edf_lag - 1) ~ scale(adcf) + (1 | dataset/person),
                    family = Gamma(link = "log"), data = d2_edf_nl)

edf_nl_ci <- confint(m_edf_nl, parm = "theta_", method = "Wald")
edf_nl_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)

# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.062137017 0.10965585 0.08254506
# Std.Dev.(Intercept)|dataset        0.009288697 0.04448894 0.02032841

summary(m_edf_nl)

# Family: Gamma  ( log )
# Formula:          I(edf_lag - 1) ~ scale(adcf) + (1 | dataset/person)
# Data: d2_edf_nl
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 8424.0    8453.8   -4207.0    8414.0      2862
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.08255  0.2873
# dataset        (Intercept) 0.02033  0.1426
# Number of obs: 2867, groups:  person:dataset, 1430; dataset, 18
#
# Dispersion estimate for Gamma family (sigma^2): 0.397
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.50899    0.03766   13.52   <2e-16 ***
#   scale(adcf)  0.23802    0.01431   16.63   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -- 2. Regime shifts: ADCF ~ number of change points
d2_cp <- d2 %>% filter(!is.na(changepoint_n), !is.na(adcf))

m_cp <- glmmTMB(changepoint_n ~ scale(adcf) + scale(logT) + (1 | dataset/person),
                family = poisson(), data = d2_cp)

cp_ci <- confint(m_cp, parm = "theta_", method = "Wald")
cp_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)

# 2.5 %    97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.24208874 0.2978737 0.26853654
# Std.Dev.(Intercept)|dataset        0.02479192 0.1083065 0.05181821

summary(m_cp)
# Family: poisson  ( log )
# Formula:          changepoint_n ~ scale(adcf) + scale(logT) + (1 | dataset/person)
# Data: d2_cp
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 43337.9   43376.7  -21663.9   43327.9     17435
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.26854  0.5182
# dataset        (Intercept) 0.05182  0.2276
# Number of obs: 17440, groups:  person:dataset, 2085; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.239120   0.056105  -4.262 2.03e-05 ***
#   scale(adcf)  0.087653   0.008457  10.365  < 2e-16 ***
#   scale(logT)  0.474946   0.031945  14.868  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# -- 3. State dependence: ADCF ~ theta_opt

# Check distributions
table(rq2_results$theta_opt, useNA = "ifany")
mean(rq2_results$theta_opt == 0, na.rm = TRUE)   # 0.3207833
mean(rq2_results$theta_opt == 8, na.rm = TRUE)   # 0.2272948

# ---> Large clusters at 0 and 8.

d2_theta <- d2 %>% filter(!is.na(theta_opt), !is.na(adcf))

# 3a) State dependence present (theta > 0)?
m_theta_bin_adcf   <- glmmTMB(I(theta_opt > 0) ~ scale(adcf) + (1 | dataset/person),
                              family = binomial, data = d2_theta)

theta_bin_ci <- confint(m_theta_bin_adcf, parm = "theta_", method = "Wald")
theta_bin_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)

# 2.5 %     97.5 %    Estimate
# Std.Dev.(Intercept)|person:dataset 0.008966941 0.10407316 0.030548616
# Std.Dev.(Intercept)|dataset        0.002571788 0.02725362 0.008372008

summary(m_theta_bin_adcf)
# Family: binomial  ( logit )
# Formula:          I(theta_opt > 0) ~ scale(adcf) + (1 | dataset/person)
# Data: d2_theta
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 20517.6   20548.4  -10254.8   20509.6     16348
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.030549 0.1748
# dataset        (Intercept) 0.008372 0.0915
# Number of obs: 16352, groups:  person:dataset, 2048; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.76481    0.02855  26.784   <2e-16 ***
#   scale(adcf)  0.02587    0.01765   1.466    0.143
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 3b) Sensitivity: theta at the upper bound of the tested range (theta = 8)
m_theta_bin_adcf.8 <- glmmTMB(I(theta_opt == 8) ~ scale(adcf)  + (1 | dataset/person),
                              family = binomial, data = d2_theta)
summary(m_theta_bin_adcf.8)
# Family: binomial  ( logit )
# Formula:          I(theta_opt == 8) ~ scale(adcf) + (1 | dataset/person)
# Data: d2_theta
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 17544.7   17575.5   -8768.3   17536.7     16348
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.037050 0.19248
# dataset        (Intercept) 0.005062 0.07115
# Number of obs: 16352, groups:  person:dataset, 2048; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.23862    0.02701  -45.86   <2e-16 ***
#   scale(adcf) -0.02507    0.01988   -1.26    0.207
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 3c) Plateau only (0 < theta < 8): degree of state dependence
d2_theta_pl <- d2_theta %>% filter(theta_opt > 0, theta_opt < 8)

m_theta_gam_adcf   <- glmmTMB(theta_opt ~ scale(adcf)  + (1 | dataset/person),
                              family = Gamma(link = "log"), data = d2_theta_pl)
summary(m_theta_gam_adcf)
# Family: Gamma  ( log )
# Formula:          theta_opt ~ scale(adcf) + (1 | dataset/person)
# Data: d2_theta_pl
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 28235.8   28270.3  -14112.9   28225.8      7375
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance  Std.Dev.
# person:dataset (Intercept) 1.888e-08 0.0001374
# dataset        (Intercept) 8.550e-04 0.0292403
# Number of obs: 7380, groups:  person:dataset, 2008; dataset, 18
#
# Dispersion estimate for Gamma family (sigma^2): 0.638
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.9650418  0.0118103   81.71   <2e-16 ***
#   scale(adcf) -0.0005833  0.0092276   -0.06     0.95
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# ==============================================================================
# 4 RQ3 (exploratory): Person-level Characteristics
# ==============================================================================

# merge dataset
d3 <- rq1_results %>%
  transmute(ID, dataset, UUID, adcf = adcf_stat, n_obs) %>%
  left_join(person_inventory %>% select(dataset, UUID, user_ID, DEP_z, NEUR_z, SWL_z),
            by = c("dataset", "UUID")) %>%
  mutate(person = interaction(dataset, user_ID, drop = TRUE))

# sanity checks before modeling
range(d3$adcf, na.rm = TRUE)          # must be within [0,1] for ordbeta (true at max_lag=1)
d3 %>% distinct(dataset, user_ID) %>% count(dataset)   # persons per dataset


# --------------------------------------
# MULTILEVEL MODELS WITH ADCF ----------

# neuroticism
d3_neur <- d3 %>%
  filter(!is.na(NEUR_z), !is.na(adcf)) %>%
  mutate(adcf = pmin(pmax(adcf, 0), 1)) # ordbeta requires exact [0;1] range

m_neur <- glmmTMB(adcf ~ NEUR_z + (1 | dataset/person),
                  family = ordbeta, data = d3_neur)

neur_ci <- confint(m_neur, parm = "theta_", method = "Wald")
neur_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)

# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.09901990 0.11969965 0.10886987
# Std.Dev.(Intercept)|dataset        0.02057832 0.08710048 0.04233652

summary(m_neur)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ NEUR_z + (1 | dataset/person)
# Data: d3_neur
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -39795.5  -39742.2   19904.7  -39809.5     14939
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.10887  0.3300
# dataset        (Intercept) 0.04234  0.2058
# Number of obs: 14946, groups:  person:dataset, 1755; dataset, 16
#
# Dispersion parameter for ordbeta family (): 19.6
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.155576   0.052603  -40.98   <2e-16 ***
#   NEUR_z      -0.011570   0.009676   -1.20    0.232
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


d3_dep <- d3 %>%
  filter(!is.na(DEP_z), !is.na(adcf)) %>%
  mutate(adcf = pmin(pmax(adcf, 0), 1))

# depression
m_dep  <- glmmTMB(adcf ~ DEP_z  + (1 | dataset/person),
                  family = ordbeta, data = d3_dep)

dep_ci <- confint(m_dep, parm = "theta_", method = "Wald")
dep_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)

# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.09515413 0.11417020 0.10422939
# Std.Dev.(Intercept)|dataset        0.01806487 0.07322457 0.03637021

summary(m_dep)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ DEP_z + (1 | dataset/person)
# Data: d3_dep
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -43793.7  -43739.7   21903.8  -43807.7     16614
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.10423  0.3228
# dataset        (Intercept) 0.03637  0.1907
# Number of obs: 16621, groups:  person:dataset, 2008; dataset, 17
#
# Dispersion parameter for ordbeta family (): 19.2
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.131319   0.047308  -45.05  < 2e-16 ***
#   DEP_z       -0.023999   0.008952   -2.68  0.00735 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# satisfaction with life
d3_swl <- d3 %>%
  filter(!is.na(SWL_z), !is.na(adcf)) %>%
  mutate(adcf = pmin(pmax(adcf, 0), 1))

m_swl  <- glmmTMB(adcf ~ SWL_z  + (1 | dataset/person),
                  family = ordbeta, data = d3_swl)

swl_ci <- confint(m_swl, parm = "theta_", method = "Wald")
swl_ci^2 # variances (Wald CI on the log-SD scale, back-transformed)

# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.09907562 0.11993065 0.10900552
# Std.Dev.(Intercept)|dataset        0.01863046 0.09108693 0.04119455

summary(m_swl)

# Family: ordbeta  ( logit )
# Formula:          adcf ~ SWL_z + (1 | dataset/person)
# Data: d3_swl
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -38746.8  -38693.8   19380.4  -38760.8     14383
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.10901  0.3302
# dataset        (Intercept) 0.04119  0.2030
# Number of obs: 14390, groups:  person:dataset, 1727; dataset, 13
#
# Dispersion parameter for ordbeta family (): 20.4
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.153362   0.057281  -37.59   <2e-16 ***
#   SWL_z        0.016116   0.009715    1.66   0.0972 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
