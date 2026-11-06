# ==============================================================================
# Supplementary analyses RQ2: state dependence (optimal S-map theta) under a
# linear AR(1) process
#
# How often does optimal_theta_from_SMap() return theta > 0 when the data are a
# linear Gaussian AR(1), i.e., when there is no state dependence? Compared with
# the empirical rates (theta > 0: 67.9%, theta = 0: 32.1%, theta = 8: 22.7%).
#
# Data: all series (dat$y, exact AR(1)) of the power simulation's Type I grid
#   (power_sim/data/ar1_data_21_08_2026_12_54.RData -> ar1_grid).
#
# Contents:
#   1. Optimal theta for every AR(1) series (or load archived results)
#   2. Rates of theta > 0, theta = 0 and theta = 8 by n_obs and phi
#   3. Comparison with the empirical series, by length bin
#   4. Figure: theta distribution, simulated AR(1) vs. empirical
#
# Run from the repository root (e.g., open the .Rproj file in RStudio).
# ==============================================================================

library(tidyverse)
library(rEDM)
library(furrr)
library(future)

source("helpers/optimal_theta_from_SMap.R")


AR1_DATA_FILE <- "power_sim/data/ar1_data_21_08_2026_12_54.RData"   # -> ar1_grid
RQ2_FILE      <- "results/rq2_results_2026-10-13_112850.Rdata"      # -> rq2_results

# TRUE: compute theta and save a new timestamped file.
# FALSE: load the archived results in THETA_AR1_FILE.
RUN_THETA_AR1  <- TRUE
THETA_AR1_FILE <- "results/supp_rq2_theta_AR1_null_<stamp>.Rdata"     # set after first run


# ------------------------------------------------------------------------------
# 1. Optimal theta for every AR(1) series --------------------------------------
# ------------------------------------------------------------------------------

load(AR1_DATA_FILE)   # -> ar1_grid

# n_obs = 22 is dropped: no empirical series is that short (n_obs >= 25), and
# with a 15-point library embeddings up to E = 15 leave
# almost no library points.
ar1_run <- ar1_grid %>%
  as_tibble() %>%
  filter(n_obs > 22) %>%
  select(sim_id, n_obs, phi, dyn_er, meas_er, rep, dat)

# A crash inside rEDM's C++ code kills the worker process ("Future ...
# interrupted"), which tryCatch() cannot catch and which aborts the whole
# future_map(). The run is therefore done in chunks (one per n_obs). If a
# chunk fails, its series are rerun one by one, each in its own worker, so a
# crashing series returns theta = NA (crashed = TRUE) instead of stopping the run.
theta_chunk <- function(d) {
  out <- tryCatch(
    future_map_dbl(d$dat, function(x) optimal_theta_from_SMap(x$y),
                   .options = furrr_options(seed = 2026, packages = "rEDM")),
    error = function(e) NULL)
  if (!is.null(out)) return(tibble(theta = out, crashed = FALSE))

  message("  chunk failed; rerunning ", nrow(d), " series one by one")
  res <- lapply(d$dat, function(x) {
    f <- future(optimal_theta_from_SMap(x$y), seed = 2026, packages = "rEDM")
    tryCatch(c(theta = value(f), crashed = 0),
             error = function(e) c(theta = NA_real_, crashed = 1))
  })
  tibble(theta   = vapply(res, `[[`, numeric(1), "theta"),
         crashed = vapply(res, `[[`, numeric(1), "crashed") == 1)
}

if (RUN_THETA_AR1) {

  plan(multisession, workers = max(1, availableCores() - 1))

  theta_ar1 <- ar1_run %>%
    group_split(n_obs) %>%
    map_dfr(function(d) {
      message("n_obs = ", d$n_obs[1], " (", nrow(d), " series)")
      bind_cols(select(d, -dat), theta_chunk(d))
    })

  plan(sequential)

  table(theta_ar1$crashed)   # number of series on which the S-map crashed

  THETA_AR1_FILE <- paste0("results/supp_rq2_theta_AR1_null_",
                           format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".Rdata")
  save(theta_ar1, file = THETA_AR1_FILE)
  message("Saved: ", THETA_AR1_FILE)

} else {
  load(THETA_AR1_FILE)   # -> theta_ar1
}


# ------------------------------------------------------------------------------
# 2. Rates under the AR(1) process ---------------------------------------------
# ------------------------------------------------------------------------------

summ_theta <- function(d) {
  d %>% summarise(
    n_series      = n(),
    na_pct        = mean(is.na(theta)) * 100,
    theta_gt0_pct = mean(theta > 0,  na.rm = TRUE) * 100,
    theta_0_pct   = mean(theta == 0, na.rm = TRUE) * 100,
    theta_8_pct   = mean(theta == 8, na.rm = TRUE) * 100,
    .groups = "drop")
}

# overall
theta_ar1 %>% summ_theta() 
# A tibble: 1 × 5
# n_series na_pct theta_gt0_pct theta_0_pct theta_8_pct
# <int>  <dbl>         <dbl>       <dbl>       <dbl>
#   1     8400  0.143          70.9        29.1        11.2

theta_ar1 %>% group_by(n_obs) %>% summ_theta() %>% print(n = Inf)
# A tibble: 7 × 6
# n_obs n_series na_pct theta_gt0_pct theta_0_pct theta_8_pct
# <dbl>    <int>  <dbl>         <dbl>       <dbl>       <dbl>
#   1    43     1200      1          73.1        26.9       18.5 
# 2    64     1200      0          70.8        29.2       14.4 
# 3    85     1200      0          73.3        26.7       13.7 
# 4   127     1200      0          73          27         10.2 
# 5   169     1200      0          71.7        28.3        8.58
# 6   190     1200      0          68.1        31.9        6.83
# 7   253     1200      0          66.8        33.2        6.58

theta_ar1 %>% group_by(n_obs, phi) %>% summ_theta() %>% print(n = Inf)
# A tibble: 14 × 7
# n_obs   phi n_series na_pct theta_gt0_pct theta_0_pct theta_8_pct
# <dbl> <dbl>    <int>  <dbl>         <dbl>       <dbl>       <dbl>
#   1    43   0.3      600  1.83           75.0        25.0       20.4 
# 2    43   0.6      600  0.167          71.1        28.9       16.7 
# 3    64   0.3      600  0              74.2        25.8       19.2 
# 4    64   0.6      600  0              67.3        32.7        9.67
# 5    85   0.3      600  0              76.5        23.5       16.8 
# 6    85   0.6      600  0              70.2        29.8       10.5 
# 7   127   0.3      600  0              76.3        23.7       15.2 
# 8   127   0.6      600  0              69.7        30.3        5.17
# 9   169   0.3      600  0              74.7        25.3       13.3 
# 10   169   0.6      600  0              68.7        31.3        3.83
# 11   190   0.3      600  0              73.2        26.8       10.3 
# 12   190   0.6      600  0              63          37          3.33
# 13   253   0.3      600  0              72.2        27.8       10.7 
# 14   253   0.6      600  0              61.3        38.7        2.5 

# ------------------------------------------------------------------------------
# 3. Comparison with the empirical series --------------------------------------
# ------------------------------------------------------------------------------
# Empirical series are binned around the simulated lengths (limits at the
# midpoints between adjacent grid values). Empirical n_obs >= 25; n_obs = 22 is not
# simulated (see 1.).

load(RQ2_FILE)   # -> rq2_results

grid_n <- sort(unique(theta_ar1$n_obs))
breaks <- c(-Inf, head(grid_n, -1) + diff(grid_n) / 2, Inf)

emp <- rq2_results %>%
  transmute(theta = theta_opt,
            n_obs = grid_n[cut(n_obs, breaks, labels = FALSE)])

bind_rows(
  theta_ar1 %>% group_by(n_obs) %>% summ_theta() %>% mutate(source = "AR(1)"),
  emp       %>% group_by(n_obs) %>% summ_theta() %>% mutate(source = "empirical")
) %>%
  arrange(n_obs, source) %>%
  print(n = Inf)

# A tibble: 13 × 7
# n_obs n_series na_pct theta_gt0_pct theta_0_pct theta_8_pct source   
# <dbl>    <int>  <dbl>         <dbl>       <dbl>       <dbl> <chr>    
#   1    43     1200  1              73.1        26.9       18.5  AR(1)    
# 2    43     4032 16.1            67.4        32.6       22.5  empirical
# 3    64     1200  0              70.8        29.2       14.4  AR(1)    
# 4    64     9848  4.69           68.9        31.1       23.5  empirical
# 5    85     1200  0              73.3        26.7       13.7  AR(1)    
# 6    85     2130  3.24           67.2        32.8       22.5  empirical
# 7   127     1200  0              73          27         10.2  AR(1)    
# 8   127      712  1.54           66.2        33.8       18.5  empirical
# 9   169     1200  0              71.7        28.3        8.58 AR(1)    
# 10   169      530  1.70           59.9        40.1       20.3  empirical
# 11   190     1200  0              68.1        31.9        6.83 AR(1)    
# 12   190      445  0.225          66.0        34.0       18.7  empirical
# 13   253     1200  0              66.8        33.2        6.58 AR(1)  

# ------------------------------------------------------------------------------
# 4. Figure: theta distribution, AR(1) vs. empirical ---------------------------
# ------------------------------------------------------------------------------

d_fig <- bind_rows(
  theta_ar1 %>% transmute(source = "Simulated AR(1)", theta),
  emp       %>% transmute(source = "Empirical", theta)
) %>%
  filter(!is.na(theta)) %>%
  count(source, theta) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_theta <- ggplot(d_fig, aes(theta, pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.45), width = 0.45) +
  scale_x_continuous(breaks = 0:8) +
  scale_fill_manual(values = c("Empirical" = "#0072B2",
                               "Simulated AR(1)" = "grey60"),
                    name = NULL) +
  labs(x = expression("Optimal " * theta), y = "Series (%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
p_theta

ggsave("figures/theta_AR1_null_vs_empirical.png", p_theta,
       width = 7, height = 4.5, dpi = 300, bg = "white")
