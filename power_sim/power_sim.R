##########################################
### Test Power Simulation (BDS / Srho / gspec / ADCF / HSIC)
##########################################
# Run from the repository root (e.g., open the .Rproj file in RStudio).
# By default (all RUN_* flags FALSE) the script loads the archived simulated
# data and results listed under FILES and reproduces the reported tables and
# figures. Set the RUN_* flags to TRUE to regenerate data and/or results;
# regenerated files are saved with a new timestamp and never overwrite the
# archived files.

# packages
library(future.apply)
library(progressr)
library(dplyr)
library(data.table)
library(ggplot2)
library(furrr)
library(energy)   # ADCF (via energy::dcor)
library(dHSIC)    # HSIC statistic
library(tidyr)    # pivot_wider() in "Tables for Word"
library(purrr)    # map_chr() in "Tables for Word"
# Also required (used inside the test helpers): tseries, tseriesEntropy

handlers(global = TRUE)

# Use data generation scripts from Failenschmid et al. (2025),
# DOI: 10.1111/bmsp.12397
source("power_sim/Failenschmid-et-al-2025/classes.R")
source("power_sim/Failenschmid-et-al-2025/functions.R")

# import tests
source("helpers/bds.R")
source("helpers/srho.R")
source("helpers/gspec.R")
source("helpers/hsic.R")
source("helpers/adcf.R")

# import refit bootstrapping algorithm
source("helpers/rb_refit.R")

# import AR generation
source("helpers/null_AR1.R")


### ==========================================================================
### FILES: archived inputs/outputs behind the reported results
### ==========================================================================
AR1_DATA_FILE  <- "power_sim/data/ar1_data_21_08_2026_12_54.RData"   # Type I grid (AR(1))
SIM_DATA_FILE  <- "power_sim/data/sim_data_08_07_2026_17_37.Rdata"   # power grid (nonlinear DGPs)
POWER_RES_FILE <- "power_sim/results/results_rb_refit_23_10_2026_19_45.rds"
TYPE1_RES_FILE <- c(
  rb_refit = "power_sim/results/ar1_results_rb_refit_23_10_2026_00_01.rds",
  rp       = "power_sim/results/ar1_results_rp_29_10_2026_23_38.rds"
)

### Regenerate simulated data? (FALSE: load the archived files above)
RUN_SIMDATA <- FALSE


### ==========================================================================
### Extraction helpers
### ==========================================================================
# Normalize any test helper's return into c(statistic, p.value).
#   - New-shape helpers (bds, adcf, hsic) return list(statistic, p.value).
#   - Legacy helpers (srho, gspec) still return a scalar p-value; for those
#     the statistic comes through as NA.

.extract_sp <- function(x) {
  if (is.list(x)) {
    stat <- x[["statistic"]]
    pv   <- x[["p.value"]]
    c(statistic = if (length(stat)) as.numeric(stat)[1] else NA_real_,
      p.value   = if (length(pv))   as.numeric(pv)[1]   else NA_real_)
  } else {
    c(statistic = NA_real_, p.value = as.numeric(x)[1])
  }
}


#   null_type = "rp"       : old null, residual permutation WITHOUT refit (helpers/*.R)
#   null_type = "rb_refit" : new null, AR(1) residual bootstrap WITH refit (helpers/rb_refit.R);
#                            all requested tests share the same B bootstrap series
#   tests                  : subset of ALL_TESTS to run
# null_type and tests are deliberately passed in explicitly (no defaults that point
# to global variables) so that `future` exports them to the workers.
ALL_TESTS <- c("bds", "srho", "gspec", "adcf", "hsic")

run_tests <- function(x, B, seed, null_type, tests) {
  null_type <- match.arg(null_type, c("rp", "rb_refit"))
  tests     <- intersect(ALL_TESTS, tests)

  raw <- if (null_type == "rb_refit") {
    tryCatch(test_rb_refit(x, tests = tests, B = B, seed = seed, estimator = "CLS"),
             error = function(e) NULL)
  } else {
    # OLD null: permute CLS residuals, no refit
    x_lag <- c(NA, x[-length(x)])
    fit   <- tryCatch(lm(x ~ x_lag, na.action = na.exclude), error = function(e) NULL)
    if (is.null(fit)) return(na_tests())
    res <- as.numeric(residuals(fit))
    specs <- list(
      bds   = function() test_bds  (res, input = "residual", m = 2, eps = 1, B = B, seed = seed),
      srho  = function() test_srho (res, input = "residual", max_lag = 1, B = B, seed = seed),
      gspec = function() test_gspec(res, input = "residual", B = B, seed = seed),
      adcf  = function() test_adcf (res, input = "residual", na_handling = "pairwise", max_lag = 1, B = B, seed = seed),
      hsic  = function() test_hsic (res, input = "residual", na_handling = "pairwise", max_lag = 1, B = B, seed = seed)
    )
    setNames(lapply(tests, function(nm) tryCatch(specs[[nm]](), error = function(e) NA_real_)), tests)
  }

  out <- setNames(rep(NA_real_, 2 * length(ALL_TESTS)),
                  as.vector(rbind(paste0("stat.", ALL_TESTS), paste0("p.", ALL_TESTS))))
  for (nm in tests) {
    sp <- .extract_sp(if (is.null(raw)) NA_real_ else raw[[nm]])
    out[paste0("stat.", nm)] <- sp[["statistic"]]
    out[paste0("p.",    nm)] <- sp[["p.value"]]
  }
  out
}

# All-NA battery result with identical names/order to run_tests() (kept for
# compatibility; the latent series is now skipped entirely when RUN_LATENT = FALSE).
na_tests <- function() {
  nm <- as.vector(rbind(paste0("stat.", ALL_TESTS), paste0("p.", ALL_TESTS)))
  setNames(rep(NA_real_, length(nm)), nm)
}


### ==========================================================================
### Simulation helpers
### ==========================================================================

### Data-only simulation function (adapted from Failenschmid et al., 2025)
simulate <- function(
    gen_model_list, conditions, repetitions,
    out_dir = getwd(), seed = 54321) {

  require(data.table)
  require(future.apply)

  set.seed(seed)

  sim_grid <- as.data.table(
    expand.grid(c(conditions, list(gen_model = gen_model_list)))
  )
  sim_grid$conditions <- apply(
    subset(sim_grid, select = -gen_model),
    1, function(x) as.list(x)
  )
  sim_grid[, gen_model := mapply(add_conditions, gen_model, conditions)]
  sim_grid <- sim_grid[rep(seq_len(nrow(sim_grid)), each = repetitions), ]

  with_progress({
    p <- progressor(along = sim_grid$gen_model)
    sim_grid[, dat := future_lapply(
      gen_model,
      function(x) {
        p()
        get_tsm_data(sim_tsm(x))
      },
      future.seed = TRUE
    )]
  })

  sim_time <- format(Sys.time(), "%d_%m_%Y_%H_%M")
  save(sim_grid, file = paste0(out_dir, "/sim_data_", sim_time, ".Rdata"))

  cat("\n\nDone. Saved to", out_dir, "\n")
  return(sim_grid)
}

### Define generative models (identical to Failenschmid et al., 2025) ----------
exp_growth <- new("gen_model",
                  model_name = "exp_growth",
                  model_type = "DE",
                  time = 200,
                  pars = list(yr = 0.02, ya = 2, dyn_er = 1, meas_er = 1),
                  delta = (50/7)/234,
                  stepsize = (50/7)/3,
                  start = list(formula(y ~ -2)),
                  state_eq = list(formula(y ~ yr * ya - yr * y)),
                  dynamic_error = list(formula(y ~ dyn_er)),
                  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er)))
)

log_growth <- new("gen_model",
                  model_name = "log_growth",
                  model_type = "DE",
                  time = 200,
                  pars = list(k = 4.3, r = 0.04, dyn_er = 1, meas_er = 1),
                  delta = (50/7)/234,
                  stepsize = (50/7)/3,
                  start = list(formula(y ~ 0.3)),
                  state_eq = list(formula(y ~ r * y * (1 - (y / k)))),
                  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er))),
                  dynamic_error = list(formula(y ~ dyn_er)),
                  boundary = function(x) ifelse(x >= 0, x, 0)
)

cusp_catastrophe <- new("gen_model",
                        model_name = "cusp_catastrophe",
                        model_type = "DE",
                        time = 200,
                        delta = (50/7)/234,
                        stepsize = (50/7)/3,
                        pars = list(a = -5, dyn_er = 1, omega = -(2*pi/50)^2, meas_er = 1),
                        start = list(formula(y ~ 1.9), formula(v ~ 0), formula(b ~ -10)),
                        state_eq = list(
                          formula(y ~ -(4*y^3 + 2*a*y + b)),
                          formula(b ~ v),
                          formula(v ~ omega * b)
                        ),
                        dynamic_error = list(
                          formula(y ~ dyn_er),
                          formula(b ~ 0),
                          formula(v ~ 0)
                        ),
                        meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er)))
)

damped_oscillator <- new("gen_model",
                         model_name = "damped_oscillator",
                         model_type = "DE",
                         time = 200,
                         delta = (50/7)/234,
                         stepsize = (50/7)/3,
                         pars = list(k = 0.01, c = 0.1, dyn_er = 1, meas_er = 1),
                         start = list(formula(y ~ 2), formula(v ~ 0)),
                         state_eq = list(
                           formula(y ~ v),
                           formula(v ~ -2*k*v - c^2*y)
                         ),
                         dynamic_error = list(
                           formula(y ~ dyn_er),
                           formula(v ~ 0)
                         ),
                         meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er)))
)

### Define conditions ----------------------------------------------------------
conditions <- list(
  time     = c(50, 100, 150, 200),             # 1, 2, 3, 4 weeks
  stepsize = c((50/7)/3, (50/7)/6, (50/7)/9),  # 3, 6, 9 obs/day
  dyn_er   = sqrt(c(0, 1)),
  meas_er  = sqrt(c(0.25, 1))
)

### Simulate -------------------------------------------------------------------
if (RUN_SIMDATA) plan(multisession, workers = availableCores() - 1)

if (RUN_SIMDATA) {
  set.seed(2310)

  # Type I grid: exact AR(1) series (y) and AR(1) + measurement error (y_obs)
  ar1_grid <- make_ar1_grid(
    n_obs   = c(22, 43, 64, 85, 127, 169, 190, 253),
    phi     = c(0.3, 0.6),
    dyn_er  = sqrt(c(0.25, 1, 2)),
    meas_er = sqrt(c(0.25, 1)),
    repetitions = 100,
    seed    = 2310
  )
  ar_sim_time <- format(Sys.time(), "%d_%m_%Y_%H_%M")
  save(ar1_grid, file = paste0("power_sim/data/ar1_data_", ar_sim_time, ".RData"))

  # Power grid: nonlinear data-generating models (saved inside simulate())
  sim <- simulate(
    gen_model_list = list(
      exp_growth, log_growth, damped_oscillator, cusp_catastrophe
    ),
    conditions  = conditions,
    repetitions = 100,
    out_dir     = "power_sim/data"
  )
} else {
  load(AR1_DATA_FILE)   # -> ar1_grid
  load(SIM_DATA_FILE)   # -> sim_grid
  sim <- sim_grid
}


### ==========================================================================
### Apply tests
### ==========================================================================

dir.create("power_sim/results", recursive = TRUE, showWarnings = FALSE)

.pkgs <- c("tseries", "tseriesEntropy", "energy", "dHSIC")


### --- RUN CONFIGURATION -------------------------------------------------------
NULL_TYPE  <- "rb_refit"
TESTS_RUN  <- c("bds", "adcf", "hsic", "gspec", "srho")
MODELS_RUN <- c("exp_growth", "log_growth", "damped_oscillator", "cusp_catastrophe")
B_PERM     <- 199
RUN_LATENT <- FALSE                       # FALSE: skip dat$y in the power grid

# Set to TRUE to regenerate results with NULL_TYPE (computationally expensive).
# With FALSE, the script loads the archived results listed under FILES.
RUN_TYPE1 <- FALSE
RUN_POWER <- FALSE
# Type I grid: restrict if needed
AR1_N_OBS  <- sort(unique(ar1_grid$n_obs))    # e.g. c(64, 127)
AR1_REPS   <- 100                         # replications per cell (<= 100)
AR1_ARMA   <- FALSE                        # also run dat$y_obs (ARMA(1,1))?


stamp <- function() format(Sys.time(), "%d_%m_%Y_%H_%M")


### --- Type I error grid (AR(1) null) ----------------------------------------
### Tested on dat$y (the clean AR(1) series) -> residuals are i.i.d. under the
### null, so empirical rejection rate here IS each test's empirical alpha.

if (RUN_TYPE1) {
  ar1_run <- ar1_grid %>% filter(n_obs %in% AR1_N_OBS, rep <= AR1_REPS)
  message("Type I grid: ", nrow(ar1_run), " series x ", length(TESTS_RUN), " tests x ",
          if (AR1_ARMA) 2 else 1, " series types; null_type = ", NULL_TYPE)
  start <- Sys.time()

  set.seed(2310)
  with_progress({
    p <- progressor(steps = nrow(ar1_run))
    ar1_results <- ar1_run %>%
      future_pmap(
        .options = furrr_options(seed = 2310, packages = .pkgs),
        function(n_obs, phi, dyn_er, meas_er, dat, gen_model, rep, sim_id, ...) {
          p()
          meta <- data.frame(null_type = NULL_TYPE, phi = phi, dyn_er = dyn_er,
                             meas_er = meas_er, n_obs = n_obs, rep = rep, sim_id = sim_id)

          # y     : latent, exactly AR(1)              -> size
          # y_obs : latent + white noise = ARMA(1,1)   -> power vs linear misspecification
          res_ar <- as.data.frame(as.list(run_tests(dat$y, B = B_PERM, seed = sim_id,
                                                    null_type = NULL_TYPE, tests = TESTS_RUN)))
          out <- cbind(gen_model = "ar1", meta, res_ar)
          if (AR1_ARMA) {
            res_arma <- as.data.frame(as.list(run_tests(dat$y_obs, B = B_PERM, seed = sim_id,
                                                        null_type = NULL_TYPE, tests = TESTS_RUN)))
            out <- rbind(out, cbind(gen_model = "arma11", meta, res_arma))
          }
          out
        }
      ) %>% bind_rows()
  })
  print(Sys.time() - start)

  arma_string <- ifelse(AR1_ARMA, "_arma_", "_")

  saveRDS(ar1_results,
          sprintf("power_sim/results/ar1%sresults_%s_%s.rds", arma_string, NULL_TYPE, stamp()))
}


### --- Power grid (nonlinear DGPs) --------------------------------------------
if (RUN_POWER) {
  # row_id = position in the FULL grid, assigned BEFORE filtering, so seeds stay
  # identical no matter which models/tests are run.
  sim$row_id <- seq_len(nrow(sim))
  sim$model  <- vapply(sim$gen_model, function(g) slot(g, "model_name"), character(1))
  sim_run    <- sim[sim$model %in% MODELS_RUN, ]
  message("Power grid: ", nrow(sim_run), " series x ", length(TESTS_RUN),
          " tests; null_type = ", NULL_TYPE, "; latent = ", RUN_LATENT)
  start <- Sys.time()

  set.seed(2311)
  with_progress({
    p <- progressor(steps = nrow(sim_run))

    results <- sim_run %>%
      future_pmap(
        .options = furrr_options(seed = TRUE, packages = .pkgs),
        function(time, stepsize, dyn_er, meas_er, gen_model, conditions, dat, row_id, model, ...) {
          p()

          # Observed series (with measurement error).
          res_obs <- run_tests(dat$y_obs, B = B_PERM, seed = row_id,
                               null_type = NULL_TYPE, tests = TESTS_RUN)
          names(res_obs) <- paste0(names(res_obs), "_obs")

          out <- cbind(
            data.frame(null_type = NULL_TYPE, model = model, time = time,
                       stepsize = stepsize, dyn_er = dyn_er, meas_er = meas_er,
                       n_obs = round(time / stepsize) + 1, row_id = row_id),
            as.data.frame(as.list(res_obs))
          )

          # Latent series (no measurement error). dat$y does not depend on meas_er,
          # so compute it once (at meas_er == 0.5). Skipped entirely unless RUN_LATENT.
          if (RUN_LATENT) {
            res_y <- if (isTRUE(all.equal(meas_er, 0.5)))
              run_tests(dat$y, B = B_PERM, seed = row_id,
                        null_type = NULL_TYPE, tests = TESTS_RUN) else na_tests()
            names(res_y) <- paste0(names(res_y), "_y")
            out <- cbind(out, as.data.frame(as.list(res_y)))
          }
          out
        }
      ) %>% bind_rows()
  })

  end <- Sys.time()
  print("Power tests:"); print(end - start)

  saveRDS(results, sprintf("power_sim/results/results_%s_%s.rds", NULL_TYPE, stamp()))
}


### ==========================================================================
### Inspect results
### ==========================================================================

# Power: refit bootstrap null only (rb_refit). Type I: both nulls are archived;
# INSPECT_NULL selects which one is inspected below. To inspect regenerated
# results, point the
# FILES entries at the new timestamped files.
INSPECT_NULL <- "rb_refit"   # "rb_refit" (main) or "rp" (residual permutation, Table S6)
results     <- readRDS(POWER_RES_FILE)
ar1_results <- readRDS(TYPE1_RES_FILE[[INSPECT_NULL]])

### --- Reshape helper ---------------------------------------------------------

to_long <- function(dt, id.vars) {
  dt <- as.data.table(copy(dt))
  if (!"row_id" %in% names(dt)) dt[, row_id := .I]
  id.vars <- union("row_id", id.vars)

  p_cols <- grep("^p\\.", names(dt), value = TRUE)
  labels <- sub("^p\\.", "", p_cols)

  long <- melt(
    dt,
    id.vars       = id.vars,
    measure.vars  = patterns(p = "^p\\.", statistic = "^stat\\."),
    variable.name = "col_idx"
  )
  long[, label  := labels[as.integer(col_idx)]]
  long[, series := fifelse(grepl("_obs$", label), "observed", "latent")]
  long[, test   := sub("_(obs|y)$", "", label)]
  long[, reject := p <= .05]
  long[]
}

### --- Type I error (AR(1) null) ---------------------------------------------
# the AR(1) vs ARMA(1,1) split is in gen_model.
ar1_long <- to_long(
  ar1_results,
  id.vars = c("gen_model", "phi", "rep", "sim_id", "dyn_er", "meas_er", "n_obs")
)

ar1_alpha <- ar1_long[!is.na(p),
                      .(alpha     = mean(reject),
                        se        = sqrt(mean(reject) * (1 - mean(reject)) / .N),
                        n_rep     = .N,
                        mean_stat = mean(statistic, na.rm = TRUE)),
                      by = .(gen_model, phi, dyn_er, meas_er, n_obs, test)]

# Quick size screen: which cells stray from nominal .05 for the new tests?
ar1_alpha[test %in% c("adcf", "hsic") & abs(alpha - .05) > .02][order(-abs(alpha - .05))]

n_obs_select <- c(22, 64, 127, 190, 253)

print(as.data.frame(ar1_alpha[ar1_alpha$n_obs %in% n_obs_select,]))

### --- Power (nonlinear DGPs) -------------------------------------------------
res_long <- to_long(
  results,
  id.vars = c("model", "time", "stepsize", "dyn_er", "meas_er", "n_obs")
)

res_power <- res_long[!is.na(p),
                      .(power     = mean(reject),
                        se        = sqrt(mean(reject) * (1 - mean(reject)) / .N),
                        n_rep     = .N,
                        n_fail    = sum(is.na(p)),
                        mean_stat = mean(statistic, na.rm = TRUE)),   # NA for srho/gspec until updated
                      by = .(model, time, stepsize, dyn_er, meas_er, n_obs, series, test)]


### --- Plot ------------------------------------------------------------------
plot_power <- function(agg, step, dyn, series_type = "observed") {
  dat <- agg[stepsize == step & dyn_er == dyn & series == series_type]
  ggplot(dat, aes(x = n_obs, y = power, colour = test, group = test)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.7) +
    geom_pointrange(aes(ymin = pmax(0, power - 1.96 * se),
                        ymax = pmin(1, power + 1.96 * se)),
                    size = 0.3) +
    facet_grid(model ~ meas_er,
               labeller = labeller(meas_er = function(x) paste("SD(meas_er) =", x))) +
    scale_colour_manual(values = c(bds   = "#0271bb", srho  = "#d55e00",
                                   gspec = "#009e73", adcf  = "#cc79a7",
                                   hsic  = "#e69f00")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Number of Observations", y = "Power", colour = "Test") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

steps   <- sort(unique(res_power$stepsize), decreasing = TRUE)
dyn_ers <- sort(unique(res_power$dyn_er),   decreasing = TRUE)

affect_models <- c("cusp_catastrophe", "damped_oscillator")

res_power.affect_models <- res_power[res_power$model %in% affect_models]

p <- plot_power(res_power.affect_models, step = steps[3], dyn = dyn_ers[1], series_type = "observed") # 9/day, dyn 1

# export
ggsave("power_sim/figures/power_9perday_dyn1.png", plot = p,
       width = 8, height = 6, units = "in", dpi = 300, bg = "white")



### --- Tables for Word --------------------------------------------------------
obs_tables <- res_long |>
  filter(series == "observed") |>
  as_tibble() |>
  summarise(power = mean(reject, na.rm = TRUE),
            .by = c(model, dyn_er, meas_er, time, stepsize, n_obs, test)) |>
  mutate(test = recode(test,
                       bds = "BDS", srho = "Srho", gspec = "Hong-Lee",
                       adcf = "ADCF", hsic = "HSIC")) |>
  pivot_wider(names_from = test, values_from = power) |>
  arrange(model, dyn_er, meas_er, time, stepsize) |>
  group_split(model)

names(obs_tables) <- map_chr(obs_tables, ~ .x$model[1])

# Optional: the same table for the observed test STATISTIC (mean across reps).
obs_stat_tables <- res_long |>
  filter(series == "observed") |>
  as_tibble() |>
  summarise(statistic = mean(statistic, na.rm = TRUE),
            .by = c(model, dyn_er, meas_er, time, stepsize, n_obs, test)) |>
  mutate(test = recode(test,
                       bds = "BDS", srho = "Srho", gspec = "Hong-Lee",
                       adcf = "ADCF", hsic = "HSIC")) |>
  pivot_wider(names_from = test, values_from = statistic) |>
  arrange(model, dyn_er, meas_er, time, stepsize) |>
  group_split(model)

names(obs_stat_tables) <- map_chr(obs_stat_tables, ~ .x$model[1])
