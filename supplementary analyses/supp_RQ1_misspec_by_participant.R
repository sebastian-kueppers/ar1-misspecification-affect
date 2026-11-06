# ==============================================================================
# Supplementary analyses RQ1: AR(1) misspecification per participant
#
# Exploratory and descriptive (not preregistered). Describes how series
# classified as AR(1)-misspecified (p <= .05) are distributed across
# participants, to check whether a subset of participants accounts for the
# overall misspecification rate.
#
# Observed: empirical cumulative distribution (ECDF) of each participant's
# proportion of misspecified series, per test (testable series of that test).
# Chance reference: the same distribution if every series of a participant had
# its dataset's misspecification rate, independently of the participant
# (binomial with the participant's number of series; R_SIM simulated samples).
# No significance test.
#
# Contents:
#   1. Saved RQ1 results
#   2. Parameters
#   3. Per-participant proportions
#   4. Figure (ECDF, observed vs. expected by chance)
#   5. Table (% of participants per range, observed vs. expected)
#   6. Key numbers for the text
# ==============================================================================

library(tidyverse)


# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ECDF of v evaluated on a fixed grid
ecdf_on_grid <- function(v, grid) ecdf(v)(grid)


# ==============================================================================
# 1 Saved RQ1 results
# ==============================================================================

# RQ1 results reported in the main text
load("results/rq1_results_2026-10-22_125950.Rdata")


# ==============================================================================
# 2 Parameters
# ==============================================================================

ALPHA <- .05                    # classification rule, as in the main analysis
R_SIM <- 500                    # simulated reference samples
SEED  <- 2026
GRID  <- seq(0, 1, by = 0.005)  # x-grid for the ECDFs

# ranges for the table (requested by Pete)
BRKS   <- c(-Inf, .1, .2, .3, .4, .5, Inf)
LABS_R <- c("0–9%", "10–19%", "20–29%", "30–39%", "40–49%", "≥ 50%")

col_ref <- c("Observed" = "#0072B2", "Expected by chance" = "grey40")

set.seed(SEED)


# ==============================================================================
# 3 Per-participant proportions
# ==============================================================================

pp <- rq1_results %>%
  dplyr::select(dataset, UUID, ADCF = adcf_p, BDS = bds_p) %>%
  pivot_longer(c(ADCF, BDS), names_to = "Test", values_to = "p") %>%
  filter(!is.na(p)) %>%                            # testable series per test
  group_by(Test, dataset) %>%
  mutate(rate_ds = mean(p <= ALPHA)) %>%           # dataset-specific rate
  group_by(Test, dataset, UUID) %>%
  summarise(k    = n(),                            # series of this participant
            x    = sum(p <= ALPHA),                # misspecified series
            prop = x / k,
            rate = first(rate_ds),
            .groups = "drop")

# participants and series per participant (for the figure note)
pp %>%
  group_by(Test) %>%
  summarise(n_persons = n(), k_mdn = median(k), k_min = min(k), k_max = max(k))


# ==============================================================================
# 4 Figure: ECDF, observed vs. expected by chance
# ==============================================================================

obs <- pp %>%
  group_by(Test) %>%
  reframe(prop_grid = GRID, F = ecdf_on_grid(prop, GRID)) %>%
  mutate(source = "Observed")

sim <- pp %>%
  group_by(Test) %>%
  reframe({
    Fs <- replicate(R_SIM, ecdf_on_grid(rbinom(n(), k, rate) / k, GRID))
    tibble(prop_grid = GRID,
           F  = rowMeans(Fs),
           lo = apply(Fs, 1, quantile, .025),
           hi = apply(Fs, 1, quantile, .975))
  }) %>%
  mutate(source = "Expected by chance")

p_pp <- ggplot(mapping = aes(x = prop_grid, y = F, colour = source)) +
  geom_ribbon(data = sim, aes(ymin = lo, ymax = hi), fill = "grey70",
              alpha = 0.4, colour = NA) +
  geom_step(data = sim, linewidth = 0.8, linetype = "dashed") +
  geom_step(data = obs, linewidth = 1) +
  facet_wrap(~ Test) +
  scale_colour_manual(values = col_ref, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x * 100, "%")) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  labs(x = "Proportion of a participant's series classified as AR(1)-misspecified",
       y = "Cumulative percentage of participants") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"))

print(p_pp)

ggsave("figures/misspecification_by_participant_ecdf_ADCF_BDS.png",
       p_pp, width = 11, height = 6, dpi = 500, bg = "white")


# ==============================================================================
# 5 Table: % of participants per range, observed vs. expected
# ==============================================================================

tab_obs <- pp %>%
  mutate(range = cut(prop, BRKS, LABS_R, right = FALSE)) %>%
  count(Test, range) %>%
  group_by(Test) %>%
  mutate(obs_pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

tab_exp <- pp %>%
  group_by(Test) %>%
  reframe({
    sims <- replicate(R_SIM,
                      table(cut(rbinom(n(), k, rate) / k, BRKS, LABS_R,
                                right = FALSE)) / n())
    tibble(range   = factor(LABS_R, LABS_R),
           exp_pct = round(100 * rowMeans(sims), 1))
  })

tab_pp <- tab_obs %>%
  dplyr::select(Test, range, n, obs_pct) %>%
  right_join(tab_exp, by = c("Test", "range")) %>%
  mutate(n = coalesce(n, 0L), obs_pct = coalesce(obs_pct, 0)) %>%
  arrange(Test, range)

print(tab_pp, n = Inf)


# ==============================================================================
# 6 Key numbers for the text
# ==============================================================================

# participants without any misspecified series (observed vs. expected),
# median proportion per participant
pp %>%
  group_by(Test) %>%
  summarise(none_obs = round(100 * mean(x == 0), 1),
            none_exp = round(100 * mean(dbinom(0, k, rate)), 1),
            half_obs = round(100 * mean(prop >= .5), 1),
            mdn_prop = round(100 * median(prop), 1))
