# ==============================================================================
# Supplementary analyses RQ2: state dependence (optimal theta from S-map) ~ ADCF
#
# Contents:
#   A. Sensitivity models on the full sample (copied from _analysis.R, 3b/3c)
#        A1. theta at the upper bound of the tested range (theta == 8), binomial
#        A2. plateau only (0 < theta < 8), Gamma(log)
#   B. Same models restricted to series with >= 40 non-missing observations
#        B1. theta > 0, binomial (model reported in Table 2)
#        B2. theta == 8, binomial
#        B3. plateau (0 < theta < 8), Gamma(log)
#
# All models: scale(adcf) + (1 | dataset/person).
# n_obs = number of non-missing values of the series (see ts_inventory in
# _analysis.R); overnight NA rows are not counted.
# ==============================================================================

library(tidyverse)
library(glmmTMB)


# load latest results (same files as in _analysis.R)
load("results/rq1_results_2026-10-22_125950.Rdata")
load("results/rq2_results_2026-10-13_112850.Rdata")

# minimum number of non-missing observations for the restricted sample
# check_RQ2_missingness.R showed 20% of S-maps failing for TS with < 50 obs.
min_obs <- 50

# ------------------------------------------------------------------------------
# Data -------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Merge rq1_results and rq2_results (as in _analysis.R)
d2 <- rq1_results %>%
  transmute(ID, adcf = adcf_stat, n_obs) %>%
  left_join(transmute(rq2_results, ID, dataset, user_ID,
                      edf_lag, theta_opt, changepoint_n),
            by = "ID") %>%
  mutate(logT    = log(n_obs),
         dataset = factor(dataset),
         person  = factor(interaction(dataset, user_ID, drop = TRUE)))

# n_obs is taken from rq1_results; check it matches rq2_results
stopifnot(all(rq1_results$n_obs == rq2_results$n_obs[match(rq1_results$ID, rq2_results$ID)],
              na.rm = TRUE))

d2_theta    <- d2 %>% filter(!is.na(theta_opt), !is.na(adcf))
d2_theta_pl <- d2_theta %>% filter(theta_opt > 0, theta_opt < 8)

# restricted samples
d2_theta.40    <- d2_theta    %>% filter(n_obs >= min_obs)
d2_theta_pl.40 <- d2_theta_pl %>% filter(n_obs >= min_obs)

# sample sizes (series / persons / datasets)
n_info <- function(d) c(series   = nrow(d),
                        persons  = n_distinct(d$person),
                        datasets = n_distinct(d$dataset))
rbind(full       = n_info(d2_theta),
      full_pl    = n_info(d2_theta_pl),
      min40      = n_info(d2_theta.40),
      min40_pl   = n_info(d2_theta_pl.40))
# series persons datasets
# full      16352    2048       18
# full_pl    7380    2008       18
# min40     14368    1800       18
# min40_pl   6493    1764       18

# theta distribution in full vs. restricted sample
table(d2_theta$theta_opt)
table(d2_theta.40$theta_opt)
mean(d2_theta.40$theta_opt == 0) # 0.3200863
mean(d2_theta.40$theta_opt == 8) # 0.2280067


# ------------------------------------------------------------------------------
# A. Sensitivity models, full sample -------------------------------------------
# ------------------------------------------------------------------------------

# A1) theta at the upper bound of the tested range (theta = 8)
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

# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_theta_bin_adcf.8, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %    Estimate
# Std.Dev.(Intercept)|person:dataset 0.010575796 0.12979652 0.037049986
# Std.Dev.(Intercept)|dataset        0.000987617 0.02594484 0.005061973


# A2) Plateau only (0 < theta < 8): degree of state dependence
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

# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_theta_gam_adcf, parm = "theta_", method = "Wald")^2
# 2.5 %      97.5 %     Estimate
# Std.Dev.(Intercept)|person:dataset 0.0000000000         Inf 1.887685e-08
# Std.Dev.(Intercept)|dataset        0.0001359145 0.005378504 8.549950e-04

# ------------------------------------------------------------------------------
# B. Restricted sample: series with >= 40 non-missing observations -------------
# ------------------------------------------------------------------------------

# B1) State dependence present (theta > 0); model of Table 2
m_theta_bin_adcf.40 <- glmmTMB(I(theta_opt > 0) ~ scale(adcf) + (1 | dataset/person),
                               family = binomial, data = d2_theta.40)
# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_theta_bin_adcf.40, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.003161178 0.13585653 0.02072358
# Std.Dev.(Intercept)|dataset        0.004344976 0.03914303 0.01304130

summary(m_theta_bin_adcf.40)
# Family: binomial  ( logit )
# Formula:          I(theta_opt > 0) ~ scale(adcf) + (1 | dataset/person)
# Data: d2_theta.40
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 18006.9   18037.1   -8999.4   17998.9     14364
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.02072  0.1440
# dataset        (Intercept) 0.01304  0.1142
# Number of obs: 14368, groups:  person:dataset, 1800; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.77640    0.03371   23.03   <2e-16 ***
#   scale(adcf)  0.01510    0.01865    0.81    0.418
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# B2) theta at the upper bound of the tested range (theta = 8)
m_theta_bin_adcf.8.40 <- glmmTMB(I(theta_opt == 8) ~ scale(adcf) + (1 | dataset/person),
                                 family = binomial, data = d2_theta.40)
# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_theta_bin_adcf.8.40, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %    Estimate
# Std.Dev.(Intercept)|person:dataset 0.004759161 0.15912952 0.027519503
# Std.Dev.(Intercept)|dataset        0.002032760 0.03073628 0.007904396

summary(m_theta_bin_adcf.8.40)
# Family: binomial  ( logit )
# Formula:          I(theta_opt == 8) ~ scale(adcf) + (1 | dataset/person)
# Data: d2_theta.40
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 15425.6   15455.9   -7708.8   15417.6     14364
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.027520 0.16589
# dataset        (Intercept) 0.007904 0.08891
# Number of obs: 14368, groups:  person:dataset, 1800; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.23081    0.03072  -40.06   <2e-16 ***
#   scale(adcf) -0.03540    0.02114   -1.67    0.094 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# B3) Plateau only (0 < theta < 8)
m_theta_gam_adcf.40 <- glmmTMB(theta_opt ~ scale(adcf) + (1 | dataset/person),
                               family = Gamma(link = "log"), data = d2_theta_pl.40)

# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_theta_gam_adcf.40, parm = "theta_", method = "Wald")^2
# 2.5 %      97.5 %     Estimate
# Std.Dev.(Intercept)|person:dataset 0.000000e+00         Inf 8.217855e-09
# Std.Dev.(Intercept)|dataset        6.995368e-05 0.006554517 6.771356e-04

summary(m_theta_gam_adcf.40)
# Family: Gamma  ( log )
# Formula:          theta_opt ~ scale(adcf) + (1 | dataset/person)
# Data: d2_theta_pl.40
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 24940.5   24974.4  -12465.2   24930.5      6488
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance  Std.Dev.
# person:dataset (Intercept) 8.218e-09 9.065e-05
# dataset        (Intercept) 6.771e-04 2.602e-02
# Number of obs: 6493, groups:  person:dataset, 1764; dataset, 18
#
# Dispersion estimate for Gamma family (sigma^2): 0.637
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  9.738e-01  1.186e-02   82.11   <2e-16 ***
#   scale(adcf) -8.277e-05  9.802e-03   -0.01    0.993
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------------------
# Overview of the ADCF effects -------------------------------------------------
# ------------------------------------------------------------------------------

theta_models <- list(
  "theta > 0, full"         = m_theta_bin_adcf.40,   # replaced below if main model is loaded
  "theta == 8, full"        = m_theta_bin_adcf.8,
  "plateau, full"           = m_theta_gam_adcf,
  "theta > 0, n_obs >= 40"  = m_theta_bin_adcf.40,
  "theta == 8, n_obs >= 40" = m_theta_bin_adcf.8.40,
  "plateau, n_obs >= 40"    = m_theta_gam_adcf.40
)
