# ==============================================================================
# Supplementary analyses RQ2: main models + time-series length as covariate
#
# Contents:
#   A. Nonlinearity (two-part), as in _analysis.R 1a/1b, + scale(logT)
#        A1. EDF > 1.001, binomial
#        A2. EDF - 1 among nonlinear series, Gamma(log)
#   B. State dependence, as in _analysis.R 3a, + scale(logT)
#        B1. theta > 0, binomial
#
# All models: scale(adcf) + scale(logT) + (1 | dataset/person).
# logT = log(n_obs); n_obs = number of non-missing values of the series (see
# ts_inventory in _analysis.R); overnight NA rows are not counted. This is the
# same covariate as in the regime-shift model (m_cp) in _analysis.R.
# ==============================================================================

library(tidyverse)
library(glmmTMB)


# load latest results
load("results/rq1_results_2026-10-22_125950.Rdata")
load("results/rq2_results_2026-10-13_112850.Rdata")

# ------------------------------------------------------------------------------
# Data -------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Merge rq1_results and rq2_results
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

d2_edf    <- d2 %>% filter(!is.na(edf_lag), !is.na(adcf))
d2_edf_nl <- d2_edf %>% filter(edf_lag > 1.001)
d2_theta  <- d2 %>% filter(!is.na(theta_opt), !is.na(adcf))

# sample sizes (series / persons / datasets); should equal the main models
n_info <- function(d) c(series   = nrow(d),
                        persons  = n_distinct(d$person),
                        datasets = n_distinct(d$dataset))
rbind(edf_bin  = n_info(d2_edf),     # main: 15814 / 2077 / 18
      edf_nl   = n_info(d2_edf_nl),  # main:  2867 / 1430 / 18
      theta    = n_info(d2_theta))   # main: 16352 / 2048 / 18

# Range of T and collinearity of ADCF with logT in each analysis sample
# (ADCF is a dependence statistic and its null distribution shifts with n)
sapply(list(edf_bin = d2_edf, edf_nl = d2_edf_nl, theta = d2_theta),
       function(d) c(T_min    = min(d$n_obs),
                     T_mdn    = median(d$n_obs),
                     T_max    = max(d$n_obs),
                     r_adcf_logT = cor(d$adcf, d$logT)))

# helper: fixed effects with Wald CIs on the link and the exponentiated scale
# (OR for binomial models, mean ratio for Gamma models)
fe_table <- function(m) {
  ci <- confint(m, parm = "beta_", method = "Wald")
  out <- cbind(b = ci[, "Estimate"], ci[, 1:2],
               exp_b = exp(ci[, "Estimate"]), exp(ci[, 1:2]))
  colnames(out) <- c("b", "b_lwr", "b_upr", "exp_b", "exp_lwr", "exp_upr")
  round(out, 3)
}


# ------------------------------------------------------------------------------
# A. Nonlinearity (two-part) ---------------------------------------------------
# ------------------------------------------------------------------------------

# A1) Is the lag-1 relation nonlinear (EDF > 1.001)?
m_edf_bin_logT <- glmmTMB(I(edf_lag > 1.001) ~ scale(adcf) + scale(logT) + (1 | dataset/person),
                          family = binomial, data = d2_edf)
summary(m_edf_bin_logT)
# Family: binomial  ( logit )
# Formula:          I(edf_lag > 1.001) ~ scale(adcf) + scale(logT) + (1 | dataset/person)
# Data: d2_edf
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 12860.3   12898.7   -6425.2   12850.3     15809
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.2829   0.5319
# dataset        (Intercept) 0.2327   0.4824
# Number of obs: 15814, groups:  person:dataset, 2077; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.80677    0.11831  -15.27   <2e-16 ***
#   scale(adcf)  0.87246    0.02615   33.36   <2e-16 ***
#   scale(logT)  0.85626    0.06629   12.92   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fe_table(m_edf_bin_logT)
# b  b_lwr  b_upr exp_b exp_lwr exp_upr
# (Intercept) -1.807 -2.039 -1.575 0.164   0.130   0.207
# scale(adcf)  0.872  0.821  0.924 2.393   2.273   2.519
# scale(logT)  0.856  0.726  0.986 2.354   2.068   2.681

# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_edf_bin_logT, parm = "theta_", method = "Wald")^2
# 2.5 %    97.5 %  Estimate
# Std.Dev.(Intercept)|person:dataset 0.2082678 0.3842261 0.2828814
# Std.Dev.(Intercept)|dataset        0.1103872 0.4905442 0.2327012


# A2) Among nonlinear series: how nonlinear (EDF - 1)?
m_edf_nl_logT <- glmmTMB(I(edf_lag - 1) ~ scale(adcf) + scale(logT) + (1 | dataset/person),
                         family = Gamma(link = "log"), data = d2_edf_nl)
summary(m_edf_nl_logT)
# Family: Gamma  ( log )
# Formula:          I(edf_lag - 1) ~ scale(adcf) + scale(logT) + (1 | dataset/person)
# Data: d2_edf_nl
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 8388.0    8423.8   -4188.0    8376.0      2861
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.078438 0.28007
# dataset        (Intercept) 0.009613 0.09805
# Number of obs: 2867, groups:  person:dataset, 1430; dataset, 18
#
# Dispersion estimate for Gamma family (sigma^2): 0.396
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.53199    0.02855  18.631  < 2e-16 ***
#   scale(adcf)  0.25650    0.01458  17.596  < 2e-16 ***
#   scale(logT)  0.16289    0.02597   6.271 3.58e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fe_table(m_edf_nl_logT)
# b b_lwr b_upr exp_b exp_lwr exp_upr
# (Intercept) 0.532 0.476 0.588 1.702   1.610   1.800
# scale(adcf) 0.256 0.228 0.285 1.292   1.256   1.330
# scale(logT) 0.163 0.112 0.214 1.177   1.118   1.238


# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_edf_nl_logT, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %    Estimate
# Std.Dev.(Intercept)|person:dataset 0.058443827 0.10527155 0.078437698
# Std.Dev.(Intercept)|dataset        0.003844974 0.02403545 0.009613307

# ------------------------------------------------------------------------------
# B. State dependence ----------------------------------------------------------
# ------------------------------------------------------------------------------

# B1) State dependence present (theta > 0)?
m_theta_bin_logT <- glmmTMB(I(theta_opt > 0) ~ scale(adcf) + scale(logT) + (1 | dataset/person),
                            family = binomial, data = d2_theta)
summary(m_theta_bin_logT)
# Family: binomial  ( logit )
# Formula:          I(theta_opt > 0) ~ scale(adcf) + scale(logT) + (1 | dataset/person)
# Data: d2_theta
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 20519.4   20557.9  -10254.7   20509.4     16347
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.031032 0.17616
# dataset        (Intercept) 0.006653 0.08156
# Number of obs: 16352, groups:  person:dataset, 2048; dataset, 18
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.76448    0.02683  28.494   <2e-16 ***
#   scale(adcf)  0.02371    0.01820   1.302    0.193
# scale(logT) -0.01453    0.02976  -0.488    0.625
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fe_table(m_theta_bin_logT)
# b  b_lwr b_upr exp_b exp_lwr exp_upr
# (Intercept)  0.764  0.712 0.817 2.148   2.038   2.264
# scale(adcf)  0.024 -0.012 0.059 1.024   0.988   1.061
# scale(logT) -0.015 -0.073 0.044 0.986   0.930   1.045

# random-effect variances (Wald CI on the log-SD scale, back-transformed)
confint(m_theta_bin_logT, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %    Estimate
# Std.Dev.(Intercept)|person:dataset 0.009244215 0.10417499 0.031032499
# Std.Dev.(Intercept)|dataset        0.001347996 0.03283297 0.006652723

