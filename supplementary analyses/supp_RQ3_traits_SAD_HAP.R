# ==============================================================================
# Supplementary analyses RQ3: Person-level characteristics and AR(1)
# misspecification for sadness (SAD) and happiness (HAP) time series only
#
# Exploratory and post hoc (not preregistered). Motivated by the higher AR(1)
# misspecification rate of negative vs. positive affect and the floor-effect
# account in the Discussion: if misspecification of negative affect series
# partly reflects floor effects, traits associated with higher negative affect
# (depression, neuroticism) should predict LESS misspecification of a negative
# affect item (SAD), and life satisfaction MORE; for a positive affect item
# (HAP), which is rarely at the floor, these associations should be weaker or
# absent. SAD and HAP are the most frequently assessed NA and PA constructs
# (17 datasets each).
#
# Models as in _analysis.R, Section 4 (RQ3): ordered beta regression (logit),
# ADCF statistic of each series predicted by one within-dataset standardized
# trait, random intercepts for datasets and persons nested in datasets.
# Each person contributes a single SAD and a single HAP series (Section 3). In
# the per-construct models the person intercept is then not identified
# separately from the ordered beta dispersion, so only a dataset intercept is
# used. The trait x construct models (Section 5) keep (1 | dataset/person).
#
# Contents:
#   1. Saved RQ1 results and person-level data
#   2. Parameters
#   3. Data
#   4. Separate models per construct and trait
#      4.1 SAD: depression, neuroticism, life satisfaction
#      4.2 HAP: depression, neuroticism, life satisfaction
#   5. Contrast: trait x construct (SAD vs. HAP)
# ==============================================================================

library(tidyverse)
library(glmmTMB)

# ==============================================================================
# 1 Saved RQ1 results and person-level data
# ==============================================================================

# RQ1 results reported in the main text (same file as _analysis.R)
load("results/rq1_results_2026-10-22_125950.Rdata")

# preprocessed data for the person-level characteristics
load("data/processed/data_processed.Rdata")

## Person inventory (identical to _analysis.R, Section 1)
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


# ==============================================================================
# 2 Parameters
# ==============================================================================

CONSTRUCTS <- c("SAD", "HAP")


# ==============================================================================
# 3 Data
# ==============================================================================

# as d3 in _analysis.R, Section 4, plus construct; restricted to SAD and HAP
d3s <- rq1_results %>%
  transmute(ID, dataset, UUID, construct, affect_type, adcf = adcf_stat, n_obs) %>%
  filter(construct %in% CONSTRUCTS, !is.na(adcf)) %>%
  left_join(person_inventory %>% dplyr::select(dataset, UUID, user_ID, DEP_z, NEUR_z, SWL_z),
            by = c("dataset", "UUID")) %>%
  mutate(person    = interaction(dataset, user_ID, drop = TRUE),
         construct = factor(construct, levels = CONSTRUCTS),
         adcf      = pmin(pmax(adcf, 0), 1))   # ordbeta requires exact [0;1] range

# sanity checks
stopifnot(all(d3s$affect_type[d3s$construct == "SAD"] == "NEG"),
          all(d3s$affect_type[d3s$construct == "HAP"] == "POS"))

# series, persons, datasets per construct
d3s %>%
  group_by(construct) %>%
  summarise(series   = n(),
            persons  = n_distinct(person),
            datasets = n_distinct(dataset))
# # Groups:   construct [2]
# construct series persons datasets
# <fct>      <int>   <int>    <int>
# 1 SAD         1968    1968       17
# 2 HAP         1927    1927       17

# series per person and construct (decides the random-effects structure)
d3s %>% count(construct, person) %>% count(construct, n, name = "persons")
# construct n persons
# 1       SAD 1    1968
# 2       HAP 1    1927

# misspecification rate (p <= .05) per construct, for orientation
rq1_results %>%
  filter(construct %in% CONSTRUCTS, !is.na(adcf_p)) %>%
  group_by(construct) %>%
  summarise(adcf_pct = round(100 * mean(adcf_p <= .05), 1), n = n())
# # A tibble: 2 × 3
# construct adcf_pct     n
# <chr>        <dbl> <int>
# 1 HAP           22.5  1927
# 2 SAD           31.5  1968

# ==============================================================================
# 4 Separate models per construct and trait
# ==============================================================================

# adcf ~ trait + (1 | dataset); one series per person and construct (Section 3)

# --------------------------------------
# 4.1 SAD -------------------------------

# depression
d3s_sad_dep <- d3s %>%
  filter(construct == "SAD", !is.na(DEP_z))

m_sad_dep <- glmmTMB(adcf ~ DEP_z + (1 | dataset),
                  family = ordbeta, data = d3s_sad_dep)

summary(m_sad_dep)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ DEP_z + (1 | dataset)
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -3764.0   -3730.7    1888.0   -3776.0      1886
#
# Random effects:
#   Groups  Name        Variance Std.Dev.
# dataset (Intercept) 0.02498  0.158
# Number of obs: 1892, groups:  dataset, 16
#
# Dispersion parameter for ordbeta family (): 9.38
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.82932    0.04516  -40.50   <2e-16 ***
#   DEP_z       -0.04480    0.01859   -2.41   0.0159 *


# neuroticism
d3s_sad_neur <- d3s %>%
  filter(construct == "SAD", !is.na(NEUR_z))

m_sad_neur <- glmmTMB(adcf ~ NEUR_z + (1 | dataset),
                  family = ordbeta, data = d3s_sad_neur)

summary(m_sad_neur)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ NEUR_z + (1 | dataset)
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -3538.1   -3505.3    1775.1   -3550.1      1745
#
# Random effects:
#   Groups  Name        Variance Std.Dev.
# dataset (Intercept) 0.02456  0.1567
# Number of obs: 1751, groups:  dataset, 16
#
# Dispersion parameter for ordbeta family (): 9.53
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.85168    0.04548  -40.72   <2e-16 ***
#   NEUR_z       0.01051    0.01940    0.54    0.588


# satisfaction with life
d3s_sad_swl <- d3s %>%
  filter(construct == "SAD", !is.na(SWL_z))

m_sad_swl <- glmmTMB(adcf ~ SWL_z + (1 | dataset),
                  family = ordbeta, data = d3s_sad_swl)

summary(m_sad_swl)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ SWL_z + (1 | dataset)
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -3068.3   -3036.0    1540.2   -3080.3      1605
#
# Random effects:
#   Groups  Name        Variance Std.Dev.
# dataset (Intercept) 0.02778  0.1667
# Number of obs: 1611, groups:  dataset, 12
#
# Dispersion parameter for ordbeta family (): 8.61
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.78786    0.05373  -33.27   <2e-16 ***
#   SWL_z        0.03996    0.02020    1.98   0.0479 *


# --------------------------------------
# 4.2 HAP -------------------------------

# depression
d3s_hap_dep <- d3s %>%
  filter(construct == "HAP", !is.na(DEP_z))

m_hap_dep <- glmmTMB(adcf ~ DEP_z + (1 | dataset),
                  family = ordbeta, data = d3s_hap_dep)

summary(m_hap_dep)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ DEP_z + (1 | dataset)
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -6307.9   -6274.8    3160.0   -6319.9      1844
#
# Random effects:
#   Groups  Name        Variance Std.Dev.
# dataset (Intercept) 0.03864  0.1966
# Number of obs: 1850, groups:  dataset, 16
#
# Dispersion parameter for ordbeta family (): 34.4
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.333093   0.051299  -45.48   <2e-16 ***
#   DEP_z       -0.004973   0.013091   -0.38    0.704


# neuroticism
d3s_hap_neur <- d3s %>%
  filter(construct == "HAP", !is.na(NEUR_z))

m_hap_neur <- glmmTMB(adcf ~ NEUR_z + (1 | dataset),
                  family = ordbeta, data = d3s_hap_neur)

summary(m_hap_neur)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ NEUR_z + (1 | dataset)
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -5495.9   -5463.7    2754.0   -5507.9      1591
#
# Random effects:
#   Groups  Name        Variance Std.Dev.
# dataset (Intercept) 0.05222  0.2285
# Number of obs: 1597, groups:  dataset, 15
#
# Dispersion parameter for ordbeta family ():   35
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.36499    0.06127  -38.60   <2e-16 ***
#   NEUR_z      -0.01336    0.01405   -0.95    0.341


# satisfaction with life
d3s_hap_swl <- d3s %>%
  filter(construct == "HAP", !is.na(SWL_z))

m_hap_swl <- glmmTMB(adcf ~ SWL_z + (1 | dataset),
                  family = ordbeta, data = d3s_hap_swl)

summary(m_hap_swl)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ SWL_z + (1 | dataset)
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -5488.3   -5456.2    2750.2   -5500.3      1563
#
# Random effects:
#   Groups  Name        Variance Std.Dev.
# dataset (Intercept) 0.04483  0.2117
# Number of obs: 1569, groups:  dataset, 12
#
# Dispersion parameter for ordbeta family ():   37
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.38192    0.06310  -37.75   <2e-16 ***
#   SWL_z        0.02057    0.01378    1.49    0.135


# ==============================================================================
# 5 Contrast: trait x construct (SAD vs. HAP)
# ==============================================================================

# Both constructs in one model; persons contribute a SAD and a HAP series, so
# the main-model random structure (1 | dataset/person) is identified.
# Reference construct = SAD: `trait` = SAD slope, `trait:constructHAP` = HAP
# slope minus SAD slope (test of a difference between SAD and HAP).

# depression
d3s_int_dep <- d3s %>%
  filter(!is.na(DEP_z))

m_int_dep <- glmmTMB(adcf ~ DEP_z * construct + (1 | dataset/person),
                  family = ordbeta, data = d3s_int_dep)

summary(m_int_dep)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ DEP_z * construct + (1 | dataset/person)
# Data: d3s_int_dep
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -10101.3  -10045.3    5059.7  -10119.3      3733
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.17180  0.4145
# dataset        (Intercept) 0.03399  0.1844
# Number of obs: 3742, groups:  person:dataset, 2006; dataset, 17
#
# Dispersion parameter for ordbeta family (): 27.1
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)        -2.03514    0.04799  -42.41  < 2e-16 ***
#   DEP_z              -0.05444    0.01614   -3.37 0.000745 ***
#   constructHAP       -0.28605    0.02033  -14.07  < 2e-16 ***
#   DEP_z:constructHAP  0.04823    0.01985    2.43 0.015139 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# neuroticism
d3s_int_neur <- d3s %>%
  filter(!is.na(NEUR_z))

m_int_neur <- glmmTMB(adcf ~ NEUR_z * construct + (1 | dataset/person),
                  family = ordbeta, data = d3s_int_neur)

summary(m_int_neur)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ NEUR_z * construct + (1 | dataset/person)
# Data: d3s_int_neur
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -9193.8   -9138.7    4605.9   -9211.8      3339
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.18448  0.4295
# dataset        (Intercept) 0.04141  0.2035
# Number of obs: 3348, groups:  person:dataset, 1753; dataset, 16
#
# Dispersion parameter for ordbeta family ():   29
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         -2.07712    0.05413  -38.37   <2e-16 ***
#   NEUR_z              -0.01800    0.01681   -1.07    0.284
# constructHAP        -0.27164    0.02081  -13.05   <2e-16 ***
#   NEUR_z:constructHAP  0.01959    0.02049    0.96    0.339
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# satisfaction with life
d3s_int_swl <- d3s %>%
  filter(!is.na(SWL_z))

m_int_swl <- glmmTMB(adcf ~ SWL_z * construct + (1 | dataset/person),
                  family = ordbeta, data = d3s_int_swl)

summary(m_int_swl)
# Family: ordbeta  ( logit )
# Formula:          adcf ~ SWL_z * construct + (1 | dataset/person)
# Data: d3s_int_swl
#
# AIC       BIC    logLik -2*log(L)  df.resid
# -8520.2   -8465.6    4269.1   -8538.2      3171
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.18541  0.4306
# dataset        (Intercept) 0.03817  0.1954
# Number of obs: 3180, groups:  person:dataset, 1725; dataset, 13
#
# Dispersion parameter for ordbeta family (): 27.1
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)        -2.021898   0.057435  -35.20   <2e-16 ***
#   SWL_z               0.031709   0.017523    1.81   0.0704 .
# constructHAP       -0.326908   0.022342  -14.63   <2e-16 ***
#   SWL_z:constructHAP -0.003791   0.021594   -0.18   0.8606
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# reference: pooled main-text estimates (Table 3; all items)
#   Neuroticism       b = -0.012, SE = 0.010, p = .232
#   Depression        b = -0.024, SE = 0.009, p = .007
#   Life satisfaction b =  0.016, SE = 0.010, p = .097


## confidence intervals:

# random-effect variances, 95% Wald CI (log-SD scale, back-transformed)
confint(m_sad_dep,  parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|dataset 0.01089454 0.05726342 0.02497716

confint(m_sad_neur, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|dataset 0.01060996 0.05682934 0.02455518

confint(m_sad_swl,  parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|dataset 0.0108727 0.07096132 0.02777663

confint(m_hap_dep,  parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|dataset 0.01831034 0.08153254 0.03863791

confint(m_hap_neur, parm = "theta_", method = "Wald")^2
# 2.5 %    97.5 %   Estimate
# Std.Dev.(Intercept)|dataset 0.02415762 0.1128788 0.05221957

confint(m_hap_swl,  parm = "theta_", method = "Wald")^2
# 2.5 %    97.5 %   Estimate
# Std.Dev.(Intercept)|dataset 0.01923407 0.1044772 0.04482768

confint(m_int_dep,  parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.15169183 0.19457310 0.17179974
# Std.Dev.(Intercept)|dataset        0.01632859 0.07074047 0.03398665

confint(m_int_neur, parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.16296645 0.20882372 0.18447563
# Std.Dev.(Intercept)|dataset        0.01945633 0.08811973 0.04140636

confint(m_int_swl,  parm = "theta_", method = "Wald")^2
# 2.5 %     97.5 %   Estimate
# Std.Dev.(Intercept)|person:dataset 0.16303132 0.21086762 0.18541312
# Std.Dev.(Intercept)|dataset        0.01667829 0.08736709 0.03817242
