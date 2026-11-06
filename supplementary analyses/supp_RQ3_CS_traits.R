# ==============================================================================
# Supplementary analyses RQ3: Person-level characteristics and
# complex-system signatures (Supplementary Materials, Section 6.1)
#
# Contents:
#   1. Data preparation and trait availability (series, persons, datasets)
#   2. Nonlinearity
#      2.1 Logistic: EDF > 1.001                         -> ORs
#      2.2 Gamma (log): EDF - 1 among EDF > 1.001        -> ratios of means
#   3. Regime shifts
#      3.1 Poisson (log): changepoint count + scale(log T) -> RRs
#   4. State dependence (optimal S-map theta)
#      4.1 Logistic: theta > 0 (main model)
#      4.2 Sensitivity: logistic, theta = 8 (upper grid bound)
#      4.3 Sensitivity: Gamma (log), 0 < theta < 8
#
# Reported in: Supplement 6.1.1-6.1.3; Table SXX (main models);
#   TableS_traits_theta_sensitivity (4.2, 4.3)
# ==============================================================================


# packages
library(dplyr)
library(glmmTMB)

# load preprocessed data
load("data/processed/data_processed.Rdata")
load("results/excluded_subjects_data.Rdata")

# load input sheet for item info
inputsESM <- readxl::read_xlsx("inputs/input_file_ESM.xlsx")

# load RQ2 results
load("results/rq2_results_2026-10-13_112850.Rdata")

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


d4 <- rq2_results %>%
  left_join(person_inventory %>% dplyr::select(dataset, UUID, DEP_z, NEUR_z, SWL_z),
            by = c("dataset", "UUID")) %>%
  mutate(
    person       = interaction(dataset, user_ID, drop = TRUE),
    dataset      = factor(dataset),
    logT         = log(n_obs),
    edf_nonlin   = edf_lag > 1.001,
    theta_nonlin = theta_opt > 0,
    theta_max    = theta_opt == 8
  )

# availability of each trait: persons and datasets
d4 %>%
  summarise(across(c(DEP_z, NEUR_z, SWL_z),
                   list(series   = ~ sum(!is.na(.x)),
                        persons  = ~ n_distinct(person[!is.na(.x)]),
                        datasets = ~ n_distinct(dataset[!is.na(.x)]))))

# --------------------------------------
# Binary outcome: nonlinear lag-1 relation (EDF > 1.001)

d4_edf <- d4 %>% filter(!is.na(edf_lag))

m_edf_bin_dep  <- glmmTMB(I(edf_lag > 1.001) ~ DEP_z + (1 | dataset/person),
                          family = binomial,
                          data = d4_edf %>% filter(!is.na(DEP_z)))

m_edf_bin_neur <- glmmTMB(I(edf_lag > 1.001) ~ NEUR_z + (1 | dataset/person),
                          family = binomial,
                          data = d4_edf %>% filter(!is.na(NEUR_z)))

m_edf_bin_swl  <- glmmTMB(I(edf_lag > 1.001) ~ SWL_z + (1 | dataset/person),
                          family = binomial,
                          data = d4_edf %>% filter(!is.na(SWL_z)))

summary(m_edf_bin_dep)
# Family: binomial  ( logit )
# Formula:          I(edf_lag > 1.001) ~ DEP_z + (1 | dataset/person)
# Data: d4_edf %>% filter(!is.na(DEP_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 13620.1   13650.6   -6806.0   13612.1     15220
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.3643   0.6036
# dataset        (Intercept) 0.2288   0.4784
# Number of obs: 15224, groups:  person:dataset, 2022; dataset, 17
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.66053    0.12019 -13.816  < 2e-16 ***
#   DEP_z       -0.10624    0.02709  -3.921 8.81e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_edf_bin_neur)
# Family: binomial  ( logit )
# Formula:          I(edf_lag > 1.001) ~ NEUR_z + (1 | dataset/person)
# Data: d4_edf %>% filter(!is.na(NEUR_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 12267.5   12297.6   -6129.8   12259.5     13738
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.3753   0.6126
# dataset        (Intercept) 0.2320   0.4817
# Number of obs: 13742, groups:  person:dataset, 1770; dataset, 16
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.65590    0.12506 -13.241   <2e-16 ***
#   NEUR_z      -0.01605    0.02855  -0.562    0.574
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_edf_bin_swl)
# Family: binomial  ( logit )
# Formula:          I(edf_lag > 1.001) ~ SWL_z + (1 | dataset/person)
# Data: d4_edf %>% filter(!is.na(SWL_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 11940.4   11970.4   -5966.2   11932.4     13316
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.3648   0.6039
# dataset        (Intercept) 0.2809   0.5300
# Number of obs: 13320, groups:  person:dataset, 1740; dataset, 13
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.64266    0.15061 -10.907  < 2e-16 ***
#   SWL_z        0.10379    0.02894   3.586 0.000335 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# --------------------------------------
# EDF ~ trait among nonlinear time series (EDF > 1.001)

d4_edf_nl <- d4 %>% filter(!is.na(edf_lag), edf_lag > 1.001)

m_edf_gam_dep  <- glmmTMB(I(edf_lag - 1) ~ DEP_z + (1 | dataset/person),
                          family = Gamma(link = "log"),
                          data = d4_edf_nl %>% filter(!is.na(DEP_z)))

m_edf_gam_neur <- glmmTMB(I(edf_lag - 1) ~ NEUR_z + (1 | dataset/person),
                          family = Gamma(link = "log"),
                          data = d4_edf_nl %>% filter(!is.na(NEUR_z)))

m_edf_gam_swl  <- glmmTMB(I(edf_lag - 1) ~ SWL_z + (1 | dataset/person),
                          family = Gamma(link = "log"),
                          data = d4_edf_nl %>% filter(!is.na(SWL_z)))

summary(m_edf_gam_dep)
# Family: Gamma  ( log )
# Formula:          I(edf_lag - 1) ~ DEP_z + (1 | dataset/person)
# Data: d4_edf_nl %>% filter(!is.na(DEP_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 8357.6    8387.2   -4173.8    8347.6      2736
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.11096  0.3331
# dataset        (Intercept) 0.01727  0.1314
# Number of obs: 2741, groups:  person:dataset, 1377; dataset, 17
#
# Dispersion estimate for Gamma family (sigma^2): 0.424
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.51653    0.03701  13.955   <2e-16 ***
#   DEP_z        0.00329    0.01652   0.199    0.842
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_edf_gam_neur)
# family: Gamma  ( log )
# Formula:          I(edf_lag - 1) ~ NEUR_z + (1 | dataset/person)
# Data: d4_edf_nl %>% filter(!is.na(NEUR_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 7454.6    7483.6   -3722.3    7444.6      2467
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.09991  0.3161
# dataset        (Intercept) 0.01690  0.1300
# Number of obs: 2472, groups:  person:dataset, 1217; dataset, 16
#
# Dispersion estimate for Gamma family (sigma^2): 0.415
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.51062    0.03794  13.460   <2e-16 ***
#   NEUR_z       0.02122    0.01730   1.226     0.22
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_edf_gam_swl)
# Family: Gamma  ( log )
# Formula:          I(edf_lag - 1) ~ SWL_z + (1 | dataset/person)
# Data: d4_edf_nl %>% filter(!is.na(SWL_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 7445.8    7474.8   -3717.9    7435.8      2423
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.10192  0.3193
# dataset        (Intercept) 0.01454  0.1206
# Number of obs: 2428, groups:  person:dataset, 1189; dataset, 13
#
# Dispersion estimate for Gamma family (sigma^2): 0.42
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.53860    0.03862  13.945   <2e-16 ***
#   SWL_z        0.01857    0.01773   1.047    0.295
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -- get odd ratios
models <- list(DEP = m_edf_bin_dep, NEUR = m_edf_bin_neur, SWL = m_edf_bin_swl)
terms  <- c(DEP = "DEP_z", NEUR = "NEUR_z", SWL = "SWL_z")

or_table <- do.call(rbind, lapply(names(models), function(k) {
  m  <- models[[k]]
  ci <- confint(m, parm = terms[k])            # Wald 95% CI on the log-odds scale
  data.frame(
    trait = k,
    OR    = exp(fixef(m)$cond[[terms[k]]]),
    lower = exp(ci[1, 1]),
    upper = exp(ci[1, 2]),
    p     = summary(m)$coefficients$cond[terms[k], "Pr(>|z|)"],
    n_series = nobs(m)
  )
}))

or_table$p_apa <- ifelse(or_table$p < .001, "< .001", sub("^0", "", sprintf("%.3f", or_table$p)))
print(or_table, digits = 3)

# trait    OR lower upper        p n_series  p_apa
# 1   DEP 0.899 0.853 0.948 8.81e-05    15224 < .001
# 2  NEUR 0.984 0.931 1.041 5.74e-01    13742   .574
# 3   SWL 1.109 1.048 1.174 3.35e-04    13320 < .001

## ----------------------------------------------------
# Change point counts

d4_cp <- d4 %>% filter(!is.na(changepoint_n))

m_cp_pois_dep  <- glmmTMB(changepoint_n ~ DEP_z + scale(logT) + (1 | dataset/person),
                          family = poisson,
                          data = d4_cp %>% filter(!is.na(DEP_z)))

m_cp_pois_neur <- glmmTMB(changepoint_n ~ NEUR_z + scale(logT) + (1 | dataset/person),
                          family = poisson,
                          data = d4_cp %>% filter(!is.na(NEUR_z)))

m_cp_pois_swl  <- glmmTMB(changepoint_n ~ SWL_z + scale(logT) + (1 | dataset/person),
                          family = poisson,
                          data = d4_cp %>% filter(!is.na(SWL_z)))


summary(m_cp_pois_dep)
# Family: poisson  ( log )
# Formula:          changepoint_n ~ DEP_z + scale(logT) + (1 | dataset/person)
# Data: d4_cp %>% filter(!is.na(DEP_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 41328.8   41367.5  -20659.4   41318.8     16870
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.27562  0.5250
# dataset        (Intercept) 0.04902  0.2214
# Number of obs: 16875, groups:  person:dataset, 2031; dataset, 17
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.26640    0.05629  -4.733 2.21e-06 ***
#   DEP_z        0.03945    0.01463   2.697    0.007 **
#   scale(logT)  0.42684    0.03065  13.927  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_cp_pois_neur)
# Family: poisson  ( log )
# Formula:          changepoint_n ~ NEUR_z + scale(logT) + (1 | dataset/person)
# Data: d4_cp %>% filter(!is.na(NEUR_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 37580.9   37619.0  -18785.4   37570.9     15189
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.28104  0.5301
# dataset        (Intercept) 0.05099  0.2258
# Number of obs: 15194, groups:  person:dataset, 1778; dataset, 16
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.25588    0.05950   -4.30  1.7e-05 ***
#   NEUR_z       0.03371    0.01575    2.14   0.0324 *
#   scale(logT)  0.45588    0.03336   13.67  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_cp_pois_swl)
# Family: poisson  ( log )
# Formula:          changepoint_n ~ SWL_z + scale(logT) + (1 | dataset/person)
# Data: d4_cp %>% filter(!is.na(SWL_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 36354.4   36392.3  -18172.2   36344.4     14617
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.27292  0.5224
# dataset        (Intercept) 0.04278  0.2068
# Number of obs: 14622, groups:  person:dataset, 1748; dataset, 13
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.24254    0.06017  -4.031 5.56e-05 ***
#   SWL_z       -0.01056    0.01576  -0.670    0.503
# scale(logT)  0.44379    0.03651  12.155  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# -- Compute rate ratios (RRs)

models <- list(DEP = m_cp_pois_dep, NEUR = m_cp_pois_neur, SWL = m_cp_pois_swl)
terms  <- c(DEP = "DEP_z", NEUR = "NEUR_z", SWL = "SWL_z")

rr_table <- do.call(rbind, lapply(names(models), function(k) {
  m  <- models[[k]]
  ci <- confint(m, parm = c(terms[k], "scale(logT)"))   # Wald 95% CI on the log scale
  data.frame(
    trait = k,
    term  = c("trait", "logT"),
    RR    = exp(c(fixef(m)$cond[[terms[k]]], fixef(m)$cond[["scale(logT)"]])),
    lower = exp(ci[, 1]),
    upper = exp(ci[, 2]),
    p     = summary(m)$coefficients$cond[c(terms[k], "scale(logT)"), "Pr(>|z|)"]
  )
}))

rr_table$p_apa <- ifelse(rr_table$p < .001, "< .001", sub("^0", "", sprintf("%.3f", rr_table$p)))
print(rr_table, digits = 3, row.names = FALSE)



## ----------------------------------------------------
# State dependence

# Logistic models for theta_opt indicating state dependence (> 0)

d4_theta <- d4 %>% filter(!is.na(theta_opt))

m_theta_bin_dep  <- glmmTMB(I(theta_opt > 0) ~ DEP_z + (1 | dataset/person),
                            family = binomial,
                            data = d4_theta %>% filter(!is.na(DEP_z)))

m_theta_bin_neur <- glmmTMB(I(theta_opt > 0) ~ NEUR_z + (1 | dataset/person),
                            family = binomial,
                            data = d4_theta %>% filter(!is.na(NEUR_z)))

m_theta_bin_swl  <- glmmTMB(I(theta_opt > 0) ~ SWL_z + (1 | dataset/person),
                            family = binomial,
                            data = d4_theta %>% filter(!is.na(SWL_z)))

summary(m_theta_bin_dep)
# Family: binomial  ( logit )
# Formula:          I(theta_opt > 0) ~ DEP_z + (1 | dataset/person)
# Data: d4_theta %>% filter(!is.na(DEP_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 19659.4   19690.0   -9825.7   19651.4     15689
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.03379  0.18383
# dataset        (Intercept) 0.00929  0.09638
# Number of obs: 15693, groups:  person:dataset, 1985; dataset, 17
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.773192   0.030328  25.494   <2e-16 ***
#   DEP_z       -0.003517   0.017708  -0.199    0.843
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_theta_bin_neur)
# Family: binomial  ( logit )
# Formula:          I(theta_opt > 0) ~ NEUR_z + (1 | dataset/person)
# Data: d4_theta %>% filter(!is.na(NEUR_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 17700.7   17730.9   -8846.4   17692.7     14119
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.01922  0.1386
# dataset        (Intercept) 0.01158  0.1076
# Number of obs: 14123, groups:  person:dataset, 1735; dataset, 16
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.76756    0.03373  22.759   <2e-16 ***
#   NEUR_z      -0.03157    0.01848  -1.708   0.0876 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_theta_bin_swl)
# Family: binomial  ( logit )
# Formula:          I(theta_opt > 0) ~ SWL_z + (1 | dataset/person)
# Data: d4_theta %>% filter(!is.na(SWL_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 17187.5   17217.6   -8589.8   17179.5     13647
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.042692 0.20662
# dataset        (Intercept) 0.007089 0.08419
# Number of obs: 13651, groups:  person:dataset, 1710; dataset, 13
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) 0.7475580  0.0309263  24.172   <2e-16 ***
#   SWL_z       0.0004907  0.0191217   0.026     0.98
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -- Sensitivity analyses for theta_opt = 8 (max of tested range)

m_theta_bin_dep.8  <- glmmTMB(I(theta_opt == 8) ~ DEP_z + (1 | dataset/person),
                              family = binomial,
                              data = d4_theta %>% filter(!is.na(DEP_z)))

m_theta_bin_neur.8 <- glmmTMB(I(theta_opt == 8) ~ NEUR_z + (1 | dataset/person),
                              family = binomial,
                              data = d4_theta %>% filter(!is.na(NEUR_z)))

m_theta_bin_swl.8  <- glmmTMB(I(theta_opt == 8) ~ SWL_z + (1 | dataset/person),
                              family = binomial,
                              data = d4_theta %>% filter(!is.na(SWL_z)))

summary(m_theta_bin_dep.8)
# Family: binomial  ( logit )
# Formula:          I(theta_opt == 8) ~ DEP_z + (1 | dataset/person)
# Data: d4_theta %>% filter(!is.na(DEP_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 16871.9   16902.5   -8431.9   16863.9     15689
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.027187 0.16488
# dataset        (Intercept) 0.003044 0.05517
# Number of obs: 15693, groups:  person:dataset, 1985; dataset, 17
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.22937    0.02518  -48.83   <2e-16 ***
#   DEP_z        0.02733    0.01935    1.41    0.158
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_theta_bin_neur.8)
# Family: binomial  ( logit )
# Formula:          I(theta_opt == 8) ~ NEUR_z + (1 | dataset/person)
# Data: d4_theta %>% filter(!is.na(NEUR_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 15147.8   15178.0   -7569.9   15139.8     14119
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.035981 0.18969
# dataset        (Intercept) 0.005105 0.07145
# Number of obs: 14123, groups:  person:dataset, 1735; dataset, 16
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.24176    0.02932  -42.35   <2e-16 ***
#   NEUR_z      -0.02477    0.02081   -1.19    0.234
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_theta_bin_swl.8)
# Family: binomial  ( logit )
# Formula:          I(theta_opt == 8) ~ SWL_z + (1 | dataset/person)
# Data: d4_theta %>% filter(!is.na(SWL_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 14641.9   14671.9   -7316.9   14633.9     13647
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 0.052232 0.22854
# dataset        (Intercept) 0.002972 0.05452
# Number of obs: 13651, groups:  person:dataset, 1710; dataset, 13
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.24249    0.02774  -44.79   <2e-16 ***
#   SWL_z       -0.02420    0.02128   -1.14    0.255
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -- Model only plateau (0 < theta < 8)
d4_theta_pl <- d4 %>% filter(!is.na(theta_opt), theta_opt > 0, theta_opt < 8)

m_theta_gam_dep  <- glmmTMB(theta_opt ~ DEP_z + (1 | dataset/person),
                            family = Gamma(link = "log"),
                            data = d4_theta_pl %>% filter(!is.na(DEP_z)))

m_theta_gam_neur <- glmmTMB(theta_opt ~ NEUR_z + (1 | dataset/person),
                            family = Gamma(link = "log"),
                            data = d4_theta_pl %>% filter(!is.na(NEUR_z)))

m_theta_gam_swl  <- glmmTMB(theta_opt ~ SWL_z + (1 | dataset/person),
                            family = Gamma(link = "log"),
                            data = d4_theta_pl %>% filter(!is.na(SWL_z)))

summary(m_theta_gam_dep)
# Family: Gamma  ( log )
# Formula:          theta_opt ~ DEP_z + (1 | dataset/person)
# Data: d4_theta_pl %>% filter(!is.na(DEP_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 27118.7   27153.1  -13554.4   27108.7      7086
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance Std.Dev.
# person:dataset (Intercept) 6.43e-09 8.019e-05
# dataset        (Intercept) 1.07e-03 3.271e-02
# Number of obs: 7091, groups:  person:dataset, 1946; dataset, 17
#
# Dispersion estimate for Gamma family (sigma^2): 0.638
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) 0.963711   0.012656   76.14   <2e-16 ***
#   DEP_z       0.004850   0.009566    0.51    0.612
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


summary(m_theta_gam_neur)
# Family: Gamma  ( log )
# Formula:          theta_opt ~ NEUR_z + (1 | dataset/person)
# Data: d4_theta_pl %>% filter(!is.na(NEUR_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 24431.6   24465.5  -12210.8   24421.6      6383
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance  Std.Dev.
# person:dataset (Intercept) 2.042e-08 0.0001429
# dataset        (Intercept) 3.092e-04 0.0175831
# Number of obs: 6388, groups:  person:dataset, 1701; dataset, 16
#
# Dispersion estimate for Gamma family (sigma^2): 0.637
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) 0.964980   0.011122   86.76   <2e-16 ***
#   NEUR_z      0.004091   0.009965    0.41    0.681
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_theta_gam_swl)
# Family: Gamma  ( log )
# Formula:          theta_opt ~ SWL_z + (1 | dataset/person)
# Data: d4_theta_pl %>% filter(!is.na(SWL_z))
#
# AIC       BIC    logLik -2*log(L)  df.resid
# 23415.1   23448.7  -11702.5   23405.1      6120
#
# Random effects:
#
#   Conditional model:
#   Groups         Name        Variance  Std.Dev.
# person:dataset (Intercept) 1.126e-08 0.0001061
# dataset        (Intercept) 1.448e-03 0.0380590
# Number of obs: 6125, groups:  person:dataset, 1675; dataset, 13
#
# Dispersion estimate for Gamma family (sigma^2): 0.637
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.962821   0.014951   64.40   <2e-16 ***
#   SWL_z       -0.005394   0.010278   -0.52      0.6
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
