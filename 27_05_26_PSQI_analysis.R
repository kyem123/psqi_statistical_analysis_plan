### PSQI DATA ANALYSES

# Authors: K.M., R.V.R.
# Date: 27/05/2026
# Purpose: Analysis of PSQI global scores (continuous), PSQI_Continuous (0–21 treated as continuous), and component scores (ordinal, 0–3) collected pre and post each session

# --------------------------------------------------------------------------------------------
# Load packages

# install.packages(c(
#   "tidyverse", "dplyr", "lubridate", "summarytools", "ggplot2",
#   "naniar", "MuMIn", "car", "lme4", "lmerTest", "emmeans",
#   "performance", "DHARMa", "broom.mixed", "tibble", "readr",
#   "forcats", "see", "nortest", "ordinal", "brant"
# ))

library(tidyverse)
library(dplyr)
library(lubridate)
library(summarytools)
library(ggplot2)
library(naniar)
library(MuMIn)
library(car)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)
library(broom.mixed)
library(tibble)
library(readr)
library(forcats)
library(see)
library(nortest)
library(ordinal)   


# --------------------------------------------------------------------------------------------
# IMPORT AND CODE DATA
# --------------------------------------------------------------------------------------------

# Dependency: requires "PSQI_scored.csv" produced by PSQI-scoring script

psqi <- read.csv("Z:/Analysis/Kye/PSQI Analysis/PSQI_scored.csv")
setwd("Z:/Analysis/Kye/PSQI Analysis/Outputs")

# Rename to match actigraphy convention
psqi <- psqi %>%
  rename(id = subject_number, treatment = ip)

# Categorical variables
psqi$id <- as.factor(psqi$id)
psqi$arm <- as.factor(psqi$arm)
psqi$treatment <- as.factor(psqi$treatment)
psqi$preorpost <- as.factor(psqi$preorpost)

# Derive treatment order (which treatment appeared first in arm 1)
psqi <- psqi %>%
  group_by(id) %>%
  mutate(order = if_else(first(treatment[arm == 1]) == "A", "A_first", "B_first")) %>%
  ungroup()
psqi$order <- as.factor(psqi$order)

# Set reference levels
psqi <- psqi %>%
  mutate(
    preorpost = fct_relevel(preorpost, "0", "1"), # preorpost coded as 0 = pre-treatment, 1 = post-treatment
    treatment = fct_relevel(treatment, "A")
  )

# Continuous global score
psqi$PSQI_Global_Score <- as.numeric(psqi$PSQI_Global_Score)

# Ordinal component scores (ordered factors, levels 0–3)
component_vars <- paste0("PSQI_Component", c(1:5, 7)) # Component 6 excluded from ordinal analysis - explained later in analysis section
for (v in component_vars) {
  psqi[[v]] <- factor(psqi[[v]], levels = 0:3, ordered = TRUE)
}

# Check structure
str(psqi)
view(dfSummary(psqi))


# --------------------------------------------------------------------------------------------
# DATA CHECKS
# --------------------------------------------------------------------------------------------

cat("\n--- Categorical variable checks ---\n")
cat("\nID (n participants):\n"); print(table(psqi$id))
cat("\nTreatment:\n"); print(table(psqi$treatment))
cat("\nArm:\n"); print(table(psqi$arm))
cat("\nPre or post:\n"); print(table(psqi$preorpost))
cat("\nOrder:\n"); print(table(psqi$order))

# Missing data
cat("\n--- Missing data ---\n")
vis_miss(psqi)
print(colSums(is.na(psqi)))
# 12.3% missing data

# Duplicate record check
cat("\n--- Duplicate record check ---\n")
dupe_check <- psqi %>%
  mutate(condition = paste0(str_to_title(as.character(preorpost)), treatment)) %>%
  group_by(id, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if (nrow(dupe_check) == 0) {
  cat("No duplicate records found.\n")
} else {
  cat("WARNING: Duplicate records detected:\n")
  print(dupe_check)
}
# No duplicate records found

# --------------------------------------------------------------------------------------------
# OUTLIER DETECTION — GLOBAL SCORE ONLY (continuous)
# --------------------------------------------------------------------------------------------

flag_mad <- function(x, threshold = 3) {
  med <- median(x, na.rm = TRUE)
  m   <- mad(x, na.rm = TRUE)
  abs(x - med) > threshold * m
}

psqi <- psqi %>%
  mutate(outlier_MAD_global = flag_mad(PSQI_Global_Score))

psqi %>%
  filter(outlier_MAD_global) %>%
  arrange(id, arm, preorpost) %>%
  select(id, arm, treatment, preorpost, PSQI_Global_Score)
# No extreme outliers detected using MAD method

# --------------------------------------------------------------------------------------------
# DESCRIPTIVE STATISTICS
# --------------------------------------------------------------------------------------------

# Global score
desc_global <- psqi %>%
  group_by(treatment, preorpost) %>%
  summarise(
    n = sum(!is.na(PSQI_Global_Score)),
    mean = mean(PSQI_Global_Score, na.rm = TRUE),
    sd = sd(PSQI_Global_Score, na.rm = TRUE),
    median = median(PSQI_Global_Score, na.rm = TRUE),
    q1 = quantile(PSQI_Global_Score, 0.25, na.rm = TRUE),
    q3 = quantile(PSQI_Global_Score, 0.75, na.rm = TRUE),
    min = min(PSQI_Global_Score, na.rm = TRUE),
    max = max(PSQI_Global_Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(outcome = "PSQI_Global_Score")

print(desc_global)
readr::write_csv(desc_global, "descriptive_statistics_PSQI_global.csv")


# Component scores (ordinal — frequency tables by treatment × preorpost)
for (v in component_vars) {
  cat("\n--- Frequency table:", v, "---\n")
  tbl <- psqi %>%
    group_by(treatment, preorpost, .data[[v]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(outcome = v)
  print(tbl)
  readr::write_csv(tbl, paste0("freq_table_", v, ".csv"))
}

# Normality tests on raw global score
sw  <- shapiro.test(psqi$PSQI_Global_Score[!is.na(psqi$PSQI_Global_Score)])
ad  <- ad.test(psqi$PSQI_Global_Score[!is.na(psqi$PSQI_Global_Score)])
norm_tbl <- data.frame(
  outcome = "PSQI_Global_Score",
  SW_W = sw$statistic, SW_p = sw$p.value,
  AD_A = ad$statistic, AD_p = ad$p.value
)
print(norm_tbl)
readr::write_csv(norm_tbl, "normality_test_PSQI_global.csv")
# Anderson-Darling and Shapiro-Wilk test show no evidence of deviations from normality

# Histogram for normality checking — PSQI Global Score
x <- psqi$PSQI_Global_Score[!is.na(psqi$PSQI_Global_Score)]

hist(x,
     breaks = 15,
     col    = "grey85",
     main   = "Distribution of PSQI Global Score",
     xlab   = "PSQI Global Score",
     freq   = FALSE)
lines(density(x), col = "steelblue", lwd = 2)
qqnorm(x, main = "QQ Plot — PSQI Global Score")
qqline(x, col = "red", lty = 2)
# Histogram and QQ-plot demonstrate that data are approximately normally distributed

# --------------------------------------------------------------------------------------------
# CARRYOVER TESTING
# --------------------------------------------------------------------------------------------

period2_pre <- psqi %>% filter(arm == "2", preorpost == "0")

carryover_results <- list()

# Global score (continuous - t-test)
tt <- t.test(period2_pre$PSQI_Global_Score ~ period2_pre$order)
cat("\nCarryover — PSQI_Global_Score\n")
cat("  t =", round(tt$statistic, 3), "| df =", round(tt$parameter, 1),
    "| p =", round(tt$p.value, 4), "\n")
if (tt$p.value > 0.05) cat("  → No evidence of carryover.\n") else
  cat("  → WARNING: Significant difference. Consider carryover in interpretation.\n")

carryover_results[["PSQI_Global_Score"]] <- data.frame(
  outcome = "PSQI_Global_Score", method = "t-test",
  statistic = round(tt$statistic, 3), df = round(tt$parameter, 1),
  p.value = round(tt$p.value, 4)
)
# No evidence of carryover (p = 0.8006)

# Component scores (ordinal — Wilcoxon rank-sum carryover test)
for (v in component_vars) {
  x <- as.numeric(as.character(period2_pre[[v]]))
  grp <- period2_pre$order
  wt <- tryCatch(wilcox.test(x ~ grp), error = function(e) NULL)
  if (!is.null(wt)) {
    cat("\nCarryover —", v, "\n")
    cat("  W =", round(wt$statistic, 3), "| p =", round(wt$p.value, 4), "\n")
    if (wt$p.value > 0.05) cat("  → No evidence of carryover.\n") else
      cat("  → WARNING: Significant difference. Consider carryover in interpretation.\n")
    carryover_results[[v]] <- data.frame(
      outcome = v, method = "Wilcoxon",
      statistic = round(wt$statistic, 3), df = NA,
      p.value = round(wt$p.value, 4)
    )
  }
}
# No evidence of carryover in any components

carryover_tbl <- bind_rows(carryover_results)
print(carryover_tbl)
readr::write_csv(carryover_tbl, "carryover_tests_PSQI_all.csv")


# ============================================================================================
# SECTION A: CONTINUOUS OUTCOMES — PSQI_Global_Score
# Pipeline: random effects selection → fixed effects AICc → assumption checks with transform fallback (raw → sqrt → log) → DID contrast → Cook's D → LOIO → DHARMa → visualisations
# ============================================================================================

# --------------------------------------------------------------------------------------------
# RANDOM EFFECTS STRUCTURE SELECTION
# --------------------------------------------------------------------------------------------

re_m1 <- lmer(PSQI_Global_Score ~ + (1 + treatment | id), data = psqi)
re_m2 <- lmer(PSQI_Global_Score ~ + (1 | id), data = psqi)

cat("\n--- Random effects structure selection (AICc) ---\n")
re_aic <- MuMIn::AICc(re_m1, re_m2)
print(re_aic)
best_re <- rownames(re_aic)[which.min(re_aic$AICc)]
cat("\nBest random effects structure:", best_re, "\n")
cat("→ (1|id) applied consistently across continuous outcomes unless otherwise noted.\n")
rm(re_m1, re_m2, re_aic, best_re)
# re_m2 appears to be the best random effects structure with the lowest AICc

# --------------------------------------------------------------------------------------------
# HELPER FUNCTIONS (continuous outcomes)
# --------------------------------------------------------------------------------------------

# Fixed-effects model selection via AICc (lmer)
select_fixed_effects <- function(outcome_var, data) {
  f_list <- list(
    a = as.formula(paste(outcome_var, "~ treatment*preorpost + (1|id)")),
    b = as.formula(paste(outcome_var, "~ treatment*preorpost + order + (1|id)")),
    c = as.formula(paste(outcome_var, "~ treatment*preorpost + arm + (1|id)")),
    d = as.formula(paste(outcome_var, "~ treatment*preorpost*order + (1|id)"))
  )
  fits <- lapply(f_list, function(f) {
    tryCatch(lmer(f, data = data, REML = FALSE, na.action = na.exclude),
             error = function(e) NULL)
  })
  valid  <- Filter(Negate(is.null), fits)
  aic_df <- MuMIn::AICc(valid[[1]], valid[[2]], valid[[3]], valid[[4]])
  best <- which.min(aic_df$AICc)
  list(best_model = valid[[best]], aicc_table = aic_df, selected = names(valid)[best])
}

# Normality check (Anderson-Darling primary, Shapiro-Wilk supplementary)
check_normality_ad <- function(model) {
  r <- resid(model); r <- r[!is.na(r)]
  ad <- ad.test(r); sw <- shapiro.test(r)
  cat("Anderson-Darling p =", round(ad$p.value, 4),
      "| Shapiro-Wilk p =", round(sw$p.value, 4), "\n")
  ad$p.value > 0.05
}

# Assumption checks with automatic transform fallback
# override: NULL (auto) | "raw" | "sqrt" | "log"
check_and_transform <- function(outcome_var, data, override = NULL) {

  plot_assumptions <- function(mod, label) {
    r <- resid(mod); f <- fitted(mod)
    dev.new(width = 10, height = 4); par(mfrow = c(1, 3))
    hist(r, breaks = 20, col = "grey85",
         main = paste0("Residuals — ", label, "\n(", outcome_var, ")"), xlab = "Residual")
    qqnorm(r, main = paste0("QQ — ", label, "\n(", outcome_var, ")")); qqline(r, col = "red", lty = 2)
    plot(f, r, main = paste0("Fitted vs Residuals — ", label, "\n(", outcome_var, ")"),
         xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2, col = "red")
    par(mfrow = c(1, 1))
  }

  transform_label <- "raw"
  model <- NULL

  # Override: return immediately if specified
  if (!is.null(override) && override == "raw") {
    sel <- select_fixed_effects(outcome_var, data)
    mod <- sel$best_model; plot_assumptions(mod, "raw (override)")
    cat("\n→ OVERRIDE: retaining raw model.\n")
    dev.new(width=12,height=10); print(check_model(mod))
    return(list(model = mod, transform = "raw", data = data))
  }
  if (!is.null(override) && override == "sqrt") {
    data <- data %>% mutate(!!paste0(outcome_var, "_sqrt") := sqrt(.data[[outcome_var]]))
    v2   <- paste0(outcome_var, "_sqrt")
    sel  <- select_fixed_effects(v2, data); mod <- sel$best_model
    plot_assumptions(mod, "sqrt (override)")
    cat("\n→ OVERRIDE: retaining sqrt model.\n")
    return(list(model = mod, transform = "sqrt", data = data))
  }
  if (!is.null(override) && override == "log") {
    data <- data %>% mutate(!!paste0(outcome_var, "_log") := log(.data[[outcome_var]]))
    v2   <- paste0(outcome_var, "_log")
    sel  <- select_fixed_effects(v2, data); mod <- sel$best_model
    plot_assumptions(mod, "log (override)")
    cat("\n→ OVERRIDE: retaining log model.\n")
    return(list(model = mod, transform = "log", data = data))
  }

  # --- Automatic selection ---
  # Raw
  cat("\n--- Testing raw model:", outcome_var, "---\n")
  sel1 <- select_fixed_effects(outcome_var, data); model1 <- sel1$best_model
  plot_assumptions(model1, "raw")
  normal1 <- check_normality_ad(model1)
  p_bp1   <- as.numeric(check_heteroscedasticity(model1))
  homosc1 <- p_bp1 > 0.05
  cat("Breusch-Pagan p =", round(p_bp1, 4), "| Homoscedastic:", homosc1, "\n")
  
  dev.new(width = 12, height = 10); print(check_model(model1))
  
  if (normal1 && homosc1) {
    model <- model1; transform_label <- "raw"
  } else {
    # SQRT
    cat("\n--- Testing sqrt transformation ---\n")
    data  <- data %>% mutate(!!paste0(outcome_var, "_sqrt") := sqrt(.data[[outcome_var]]))
    v2    <- paste0(outcome_var, "_sqrt")
    sel2  <- select_fixed_effects(v2, data); model2 <- sel2$best_model
    plot_assumptions(model2, "sqrt")
    normal2 <- check_normality_ad(model2)
    p_bp2   <- as.numeric(check_heteroscedasticity(model2))
    homosc2 <- p_bp2 > 0.05
    cat("Breusch-Pagan p =", round(p_bp2, 4), "| Homoscedastic:", homosc2, "\n")
    
    if (normal2 && homosc2) {
      model <- model2; transform_label <- "sqrt"
    } else {
      # LOG
      cat("\n--- Testing log transformation ---\n")
      data  <- data %>% mutate(!!paste0(outcome_var, "_log") := log(.data[[outcome_var]]))
      v3    <- paste0(outcome_var, "_log")
      sel3  <- select_fixed_effects(v3, data); model3 <- sel3$best_model
      plot_assumptions(model3, "log")
      normal3 <- check_normality_ad(model3)
      p_bp3   <- as.numeric(check_heteroscedasticity(model3))
      homosc3 <- p_bp3 > 0.05
      cat("Breusch-Pagan p =", round(p_bp3, 4), "| Homoscedastic:", homosc3, "\n")
      
      if (normal3 && homosc3) {
        model <- model3; transform_label <- "log"
      } else {
        message("  → LOG failed. No further transforms available. Reporting log model results with caution.")
        model <- model3; transform_label <- "log (non-normal)"
      }
    }
  }

  # Collinearity check
  if (!is.null(model)) {
    cat("\n--- Collinearity check (VIF) —", outcome_var, "---\n")
    dev.new(width=12, height=10); print(check_model(model))
  }

  cat("\n→ Final model transform:", transform_label, "\n\n")
  list(model = model, transform = transform_label, data = data)
}

# Export results (fixed effects, Type-III ANOVA, DID contrast, emmeans)
export_results_continuous <- function(outcome_var, model, transform, data, prefix) {
  # Fixed effects
  coef_tbl <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(outcome = outcome_var, transform = transform) %>%
    select(outcome, transform, everything())
  readr::write_csv(coef_tbl, paste0(prefix, "_fixed_effects.csv"))

  # Type-III ANOVA
  anova_tbl <- car::Anova(model, type = "III") %>%
    as.data.frame() %>% rownames_to_column("Effect") %>%
    rename_with(~ sub("^Pr\\(>F\\)$", "p.value", .x)) %>%
    rename_with(~ sub("^Pr\\(>Chisq\\)$", "p.value", .x)) %>%
    mutate(outcome = outcome_var, transform = transform)
  readr::write_csv(anova_tbl, paste0(prefix, "_typeIII.csv"))
  cat("\n--- Type III ANOVA:", outcome_var, "---\n"); print(anova_tbl)

  # DID contrast: (B_post − B_pre) − (A_post − A_pre)
  emm <- emmeans(model, ~ treatment * preorpost)
  cat("\nemmeans cell order check:\n"); print(emm)

  did_tbl <- tryCatch({
    did_raw <- contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>% as.data.frame()
    if ("asymp.LCL" %in% names(did_raw))
      did_raw <- did_raw %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    did_raw %>% transmute(
      outcome = outcome_var, transform = transform,
      contrast = "DID: (B_post−B_pre) − (A_post−A_pre)",
      estimate, SE, df, lower.CL, upper.CL, p.value
    )
  }, error = function(e) { message("DID contrast failed: ", e$message); NULL })

  if (!is.null(did_tbl)) {
    readr::write_csv(did_tbl, paste0(prefix, "_DID.csv"))
    cat("\n--- DID:", outcome_var, "---\n"); print(did_tbl)
  }

  # Marginal means
  emm_summary <- emmeans(model, ~ treatment * preorpost) %>%
    as.data.frame() %>% mutate(outcome = outcome_var, transform = transform)
  readr::write_csv(emm_summary, paste0(prefix, "_emmeans.csv"))

  invisible(list(coef = coef_tbl, anova = anova_tbl, did = did_tbl, emm = emm_summary))
}

# Cook's D influential case check (lmer only)
check_influential <- function(model, outcome_var, data) {
  explicit_formula <- as.formula(paste(outcome_var, "~ treatment * preorpost + (1|id)"))
  m_full <- tryCatch(lmer(explicit_formula, data = data, REML = FALSE),
                     error = function(e) { message("Full model refit failed: ", e$message); NULL })
  if (is.null(m_full)) return(invisible(NULL))

  beta_full <- fixef(m_full); vcov_full <- as.matrix(vcov(m_full))
  ids <- unique(data$id[!is.na(data[[outcome_var]])]); p <- length(beta_full)

  cd <- sapply(ids, function(i) {
    dat_i <- data[data$id != i, ]
    m_i <- tryCatch(lmer(explicit_formula, data = dat_i, REML = FALSE),
                      error = function(e) NULL)
    if (is.null(m_i)) return(NA_real_)
    diff <- fixef(m_i) - beta_full
    tryCatch(as.numeric(t(diff) %*% solve(vcov_full) %*% diff) / p,
             error = function(e) NA_real_)
  })

  names(cd) <- as.character(ids); cd <- cd[!is.na(cd)]
  cut <- 4 / length(cd)
  cat("\n--- Influential cases (Cook's D > 4/n =", round(cut, 3), ") for", outcome_var, "---\n")
  flagged <- cd[cd > cut]
  if (length(flagged) == 0) cat("  None flagged.\n") else print(sort(flagged, decreasing = TRUE))
  invisible(list(cooks_d = cd, cutoff = cut, flagged = flagged))
}

# LOIO sensitivity analysis
loio_sensitivity <- function(outcome_var, influential_result, data, prefix) {
  if (is.null(influential_result)) return(invisible(NULL))
  flagged_ids <- names(influential_result$flagged)
  if (length(flagged_ids) == 0) { cat("No influential cases to drop for LOIO.\n"); return(invisible(NULL)) }

  top_id <- names(sort(influential_result$flagged, decreasing = TRUE))[1]
  cat("\nLOIO: dropping", top_id, "for", outcome_var, "\n")

  dat_drop <- filter(data, id != top_id)
  sel_drop <- select_fixed_effects(outcome_var, dat_drop)
  m_drop <- sel_drop$best_model

  did_from <- function(fit) {
    emm <- emmeans(fit, ~ treatment * preorpost)
    contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>% as.data.frame() %>%
      transmute(estimate, SE, df, lower.CL, upper.CL, p.value)
  }

  m_all_fit <- select_fixed_effects(outcome_var, data)$best_model
  did_all <- did_from(m_all_fit)  %>% mutate(spec = "All participants")
  did_drop <- did_from(m_drop)     %>% mutate(spec = paste0("Drop ", top_id))

  did_comp <- bind_rows(did_all, did_drop) %>%
    mutate(outcome = outcome_var, delta_est = estimate - first(estimate)) %>%
    select(outcome, spec, estimate, SE, df, lower.CL, upper.CL, p.value, delta_est)

  readr::write_csv(did_comp, paste0(prefix, "_LOIO.csv"))
  cat("\n--- LOIO comparison:", outcome_var, "---\n"); print(did_comp)
  invisible(did_comp)
}

# DHARMa diagnostics (lmer models on raw outcome)
dharma_diagnostics <- function(outcome_var, data) {
  mf <- model.frame(
    as.formula(paste(outcome_var, "~ treatment * preorpost")),
    data = data, na.action = na.omit
  )
  dat_cc <- data[as.integer(rownames(mf)), ]
  m_cc <- tryCatch(
    lmer(as.formula(paste(outcome_var, "~ treatment * preorpost + (1|id)")),
         data = dat_cc, REML = FALSE, na.action = na.omit),
    error = function(e) NULL
  )
  if (is.null(m_cc)) { message("DHARMa model failed to fit."); return(invisible(NULL)) }
  set.seed(123)
  simres <- simulateResiduals(fittedModel = m_cc, n = 1000, refit = FALSE)
  plot(simres, main = paste("DHARMa —", outcome_var))
  print(testUniformity(simres))
  print(testOutliers(simres))
  invisible(simres)
}


# --------------------------------------------------------------------------------------------
# RUN CONTINUOUS PIPELINE — PSQI_Global_Score
# --------------------------------------------------------------------------------------------

all_did_results <- list()
all_anova_results <- list()

continuous_outcomes <- list(
  list(var = "PSQI_Global_Score", label = "PSQI Global Score (0–21)")
)

for (o in continuous_outcomes) {
  v <- o$var
  lbl <- o$label
  pfx <- paste0("results_", v)

  cat("\n\n##############################################################\n")
  cat("OUTCOME:", lbl, "\n")
  cat("##############################################################\n")

  hist(psqi[[v]], breaks = 20, col = "grey85",
       main = paste("Raw distribution —", lbl), xlab = lbl)

  result <- check_and_transform(v, psqi)
  model <- result$model
  transform <- result$transform
  data_used <- result$data

  if (is.null(model)) { warning("Model fitting failed for: ", v, ". Skipping."); next }

  cat("\n--- Model summary:", v, "---\n"); print(summary(model))

  par(mfrow = c(1, 2))
  hist(resid(model), main = paste("Residuals —", v), xlab = "Residual")
  plot(fitted(model), resid(model), main = paste("Fitted vs Residuals —", v))
  abline(h = 0, lty = 2); par(mfrow = c(1, 1))
  dharma_diagnostics(v, psqi)

  res <- export_results_continuous(v, model, transform, data_used, pfx)
  if (!is.null(res$did))   all_did_results[[v]]   <- res$did
  if (!is.null(res$anova)) all_anova_results[[v]] <- res$anova

  infl_res <- check_influential(model, v, data_used)
  loio_sensitivity(v, infl_res, psqi, pfx)
}

# Raw model selected 
  # No evidence of non-normality of distribution of residuals on Anderson-Darling and Shapiro-Wilk test
  # No heteroscedasticity on visual inspection
  # No collinearity issues
  # No outliers flagged on DHARMa outlier test
# No evidence of a main or interaction effect on Global PSQI score
# No values flagged for LOIO analysis


# ============================================================================================
# SECTION B: ORDINAL OUTCOMES — PSQI Components 1–7 (0–3 ordered)
# Model: cumulative link mixed model (clmm, ordinal package)
# Fixed-effects selection via AICc; DID-equivalent contrast via emmeans; LOO influence screening → LOIO; combined summary export.
# DHARMa does not support clmm
# ============================================================================================

# NOTE ON COMPONENT 6 (regularity of sleep medication use)
  # Component 6 has a bimodal score distribution with majority of participants scoring 0 or 3
  # This causes Hessian non-positive definitive convergence failure in clmm
  # To get around this, component 6 has been collapsed into a binary variable (0 = never use sleep medication, 1 = any sleep medication use)
  # Component 6 will then be analysed using a logistic mixed model for binary data (glmer)

# --------------------------------------------------------------------------------------------
# ORDINAL HELPER FUNCTIONS
# --------------------------------------------------------------------------------------------

# Fixed-effects model selection for clmm via manual AICc
select_fixed_effects_ordinal <- function(outcome_var, data) {
  f_list <- list(
    a = as.formula(paste(outcome_var, "~ treatment*preorpost + (1|id)")),
    b = as.formula(paste(outcome_var, "~ treatment*preorpost + order + (1|id)")),
    c = as.formula(paste(outcome_var, "~ treatment*preorpost + arm + (1|id)")),
    d = as.formula(paste(outcome_var, "~ treatment*preorpost*order + (1|id)"))
  )
  fits <- lapply(f_list, function(f) {
    tryCatch(clmm(f, data = data, na.action = na.exclude, Hess = TRUE),
             error = function(e) NULL)
  })
  valid <- Filter(Negate(is.null), fits)
  if (length(valid) == 0) { message("All clmm models failed for ", outcome_var); return(NULL) }

  n_obs <- nrow(na.omit(data[, c(outcome_var, "treatment", "preorpost", "id")]))
  aic_df <- do.call(rbind, lapply(seq_along(valid), function(i) {
    m  <- valid[[i]]
    k <- length(coef(m)) + length(m$ST) 
    aic <- AIC(m)
    aicc <- aic + (2*k^2 + 2*k) / max(n_obs - k - 1, 1)
    data.frame(row.names = names(valid)[i], df = k, AICc = round(aicc, 2))
  }))

  best <- which.min(aic_df$AICc)
  cat("Selected clmm fixed-effects structure:", names(valid)[best], "\n")
  print(aic_df)
  list(best_model = valid[[best]], aicc_table = aic_df, selected = names(valid)[best])
}

# Export results for clmm
export_results_ordinal <- function(outcome_var, model, prefix) {

  # Coefficients
  coef_tbl <- broom.mixed::tidy(model, conf.int = TRUE) %>%
    mutate(outcome = outcome_var, model_type = "clmm") %>%
    select(outcome, model_type, everything())
  readr::write_csv(coef_tbl, paste0(prefix, "_fixed_effects.csv"))

  # Likelihood ratio test (Type II, most appropriate for clmm)
  anova_tbl <- tryCatch(
    car::Anova(model, type = "II") %>%
      as.data.frame() %>% rownames_to_column("Effect") %>%
      rename_with(~ sub("^Pr\\(>Chisq\\)$", "p.value", .x)) %>%
      mutate(outcome = outcome_var, model_type = "clmm"),
    error = function(e) {
      message("Anova() failed for ", outcome_var, ": ", e$message); NULL
    }
  )
  if (!is.null(anova_tbl)) {
    readr::write_csv(anova_tbl, paste0(prefix, "_LRT.csv"))
    cat("\n--- LRT:", outcome_var, "---\n"); print(anova_tbl)
  }

  # Marginal means and DID-equivalent contrast on the latent scale
  emm <- tryCatch(
    emmeans(model, ~ treatment * preorpost, mode = "latent"),
    error = function(e) { message("emmeans failed for ", outcome_var); NULL }
  )

  did_tbl <- NULL
  if (!is.null(emm)) {
    cat("\nemmeans cell order check (latent scale):\n"); print(emm)
    did_tbl <- tryCatch({
      did_raw <- contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
        summary(infer = TRUE, adjust = "none") %>% as.data.frame()
      if ("asymp.LCL" %in% names(did_raw))
        did_raw <- did_raw %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
      did_raw %>% transmute(
        outcome = outcome_var, model_type = "clmm",
        contrast = "DID (latent): (B_post−B_pre) − (A_post−A_pre)",
        estimate, SE, df, lower.CL, upper.CL, p.value
      )
    }, error = function(e) { message("DID contrast failed for ", outcome_var); NULL })

    if (!is.null(did_tbl)) {
      readr::write_csv(did_tbl, paste0(prefix, "_DID.csv"))
      cat("\n--- DID (latent scale):", outcome_var, "---\n"); print(did_tbl)
    }

    emm_summary <- as.data.frame(emm) %>% mutate(outcome = outcome_var, model_type = "clmm")
    readr::write_csv(emm_summary, paste0(prefix, "_emmeans.csv"))
  }

  invisible(list(coef = coef_tbl, anova = anova_tbl, did = did_tbl))
}

# Predicted probabilities for each ordinal category by treatment × preorpost cell computed manually from threshold and fixed-effect coefficients because ordinal::predict.clmm() does not support newdata in all package versions
predicted_probs_ordinal <- function(outcome_var, model, data, prefix) {

  cf     <- coef(model)
  th_nms <- names(cf)[grepl("\\|", names(cf))]   # threshold names (e.g. "0|1")
  fe_nms <- names(cf)[!grepl("\\|", names(cf))]  # fixed-effect names
  thetas <- cf[th_nms]
  betas  <- cf[fe_nms]

  # Build the four treatment × preorpost cells and their linear predictors
  # Reference cell: treatment = A, preorpost = 0 → eta = 0
  cells <- expand.grid(
    treatment = levels(data$treatment),
    preorpost = levels(data$preorpost)
  )

  eta <- numeric(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    trt <- as.character(cells$treatment[i])
    pop <- as.character(cells$preorpost[i])

    b_trt <- if (trt == "B" && "treatmentB" %in% fe_nms) betas["treatmentB"] else 0
    b_pop <- if (pop == "1" && "preorpost1" %in% fe_nms) betas["preorpost1"] else 0
    b_int <- if (trt == "B" && pop == "1" && "treatmentB:preorpost1" %in% fe_nms)
               betas["treatmentB:preorpost1"] else 0

    eta[i] <- b_trt + b_pop + b_int
  }

  # Cumulative probabilities via logistic CDF: P(Y <= k) = plogis(theta_k - eta)
  cum_probs <- outer(thetas, eta, function(th, e) plogis(th - e))  # n_thresholds × n_cells
  # Category probabilities: P(Y = k) = P(Y <= k) - P(Y <= k-1)
  n_cats <- length(thetas) + 1
  cat_probs <- matrix(NA_real_, nrow = nrow(cells), ncol = n_cats)
  for (k in seq_len(n_cats)) {
    p_le_k <- if (k <= length(thetas)) cum_probs[k, ]    else rep(1, nrow(cells))
    p_le_km1 <- if (k > 1)              cum_probs[k - 1, ] else rep(0, nrow(cells))
    cat_probs[, k] <- p_le_k - p_le_km1
  }

  colnames(cat_probs) <- paste0("P(score=", 0:(n_cats - 1), ")")

  prob_tbl <- cbind(cells, cat_probs) %>%
    mutate(
      outcome = outcome_var,
      `P(score<=1)` = `P(score=0)` + `P(score=1)`  # Probability of being in the better half of the scale
    )

  cat("\n--- Predicted probabilities (population-level):", outcome_var, "---\n")
  print(prob_tbl)
  readr::write_csv(prob_tbl, paste0(prefix, "_predicted_probs.csv"))
  invisible(prob_tbl)
}

# LOO influence check for clmm
# Computes a Wald-based coefficient perturbation statistic by refitting the model with each participant removed
# The 4/n threshold is a heuristic adapted from OLS Cook's D; treat flagged cases as candidates for LOIO, not definitive outliers.
check_loo_influence_ordinal <- function(model, outcome_var, data) {
  dat_cc <- na.omit(data[, c(outcome_var, "treatment", "preorpost", "id", "order", "arm")])
  f <- as.formula(paste(outcome_var, "~ treatment * preorpost + (1|id)"))
  m_full <- tryCatch(clmm(f, data = dat_cc, Hess = TRUE),
                     error = function(e) { message("Full clmm refit failed"); NULL })
  if (is.null(m_full)) return(invisible(NULL))

  beta_full <- coef(m_full)
  vcov_full <- tryCatch(as.matrix(vcov(m_full)), error = function(e) NULL)
  if (is.null(vcov_full)) return(invisible(NULL))

  ids <- unique(dat_cc$id); p <- length(beta_full)

  loo_stat <- sapply(ids, function(i) {
    dat_i <- dat_cc[dat_cc$id != i, ]
    m_i <- tryCatch(clmm(f, data = dat_i, Hess = TRUE), error = function(e) NULL)
    if (is.null(m_i)) return(NA_real_)
    # Align coefficient names (threshold params may differ slightly)
    b_i <- coef(m_i)
    nms <- intersect(names(beta_full), names(b_i))
    diff <- b_i[nms] - beta_full[nms]
    vc <- vcov_full[nms, nms]
    tryCatch(as.numeric(t(diff) %*% solve(vc) %*% diff) / length(nms),
             error = function(e) NA_real_)
  })

  names(loo_stat) <- as.character(ids); loo_stat <- loo_stat[!is.na(loo_stat)]
  cut <- 4 / length(loo_stat)  # Heuristic threshold adapted from OLS Cook's D; treat flagged cases as candidates for LOIO, not definitive outliers
  cat("\n--- LOO influence (Wald-based coefficient perturbation > 4/n =", round(cut, 3), ") for", outcome_var, "---\n")
  flagged <- loo_stat[loo_stat > cut]
  if (length(flagged) == 0) cat("  None flagged.\n") else print(sort(flagged, decreasing = TRUE))
  invisible(list(loo_stat = loo_stat, cutoff = cut, flagged = flagged))
}

loio_sensitivity_ordinal <- function(outcome_var, influential_result, data, prefix) {
  if (is.null(influential_result) || length(influential_result$flagged) == 0) {
    cat("No influential cases for LOIO (", outcome_var, ").\n"); return(invisible(NULL))
  }
  top_id <- names(sort(influential_result$flagged, decreasing = TRUE))[1]
  cat("\nLOIO (ordinal): dropping", top_id, "for", outcome_var, "\n")

  f <- as.formula(paste(outcome_var, "~ treatment * preorpost + (1|id)"))
  dat_all  <- na.omit(data[, c(outcome_var, "treatment", "preorpost", "id")])
  dat_drop <- filter(dat_all, id != top_id)

  did_from_clmm <- function(dat) {
    m   <- tryCatch(clmm(f, data = dat, Hess = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    emm <- tryCatch(emmeans(m, ~ treatment * preorpost, mode = "latent"),
                    error = function(e) NULL)
    if (is.null(emm)) return(NULL)
    did_raw <- contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>% as.data.frame()
    if ("asymp.LCL" %in% names(did_raw))
      did_raw <- did_raw %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    did_raw %>% transmute(estimate, SE, df, lower.CL, upper.CL, p.value)
  }

  did_all  <- did_from_clmm(dat_all)  %>% mutate(spec = "All participants")
  did_drop <- did_from_clmm(dat_drop) %>% mutate(spec = paste0("Drop ", top_id))
  if (is.null(did_all) || is.null(did_drop)) return(invisible(NULL))

  did_comp <- bind_rows(did_all, did_drop) %>%
    mutate(outcome = outcome_var, delta_est = estimate - first(estimate))
  readr::write_csv(did_comp, paste0(prefix, "_LOIO.csv"))
  cat("\n--- LOIO (ordinal):", outcome_var, "---\n"); print(did_comp)
  invisible(did_comp)
}


# --------------------------------------------------------------------------------------------
# RUN ORDINAL PIPELINE — Components 1–5 + 7 (Component 6 not analysed using clmm - see above)
# --------------------------------------------------------------------------------------------

all_did_results_ordinal   <- list()
all_anova_results_ordinal <- list()

for (v in component_vars) {
  lbl <- gsub("PSQI_Component", "Component ", v)
  pfx <- paste0("results_", v)

  cat("\n\n##############################################################\n")
  cat("OUTCOME (ordinal):", lbl, "\n")
  cat("##############################################################\n")

  # Frequency distribution
  print(table(psqi[[v]], useNA = "ifany"))

  sel <- select_fixed_effects_ordinal(v, psqi)
  if (is.null(sel)) next

  model <- sel$best_model
  cat("\n--- Model summary:", v, "---\n"); print(summary(model))

  # Visual diagnostics
    # Random effects inspected visually to characterise between-participant heterogeneity
    # clmm is robust to violations of normality of random effects - results are not conditional to these checks
  dat_cc <- na.omit(psqi[, c(v, "treatment", "preorpost", "id")])
  f_diag <- as.formula(paste(v, "~ treatment * preorpost + (1|id)"))
  m_diag <- tryCatch(clmm(f_diag, data = dat_cc, Hess = TRUE), error = function(e) NULL)
  
  if (!is.null(m_diag)) {
    dev.new(width = 8, height = 4); par(mfrow = c(1, 2))
    
    # Random effects QQ plot
    re <- ranef(m_diag)$id[, 1]
    qqnorm(re, main = paste0("Random effects QQ\n(", v, ")"))
    qqline(re, col = "red", lty = 2)
    
    par(mfrow = c(1, 1))
  }

  # Collinearity check
  tryCatch({
    cat("\n--- Collinearity check (VIF) —", v, "---\n")
    print(check_collinearity(model))
  }, error = function(e) {
    message("Collinearity check failed for ", v, ": ", e$message)
  })
  
  # Export
  res <- export_results_ordinal(v, model, pfx)
  if (!is.null(res$did))   all_did_results_ordinal[[v]]   <- res$did
  if (!is.null(res$anova)) all_anova_results_ordinal[[v]] <- res$anova

  # Predicted probabilities
  predicted_probs_ordinal(v, model, psqi, pfx)

  # Influential cases and LOIO
  infl_res <- check_loo_influence_ordinal(model, v, psqi)
  loio_sensitivity_ordinal(v, infl_res, psqi, pfx)
}

# Model a was chosen as the best fit via AICc for every component
# No issues with multicollinearity found for any component

# INTERPRETATION OF VISUAL DIAGNOSTICS
# Component 1: mild heavy-tailedness
# Component 2: heavy-tailedness and possible skew
# Component 3: heavy left tail (floor effect) - reflects heavily skewed distribution of scores 
# Component 4: minor deviation at lower tail
# Component 5: S-shaped pattern with several extremes at each end - high missingness of data for component 5
# Component 7: mild heavy-tailedness
# clmm is robust to these violations but caution should be taken when interpreting components 3 and 5 due to floor effect and high missingness, respectively

# INFLUENTIAL CASES AND LOIO
# Component 1:
  # No participants flagged for LOIO analysis
# Component 2:
  # R27 flagged for LOIO analysis
  # Both CIs overlap and remain non-significant, CIs both cross zero
  # The results are robust to exclusion of R27
# Component 3:
  # R02 flagged for LOIO analysis
  # Both CIs overlap and remain non-significant, CIs both cross zero
  # The results are robust to exclusion of R02
# Component 4:
  # No participants flagged for LOIO analysis
# Component 5:
  # R03 flagged for LOIO analysis
  # Both CIs overlap and remain non-significant, CIs both cross zero
  # The results are robust to exclusion of R03
# Component 7:
  # R32 flagged for LOIO analysis
  # Both CIs overlap and remain non-significant, CIs both cross zero
  # The results are robust to exclusion of R32


# --------------------------------------------------------------------------------------------
# COMPONENT 6 — BINARY LOGISTIC MIXED MODEL
# --------------------------------------------------------------------------------------------
# Component 6 has been recoded as a binary variable (0 = never use sleep medication vs 1 = any use) and analysed using a logistic mixed model (glmer, binomial family).

# Recode Component 6 as binary
psqi$PSQI_Component6_binary <- ifelse(psqi$PSQI_Component6 == 0, 0, 1)
psqi$PSQI_Component6_binary <- as.integer(psqi$PSQI_Component6_binary)

cat("\n\n##############################################################\n")
cat("OUTCOME (binary logistic): Component 6 — Sleep Medication Use\n")
cat("##############################################################\n")

# Frequency distribution
cat("\nOriginal Component 6 distribution:\n")
print(table(psqi$PSQI_Component6, useNA = "ifany"))
cat("\nBinary Component 6 distribution (0 = never, 1 = any use):\n")
print(table(psqi$PSQI_Component6_binary, useNA = "ifany"))

# --- Model selection via AICc ---
f6_list <- list(
  a = PSQI_Component6_binary ~ treatment * preorpost + (1|id),
  b = PSQI_Component6_binary ~ treatment * preorpost + order + (1|id),
  c = PSQI_Component6_binary ~ treatment * preorpost + arm + (1|id),
  d = PSQI_Component6_binary ~ treatment * preorpost * order + (1|id)
)

fits6 <- lapply(f6_list, function(f) {
  tryCatch(
    glmer(f, data = psqi, family = binomial(link = "logit"),
          na.action = na.exclude,
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL
  )
})

valid6 <- Filter(Negate(is.null), fits6)
n_obs6 <- nrow(na.omit(psqi[, c("PSQI_Component6_binary", "treatment", "preorpost", "id")]))
aic6 <- do.call(rbind, lapply(seq_along(valid6), function(i) {
  m <- valid6[[i]]
  k <- attr(logLik(m), "df")
  data.frame(row.names = names(valid6)[i], df = k,
             AICc = AIC(m) + (2*k^2 + 2*k) / max(n_obs6 - k - 1, 1))
}))
best6 <- which.min(aic6$AICc)

cat("\nComponent 6 binary — model selection:\n")
print(aic6)
cat("Selected fixed-effects structure:", names(valid6)[best6], "\n")

model6 <- valid6[[best6]]
cat("\n--- Model summary: PSQI_Component6_binary ---\n")
print(summary(model6))

# --- Assumption checks ---
# Collinearity check
tryCatch({
  cat("\n--- Collinearity check (VIF) — PSQI_Component6_binary ---\n")
  print(check_collinearity(model6))
}, error = function(e) {
  message("Collinearity check failed for Component 6 binary: ", e$message)
})
# No multicollinearity issues (VIF < 5)

# Normality of random effects (visual only — glmer is robust to violations)
re6 <- ranef(model6)$id[, 1]
dev.new(width = 5, height = 5)
qqnorm(re6, main = "Random effects QQ\n(PSQI_Component6_binary)")
qqline(re6, col = "red", lty = 2)
# S-shaped distribution with heavy tails

# DHARMa diagnostics on complete cases to avoid NA dimension mismatch
set.seed(123)
simres6 <- tryCatch({
  dat_dharma6 <- na.omit(psqi[, c("PSQI_Component6_binary", "treatment", "preorpost", "id")])
  m_dharma6 <- glmer(PSQI_Component6_binary ~ treatment * preorpost + (1|id),
                       data    = dat_dharma6,
                       family  = binomial(link = "logit"),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  sim <- simulateResiduals(fittedModel = m_dharma6, n = 1000, refit = FALSE)
  plot(sim, main = "DHARMa — PSQI_Component6_binary")
  print(testUniformity(sim))
  print(testOutliers(sim))
  sim
}, error = function(e) {
  message("DHARMa failed for Component 6 binary: ", e$message)
  NULL
})
# No outliers flagged

# --- Type-III Anova ---
anova6 <- tryCatch(
  car::Anova(model6, type = "III") %>%
    as.data.frame() %>%
    rownames_to_column("Effect") %>%
    rename_with(~ sub("^Pr\\(>Chisq\\)$", "p.value", .x)) %>%
    mutate(outcome = "PSQI_Component6_binary", model_type = "glmer_binary"),
  error = function(e) { message("Anova() failed for Component 6 binary: ", e$message); NULL }
)
if (!is.null(anova6)) {
  readr::write_csv(anova6, "results_PSQI_Component6_binary_typeIII.csv")
  cat("\n--- Type III ANOVA: PSQI_Component6_binary ---\n"); print(anova6)
}

# --- Fixed effects table ---
coef6 <- tryCatch(
  broom.mixed::tidy(model6, effects = "fixed", conf.int = TRUE) %>%
    mutate(outcome = "PSQI_Component6_binary", model_type = "glmer_binary") %>%
    select(outcome, model_type, everything()),
  error = function(e) {
    message("Fixed effects table failed for Component 6 binary: ", e$message)
    NULL
  }
)
if (!is.null(coef6)) readr::write_csv(coef6, "results_PSQI_Component6_binary_fixed_effects.csv")

# --- Emmeans and DID contrast (log-odds scale) ---
emm6 <- tryCatch(
  emmeans(model6, ~ treatment * preorpost, type = "link"),
  error = function(e) { message("emmeans failed for Component 6 binary"); NULL }
)

did6_tbl <- NULL
if (!is.null(emm6)) {
  cat("\nemmeans cell order check (log-odds scale):\n"); print(emm6)

  did6_tbl <- tryCatch({
    did6_raw <- contrast(emm6, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>%
      as.data.frame()
    if ("asymp.LCL" %in% names(did6_raw))
      did6_raw <- did6_raw %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    did6_raw %>% transmute(
      outcome    = "PSQI_Component6_binary",
      model_type = "glmer_binary",
      contrast   = "DID (log-odds): (B_post-B_pre) - (A_post-A_pre)",
      estimate, SE, df, lower.CL, upper.CL, p.value
    )
  }, error = function(e) { message("DID contrast failed for Component 6 binary"); NULL })

  if (!is.null(did6_tbl)) {
    readr::write_csv(did6_tbl, "results_PSQI_Component6_binary_DID.csv")
    cat("\n--- DID (log-odds scale): PSQI_Component6_binary ---\n"); print(did6_tbl)
  }

  emm6_summary <- as.data.frame(emm6) %>%
    mutate(outcome = "PSQI_Component6_binary", model_type = "glmer_binary")
  readr::write_csv(emm6_summary, "results_PSQI_Component6_binary_emmeans.csv")
}

# --- Influential cases (LOO coefficient perturbation — consistent with ordinal pipeline) ---
cat("\n--- Influential cases (LOO perturbation) for PSQI_Component6_binary ---\n")

dat_cc6  <- na.omit(psqi[, c("PSQI_Component6_binary", "treatment", "preorpost", "id")])
f6_infl  <- PSQI_Component6_binary ~ treatment * preorpost + (1|id)
m6_full  <- tryCatch(
  glmer(f6_infl, data = dat_cc6, family = binomial(link = "logit"),
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
  error = function(e) NULL
)

if (!is.null(m6_full)) {
  beta_full <- fixef(m6_full)
  vcov_full <- tryCatch(as.matrix(vcov(m6_full)), error = function(e) NULL)
  
  if (!is.null(vcov_full)) {
    ids6 <- unique(dat_cc6$id)
    
    cd6 <- sapply(ids6, function(i) {
      dat_i <- dat_cc6[dat_cc6$id != i, ]
      m_i <- tryCatch(
        glmer(f6_infl, data = dat_i, family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
        error = function(e) NULL
      )
      if (is.null(m_i)) return(NA_real_)
      b_i <- fixef(m_i)
      nms <- intersect(names(beta_full), names(b_i))
      diff <- b_i[nms] - beta_full[nms]
      vc <- vcov_full[nms, nms]
      tryCatch(as.numeric(t(diff) %*% solve(vc) %*% diff) / length(nms),
               error = function(e) NA_real_)
    })
    
    names(cd6) <- as.character(ids6)
    cd6 <- cd6[!is.na(cd6)]
    cut6 <- 4 / length(cd6)
    
    cat("Cook's D cutoff (4/n =", round(cut6, 3), "):\n")
    flagged6 <- cd6[cd6 > cut6]
    
    if (length(flagged6) == 0) {
      cat("  None flagged.\n")
    } else {
      print(sort(flagged6, decreasing = TRUE))
      
      # LOIO sensitivity
      top6 <- names(sort(flagged6, decreasing = TRUE))[1]
      cat("\nLOIO: dropping", top6, "for PSQI_Component6_binary\n")
      dat_drop6 <- filter(dat_cc6, id != top6)
      
      did_from_glmer <- function(dat) {
        m <- tryCatch(
          glmer(f6_infl, data = dat, family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
          error = function(e) NULL
        )
        if (is.null(m)) return(NULL)
        emm <- tryCatch(
          emmeans(m, ~ treatment * preorpost, type = "link"),
          error = function(e) NULL
        )
        if (is.null(emm)) return(NULL)
        did_raw <- contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
          summary(infer = TRUE, adjust = "none") %>% as.data.frame()
        if ("asymp.LCL" %in% names(did_raw))
          did_raw <- did_raw %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
        did_raw %>% transmute(estimate, SE, df, lower.CL, upper.CL, p.value)
      }
      
      did6_all <- did_from_glmer(dat_cc6)  %>% mutate(spec = "All participants")
      did6_drop <- did_from_glmer(dat_drop6) %>% mutate(spec = paste0("Drop ", top6))
      
      if (!is.null(did6_all) && !is.null(did6_drop)) {
        did6_comp <- bind_rows(did6_all, did6_drop) %>%
          mutate(outcome = "PSQI_Component6_binary",
                 delta_est = estimate - first(estimate))
        readr::write_csv(did6_comp, "results_PSQI_Component6_binary_LOIO.csv")
        cat("\n--- LOIO: PSQI_Component6_binary ---\n"); print(did6_comp)
      }
    }
  }
}
# R31 flagged for LOIO analysis
# Both CIs overlap and remain non-significant, CIs both cross zero
# The results are robust to exclusion of R31

# Store results for combined summary
all_did_results_binary   <- if (!is.null(did6_tbl)) list(PSQI_Component6_binary = did6_tbl) else list()
all_anova_results_binary <- if (!is.null(anova6))   list(PSQI_Component6_binary = anova6)   else list()

# Stacked bar chart for Component 6 binary
comp6_summ <- psqi %>%
  filter(!is.na(PSQI_Component6_binary)) %>%
  mutate(score = factor(PSQI_Component6_binary, labels = c("Never", "Any use"))) %>%
  group_by(treatment, preorpost, score) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(treatment, preorpost) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp6 <- ggplot(comp6_summ, aes(x = preorpost, y = pct, fill = score)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.3) +
  facet_wrap(~ treatment,
             labeller = labeller(treatment = c("A" = "Treatment A", "B" = "Treatment B"))) +
  scale_fill_brewer(palette = "Blues", name = "Medication use") +
  labs(title = "Score distribution — Component 6 (binary)",
       x = "Period (pre/post)", y = "Percentage (%)") +
  theme_bw()
print(p_comp6)
ggsave("stacked_bar_PSQI_Component6_binary.png", plot = p_comp6, width = 7, height = 4, dpi = 150)


# --------------------------------------------------------------------------------------------
# COMBINED SUMMARY TABLES
# --------------------------------------------------------------------------------------------

# Continuous DID summary
did_combined_continuous <- bind_rows(all_did_results)
readr::write_csv(did_combined_continuous, "SUMMARY_DID_PSQI_continuous.csv")
cat("\n\n=== COMBINED DID SUMMARY (continuous) ===\n")
print(did_combined_continuous %>% select(outcome, transform, estimate, SE, lower.CL, upper.CL, p.value))

# Continuous Type-III ANOVA — interaction row
anova_combined_continuous <- bind_rows(all_anova_results) %>%
  filter(grepl("treatment:preorpost|treatment.preorpost", Effect, ignore.case = TRUE))
readr::write_csv(anova_combined_continuous, "SUMMARY_typeIII_interaction_PSQI_continuous.csv")
cat("\n=== TYPE-III INTERACTION TERMS (continuous) ===\n")
print(anova_combined_continuous)

# Ordinal DID summary (latent scale)
did_combined_ordinal <- bind_rows(all_did_results_ordinal)
readr::write_csv(did_combined_ordinal, "SUMMARY_DID_PSQI_ordinal.csv")
cat("\n\n=== COMBINED DID SUMMARY (ordinal, latent scale) ===\n")
print(did_combined_ordinal %>% select(outcome, estimate, SE, lower.CL, upper.CL, p.value))

# Ordinal LRT — interaction row
anova_combined_ordinal <- bind_rows(all_anova_results_ordinal) %>%
  filter(grepl("treatment:preorpost|treatment.preorpost", Effect, ignore.case = TRUE))
readr::write_csv(anova_combined_ordinal, "SUMMARY_LRT_interaction_PSQI_ordinal.csv")
cat("\n=== LRT INTERACTION TERMS (ordinal) ===\n")
print(anova_combined_ordinal)

# Component 6 binary DID summary
if (length(all_did_results_binary) > 0) {
  did_combined_binary <- bind_rows(all_did_results_binary)
  readr::write_csv(did_combined_binary, "SUMMARY_DID_PSQI_Component6_binary.csv")
  cat("\n\n=== DID SUMMARY (Component 6 binary, log-odds scale) ===\n")
  print(did_combined_binary %>% select(outcome, estimate, SE, lower.CL, upper.CL, p.value))
}

# Component 6 binary ANOVA — interaction row
if (length(all_anova_results_binary) > 0) {
  anova_combined_binary <- bind_rows(all_anova_results_binary) %>%
    filter(grepl("treatment:preorpost|treatment.preorpost", Effect, ignore.case = TRUE))
  readr::write_csv(anova_combined_binary, "SUMMARY_typeIII_interaction_PSQI_Component6_binary.csv")
  cat("\n=== TYPE-III INTERACTION TERM (Component 6 binary) ===\n")
  print(anova_combined_binary)
}


# --------------------------------------------------------------------------------------------
# VISUALISATIONS
# --------------------------------------------------------------------------------------------

# Spaghetti plots — global score
p <- ggplot(psqi, aes(x = preorpost, y = PSQI_Global_Score,
                             group = id, colour = id)) +
  geom_line(alpha = 0.5) + geom_point(size = 1.8) +
  facet_wrap(~ treatment,
             labeller = labeller(treatment = c("A" = "Treatment A", "B" = "Treatment B"))) +
  labs(title = "Individual trajectories — PSQI Global Score",
       x = "Period (pre/post)", y = "PSQI Global Score") +
  theme_bw() + theme(legend.position = "none")
print(p)
ggsave("spaghetti_PSQI_Global_Score.png", plot = p, width = 8, height = 5, dpi = 150)

# Mean ± SE — global score
summ_global <- psqi %>%
  group_by(treatment, preorpost) %>%
  summarise(
    mean_val = mean(PSQI_Global_Score, na.rm = TRUE),
    se_val = sd(PSQI_Global_Score, na.rm = TRUE) / sqrt(sum(!is.na(PSQI_Global_Score))),
    .groups = "drop"
  )

p2 <- ggplot(summ_global, aes(x = preorpost, y = mean_val,
                               colour = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                width = 0.1, position = position_dodge(0.2)) +
  geom_line(position = position_dodge(0.2), linewidth = 0.9) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  scale_colour_manual(values = c("A" = "#2196F3", "B" = "#F44336"),
                      labels = c("A" = "Treatment A", "B" = "Treatment B")) +
  labs(title = "Mean ± SE — PSQI Global Score",
       x = "Period (pre/post)", y = "PSQI Global Score", colour = "") +
  theme_bw()
print(p2)
ggsave("mean_SE_PSQI_Global_Score.png", plot = p2, width = 6, height = 4, dpi = 150)

# Stacked bar charts — component scores (ordinal)
for (v in component_vars) {
  lbl <- gsub("PSQI_Component", "Component ", v)
  comp_summ <- psqi %>%
    filter(!is.na(.data[[v]])) %>%
    group_by(treatment, preorpost, score = .data[[v]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(treatment, preorpost) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  p3 <- ggplot(comp_summ, aes(x = preorpost, y = pct, fill = score)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.3) +
    facet_wrap(~ treatment,
               labeller = labeller(treatment = c("A" = "Treatment A", "B" = "Treatment B"))) +
    scale_fill_brewer(palette = "Blues", name = "Score") +
    labs(title = paste("Score distribution —", lbl),
         x = "Period (pre/post)", y = "Percentage (%)") +
    theme_bw()
  print(p3)
  ggsave(paste0("stacked_bar_", v, ".png"), plot = p3, width = 7, height = 4, dpi = 150)
}

cat("\n\nAnalysis complete. All CSV and PNG files written to the working directory.\n")
