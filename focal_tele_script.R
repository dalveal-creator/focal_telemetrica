#Load dataset daniel
source("config_local.R")
dataset <- read_csv(file.path(data_path, "mediciones_daniel.csv"))

#rename vars
library(dplyr)
library(stringr)

names(dataset) <- str_to_lower(names(dataset)) %>% str_replace_all(" ", "_")

names(dataset)

#filter by complete cases
dataset <- dataset %>%
  # Conditions 1 & 2: Remove rows with NA in either variable
  filter(!is.na(pt_tele), !is.na(pi_tele),!is.na(ss_tele),!is.na(pt_focal),
         !is.na(pi_focal),!is.na(ss_focal)) 

dataset <- dataset |>
  mutate(
    nombre_paciente = str_replace_all(nombre_paciente, " ", "_")          # spaces -> underscore
  )

#Data-integrity check
dataset %>%
  mutate(pi_resid = pi_tele - (pt_tele + ss_tele)) %>%
  summarise(mean_resid = mean(pi_resid, na.rm = TRUE),
            max_abs    = max(abs(pi_resid), na.rm = TRUE),
            n_off_3deg = sum(abs(pi_resid) > 3, na.rm = TRUE))

dataset %>%
  mutate(pi_resid = pi_tele - (pt_tele + ss_tele)) %>%
  filter(abs(pi_resid) > 3) %>%
  mutate(
    pi_recalc = pt_tele + ss_tele,
    ratio_pt  = pi_resid / pt_tele,   # ≈ 2 → PT sign flip
    ratio_ss  = pi_resid / ss_tele      # ≈ 2 → SS sign flip
  ) %>%
  select(nombre_paciente, pi_tele, pt_tele, ss_tele,
         pi_recalc, pi_resid, ratio_pt, ratio_ss) %>%
  arrange(desc(abs(pi_resid))) #OK

str(dataset)

dataset %>%
  mutate(pi_resid = pi_focal - (pt_focal + ss_focal)) %>%
  summarise(mean_resid = mean(pi_resid, na.rm = TRUE),
            max_abs    = max(abs(pi_resid), na.rm = TRUE),
            n_off_3deg = sum(abs(pi_resid) > 3, na.rm = TRUE))

dataset %>%
  mutate(pi_resid = pi_focal - (pt_focal + ss_focal)) %>%
  filter(abs(pi_resid) > 3) %>%
  mutate(
    pi_recalc = pt_focal + ss_focal,
    ratio_pt  = pi_resid / pt_focal,   # ≈ 2 → PT sign flip
    ratio_ss  = pi_resid / ss_focal      # ≈ 2 → SS sign flip
  ) %>%
  select(nombre_paciente, pi_focal, pt_focal, ss_focal,
         pi_recalc, pi_resid, ratio_pt, ratio_ss) %>%
  arrange(desc(abs(pi_resid)))

# ==============================================================================
# Pelvic Incidence Agreement Analysis: Tele-RX vs Focalised RX
# ==============================================================================

# 1. Load Libraries ------------------------------------------------------------
library(tidyverse)
library(irr)        # For ICC calculation
library(BlandAltmanLeh) # For Bland-Altman plots

# 2. Data Cleaning & Geometric QC ---------------------------------------------
# Verify identity PI = PT + SS (tolerance +/- 2 degrees for manual rounding)
dataset_clean <- dataset %>%
  mutate(
    qc_tele_diff  = abs(pi_tele - (pt_tele + ss_tele)),
    qc_focal_diff = abs(pi_focal - (pt_focal + ss_focal)),
    pi_diff       = pi_tele - pi_focal
  )

# Flag potential measurement inconsistencies
qc_flags <- dataset_clean %>%
  filter(qc_tele_diff > 2 | qc_focal_diff > 2)

if (nrow(qc_flags) > 0) {
  message("Warning: Found ", nrow(qc_flags), " records with PI != PT + SS discrepancies.")
  print(qc_flags %>% select(nombre_paciente, pi_tele, pt_tele, ss_tele, pi_focal, pt_focal, ss_focal))
}

# 3. Descriptive Statistics ----------------------------------------------------
summary_stats <- dataset_clean %>%
  summarise(
    n             = n(),
    mean_tele     = mean(pi_tele, na.rm = TRUE),
    sd_tele       = sd(pi_tele, na.rm = TRUE),
    mean_focal    = mean(pi_focal, na.rm = TRUE),
    sd_focal      = sd(pi_focal, na.rm = TRUE),
    mean_diff     = mean(pi_diff, na.rm = TRUE),
    sd_diff       = sd(pi_diff, na.rm = TRUE),
    median_diff   = median(pi_diff, na.rm = TRUE),
    iqr_diff      = IQR(pi_diff, na.rm = TRUE)
  )
print(summary_stats)

# 4. Hypothesis Testing for Mean/Median Bias -----------------------------------
# Test normality of differences
shapiro_res <- shapiro.test(dataset_clean$pi_diff)
print(shapiro_res)

if (shapiro_res$p.value > 0.05) {
  message("Differences are normally distributed. Running Paired t-test:")
  test_res <- t.test(dataset_clean$pi_tele, dataset_clean$pi_focal, paired = TRUE)
} else {
  message("Differences deviate from normality. Running Wilcoxon signed-rank test:")
  test_res <- wilcox.test(dataset_clean$pi_tele, dataset_clean$pi_focal, paired = TRUE)
}
print(test_res)

# 5. Reliability: Intraclass Correlation Coefficient (ICC) --------------------
# Model: Two-way random/mixed, single score, absolute agreement: ICC(2,1)
icc_data <- dataset_clean %>% select(pi_tele, pi_focal)
icc_res  <- icc(icc_data, model = "twoway", type = "agreement", unit = "single")
print(icc_res)

# 6. Agreement: Bland-Altman Analysis -----------------------------------------
diff_mean <- mean(dataset_clean$pi_diff, na.rm = TRUE)
diff_sd   <- sd(dataset_clean$pi_diff, na.rm = TRUE)
upper_loa <- diff_mean + 1.96 * diff_sd
lower_loa <- diff_mean - 1.96 * diff_sd

cat("Mean Bias:", round(diff_mean, 2), "degrees\n")
cat("95% Limits of Agreement:", round(lower_loa, 2), "to", round(upper_loa, 2), "degrees\n")

# Bland-Altman ggplot
dataset_clean <- dataset_clean %>%
  mutate(pi_mean = (pi_tele + pi_focal) / 2)

p_ba <- ggplot(dataset_clean, aes(x = pi_mean, y = pi_diff)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = diff_mean, color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = c(upper_loa, lower_loa), color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Bland-Altman Plot: Pelvic Incidence (Tele vs Focal)",
    x = "Mean of PI Tele and PI Focal (°)",
    y = "Difference: PI Tele - PI Focal (°)"
  ) +
  annotate("text", x = max(dataset_clean$pi_mean), y = diff_mean, 
           label = paste("Bias:", round(diff_mean, 1)), vjust = -1, hjust = 1, color = "blue") +
  annotate("text", x = max(dataset_clean$pi_mean), y = upper_loa, 
           label = paste("+1.96 SD:", round(upper_loa, 1)), vjust = -1, hjust = 1, color = "red") +
  annotate("text", x = max(dataset_clean$pi_mean), y = lower_loa, 
           label = paste("-1.96 SD:", round(lower_loa, 1)), vjust = 1.5, hjust = 1, color = "red") +
  theme_minimal()

print(p_ba)

library(psych)

icc_psych <- ICC(dataset_clean[, c("pi_tele", "pi_focal")])
# Inspect ICC2 (Two-way random single rater, identical mathematically to ICC(A,1))
print(icc_psych$results)

#XAVI----
#Load dataset_xavi daniel
dataset_xavi <- read_csv(file.path(data_path, "mediciones_xavi.csv"))

#rename vars
library(dplyr)
library(stringr)

names(dataset_xavi) <- str_to_lower(names(dataset_xavi)) %>% str_replace_all(" ", "_")

names(dataset_xavi)

#filter by complete cases
dataset_xavi <- dataset_xavi %>%
  # Conditions 1 & 2: Remove rows with NA in either variable
  filter(!is.na(pt_tele), !is.na(pi_tele),!is.na(ss_tele),!is.na(pt_focal),
         !is.na(pi_focal),!is.na(ss_focal)) 

dataset_xavi <- dataset_xavi |>
  mutate(
    nombre_paciente = str_replace_all(nombre_paciente, " ", "_")          # spaces -> underscore
  )

#Data-integrity check
dataset_xavi %>%
  mutate(pi_resid = pi_tele - (pt_tele + ss_tele)) %>%
  summarise(mean_resid = mean(pi_resid, na.rm = TRUE),
            max_abs    = max(abs(pi_resid), na.rm = TRUE),
            n_off_3deg = sum(abs(pi_resid) > 3, na.rm = TRUE))

dataset_xavi %>%
  mutate(pi_resid = pi_tele - (pt_tele + ss_tele)) %>%
  filter(abs(pi_resid) > 3) %>%
  mutate(
    pi_recalc = pt_tele + ss_tele,
    ratio_pt  = pi_resid / pt_tele,   # ≈ 2 → PT sign flip
    ratio_ss  = pi_resid / ss_tele      # ≈ 2 → SS sign flip
  ) %>%
  select(nombre_paciente, pi_tele, pt_tele, ss_tele,
         pi_recalc, pi_resid, ratio_pt, ratio_ss) %>%
  arrange(desc(abs(pi_resid))) #OK

str(dataset_xavi)

dataset_xavi %>%
  mutate(pi_resid = pi_focal - (pt_focal + ss_focal)) %>%
  summarise(mean_resid = mean(pi_resid, na.rm = TRUE),
            max_abs    = max(abs(pi_resid), na.rm = TRUE),
            n_off_3deg = sum(abs(pi_resid) > 3, na.rm = TRUE))

dataset_xavi %>%
  mutate(pi_resid = pi_focal - (pt_focal + ss_focal)) %>%
  filter(abs(pi_resid) > 3) %>%
  mutate(
    pi_recalc = pt_focal + ss_focal,
    ratio_pt  = pi_resid / pt_focal,   # ≈ 2 → PT sign flip
    ratio_ss  = pi_resid / ss_focal      # ≈ 2 → SS sign flip
  ) %>%
  select(nombre_paciente, pi_focal, pt_focal, ss_focal,
         pi_recalc, pi_resid, ratio_pt, ratio_ss) %>%
  arrange(desc(abs(pi_resid)))

# 2. Data Cleaning & Geometric QC ---------------------------------------------
# Verify identity PI = PT + SS (tolerance +/- 2 degrees for manual rounding)
dataset_clean <- dataset_xavi %>%
  mutate(
    qc_tele_diff  = abs(pi_tele - (pt_tele + ss_tele)),
    qc_focal_diff = abs(pi_focal - (pt_focal + ss_focal)),
    pi_diff       = pi_tele - pi_focal
  )

# Flag potential measurement inconsistencies
qc_flags <- dataset_clean %>%
  filter(qc_tele_diff > 2 | qc_focal_diff > 2)

if (nrow(qc_flags) > 0) {
  message("Warning: Found ", nrow(qc_flags), " records with PI != PT + SS discrepancies.")
  print(qc_flags %>% select(nombre_paciente, pi_tele, pt_tele, ss_tele, pi_focal, pt_focal, ss_focal))
}

# 3. Descriptive Statistics ----------------------------------------------------
summary_stats <- dataset_clean %>%
  summarise(
    n             = n(),
    mean_tele     = mean(pi_tele, na.rm = TRUE),
    sd_tele       = sd(pi_tele, na.rm = TRUE),
    mean_focal    = mean(pi_focal, na.rm = TRUE),
    sd_focal      = sd(pi_focal, na.rm = TRUE),
    mean_diff     = mean(pi_diff, na.rm = TRUE),
    sd_diff       = sd(pi_diff, na.rm = TRUE),
    median_diff   = median(pi_diff, na.rm = TRUE),
    iqr_diff      = IQR(pi_diff, na.rm = TRUE)
  )
print(summary_stats)

# 4. Hypothesis Testing for Mean/Median Bias -----------------------------------
# Test normality of differences
shapiro_res <- shapiro.test(dataset_clean$pi_diff)
print(shapiro_res)

if (shapiro_res$p.value > 0.05) {
  message("Differences are normally distributed. Running Paired t-test:")
  test_res <- t.test(dataset_clean$pi_tele, dataset_clean$pi_focal, paired = TRUE)
} else {
  message("Differences deviate from normality. Running Wilcoxon signed-rank test:")
  test_res <- wilcox.test(dataset_clean$pi_tele, dataset_clean$pi_focal, paired = TRUE)
}
print(test_res)

# 5. Reliability: Intraclass Correlation Coefficient (ICC) --------------------
# Model: Two-way random/mixed, single score, absolute agreement: ICC(2,1)
icc_data <- dataset_clean %>% select(pi_tele, pi_focal)
icc_res  <- icc(icc_data, model = "twoway", type = "agreement", unit = "single")
print(icc_res)

# 6. Agreement: Bland-Altman Analysis -----------------------------------------
diff_mean <- mean(dataset_clean$pi_diff, na.rm = TRUE)
diff_sd   <- sd(dataset_clean$pi_diff, na.rm = TRUE)
upper_loa <- diff_mean + 1.96 * diff_sd
lower_loa <- diff_mean - 1.96 * diff_sd

cat("Mean Bias:", round(diff_mean, 2), "degrees\n")
cat("95% Limits of Agreement:", round(lower_loa, 2), "to", round(upper_loa, 2), "degrees\n")

# Bland-Altman ggplot
dataset_clean <- dataset_clean %>%
  mutate(pi_mean = (pi_tele + pi_focal) / 2)

p_ba <- ggplot(dataset_clean, aes(x = pi_mean, y = pi_diff)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = diff_mean, color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = c(upper_loa, lower_loa), color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Bland-Altman Plot: Pelvic Incidence (Tele vs Focal)",
    x = "Mean of PI Tele and PI Focal (°)",
    y = "Difference: PI Tele - PI Focal (°)"
  ) +
  annotate("text", x = max(dataset_clean$pi_mean), y = diff_mean, 
           label = paste("Bias:", round(diff_mean, 1)), vjust = -1, hjust = 1, color = "blue") +
  annotate("text", x = max(dataset_clean$pi_mean), y = upper_loa, 
           label = paste("+1.96 SD:", round(upper_loa, 1)), vjust = -1, hjust = 1, color = "red") +
  annotate("text", x = max(dataset_clean$pi_mean), y = lower_loa, 
           label = paste("-1.96 SD:", round(lower_loa, 1)), vjust = 1.5, hjust = 1, color = "red") +
  theme_minimal()

print(p_ba)

library(psych)

icc_psych <- ICC(dataset_clean[, c("pi_tele", "pi_focal")])
# Inspect ICC2 (Two-way random single rater, identical mathematically to ICC(A,1))
print(icc_psych$results)

# ==============================================================================
#Daniel vs Xavi (Tele vs Focal)----
# ==============================================================================

library(tidyverse)
library(irr)
library(lme4)
library(lmerTest)

# 1. Add rater identifiers -----------------------------------------------------
dataset      <- dataset %>% mutate(rater = "daniel")
dataset_xavi <- dataset_xavi %>% mutate(rater = "xavi")

# 2 & 3. Merge, filter complete cases, and keep shared patients (n = 18) --------
merged_long <- bind_rows(dataset, dataset_xavi) %>%
  # Filter complete cases for the primary angles
  filter(!is.na(pi_tele) & !is.na(pi_focal)) %>%
  # Keep only patients evaluated by both raters
  group_by(nombre_paciente) %>%
  filter(n_distinct(rater) == 2) %>%
  ungroup()

# Reshape to wide format for pairwise inter-rater analyses
df_paired <- merged_long %>%
  select(nombre_paciente, rater, pi_tele, pi_focal, pt_tele, pt_focal, ss_tele, ss_focal) %>%
  pivot_wider(
    id_cols = nombre_paciente,
    names_from = rater,
    values_from = c(pi_tele, pi_focal, pt_tele, pt_focal, ss_tele, ss_focal)
  )

cat("Number of complete paired patients:", nrow(df_paired), "\n")

# Inter-rater ICC for Full-Standing (Tele)
icc_tele <- icc(
  df_paired %>% select(pi_tele_daniel, pi_tele_xavi),
  model = "twoway", type = "agreement", unit = "single"
)
cat("\n--- Inter-Rater ICC: PI Full-Standing (Tele) ---\n")
print(icc_tele)

# Inter-rater ICC for Focalised
icc_focal <- icc(
  df_paired %>% select(pi_focal_daniel, pi_focal_xavi),
  model = "twoway", type = "agreement", unit = "single"
)
cat("\n--- Inter-Rater ICC: PI Focalised ---\n")
print(icc_focal)

# Normality test on inter-rater differences
shapiro_tele_diff  <- shapiro.test(df_paired$pi_tele_daniel - df_paired$pi_tele_xavi)
shapiro_focal_diff <- shapiro.test(df_paired$pi_focal_daniel - df_paired$pi_focal_xavi)

# Paired test on Tele
if (shapiro_tele_diff$p.value > 0.05) {
  test_rater_tele <- t.test(df_paired$pi_tele_daniel, df_paired$pi_tele_xavi, paired = TRUE)
} else {
  test_rater_tele <- wilcox.test(df_paired$pi_tele_daniel, df_paired$pi_tele_xavi, paired = TRUE)
}
cat("\n--- Inter-Rater Bias Test (Tele) ---\n")
print(test_rater_tele)

# Paired test on Focal
if (shapiro_focal_diff$p.value > 0.05) {
  test_rater_focal <- t.test(df_paired$pi_focal_daniel, df_paired$pi_focal_xavi, paired = TRUE)
} else {
  test_rater_focal <- wilcox.test(df_paired$pi_focal_daniel, df_paired$pi_focal_xavi, paired = TRUE)
}
cat("\n--- Inter-Rater Bias Test (Focal) ---\n")
print(test_rater_focal)

#Examine whether Xavi independently reproduces the $+6.5^\circ$ tele-vs-focal bias that Daniel observed.
df_paired <- df_paired %>%
  mutate(
    diff_modality_daniel = pi_tele_daniel - pi_focal_daniel,
    diff_modality_xavi   = pi_tele_xavi - pi_focal_xavi
  )

cat("\n--- Modality Differences (Tele - Focal) ---\n")
cat("Daniel Mean Bias (Tele - Focal):", round(mean(df_paired$diff_modality_daniel), 2), "°\n")
cat("Xavi Mean Bias (Tele - Focal):  ", round(mean(df_paired$diff_modality_xavi), 2), "°\n")

# Paired t-test for Xavi (Tele vs Focal)
t_test_xavi_modality <- t.test(df_paired$pi_tele_xavi, df_paired$pi_focal_xavi, paired = TRUE)
print(t_test_xavi_modality)

#This is the definitive test. 
#It models all four measurements per patient simultaneously and separates variance attributable to the patient, 
#the rater, the imaging modality, and their interaction.

# Prepare long dataset with standard factors
df_lmm <- merged_long %>%
  select(nombre_paciente, rater, pi_tele, pi_focal) %>%
  pivot_longer(
    cols = c(pi_tele, pi_focal),
    names_to = "modality",
    names_prefix = "pi_",
    values_to = "pi"
  )

# Model: Fixed effects for modality and rater; random intercept per patient
lmm_fit <- lmer(pi ~ modality * rater + (1 | nombre_paciente), data = df_lmm)

cat("\n--- Linear Mixed Model ANOVA (Type III) ---\n")
print(anova(lmm_fit))

cat("\n--- Model Coefficients & Fixed Effects ---\n")
print(summary(lmm_fit))