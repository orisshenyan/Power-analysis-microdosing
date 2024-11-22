# clear workspace
rm(list=ls())

# set wd
hablar::set_wd_to_script_path()

# libraries
library(glmmTMB)
library(tidyverse)

# ------------------------------------------------------------------------------
# data
Count <- read.csv("Count.csv")

Count$X <- NULL
Count <- Count %>%
  dplyr::select(-Proportional_count)
str(Count)

#hallucination count from original study
summary_stats <- Count %>%
  group_by(Simple_Complex, Condition) %>%
  summarise(
    mean_hallucination = mean(Hallucination_Count, na.rm = TRUE),
    sd_hallucination = sd(Hallucination_Count, na.rm = TRUE)
  )

summary_stats

# ------------------------------------------------------------------------------
# estimated desired effect size as rate-ratio

ganzflicker_simple_mu_count <- summary_stats %>%
  filter(Simple_Complex == "Simple" & Condition == "Flicker") %>%
  pull(mean_hallucination)

ganzflicker_simple_sd_count <- summary_stats %>%
  filter(Simple_Complex == "Simple" & Condition == "Flicker") %>%
  pull(sd_hallucination)

ganzflicker_complex_mu_count <- summary_stats %>%
  filter(Simple_Complex == "Complex" & Condition == "Flicker") %>%
  pull(mean_hallucination)

ganzflicker_complex_sd_count <- summary_stats %>%
  filter(Simple_Complex == "Complex" & Condition == "Flicker") %>%
  pull(sd_hallucination)

ganzfeld_simple_mu_count <- summary_stats %>%
  filter(Simple_Complex == "Simple" & Condition == "Ganzfeld") %>%
  pull(mean_hallucination)

ganzfeld_simple_sd_count <- summary_stats %>%
  filter(Simple_Complex == "Simple" & Condition == "Ganzfeld") %>%
  pull(sd_hallucination)

ganzfeld_complex_mu_count <- summary_stats %>%
  filter(Simple_Complex == "Complex" & Condition == "Ganzfeld") %>%
  pull(mean_hallucination)

ganzfeld_complex_sd_count <- summary_stats %>%
  filter(Simple_Complex == "Complex" & Condition == "Ganzfeld") %>%
  pull(sd_hallucination)

effect_size_d <- 0.8

#ganzflicker simple
est_mean_psilo_ganzflicker_simple_count <- effect_size_d *ganzflicker_simple_sd_count + ganzflicker_simple_mu_count

#ganzflicker complex
est_mean_psilo_ganzflicker_complex_count <- effect_size_d *ganzflicker_complex_sd_count + ganzflicker_complex_mu_count

#ganzfeld simple
est_mean_psilo_ganzfeld_simple_count <- effect_size_d *ganzfeld_simple_sd_count + ganzfeld_simple_mu_count

#ganzfeld complex
est_mean_psilo_ganzfeld_complex_count <- 0.8*ganzfeld_complex_sd_count + ganzfeld_complex_mu_count


# 
rate_ratios <- c(est_mean_psilo_ganzflicker_simple_count/ganzflicker_simple_mu_count,
                 est_mean_psilo_ganzflicker_complex_count/ganzflicker_complex_mu_count,
                 est_mean_psilo_ganzfeld_simple_count /ganzfeld_simple_mu_count,
                 est_mean_psilo_ganzfeld_complex_count /ganzfeld_complex_mu_count )

# effect size in rate-ratio
rate_ratio_drug <- mean(rate_ratios)


# ------------------------------------------------------------------------------
Count$Duration <- ifelse(Count$Condition == 'Ganzfeld', 25, 15)
mm_nb <- glmmTMB(Hallucination_Count ~ Condition*Simple_Complex + offset(log(Duration)) + (1|Participant_Number),
                               family = "nbinom2",
                               data=Count)
summary(mm_nb)

# extract parameters from fitted model (these will be used in the function below)
fx <- fixef(mm_nb)$cond
dispersion_theta <- sigma(mm_nb) # this is the estimate dispersion parameters
random_intercept_sd <- unname(sqrt(VarCorr(mm_nb)[[1]]$Participant_Number)[1,1])

# ------------------------------------------------------------------------------
# Function that creates simulated count data
simulate_count_data <- function(n) {
  
  # Fixed effects coefficients (based on original model and desired effect sizes)
  beta0 <- fx[1]      # Adjusted Intercept to control for baseline counts
  beta1 <- fx[2]      # ConditionGanzfeld
  beta2 <- fx[3]      # Simple_ComplexSimple
  beta3 <- fx[4]      # ConditionGanzfeld:Simple_ComplexSimple
  
  # Effect size for 'Drug' in terms of rate ratio
  beta4 <- log(rate_ratio_drug)  # DrugPsilocybin effect
  
  # Interaction effects involving Drug (set to zero or specify as needed)
  beta5 <- 0          # ConditionGanzfeld:DrugPsilocybin
  beta6 <- 0          # Simple_ComplexSimple:DrugPsilocybin
  beta7 <- 0          # ConditionGanzfeld:Simple_ComplexSimple:DrugPsilocybin
  
  # beta <- c(beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7)
  # names(beta) <- c("(Intercept)", "ConditionGanzfeld", "Simple_ComplexSimple", 
  #                  "ConditionGanzfeld:Simple_ComplexSimple", "DrugPsilocybin", 
  #                  "ConditionGanzfeld:DrugPsilocybin", "Simple_ComplexSimple:DrugPsilocybin", 
  #                  "ConditionGanzfeld:Simple_ComplexSimple:DrugPsilocybin")
  
  
  # Random effects variance components
  Var_b <- c(random_intercept_sd^2,   # Random intercept
             (random_intercept_sd^2)/2, # Random slope for main effect of Drug
             0,                          # Random slope for Ganzfeld:Psilocybin
             0,                          # Random slope for Simple:Psilocybin
             0)                          # Random slope for Ganzfeld:Simple:Psilocybin
  
  # Random effects correlation
  correlation <- 0.5
  
  # Construct variance-covariance matrix for random effects
  Sigma <- matrix(0, nrow = 5, ncol = 5)
  for (i in 1:5) {
    for (j in 1:5) {
      if (i == j) {
        Sigma[i, j] <- Var_b[i]
      } else {
        Sigma[i, j] <- correlation * sqrt(Var_b[i] * Var_b[j])
      }
    }
  }
  
  # Generate random effects for each participant
  random_effects <- MASS::mvrnorm(n = n, mu = rep(0, 5), Sigma = Sigma)
  colnames(random_effects) <- c("b0i", "b1i", "b2i", "b3i", "b4i")
  
  # Create data frame with all combinations of factors for each participant
  participants <- 1:n
  conditions <- c("Flicker", "Ganzfeld")
  complexity <- c("Complex", "Simple")
  drugs <- c("Placebo", "Psilocybin")
  
  df <- expand_grid(
    Participant_Number = participants,
    Condition = conditions,
    Simple_Complex = complexity,
    Drug = drugs
  )
  
  # Merge random effects with data frame
  random_effects_df <- data.frame(Participant_Number = participants, random_effects)
  df <- merge(df, random_effects_df, by = "Participant_Number")
  
  # Create indicator variables for fixed effects
  df$ConditionGanzfeld <- ifelse(df$Condition == "Ganzfeld", 1, 0)
  df$Simple_ComplexSimple <- ifelse(df$Simple_Complex == "Simple", 1, 0)
  df$DrugPsilocybin <- ifelse(df$Drug == "Psilocybin", 1, 0)
  
  # Interaction terms
  df$ConditionGanzfeld_Simple_ComplexSimple <- df$ConditionGanzfeld * df$Simple_ComplexSimple
  df$ConditionGanzfeld_DrugPsilocybin <- df$ConditionGanzfeld * df$DrugPsilocybin
  df$Simple_ComplexSimple_DrugPsilocybin <- df$Simple_ComplexSimple * df$DrugPsilocybin
  df$ConditionGanzfeld_Simple_ComplexSimple_DrugPsilocybin <- df$ConditionGanzfeld * df$Simple_ComplexSimple * df$DrugPsilocybin
  
  # Create random effects design variables (indexes)
  df$z1 <- df$DrugPsilocybin
  df$z2 <- df$DrugPsilocybin * df$ConditionGanzfeld
  df$z3 <- df$DrugPsilocybin * df$Simple_ComplexSimple
  df$z4 <- df$DrugPsilocybin * df$ConditionGanzfeld * df$Simple_ComplexSimple
  
  # Compute random effects contribution
  df$random_effects_contribution <- df$b0i +
    df$b1i * df$z1 +
    df$b2i * df$z2 +
    df$b3i * df$z3 +
    df$b4i * df$z4
  
  # Compute linear predictor (eta) with fixed effects and random effects
  df$eta_fixed <- beta0 +
    beta1 * df$ConditionGanzfeld +
    beta2 * df$Simple_ComplexSimple +
    beta3 * df$ConditionGanzfeld_Simple_ComplexSimple +
    beta4 * df$DrugPsilocybin +
    beta5 * df$ConditionGanzfeld_DrugPsilocybin +
    beta6 * df$Simple_ComplexSimple_DrugPsilocybin +
    beta7 * df$ConditionGanzfeld_Simple_ComplexSimple_DrugPsilocybin
  
  df$eta <- df$eta_fixed + df$random_effects_contribution
  
  # Offset term (Duration depends on Condition)
  df$Duration <- ifelse(df$Condition == "Ganzfeld", 25, 15)
  df$offset <- log(df$Duration)
  
  # Adjust eta with offset
  df$eta_offset <- df$eta + df$offset
  
  # Compute expected counts (mu)
  df$mu <- exp(df$eta_offset)
  
  # Generate counts from negative binomial distribution
  df$Hallucination_Count <- rnbinom(nrow(df), size = dispersion_theta, mu = df$mu)
  
  # Return the data frame with relevant columns
  return(df[, c("Participant_Number", "Condition", "Simple_Complex", "Drug", "Duration", "Hallucination_Count")])
}

# ------------------------------------------------------------------------------
# test
df <- simulate_count_data(30)

# some sanity checks

# counts of original data
with(Count, tapply(Hallucination_Count, list(Condition, Simple_Complex), mean))
with(Count, tapply(Hallucination_Count, list(Condition, Simple_Complex), sd))

# placebo should be ~ the same (on average)
with(df[df$Drug=="Placebo",], tapply(Hallucination_Count, list(Condition, Simple_Complex), mean))
with(df[df$Drug=="Placebo",], tapply(Hallucination_Count, list(Condition, Simple_Complex), sd))

# here the mean cound should be ~twice as large (on average)
with(df[df$Drug=="Psilocybin",], tapply(Hallucination_Count, list(Condition, Simple_Complex), mean))
with(df[df$Drug=="Psilocybin",], tapply(Hallucination_Count, list(Condition, Simple_Complex), sd))

# fit model on simulated data
mm_nb_sim <- glmmTMB(Hallucination_Count ~ Condition*Simple_Complex * Drug + offset(log(Duration)) + (Drug|Participant_Number),
                 family = "nbinom2",
                 data=df)
summary(mm_nb_sim)


