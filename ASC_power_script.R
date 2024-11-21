# clear workspace
rm(list=ls())

# set wd
hablar::set_wd_to_script_path()

# libraries
library(truncnorm)
library(lme4)
library(lmerTest)
library(tidyverse)

# ------------------------------------------------------------------------------
# data
ASC <- read.csv("ASC.csv")

# ------------------------------------------------------------------------------
# tydying
ASC <- ASC%>%
  rename(Participant_Number = Participant.number,
         Simple_Complex = ASC.dimension,
         Average_Score = Average.score) %>%
  mutate(Simple_Complex = case_when(
    Simple_Complex == "Elementary imagery" ~ "Simple",
    Simple_Complex == "Complex imagery" ~ "Complex",
    TRUE ~ Simple_Complex)) %>%
  dplyr::select(Participant_Number, Simple_Complex, Average_Score, Condition)

# ------------------------------------------------------------------------------
# estimate original model
ASCmod <- lmer(Average_Score ~ Condition*Simple_Complex + (1|Participant_Number),
               data=ASC)
summary(ASCmod)

# extract parameters from fitted model (these will be used in the function below)
fx <- fixef(ASCmod)
residual_sd <- sigma(ASCmod)
random_intercept_sd <- as.data.frame(VarCorr(ASCmod))$sdcor[1]

# ------------------------------------------------------------------------------
# Function that creates simulated data
simulate_data <- function(n) {
  # Load necessary package
  library(MASS)
  
  # additional assumptions
  random_correlation <- 0.5
  effect_size <- 0.8 # standardized
  
  # Random effects variance components
  # assume smaller variances in the effect of drugs than the variance at baseline
  Var_b <- c(random_intercept_sd^2,   # random intercept
             random_intercept_sd^2/2, # random slope for main effect of Drug
             0                      , # random slope for Ganzfeld:Psilocybin
             0                      , # random slope for Simple:Psilocybin
             0                      ) # random slope for Ganzfeld:Simple:Psilocybin
  
  # effect size of psylocibin (main effect) in units of the dependent variable
  # the effect size is defined as cohen d so we need to pool together all sources
  # of variability (residuals & random effects)
  # see example here: https://mlisi.xyz/files/workshops/LMM101/LMM_part2.html#42
  drug_main_effect <- effect_size * sqrt(residual_sd^2 + 
                                           Var_b[1] +
                                           0.5*Var_b[2] +
                                           random_correlation * sqrt(Var_b[1]) * sqrt(Var_b[2]))
  
  # Fixed parameters 
  beta0 <-fx[1]     # Intercept
  beta1 <- fx[2]    # ConditionGanzfeld
  beta2 <- fx[3]    # Simple_ComplexSimple
  beta3 <- drug_main_effect  # DrugPsilocybin 
  beta4 <- fx[4]    # ConditionGanzfeld:Simple_ComplexSimple
  beta5 <- 0        # ConditionGanzfeld:DrugPsilocybin
  beta6 <- 0        # Simple_ComplexSimple:DrugPsilocybin
  beta7 <- 0        # ConditionGanzfeld:Simple_ComplexSimple:DrugPsilocybin
  
  beta <- c(beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7)
  
  # assigne names (to avoid errors below)
  names(beta) <- c("Intercept","Ganzfeld","Simple","Psilocybin",
                   "Ganzfeld:Simple","Ganzfeld:Psilocybin","Simple:Psilocybin",
                   "Ganzfeld:Simple:Psilocybin")
  
  # note that currently I set some of the variances to zero, which 
  # effectively means we are simulating only random slopes for main effect of Drug
  # (Not all of these might be estimable from the data in practice)
  # if you want to test for an interaction, you should add some heterogeneity in the 
  # interaction parameter of interest by simulating random slopes
  
  correlation <- 0.5  # Correlation among random effects
  
  # Construct variance-covariance matrix for random effects
  Sigma <- matrix(0, nrow=5, ncol=5)
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
  random_effects <- mvrnorm(n = n, mu = rep(0, 5), Sigma = Sigma)
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

  # Create design matrix for fixed effects
  X <- model.matrix(~ Condition * Simple_Complex * Drug, data = df)
  
  # Map random effects to observations
  random_effects_df <- data.frame(Participant_Number = 1:n, random_effects)
  df <- merge(df, random_effects_df, by = "Participant_Number")
  
  # Create random effects design variables (indexes)
  df$z1 <- ifelse(df$Drug == "Psilocybin", 1, 0)
  df$z2 <- ifelse(df$Drug == "Psilocybin" & df$Condition == "Ganzfeld", 1, 0)
  df$z3 <- ifelse(df$Drug == "Psilocybin" & df$Simple_Complex == "Simple", 1, 0)
  df$z4 <- ifelse(df$Drug == "Psilocybin" & df$Condition == "Ganzfeld" & df$Simple_Complex == "Simple", 1, 0)
  
  # Compute random effects contribution
  df$random_effects_contribution <- df$b0i +
    df$b1i * df$z1 +
    df$b2i * df$z2 +
    df$b3i * df$z3 +
    df$b4i * df$z4
  
  # Compute linear predictor (eta)
  eta_fixed <- X %*% beta
  eta_fixed <- as.numeric(eta_fixed)  # Convert matrix to vector
  df$eta <- eta_fixed + df$random_effects_contribution
  
  # Add residual error (bounded)
  df$ASC_Score <- rtruncnorm(nrow(df), mean = df$eta, sd = residual_sd, a=0, b=100)
  
  # Return the data frame with relevant columns
  return(df[, c("Participant_Number", "Condition", "Simple_Complex", "Drug", "ASC_Score")])
}

# test
simulated_data <- simulate_data(n = 30)

# fit model
sim_m <- lmer(ASC_Score ~ Condition * Simple_Complex * Drug + ( Drug | Participant_Number), data=simulated_data)
summary(sim_m)

