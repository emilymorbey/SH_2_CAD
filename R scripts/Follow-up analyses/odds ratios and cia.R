# betas and ses to ORs and CIs


# women 

betas_new <- c(0.237, 0.081,
               0.003,
               0.061,
               -0.165,
               -0.006,
               -0.004,
               -0.179,
               0.334,
               0.128,
               0.004,
               -0.019,
               -0.010,
               -0.058,
               0.002,
               -0.080)

# New Standard errors
se_betas_new <- c(0.061,
                  0.092,
                  0.001,
                  0.073,
                  0.040,
                  0.073,
                  0.002,
                  0.056,
                  0.073,
                  0.118,
                  0.002,
                  0.099,
                  0.039,
                  0.071,
                  0.002,
                  0.059)

# Function to convert beta and standard error to odds ratio and confidence interval
convert_to_odds_ratio <- function(beta, se_beta) {
  # Calculate odds ratio
  odds_ratio <- exp(beta)
  
  # Calculate standard error of log odds ratio
  se_log_odds_ratio <- se_beta
  
  # Calculate 95% confidence interval for log odds ratio
  ci_log_odds_ratio <- c(beta - 1.96 * se_beta, beta + 1.96 * se_beta)
  
  # Calculate confidence interval for odds ratio
  ci_odds_ratio <- exp(ci_log_odds_ratio)
  
  result <- list(odds_ratio = odds_ratio,
                 ci_odds_ratio_lower = ci_odds_ratio[1],
                 ci_odds_ratio_upper = ci_odds_ratio[2])
  
  return(result)
}

results_new <- mapply(convert_to_odds_ratio, betas_new, se_betas_new)

results_new
