# Function to run Gibbs sampler for different prior settings
source("gibbsSamplers.R")
run_gibbs_sampler <- function(prior_setting) {
  intercept_prior_prec <- 1 / (prior_setting$intercept_prior_sd^2)
  slope1_prior_prec <- 1 / (prior_setting$slope1_prior_sd^2)
  slope2_prior_prec <- 1 / (prior_setting$slope2_prior_sd^2)
  
  sim_results <- gibbs_sampler(
    seed = 123, 
    x1 = as.numeric(as.character(df$x1)), 
    x2 = df$x2,
    y = df$y, 
    intercept_prior_mu = prior_setting$intercept_prior_mu, 
    intercept_prior_prec = intercept_prior_prec,
    slope1_prior_mu = prior_setting$slope1_prior_mu, 
    slope1_prior_prec = slope1_prior_prec, 
    slope2_prior_mu = prior_setting$slope2_prior_mu, 
    slope2_prior_prec = slope2_prior_prec,
    tau_shape_prior = prior_setting$tau_shape_prior, 
    tau_rate_prior = prior_setting$tau_rate_prior,
    intercept_start_value = 0, 
    slope1_start_value = 0,
    slope2_start_value = 0,
    n_sims = 6002,
    burn_in = 1002
  )
  return(sim_results)
}



# Function to summarize simulation results
summarize_simulation_results <- function(sim_results, prior_type) {
  intercept_posterior_mean <- round(mean(sim_results$intercept, na.rm = TRUE), 2)
  intercept_posterior_sd <- round(sd(sim_results$intercept, na.rm = TRUE), 2)
  intercept_posterior_cred_int <- round(quantile(sim_results$intercept, probs = c(0.025, 0.975), na.rm = TRUE), 2)
  
  slope1_posterior_mean <- round(mean(sim_results$slope1, na.rm = TRUE), 2)
  slope1_posterior_sd <- round(sd(sim_results$slope1, na.rm = TRUE), 2)
  slope1_posterior_cred_int <- round(quantile(sim_results$slope1, probs = c(0.025, 0.975), na.rm = TRUE), 2)
  
  slope2_posterior_mean <- round(mean(sim_results$slope2, na.rm = TRUE), 2)
  slope2_posterior_sd <- round(sd(sim_results$slope2, na.rm = TRUE), 2)
  slope2_posterior_cred_int <- round(quantile(sim_results$slope2, probs = c(0.025, 0.975), na.rm = TRUE), 2)
  
  sigma_posterior_mean <- round(mean(sim_results$sigma, na.rm = TRUE), 2)
  sigma_posterior_sd <- round(sd(sim_results$sigma, na.rm = TRUE), 2)
  sigma_posterior_cred_int <- round(quantile(sim_results$sigma, probs = c(0.025, 0.975), na.rm = TRUE), 2)
  
  summary_stats <- data.frame(
    Prior_Type = prior_type,
    Parameter = c("Intercept", "Slope1", "Slope2", "Sigma"),
    Mean = c(intercept_posterior_mean, slope1_posterior_mean, slope2_posterior_mean, sigma_posterior_mean),
    SD = c(intercept_posterior_sd, slope1_posterior_sd, slope2_posterior_sd, sigma_posterior_sd),
    Credible_Intervals = c(
      paste(intercept_posterior_cred_int, collapse = " - "), 
      paste(slope1_posterior_cred_int, collapse = " - "), 
      paste(slope2_posterior_cred_int, collapse = " - "),
      paste(sigma_posterior_cred_int, collapse = " - ")
    )
  )
  
  return(summary_stats)
}

# Function to summarize simulation results
summarize_simulation_results_v2 <- function(sim_results) {
  ## Extract summary stats
  intercept_posterior_mean <- round(mean(sim_results$intercept, na.rm = TRUE), 2)
  intercept_posterior_sd <- round(sd(sim_results$intercept, na.rm = TRUE), 2)
  intercept_posterior_cred_int <- round(quantile(sim_results$intercept, c(0.025, 0.975)), 2)
  
  slope1_posterior_mean <- round(mean(sim_results$slope1, na.rm = TRUE), 2)
  slope1_posterior_sd <- round(sd(sim_results$slope1, na.rm = TRUE), 2)
  slope1_posterior_cred_int <- round(quantile(sim_results$slope1, c(0.025, 0.975)), 2)
  
  slope2_posterior_mean <- round(mean(sim_results$slope2, na.rm = TRUE), 2)
  slope2_posterior_sd <- round(sd(sim_results$slope2, na.rm = TRUE), 2)
  slope2_posterior_cred_int <- round(quantile(sim_results$slope2, c(0.025, 0.975)), 2)
  
  sigma_posterior_mean <- round(mean(sim_results$sigma, na.rm = TRUE), 2)
  sigma_posterior_sd <- round(sd(sim_results$sigma, na.rm = TRUE), 2)
  sigma_posterior_cred_int <- round(quantile(sim_results$sigma, c(0.025, 0.975)), 2)
  
  summary_stats <- data.frame(
    Parameter = c("$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\sigma$"),
    Mean = c(intercept_posterior_mean, slope1_posterior_mean, slope2_posterior_mean, sigma_posterior_mean),
    SD = c(intercept_posterior_sd, slope1_posterior_sd, slope2_posterior_sd, sigma_posterior_sd),
    `95% Credible Interval` = c(
      paste(intercept_posterior_cred_int, collapse = " - "), 
      paste(slope1_posterior_cred_int, collapse = " - "), 
      paste(slope2_posterior_cred_int, collapse = " - "),
      paste(sigma_posterior_cred_int, collapse = " - ")
    )
  )
  
  return(summary_stats)
}


# Function to extract summary statistics
extract_summary <- function(sim_results) {
  data.frame(
    Parameter = c("Intercept", "Slope1", "Slope2", "Sigma"),
    Mean = c(mean(sim_results$intercept), mean(sim_results$slope1), mean(sim_results$slope2), mean(sim_results$sigma)),
    SD = c(sd(sim_results$intercept), sd(sim_results$slope1), sd(sim_results$slope2), sd(sim_results$sigma)),
    CI_Lower = c(quantile(sim_results$intercept, 0.025), quantile(sim_results$slope1, 0.025), quantile(sim_results$slope2, 0.025), quantile(sim_results$sigma, 0.025)),
    CI_Upper = c(quantile(sim_results$intercept, 0.975), quantile(sim_results$slope1, 0.975), quantile(sim_results$slope2, 0.975), quantile(sim_results$sigma, 0.975))
  )
}
