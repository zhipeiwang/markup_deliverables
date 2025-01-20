


# gibbs sampler

gibbs_sampler <- function(seed, x1, x2, y, intercept_prior_mu, intercept_prior_prec, 
                          slope1_prior_mu, slope1_prior_prec,
                          slope2_prior_mu, slope2_prior_prec,
                          tau_shape_prior, tau_rate_prior, 
                          intercept_start_value, slope1_start_value, slope2_start_value,
                          n_sims, burn_in) {
  
  ## set seed to ensure reproducibility
  set.seed(seed)
  ## get sample size
  n_obs <- length(y)
  
  ## initial predictions with starting values
  preds1 <- intercept_start_value + slope1_start_value * x1 + slope2_start_value * x2
  sse1 <- sum((y - preds1)^2)
  tau_shape <- tau_shape_prior + n_obs / 2
  
  ## vectors to store values
  sse <- c(sse1, rep(NA, n_sims))
  
  intercept <- c(intercept_start_value, rep(NA, n_sims))
  slope1 <- c(slope1_start_value, rep(NA, n_sims))
  slope2 <- c(slope2_start_value, rep(NA, n_sims))
  tau_rate <- c(NA, rep(NA, n_sims))
  tau <- c(NA, rep(NA, n_sims))
  
  for(i in 2:n_sims){
    
    # Tau Values
    tau_rate[i] <- tau_rate_prior + sse[i - 1]/2
    tau[i] <- rgamma(n = 1, shape = tau_shape, rate = tau_rate[i])
    
    # Intercept Values
    intercept_mu <- (intercept_prior_prec*intercept_prior_mu + 
                       tau[i] * sum(y - slope1[i - 1]*x1 - slope2[i - 1]*x2)) / 
      (intercept_prior_prec + n_obs*tau[i])
    intercept_prec <- intercept_prior_prec + n_obs*tau[i]
    intercept[i] <- rnorm(n = 1, mean = intercept_mu, sd = sqrt(1 / intercept_prec))
    
    # Slope Values
    slope1_mu <- (slope1_prior_prec * slope1_prior_mu + 
                    tau[i] * sum(x1 * (y - intercept[i] - slope2[i - 1] * x2))) / 
      (slope1_prior_prec + tau[i] * sum(x1^2))
    slope1_prec <- slope1_prior_prec + tau[i] * sum(x1^2)
    slope1[i] <- rnorm(n = 1, mean = slope1_mu, sd = sqrt(1 / slope1_prec))
    
    slope2_mu <- (slope2_prior_prec * slope2_prior_mu + 
                    tau[i] * sum(x2 * (y - intercept[i] - slope1[i - 1] * x1))) / 
      (slope2_prior_prec + tau[i] * sum(x2^2))
    slope2_prec <- slope2_prior_prec + tau[i] * sum(x2^2)
    slope2[i] <- rnorm(n = 1, mean = slope2_mu, sd = sqrt(1 / slope2_prec))
    
    preds <- intercept[i] + slope1[i] * x1 + slope2[i] * x2
    sse[i] <- sum((y - preds)^2)
    
  }
  
  sigma <- sqrt(1 / na.omit(tau[-1:-burn_in]))
  
  list(
    intercept = na.omit(intercept[-1:-burn_in]), 
    slope1 = na.omit(slope1[-1:-burn_in]), 
    slope2 = na.omit(slope2[-1:-burn_in]), 
    sigma = sigma
  )
}



#### gibbs sampler that uses historic data ------------------------------------------------




gibbs_sampler_hist <- function(seed, x1, x2, y, x1_hist, x2_hist, y_hist, 
                               intercept_prior_mu, intercept_prior_prec, 
                               slope1_prior_mu, slope1_prior_prec,
                               slope2_prior_mu, slope2_prior_prec,
                               tau_shape_prior, tau_rate_prior, 
                               intercept_start_value, slope1_start_value, slope2_start_value,
                               n_sims, burn_in, a0) {
  
  ## set seed to ensure reproducibility
  set.seed(seed)
  ## get sample sizes
  n_obs <- length(y)
  n_hist <- length(y_hist)
  
  ## initial predictions with starting values
  preds1 <- intercept_start_value + slope1_start_value * x1 + slope2_start_value * x2
  preds1_hist <- intercept_start_value + slope1_start_value * x1_hist + slope2_start_value * x2_hist
  sse1 <- sum((y - preds1)^2) + a0 * sum((y_hist - preds1_hist)^2)
  tau_shape <- tau_shape_prior + (n_obs + a0 * n_hist) / 2
  
  ## vectors to store values
  sse <- c(sse1, rep(NA, n_sims))
  
  intercept <- c(intercept_start_value, rep(NA, n_sims))
  slope1 <- c(slope1_start_value, rep(NA, n_sims))
  slope2 <- c(slope2_start_value, rep(NA, n_sims))
  tau_rate <- c(NA, rep(NA, n_sims))
  tau <- c(NA, rep(NA, n_sims))
  
  for(i in 2:n_sims){
    
    # Tau Values
    tau_rate[i] <- tau_rate_prior + sse[i - 1]/2
    tau[i] <- rgamma(n = 1, shape = tau_shape, rate = tau_rate[i])
    
    # Intercept Values
    intercept_mu <- (intercept_prior_prec * intercept_prior_mu + 
                       tau[i] * (sum(y - slope1[i - 1] * x1 - slope2[i - 1] * x2) + 
                                   a0 * sum(y_hist - slope1[i - 1] * x1_hist - slope2[i - 1] * x2_hist))) / 
      (intercept_prior_prec + tau[i] * (n_obs + a0 * n_hist))
    intercept_prec <- intercept_prior_prec + tau[i] * (n_obs + a0 * n_hist)
    intercept[i] <- rnorm(n = 1, mean = intercept_mu, sd = sqrt(1 / intercept_prec))
    
    # Slope1 Values
    slope1_mu <- (slope1_prior_prec * slope1_prior_mu + 
                    tau[i] * (sum(x1 * (y - intercept[i] - slope2[i - 1] * x2)) + 
                                a0 * sum(x1_hist * (y_hist - intercept[i] - slope2[i - 1] * x2_hist)))) / 
      (slope1_prior_prec + tau[i] * (sum(x1^2) + a0 * sum(x1_hist^2)))
    slope1_prec <- slope1_prior_prec + tau[i] * (sum(x1^2) + a0 * sum(x1_hist^2))
    slope1[i] <- rnorm(n = 1, mean = slope1_mu, sd = sqrt(1 / slope1_prec))
    
    # Slope2 Values
    slope2_mu <- (slope2_prior_prec * slope2_prior_mu + 
                    tau[i] * (sum(x2 * (y - intercept[i] - slope1[i - 1] * x1)) + 
                                a0 * sum(x2_hist * (y_hist - intercept[i] - slope1[i - 1] * x1_hist)))) / 
      (slope2_prior_prec + tau[i] * (sum(x2^2) + a0 * sum(x2_hist^2)))
    slope2_prec <- slope2_prior_prec + tau[i] * (sum(x2^2) + a0 * sum(x2_hist^2))
    slope2[i] <- rnorm(n = 1, mean = slope2_mu, sd = sqrt(1 / slope2_prec))
    
    preds <- intercept[i] + slope1[i] * x1 + slope2[i] * x2
    preds_hist <- intercept[i] + slope1[i] * x1_hist + slope2[i] * x2_hist
    sse[i] <- sum((y - preds)^2) + a0 * sum((y_hist - preds_hist)^2)
    
  }
  
  sigma <- sqrt(1 / na.omit(tau[-1:-burn_in]))
  
  list(
    intercept = na.omit(intercept[-1:-burn_in]), 
    slope1 = na.omit(slope1[-1:-burn_in]), 
    slope2 = na.omit(slope2[-1:-burn_in]), 
    sigma = sigma
  )
}