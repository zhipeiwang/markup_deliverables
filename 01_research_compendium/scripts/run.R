
setwd("scripts")

# Load the RData file from the data folder
load("../data/current_data.RData")

########### Without historic data ################################################

# four types of priors ("correct" ones) ----------------------------------------
# Priors specifications
prior_settings <- list(
  uninformative = list(intercept_prior_mu = 0, intercept_prior_sd = 1000, slope1_prior_mu = 0, slope1_prior_sd = 1000, slope2_prior_mu = 0, slope2_prior_sd = 1000, tau_shape_prior = 0.001, tau_rate_prior = 0.001),
  regularizing = list(intercept_prior_mu = 0, intercept_prior_sd = 10, slope1_prior_mu = 0, slope1_prior_sd = 10, slope2_prior_mu = 0, slope2_prior_sd = 10, tau_shape_prior = 2, tau_rate_prior = 2),
  principled = list(intercept_prior_mu = 30, intercept_prior_sd = 5, slope1_prior_mu = 15, slope1_prior_sd = 5, slope2_prior_mu = 1, slope2_prior_sd = 5, tau_shape_prior = 3, tau_rate_prior = 3),
  informative = list(intercept_prior_mu = 30, intercept_prior_sd = 1, slope1_prior_mu = 15, slope1_prior_sd = 1, slope2_prior_mu = 1, slope2_prior_sd = 0.1, tau_shape_prior = 5, tau_rate_prior = 150)
)

# Run Gibbs sampler for each prior setting
source("workflowFunctions.R")
sim_results_uninformative <- run_gibbs_sampler(prior_settings$uninformative)
sim_results_regularizing <- run_gibbs_sampler(prior_settings$regularizing)
sim_results_principled <- run_gibbs_sampler(prior_settings$principled)
sim_results_informative <- run_gibbs_sampler(prior_settings$informative)

# Summarize results for each set of priors
summary_uninformative <- summarize_simulation_results(sim_results_uninformative, "Uninformative")
summary_regularizing <- summarize_simulation_results(sim_results_regularizing, "Regularizing")
summary_principled <- summarize_simulation_results(sim_results_principled, "Principled")
summary_informative <- summarize_simulation_results(sim_results_informative, "Informative")

# Combine summaries into a single data frame
summary_all <- rbind(summary_uninformative, summary_regularizing, summary_principled, summary_informative)


# Save the summary results
save(summary_all, file = "../output/summary_results.RData")



## false informative priors ----------------------------------------------------
# New informative priors (potentially "false" beliefs)
prior_settings_false <- list(
  informative_false = list(
    intercept_prior_mu = 50, intercept_prior_sd = 1, 
    slope1_prior_mu = 10, slope1_prior_sd = 1, 
    slope2_prior_mu = 2, slope2_prior_sd = 0.1, 
    tau_shape_prior = 20, tau_rate_prior = 20
  )
)

# Run the Gibbs sampler with the new informative priors
source("gibbsSamplers.R")
sim_results_informative_false <- gibbs_sampler(
  seed = 123, 
  x1 = as.numeric(as.character(df$x1)), 
  x2 = df$x2,
  y = df$y, 
  intercept_prior_mu = prior_settings_false$informative_false$intercept_prior_mu, 
  intercept_prior_prec = 1 / (prior_settings_false$informative_false$intercept_prior_sd^2),
  slope1_prior_mu = prior_settings_false$informative_false$slope1_prior_mu, 
  slope1_prior_prec = 1 / (prior_settings_false$informative_false$slope1_prior_sd^2), 
  slope2_prior_mu = prior_settings_false$informative_false$slope2_prior_mu, 
  slope2_prior_prec = 1 / (prior_settings_false$informative_false$slope2_prior_sd^2),
  tau_shape_prior = prior_settings_false$informative_false$tau_shape_prior, 
  tau_rate_prior = prior_settings_false$informative_false$tau_rate_prior,
  intercept_start_value = 0, 
  slope1_start_value = 0,
  slope2_start_value = 0,
  n_sims = 6002,
  burn_in = 1002
)


# Summarize and output results
source("workflowFunctions.R")
summary_false <- summarize_simulation_results_v2(sim_results_informative_false)

# Save the summary results
save(summary_false, file = "../output/Summary_False_Informative_Priors.RData")


## false principled priors -----------------------------------------------------

# Generate False Principled Priors
prior_settings_false_principled <- list(
  intercept_prior_mu = 50, intercept_prior_sd = 20, 
  slope1_prior_mu = 10, slope1_prior_sd = 10, 
  slope2_prior_mu = 2, slope2_prior_sd = 2, 
  tau_shape_prior = 3, tau_rate_prior = 3
)

# Run the Gibbs sampler for False Principled Priors
source("gibbsSamplers.R")
sim_results_false_principled <- gibbs_sampler(
  seed = 123, 
  x1 = as.numeric(as.character(df$x1)), 
  x2 = df$x2,
  y = df$y, 
  intercept_prior_mu = prior_settings_false_principled$intercept_prior_mu, 
  intercept_prior_prec = 1 / (prior_settings_false_principled$intercept_prior_sd^2),
  slope1_prior_mu = prior_settings_false_principled$slope1_prior_mu, 
  slope1_prior_prec = 1 / (prior_settings_false_principled$slope1_prior_sd^2), 
  slope2_prior_mu = prior_settings_false_principled$slope2_prior_mu, 
  slope2_prior_prec = 1 / (prior_settings_false_principled$slope2_prior_sd^2),
  tau_shape_prior = prior_settings_false_principled$tau_shape_prior, 
  tau_rate_prior = prior_settings_false_principled$tau_rate_prior,
  intercept_start_value = 0, 
  slope1_start_value = 0,
  slope2_start_value = 0,
  n_sims = 6002,
  burn_in = 1002
)

# Summarize and output results
source("workflowFunctions.R")
summary_false_principled <- summarize_simulation_results_v2(sim_results_false_principled)

# Save the summary results
save(summary_false_principled, file = "../output/Summary_False_Principled_Priors.RData")


########### With historic data ################################################

# Load the RData file from the data folder
load("../data/historic_data.RData")

# Set priors
intercept_prior_mu <- 0
intercept_prior_prec <- 0.001
slope1_prior_mu <- 0
slope1_prior_prec <- 0.001
slope2_prior_mu <- 0
slope2_prior_prec <- 0.001
tau_shape_prior <- 0.001
tau_rate_prior <- 0.001

intercept_start_value <- 0
slope1_start_value <- 0
slope2_start_value <- 0

# Run the Gibbs sampler for different values of a0 (the weight given to the historic data)
a0_values <- seq(0, 1, by = 0.2)
results_list <- list()

source("gibbsSamplers.R")
for (a0 in a0_values) {
  results_list[[paste0("a0_", a0)]] <- gibbs_sampler_hist(
    seed = 123, 
    x1 = as.numeric(as.character(df$x1)), 
    x2 = df$x2, 
    y = df$y, 
    x1_hist = as.numeric(as.character(df_hist1$x1)), 
    x2_hist = df_hist1$x2, 
    y_hist = df_hist1$y, 
    intercept_prior_mu = intercept_prior_mu, 
    intercept_prior_prec = intercept_prior_prec, 
    slope1_prior_mu = slope1_prior_mu, 
    slope1_prior_prec = slope1_prior_prec, 
    slope2_prior_mu = slope2_prior_mu, 
    slope2_prior_prec = slope2_prior_prec, 
    tau_shape_prior = tau_shape_prior, 
    tau_rate_prior = tau_rate_prior, 
    intercept_start_value = intercept_start_value, 
    slope1_start_value = slope1_start_value, 
    slope2_start_value = slope2_start_value, 
    n_sims = 6002, 
    burn_in = 1002, 
    a0 = a0
  )
}

# Extract summaries for all results
source("workflowFunctions.R")
summaries <- lapply(results_list, extract_summary)

# Combine summaries into one data frame
summaries_df <- do.call(rbind, summaries)
summaries_df$a0 <- rep(a0_values, each = 4)

# Plotting the results
library(ggplot2)

# Define the real values for each parameter
real_values <- data.frame(Parameter = c("Intercept", "Slope1", "Slope2", "Sigma"),
                          Value = c(30, 15, 1, 5))

# Define y-axis limits for each parameter
y_limits <- list(
  Intercept = c(20, 50),
  Slope1 = c(10, 20),
  Slope2 = c(0.5, 1.5),
  Sigma = c(2, 8)
)

# Plotting
ggplot(summaries_df, aes(x = a0, y = Mean, color = Parameter)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1) +
  geom_hline(data = real_values, aes(yintercept = Value), linetype = "dashed", color = "black") +
  facet_wrap(~ Parameter, scales = "free_y", labeller = label_both) +
  theme_minimal() +
  labs(x = "a0",
       y = "Value", title = "Posterior Means and Credible Intervals for Different Values of a0", subtitle = "where a0 is the weight of the historic data") +
  scale_color_manual(values = c("Intercept" = "blue", "Slope1" = "red", "Slope2" = "green", "Sigma" = "purple")) +
  theme(legend.position = "none")

# Save the plot
ggsave("../output/a0_results_plot.png", width = 10, height = 8, units = "in")




