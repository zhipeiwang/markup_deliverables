
#### generate the fake current data ------------------------------------------------

setwd("scripts")
# Set seed for reproducibility
set.seed(123)

# Simulate the independent binary variable x1 as a factor
x1 <- factor(c(rep(1, 8), rep(0, 12)))

# Simulate the independent continuous variable x2 as normally distributed with mean 40 and sd 10
x2 <- rnorm(20, mean = 40, sd = 10)

# Randomly select indices to adjust
indices_to_adjust <- sample(1:20, size = 10)

# Adjust selected x2 values to have smaller values
x2[indices_to_adjust] <- x2[indices_to_adjust] - 10

# Calculate the dependent variable y
y <- 30 + 15 * as.numeric(as.character(x1)) + x2 + rnorm(20, mean = 0, sd = 5)

# Create the dataframe
df <- data.frame(x1 = x1, x2 = x2, y = y)

# Save the dataframe
save(df, file = "../data/current_data.RData")


#### generate the fake historic data ------------------------------------------------

# Set seed for reproducibility
set.seed(456)

# Simulate the independent binary variable x1 as a factor
x1_hist1 <- factor(c(rep(1, 10), rep(0, 10)))

# Simulate the independent continuous variable x2 as normally distributed with mean 40 and sd 10
x2_hist1 <- rnorm(20, mean = 40, sd = 10)

# Calculate the dependent variable y
y_hist1 <- 30 + 15 * as.numeric(as.character(x1_hist1)) + x2_hist1 + rnorm(20, mean = 0, sd = 5)

# Create the dataframe
df_hist1 <- data.frame(x1 = x1_hist1, x2 = x2_hist1, y = y_hist1)

# Save the dataframe
save(df_hist1, file = "../data/historic_data.RData")

