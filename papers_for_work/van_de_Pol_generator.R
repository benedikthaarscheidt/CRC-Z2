# Load required libraries
library(lme4) # For fitting mixed models
library(MASS) # For generating multivariate normal distributions
library(ggplot2) # For plotting

# 1. Function to simulate data for a single dataset
simulate_random_regression_data <- function(I = 5, J = 10, 
                                            VE = 0.4, VR = 0.6, VS = 0.2, rES = 0.25, plot = TRUE) {
  
  # Total number of observations
  N <- I * J 
  
  # Generate a shared environmental spectrum for all individuals
  X <- rep(seq(-3, 3, length.out = J), I) 
  
  # Generate the covariance matrix for intercepts and slopes
  cov_matrix <- matrix(c(VE, rES * sqrt(VE * VS), 
                         rES * sqrt(VE * VS), VS), nrow = 2)
  
  # Generate random intercepts (u0) and slopes (u1) for each individual (multivariate normal distribution)
  random_effects <- mvrnorm(I, mu = c(0, 0), Sigma = cov_matrix)
  u0 <- rep(random_effects[, 1], each = J) # Intercept for each observation
  u1 <- rep(random_effects[, 2], each = J) # Slope for each observation
  
  # Generate the population mean intercept (b0) and slope (b1) (both set to 0 as per instructions)
  b0 <- 0
  b1 <- 0
  
  # Generate the residuals for each observation
  residuals <- rnorm(N, mean = 0, sd = sqrt(VR))
  
  # Calculate Y according to the model Y_ij = b0 + u0_i + (b1 + u1_i) * X_ij + e_ij
  Y <- b0 + u0 + (b1 + u1) * X + residuals
  
  # Create a data frame for the generated data
  data <- data.frame(
    Individual = factor(rep(1:I, each = J)), # Factor to visualize groups
    X = X,
    Y = Y
  )
  
  # Plot the reaction norms for each individual (optional)
  if (plot) {
    p <- ggplot(data, aes(x = X, y = Y, color = Individual, group = Individual)) +
      geom_line() +
      labs(title = 'Random Regression for Individuals', 
           x = 'X (Environmental Gradient)', 
           y = 'Y (Phenotypic Trait)') +
      theme_minimal() +
      theme(legend.position = 'none') # Remove legend for clarity
    print(p) # This ensures the plot is displayed
  }
  
  return(data)
}

# 2. Function to fit a random regression model to the data
fit_random_regression_model <- function(data) {
  
  # Fit the mixed model with random intercepts and random slopes
  model <- lmer(Y ~ X + (1 + X | Individual), data = data, REML = TRUE)
  
  var_components <- as.data.frame(VarCorr(model))
  ve <- var_components[var_components$grp == "Individual" & var_components$var1 == "(Intercept)", "vcov"]
  vs <- var_components[var_components$grp == "Individual" & var_components$var1 == "X", "vcov"]
  cov_es <- var_components[var_components$grp == "Individual" & var_components$var1 == "(Intercept)" & var_components$var2 == "X", "vcov"]
  vr <- var_components[var_components$grp == "Residual", "vcov"]
  
  # Calculate the correlation between intercepts and slopes
  correlation_es <- cov_es / (sqrt(ve) * sqrt(vs))
  
  return(list(
    model = model,
    VE = ve,
    VS = vs,
    Cov_ES = cov_es,
    Cor_ES = correlation_es,
    VR = vr
  ))
}


set.seed(42) 

# Simulate one dataset and ensure the plot is displayed
simulated_data <- simulate_random_regression_data(I = 3, J = 10, VE = 0.5, VR = 0.5, VS = 0, rES = 1, plot = TRUE)

# Fit a random regression model to the simulated data
fit <- fit_random_regression_model(simulated_data)

# View the variance components from one model
print(fit)
