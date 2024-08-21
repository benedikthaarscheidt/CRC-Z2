#' Generate Synthetic Plant Plasticity Data with Sequential Environments
#'
#' This function generates a synthetic dataset for plant plasticity measurements across multiple environments. 
#' It simulates traits for a specified number of plants in different environments, with trait values normally 
#' distributed around given baseline values. The generated dataset is structured such that the first set of rows 
#' corresponds to all plants in the first environment, followed by all plants in the second environment, and so on. 
#'
#' @param n_plants An integer specifying the number of plants (rows) to generate in each environment.
#' @param baseline_values A matrix or data frame where each row represents an environment and each column 
#' represents a trait. The values in this matrix provide the baseline (mean) values for each trait in each environment.
#' @param within_variance A matrix or data frame where each row represents an environment and each column 
#' represents a trait. The values in this matrix specify the standard deviation for each trait in each environment, 
#' determining the variability around the baseline values.
#' 
#' @return A data frame containing the synthetic dataset. The data frame will have `n_plants * n_environments` rows 
#' and `n_traits` columns. The row names indicate the plant number and environment (e.g., 1.1 for plant 1 in environment 1).
#' The column names represent the traits (e.g., `Trait_1`, `Trait_2`, etc.).
#' 
#' @details This function can be used for simulating data where different environments may have different 
#' levels of variability for the same trait. The data is structured sequentially by environment, making it easy to analyze 
#' the effects of different environmental conditions on plant traits.
#' 
#' @examples
#' # Define the parameters
#' n_plants = 100
#' 
#' # Define baseline values for each trait in each environment
#' # Rows are environments, columns are traits
#' baseline_values = matrix(
#'   c(100, 110, 120,  # Trait 1 baseline values in Env 1, 2, 3
#'     200, 195, 205,  # Trait 2 baseline values in Env 1, 2, 3
#'     300, 310, 290), # Trait 3 baseline values in Env 1, 2, 3
#'   nrow = 3, ncol = 3, byrow = TRUE
#' )
#' 
#' # Define within-environment variance for each trait in each environment
#' within_variance = matrix(
#'   c(10, 15, 20,   # Trait 1 variance in Env 1, 2, 3
#'     8, 12, 10,    # Trait 2 variance in Env 1, 2, 3
#'     5, 7, 6),     # Trait 3 variance in Env 1, 2, 3
#'   nrow = 3, ncol = 3, byrow = TRUE
#' )
#' 
#' # Generate the synthetic dataset
#' synthetic_data = generate_synthetic_data(n_plants, baseline_values, within_variance)
#' 
#' # View the first few rows of the generated dataset
#' head(synthetic_data)
#' 
#' @export
generate_synthetic_data = function(n_plants, baseline_values, within_variance, environmental_factor) {
  n_traits = ncol(baseline_values)
  n_environments = nrow(baseline_values)
  
  # Initialize empty data frame to store the results, with an additional column for environmental factor
  synthetic_data = data.frame(matrix(nrow = n_plants * n_environments, ncol = n_traits + 2))
  
  # Generate row names where the integer part is the plant id and the fractional part is the environment indicator 
  rownames(synthetic_data) = rep(1:n_plants, times = n_environments) + 
    rep(seq(0.1, (n_environments - 1) / 10 + 0.1, by = 0.1), each = n_plants)
  
  # Set column names for the traits
  colnames(synthetic_data) = c("Species identificator", "Environmental Factor", paste("Trait", 1:n_traits, sep = "_"))
  
  for (env in 1:n_environments) {
    # Calculate the row indices for the current environment
    start_idx = (env - 1) * n_plants + 1
    end_idx = env * n_plants
    
    # Assign environmental factor and environment indicator
    synthetic_data[start_idx:end_idx, 1] = env
    synthetic_data[start_idx:end_idx, 2] = environmental_factor[env]
    
    for (i in 1:n_traits) {
      # Generate random data for the current trait and environment with specific variance
      synthetic_data[start_idx:end_idx, i + 2] = rnorm(n_plants, mean = baseline_values[env, i], sd = within_variance[env, i])
    }
  }
  
  return(synthetic_data)
}

n_plants = 100

# Define baseline values for each trait in each environment
baseline_values = matrix(
  c(100, 11, 1,  
    100, 52, 500,  
    100, 101, 1020), 
  nrow = 3, ncol = 3, byrow = TRUE
)

# Define within-environment variance for each trait in each environment
within_variance = matrix(
  c(10, 5, 5,   
    10, 12, 10,    
    10, 7, 6),     
  nrow = 3, ncol = 3, byrow = TRUE
)
environmental_factor=c(0.2,0.5,0.9)

set.seed(12345)
synthetic_data1 = generate_synthetic_data(n_plants, baseline_values, within_variance,environmental_factor)
print(class(synthetic_data1))

########################################

# Function to check and install missing packages
check_and_install_packages = function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}

########################################

#' Calculate the Phenotypic Stability Index (PSI) Using Linear Regression
#'
#' This function calculates the Phenotypic Stability Index (PSI) for a specific trait across multiple environments.
#' The PSI is derived from the regression coefficient obtained by regressing the trait values against the environmental factor.
#' A lower regression coefficient indicates higher phenotypic stability.
#'
#' The environmental factors can be provided either as a column in the data frame or as an external vector.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number or name representing the environmental conditions within the data frame, or NULL if using an external vector.
#' @param external_env_factor (Optional) A numeric vector representing external environmental factors. If provided, this will override `env_col`.
#' @return The Phenotypic Stability Index (PSI), with higher values indicating greater stability.
#' @examples
#' df = data.frame(
#'   Environment = rep(1:3, each = 10),
#'   Trait_1 = c(rnorm(10, 100, 10), rnorm(10, 110, 15), rnorm(10, 120, 20)),
#'   Trait_2 = c(rnorm(10, 200, 8), rnorm(10, 195, 12), rnorm(10, 205, 10))
#' )
#' 
#' # Calculate PSI for Trait 1 using the internal environment column
#' PSI_internal = calculate_PSI(df, trait_col = "Trait_1", env_col = "Environment")
#' print(PSI_internal)
#' 
#' # Calculate PSI for Trait 1 using an external environmental factor vector
#' external_env = rep(1:3, each = 10)
#' PSI_external = calculate_PSI(df, trait_col = "Trait_1", env_col = NULL, external_env_factor = external_env)
#' print(PSI_external)
#' @export
calculate_PSI = function(data, trait_col, env_col = NULL, external_env_factor = NULL) {
  # Extract the trait values
  trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Determine the environmental values to use
  if (!is.null(external_env_factor)) {
    environment_values = external_env_factor
  } else if (!is.null(env_col)) {
    environment_values = if (is.numeric(env_col)) data[[env_col]] else data[[env_col]]
  } else {
    stop("Either 'env_col' or 'external_env_factor' must be provided.")
  }
  
  # Ensure environment values are treated as numeric
  environment_values = as.numeric(environment_values)
  
  # Check for near-constant environmental values
  if (sd(environment_values) < .Machine$double.eps) {
    warning("Environmental factor has too little variation, making the PSI calculation unreliable.")
  }
  
  # Perform a linear regression of trait values on environment values
  model = lm(trait_values ~ environment_values)
  
  # Extract the regression coefficient (slope)
  stability_coefficient = coef(model)[2]
  
  # Calculate the stability score
  stability_score = 1 / (1 + abs(stability_coefficient))
  

  return(stability_score)
}


external_light = rep(c(0.4, 0.6, 0.8), each = 100)
PSI = calculate_PSI(synthetic_data1, trait_col = 5, external_env_factor = external_light)
print(PSI)

##################################







