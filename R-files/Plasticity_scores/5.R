#this files contains the following indices/functions: 
#Coefficient sum (calculate_plasticity),
#Relative Plasticity Index (calculate_RPI) - tested,
#Plasticity Quotient (calculate_PQ) - tested,
#Phenotypic Range (PR) (calculate_PR) - tested,
#Norm of reaction width (calculate_NRW) - tested,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
#Calculate Plasticity Differential (PD) (calculate_PD) - tested,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,



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



############################################


#' Calculate the Plasticity Score from a Polynomial Reaction Norm
#'
#' This function fits polynomial regression models of varying degrees to describe 
#' the relationship between a phenotype and an environmental gradient.
#' The best-fitting polynomial is selected using both AIC and BIC criteria, with 
#' a preference for lower-degree models when performance differences are small.
#'
#' The plasticity score is calculated as the sum of the absolute values of the polynomial coefficients.
#'
#' @param trait_values A numeric vector representing phenotype values.
#' @param env_values (Optional) A numeric vector representing environmental values. 
#' If NULL, equidistant values will be generated.
#' @param max_degree The maximum polynomial degree to consider (default is 3).
#' @param criterion Selection criterion for model complexity. Options are "AIC" (default) or "BIC".
#' @return A list containing:
#'   - `best_degree`: The best-fitting polynomial degree.
#'   - `plasticity_score`: The sum of absolute values of the polynomial coefficients.
#'   - `coefficients`: The estimated polynomial coefficients for the best model.
#'
#' @examples
#' trait_values = c(100, 110, 120, 105, 115, 125)
#' env_values = c(1, 2, 3, 1, 2, 3)
#'
#' # Calculate Plasticity Score with given environmental values
#' plasticity_result = calculate_plasticity(trait_values, env_values)
#' print(plasticity_result)
#'
#' # Calculate Plasticity Score assuming equidistant environments
#' plasticity_result_no_env = calculate_plasticity(trait_values)
#' print(plasticity_result_no_env)
#'@references
#' GENETICS AND EVOLUTION OF PHENOTYPIC PLASTICITY - Samuel M. Scheiner
#' @export
calculate_plasticity = function(trait_values, env_values = NULL, max_degree = 3, criterion = "BIC") {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }
  
  n_values = length(trait_values)
  
  # If no environmental values are provided, create equidistant values
  if (is.null(env_values)) {
    env_values = seq_len(n_values)
  }
  
  # Ensure environmental values are numeric
  if (!is.numeric(env_values)) {
    stop("env_values must be numeric")
  }
  
  # Check for length mismatch
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }
  
  # Store best model information
  best_degree = 1
  best_criterion_value = Inf
  best_model = NULL
  
  # Try polynomial degrees from 1 to max_degree
  for (degree in 1:max_degree) {
    formula = as.formula(paste("trait_values ~ poly(env_values, ", degree, ", raw = TRUE)", sep = ""))
    model = lm(formula)
    
    # Use AIC or BIC for model selection
    model_criterion = ifelse(criterion == "AIC", AIC(model), BIC(model))
    
    # Penalize higher-degree models if difference is small (<2)
    if ((model_criterion) < best_criterion_value) {
      best_criterion_value = model_criterion
      best_degree = degree
      best_model = model
    }
  }
  
  
  coefficients = coef(best_model)[-1]
  
  # Compute plasticity score as sum of absolute values of coefficients
  plasticity_score = sum(abs(coefficients))
  
  return(list(
    best_degree = best_degree,
    plasticity_score = plasticity_score,
    coefficients = coefficients
  ))
}

##############################################



#' Calculate Cross-Environment Covariance and Correlation for Linear Reaction Norms
#'
#' This function computes the covariance (and optionally correlation) of trait values across different environments.
#' If no environment vector is provided, it assumes equidistant environments.
#'
#' @param trait_values A numeric vector representing trait values measured across different environments.
#' @param env_values (Optional) A numeric vector specifying the environment associated with each trait measurement.
#'                   If NULL, environments are assumed to be equidistant.
#' @param return_correlation Logical; if TRUE, returns both covariance and correlation (default: FALSE).
#'
#' @return A named list containing:
#'   - `covariance`: The covariance between trait values across environments.
#'   - `correlation` (if `return_correlation = TRUE`): The Pearson correlation coefficient.
#'
#' @examples
#' trait_values = c(10, 12, 15, 17, 20, 22, 25, 27, 30, 32)
#'
#' # Calculate covariance with equidistant environments
#' result = cross_env_cov(trait_values)
#' print(result)
#'
#' # Provide an explicit environment vector
#' env_values = seq_along(trait_values)  # Example of environmental values
#' result_with_env = cross_env_cov(trait_values, env_values, return_correlation = TRUE)
#' print(result_with_env)
#'
#' @export
cross_env_cov = function(trait_values, env_values = NULL, return_correlation = FALSE) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }
  
  n = length(trait_values)
  
  # If no env_values are provided, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)  # Equidistant environments
  }
  
  # Ensure correct length
  if (length(env_values) != n) {
    stop("trait_values and env_values must have the same length")
  }
  
  # Calculate covariance
  covariance = cov(trait_values, env_values)
  
  # Calculate correlation if requested
  if (return_correlation) {
    correlation = cor(trait_values, env_values)
    return(list(covariance = covariance, correlation = correlation))
  } else {
    return(list(covariance = covariance))
  }
}

