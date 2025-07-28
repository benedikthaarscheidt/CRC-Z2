#Calculate Developmental Plasticity Index (DPI)(calculate_DPI) - tested,
#Calculate Coefficient of Environmental Variation (CEV)(calculate_CEV) - tested,
#Calculate Plasticity Response Index (PRI)(calculate_PRI) - tested,
#Calculate Phenotypic Flexibility Index (PFI)(calculate_PFI) - tested,
#Calculate Standardized Plasticity Index (SPI)(calculate_SPI) - tested,
#Calculate Absolute Plasticity Coefficient (APC)(calculate_APC) - tested,
#Calculate Stability Index (SI)(calculate_SI) - tested,
#Calculate Relative Stability Index (RSI)(calculate_RSI),
#Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
#Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
#Calculate Multivariate Plasticity Index (MVPi)(calculate_MVPi) NEEDS TO BE TESTED WITH THE REQUESTED DATASET FROM PROF. BARBOSA,
#Calculate Standardized Plasticity Metric (SPM)(calculate_SPM) - tested,
#Calculate SSpop/SStotal Plasticity Ratio(calculate_Plasticity_Ratio) - tested


source("~/CRC_1644_Z2/R-files/Plasticity_scores/2.R")

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

##################### datasets for testing

df_test1 = data.frame(Column1 = c(rep(4, 10), rep(2, 10)), 
                      Column2 = c(rep(10, 10), rep(1, 10)))

df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),
                      Column1 = c(rep(2, 10), rep(1, 10)), 
                      Column2 = c(rep(2, 10), rep(4, 10)))

df_test2.2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),
                      Column1 = c(rep(2, 10), rep(1, 5),rep(2,5)), 
                      Column2 = c(rep(2, 10), rep(4, 10)))



df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),
                    Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), 
                    Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))

df_test4 = data.frame(
  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   # Response variable
  Column1 = c(rep(2, 10), rep(3, 10), rep(1, 10)),   # Control (2) and Treatment (3)
  Column2 = c(rep(1, 10), rep(1, 10), rep(1, 10))    # Covariate (matches values of Column0)
)

df_test5 = data.frame(
  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   # Response variable
  Column1 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   # Control (2) and Treatment (3)
  Column2 = c(rep(3, 10), rep(2, 10), rep(1, 10))    # Covariate (matches values of Column0)
)

df_test_simple = data.frame(
  Column1 = c(2, 2, 3, 3),      # Control (2) and Treatment (3)
  Column0 = c(10, 12, 20, 22),  # Response variable (trait)
  Column2 = c(1, 1, 2, 2)       # Covariate (e.g., biomass)
)



#########################################

#' @title Calculate Developmental Plasticity Index (DPI)
#' @description This function calculates the Developmental Plasticity Index (DPI) by quantifying how a trait value 
#' changes between two time points during development.
#'
#' The formula used is:
#' \deqn{DPI = \frac{Trait Value at Time2 - Trait Value at Time1}{Time Interval}}
#'
#' @param trait_values_time1 A numeric vector containing trait values at Time 1.
#' @param trait_values_time2 A numeric vector containing trait values at Time 2.
#' @param time_interval Optional numeric vector of the time intervals between the two time points for each sample.
#' If not provided, a vector of ones is assumed.
#'
#' @return A numeric vector of DPI values.
#' 
#' @examples
#' # Example usage with synthetic data and explicit time intervals
#' trait_time1 = c(10, 12, 14)
#' trait_time2 = c(15, 18, 22)
#' time_intervals = c(5, 5, 5)
#' dpi_values = calculate_DPI(trait_time1, trait_time2, time_intervals)
#' print(dpi_values)
#'
#' # Example usage with default time interval (all ones)
#' dpi_values_default = calculate_DPI(trait_time1, trait_time2)
#' print(dpi_values_default)
#'
#' @export
calculate_DPI = function(trait_values_time1, trait_values_time2, time_interval = NULL) {
  # Ensure trait values are numeric vectors
  if (!is.numeric(trait_values_time1) || !is.numeric(trait_values_time2)) {
    stop("trait_values_time1 and trait_values_time2 must be numeric vectors.")
  }
  
  # If no time_interval vector is provided, assume a vector of ones
  if (is.null(time_interval)) {
    time_interval = rep(1, length(trait_values_time1))
  } else if (!is.numeric(time_interval)) {
    stop("time_interval must be a numeric vector.")
  }
  
  # Ensure all vectors have the same length
  if (length(trait_values_time1) != length(trait_values_time2) || length(trait_values_time1) != length(time_interval)) {
    stop("trait_values_time1, trait_values_time2, and time_interval must all have the same length.")
  }
  
  # Ensure all elements of time_interval are positive numeric values
  if (any(time_interval <= 0)) {
    stop("All elements of time_interval must be positive numeric values.")
  }
  
  # Calculate DPI for each sample
  dpi = (trait_values_time2 - trait_values_time1) / time_interval
  
  return(dpi)
}




## test - passed on synthetic dataset - however it is not clear how one would structure a dataset with time resolved data (will the measurements at time x for one sample be in different columns all with equally separated time intervals or in rows?)
# mathematically the function does the wanted if one knows how to structure the data

#calculate_DPI(df_test2,2,3,3)



###############################



#' Calculate Coefficient of Environmental Variation (CEV) for a Single Trait
#'
#' This function calculates the Coefficient of Environmental Variation (CEV) for a single trait,
#' given a numeric vector of trait values measured across environments. If the environmental values
#' are not provided, they are assumed to be equidistant.
#'
#' The formula used is:
#' \deqn{CEV = \frac{Standard\ Deviation}{Mean} \times 100}
#'
#' @param trait_values A numeric vector of trait values for a single trait (each value measured in an
#' environment, assumed equidistant).
#'
#' @return A list containing:
#'   - `CEV`: The calculated Coefficient of Environmental Variation.
#'   - `Mean`: The mean of the trait values.
#'   - `SD`: The standard deviation of the trait values.
#'   - `Valid`: Logical indicating whether the calculation was successful.
#'
#' @examples
#' trait_values = c(10, 15, 20, 25, 30)
#' result = calculate_CEV(trait_values)
#' print(result)
#'
#' @export
calculate_CEV = function(trait_values) {
  # Ensure input is a numeric vector
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  # Check that there are enough data points to compute a standard deviation
  if (length(trait_values) < 2) {
    return(list(CEV = NA, Mean = NA, SD = NA, Valid = FALSE))
  }
  
  # Calculate mean and standard deviation (ignoring NA values)
  mean_val = mean(trait_values, na.rm = TRUE)
  sd_val = sd(trait_values, na.rm = TRUE)
  
  # Avoid division by zero if the mean is zero
  if (mean_val == 0) {
    return(list(CEV = NA, Mean = mean_val, SD = sd_val, Valid = FALSE))
  }
  
  # Calculate the Coefficient of Environmental Variation (CEV)
  cev = (sd_val / mean_val) * 100
  
  return(cev)
}



## tested - passed on synthetic dataset 

#calculate_CEV(df_test2,c(2,3))
#sd(df_test2[,3])/mean(df_test2[,3])*100

################################


#' Calculate Plasticity Response Index (PRI)
#'
#' This function calculates the Plasticity Response Index (PRI) for a single trait. It requires a
#' numeric vector of trait values and a corresponding binary vector of the same length indicating
#' the environment for each measurement (0 for control, 1 for extreme). The PRI is calculated as:
#'
#' \deqn{PRI = \frac{\text{Mean(Extreme)} - \text{Mean(Control)}}{\text{Mean(Overall)}}}
#'
#' @param trait_values A numeric vector of trait values.
#' @param env_indicator A binary vector of the same length as trait_values (0 for control, 1 for extreme).
#'
#' @return A numeric value representing the Plasticity Response Index (PRI).
#'
#' @examples
#' trait_values <- c(10, 15, 20, 25, 30)
#' env_indicator <- c(0, 0, 1, 1, 1)
#' pri_value <- calculate_PRI(trait_values, env_indicator)
#' print(pri_value)
#'
#' @export
calculate_PRI <- function(trait_values, env_indicator) {
  # Validate input
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  if (length(trait_values) != length(env_indicator)) {
    stop("trait_values and env_indicator must have the same length.")
  }
  if (!all(env_indicator %in% c(0, 1))) {
    stop("env_indicator must be a binary vector with values 0 (control) and 1 (extreme).")
  }
  
  # Calculate overall mean
  overall_mean <- mean(trait_values, na.rm = TRUE)
  if (overall_mean == 0) {
    warning("Overall mean of trait values is 0; PRI is undefined. Returning NA.")
    return(NA)
  }
  
  # Calculate means for extreme and control environments
  mean_extreme <- mean(trait_values[env_indicator == 1], na.rm = TRUE)
  mean_control <- mean(trait_values[env_indicator == 0], na.rm = TRUE)
  
  # Calculate PRI
  pri_value <- (mean_extreme - mean_control) / overall_mean
  
  return(pri_value)
}


## test - passed on synthetic dataset

#calculate_PRI(df_test2,1,2,1,2)

############################

#' Calculate Phenotypic Flexibility Index (PFI) for a Single Trait
#'
#' This function calculates the Phenotypic Flexibility Index (PFI) for a single trait by finding
#' the sample with the maximum absolute deviation from its baseline value, then dividing this deviation
#' by that sample's baseline value. If a baseline is provided, it is used; otherwise, the overall mean
#' of the trait values is used as the baseline for all samples.
#'
#' @param trait_values A numeric vector of trait measurements.
#' @param baseline_values Optional numeric vector of baseline values of the same length as `trait_values`.
#' If not provided, the overall mean of `trait_values` is used as the baseline for every sample.
#'
#' @return A numeric value representing the Phenotypic Flexibility Index (PFI).
#'
#' @examples
#' # Example with a provided baseline:
#' trait_values <- c(18, 22, 25, 20, 19)
#' baseline_values <- c(20, 20, 20, 20, 20)
#' pfi <- calculate_PFI(trait_values, baseline_values)
#' print(pfi)
#'
#' # Example without a provided baseline (using the mean as baseline):
#' trait_values <- c(18, 22, 25, 20, 19)
#' pfi <- calculate_PFI(trait_values)
#' print(pfi)
#'
#' @export
calculate_PFI <- function(trait_values, baseline_values = NULL) {
  # Validate input for trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  n <- length(trait_values)
  
  # If baseline_values is not provided, use the overall mean for all samples.
  if (is.null(baseline_values)) {
    baseline_values <- rep(mean(trait_values, na.rm = TRUE), n)
  } else {
    # Validate baseline_values
    if (!is.numeric(baseline_values)) {
      stop("baseline_values must be a numeric vector.")
    }
    if (length(baseline_values) != n) {
      stop("baseline_values must have the same length as trait_values.")
    }
  }
  
  # Calculate the absolute deviations from the baseline for each sample
  deviations <- abs(trait_values - baseline_values)
  
  # Identify the sample with the maximum deviation
  max_deviation_index <- which.max(deviations)
  max_deviation <- deviations[max_deviation_index]
  
  # Retrieve the baseline corresponding to the maximum deviation sample
  baseline_for_max <- baseline_values[max_deviation_index]
  
  if (baseline_for_max == 0) {
    warning("The baseline value corresponding to the maximum deviation is 0; PFI is undefined. Returning NA.")
    return(NA)
  }
  
  # Calculate the PFI as the ratio of the maximum deviation to its corresponding baseline
  pfi_value <- max_deviation / baseline_for_max
  
  return(pfi_value)
}

## test - passed on synthetic dataset - however I am unsure wether this is supposed to be one score or a separate score for each timepoint/sample

#calculate_PFI(df_test2,baseline_cols=c(1,1),trait_cols=c(2,3))


##########################

#' Calculate Standardized Plasticity Index (SPI)
#'
#' This function calculates the Standardized Plasticity Index (SPI) for multiple traits by comparing trait values 
#' across environments relative to a reference environment.
#'
#' The formula used is:
#' \deqn{SPI = \frac{Trait Value Env2 - Trait Value Env1}{Standard Deviation of Reference Trait Values}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A vector of strings or numeric values indicating the columns in the dataset containing the trait values.
#' @param trait_env1 A numeric or string value representing the identifier for environment 1 in `env_col`.
#' @param trait_env2 A numeric or string value representing the identifier for environment 2 in `env_col`.
#' @param reference_env A numeric or string value representing the identifier for the reference environment in `env_col`.
#' @return A named vector containing the SPI for each trait.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1_Values = c(10, 12, 15, 18),
#'   Trait2_Values = c(5, 6, 7, 8),
#'   Environment = factor(c(1, 1, 2, 2))  # Env1 = 1, Env2 = 2, Reference = 1
#' )
#' spi_values = calculate_SPI(synthetic_data, env_col = "Environment", trait_cols = c("Trait1_Values", "Trait2_Values"), trait_env1 = 1, trait_env2 = 2, reference_env = 1)
#' print(spi_values)
#'
#' @export
calculate_SPI = function(data, env_col, trait_cols, env1, env2, reference_env) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a list to store SPI values for each trait
  spi_results = numeric(length(trait_cols))
  names(spi_results) = trait_cols
  
  i=0
  # Loop through each trait column and calculate SPI
  for (trait_col in trait_cols) {
    i=i+1
    # Extract trait data
    trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Calculate mean trait values for env1 and env2
    trait_env1_value = mean(trait_values[env_col == env1], na.rm = TRUE)
    trait_env2_value = mean(trait_values[env_col == env2], na.rm = TRUE)
    
    # Calculate the standard deviation of the reference environment's trait values
    sd_reference = sd(trait_values[env_col == reference_env], na.rm = TRUE)
    
    # Calculate the Standardized Plasticity Index (SPI) for this trait
    spi = (trait_env2_value - trait_env1_value) / sd_reference
    
    # Store the result
    spi_results[i] = spi
  }
  
  return(spi_results)
}

## test - passed on synthetic dataset 

#calculate_SPI(df_test2.2,1,2,1,2,2)
#-0.5/sd(c(rep(1, 5),rep(2,5)))




########################

#' Calculate Absolute Plasticity Coefficient (APC) for a Single Trait
#'
#' This function calculates the Absolute Plasticity Coefficient (APC) for a single trait by
#' computing the average absolute difference between the mean trait values of consecutive environments.
#' The trait measurements are grouped by the provided environment identifiers. If no environment vector
#' is given, the function assumes that the environments are equidistant and uses \code{seq_along(trait_values)}.
#'
#' The formula used is:
#' \deqn{APC = \frac{1}{n-1} \sum_{i=1}^{n-1} \left| \mu_{i+1} - \mu_i \right|}
#'
#' where \eqn{\mu_i} is the mean trait value in environment \eqn{i} and \eqn{n} is the number of unique environments.
#'
#' @param trait_values A numeric vector of trait measurements.
#' @param env_labels Optional vector of environment identifiers (numeric, character, or factor).
#'                   Must be the same length as \code{trait_values}. If not provided, equidistant environments
#'                   are assumed.
#' @param sequential_env Logical indicating whether the environments should be treated as sequentially ordered.
#'                       Defaults to \code{TRUE}. If \code{TRUE}, the function orders the environments in increasing order.
#'
#' @return A numeric value representing the Absolute Plasticity Coefficient (APC).
#'
#' @examples
#' # Example with explicit environment labels:
#' trait_values <- c(10, 12, 15, 20, 22, 25)
#' env_labels <- c(1, 1, 2, 2, 3, 3)
#' apc_value <- calculate_APC(trait_values, env_labels, sequential_env = TRUE)
#' print(apc_value)
#'
#' # Example without providing env_labels (assumes equidistant environments):
#' trait_values <- c(10, 12, 15, 20, 22, 25)
#' apc_value <- calculate_APC(trait_values)
#' print(apc_value)
#'
#' @export
calculate_APC <- function(trait_values, env_labels = NULL, sequential_env = TRUE) {
  # Input validation for trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  # If no environment labels are provided, assume equidistant environments.
  if (is.null(env_labels)) {
    env_labels <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env_labels)) {
      stop("trait_values and env_labels must have the same length.")
    }
  }
  
  # Order the data based on the environment labels.
  if (sequential_env) {
    if (is.numeric(env_labels)) {
      ordering <- order(env_labels)
    } else if (is.factor(env_labels)) {
      if (is.ordered(env_labels)) {
        ordering <- order(as.numeric(env_labels))
      } else {
        ordering <- seq_along(env_labels)
      }
    } else {  # For character vectors
      ordering <- order(env_labels)
    }
  } else {
    # Use the order of first appearance
    uniq_env <- unique(env_labels)
    ordering <- order(match(env_labels, uniq_env))
  }
  
  trait_values <- trait_values[ordering]
  env_labels <- env_labels[ordering]
  
  # Calculate mean trait value for each unique environment
  env_means <- tapply(trait_values, env_labels, mean, na.rm = TRUE)
  if (length(env_means) < 2) {
    warning("Less than two unique environments found; APC is undefined. Returning NA.")
    return(NA)
  }
  
  # Compute absolute differences between consecutive environment means
  abs_diff <- abs(diff(env_means))
  
  # Calculate APC as the mean of these absolute differences
  apc_value <- mean(abs_diff, na.rm = TRUE)
  
  return(apc_value)
}

## test - passed on synthetic dataset 
#calculate_APC(df_test2,1,2)

#############################


#' Calculate Stability Index (SI) for a Single Trait
#'
#' This function calculates the Stability Index (SI) for a single trait, quantifying the stability
#' of trait expression across different environments. It groups the trait measurements by the provided
#' environment vector, computes the mean trait value for each environment, and then calculates SI as the
#' ratio of the variance among these environment means to their overall mean.
#'
#' The formula used is:
#' \deqn{SI = \frac{\text{Variance of Environment Means}}{\text{Mean of Environment Means}}}
#'
#' @param trait_values A numeric vector of trait measurements.
#' @param env Optional vector of environments (numeric, character, or factor) corresponding to each measurement.
#'            If not provided, the function assumes equidistant environments using \code{seq_along(trait_values)}.
#'
#' @return A numeric value representing the Stability Index (SI).
#'
#' @examples
#' # Example with explicit environments:
#' trait_values <- c(10, 12, 15, 20)
#' env <- c(1, 2, 3, 4)
#' si_value <- calculate_SI(trait_values, env)
#' print(si_value)
#'
#' # Example without providing env (assumes equidistant environments):
#' trait_values <- c(10, 12, 15, 20)
#' si_value <- calculate_SI(trait_values)
#' print(si_value)
#'
#' @export
calculate_SI <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }
  
  # Compute the mean trait value for each unique environment.
  env_means <- tapply(trait_values, env, mean, na.rm = TRUE)
  
  # Calculate the variance among environment means and their overall mean.
  variance_env <- var(env_means, na.rm = TRUE)
  overall_mean <- mean(env_means, na.rm = TRUE)
  
  if (overall_mean == 0) {
    warning("The overall mean of environment means is 0; SI is undefined. Returning NA.")
    return(NA)
  }
  
  # Compute the Stability Index (SI)
  SI_value <- variance_env / overall_mean
  
  return(SI_value)
}


### test - passed on synthetic dataset 
#var(df_test2[,2])/mean(df_test2[,2])
#calculate_SI(df_test2,1,2)
#######################


#' Calculate Relative Stability Index (RSI) for a Single Trait
#'
#' This function calculates the Relative Stability Index (RSI) for a single trait, measuring the relative
#' stability of trait expression across different environments. It groups the trait measurements by the provided
#' environment vector, computes the mean trait value for each environment, and then calculates RSI as:
#'
#' \deqn{RSI = 1 - \frac{SD(\text{Environment Means})}{Mean(\text{Environment Means})}}
#'
#' If no environment vector is provided, the function assumes equidistant environments using \code{seq_along(trait_values)}.
#'
#' @param trait_values A numeric vector of trait measurements.
#' @param env Optional vector of environments (numeric, character, or factor) corresponding to each measurement.
#'            If not provided, the function assumes equidistant environments.
#'
#' @return A numeric value representing the Relative Stability Index (RSI).
#'
#' @examples
#' # Example with explicit environments:
#' trait_values <- c(10, 12, 15, 20)
#' env <- c(1, 2, 3, 4)
#' rsi_value <- calculate_RSI(trait_values, env)
#' print(rsi_value)
#'
#' # Example without providing env (assumes equidistant environments):
#' trait_values <- c(10, 12, 15, 20)
#' rsi_value <- calculate_RSI(trait_values)
#' print(rsi_value)
#'
#' @export
calculate_RSI <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }
  
  # Compute the mean trait value for each unique environment.
  env_means <- tapply(trait_values, env, mean, na.rm = TRUE)
  
  # Calculate the standard deviation and mean of the environment means.
  overall_sd <- sd(env_means, na.rm = TRUE)
  overall_mean <- mean(env_means, na.rm = TRUE)
  
  if (overall_mean == 0) {
    warning("The overall mean of environment means is 0; RSI is undefined. Returning NA.")
    return(NA)
  }
  
  # Compute the Relative Stability Index (RSI)
  RSI_value <- 1 - (overall_sd / overall_mean)
  
  return(RSI_value)
}


### test - passed on synthetic dataset
#1-sd(df_test2[,2])/mean(df_test2[,2])
#calculate_RSI(df_test2,1,2)
#
###test2
#1-sd(synthetic_data1[,3])/mean(synthetic_data1[,3])
#calculate_RSI(synthetic_data1,1,3)
#
###test3
#1-sd(df_test5[1:20,2])/mean(df_test5[1:20,2])
#calculate_RSI(df_test5,1,2,specific_envs=c(3,2))



##########################

#' Calculate Environmental Variance Sensitivity (EVS) for a Single Trait
#'
#' This function calculates the Environmental Variance Sensitivity (EVS) for a single trait,
#' quantifying how sensitive the trait is to changes in the environment. EVS is computed as:
#'
#' \deqn{EVS = \frac{Var(Trait\ Values)}{Var(Environment)}}
#'
#' If no environment vector is provided, the function assumes equidistant environments using
#' \code{seq_along(trait_values)}.
#'
#' @param trait_values A numeric vector of trait measurements.
#' @param env Optional vector of environments (numeric, character, or factor) corresponding to each measurement.
#'            If not provided, equidistant environments are assumed.
#'
#' @return A numeric value representing the Environmental Variance Sensitivity (EVS).
#'
#' @examples
#' # Example with explicit environments:
#' trait_values <- c(10, 12, 15, 18)
#' env <- c(1, 2, 3, 4)
#' evs_value <- calculate_EVS(trait_values, env)
#' print(evs_value)
#'
#' # Example without providing env (assumes equidistant environments):
#' trait_values <- c(10, 12, 15, 18)
#' evs_value <- calculate_EVS(trait_values)
#' print(evs_value)
#'
#' @export
calculate_EVS <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }
  
  # Compute the variance in trait values.
  trait_variance <- var(trait_values, na.rm = TRUE)
  
  # If env is not numeric, convert it to numeric via as.factor.
  if (!is.numeric(env)) {
    env_numeric <- as.numeric(as.factor(env))
  } else {
    env_numeric <- env
  }
  
  # Compute the variance in the environment.
  env_variance <- var(env_numeric, na.rm = TRUE)
  
  if (env_variance == 0) {
    warning("Environmental variance is zero; EVS is undefined. Returning NA.")
    return(NA)
  }
  
  # Calculate EVS as the ratio of the two variances.
  EVS_value <- trait_variance / env_variance
  
  return(EVS_value)
}

## test - passed on synthetic dataset

#var(df_test2[,2])/var(df_test2[,1])
#calculate_EVS(df_test2,1,2)
#
###test2 - passed 
#var(synthetic_data1[,3])/var(synthetic_data1[,2])
#calculate_EVS(synthetic_data1,2,3)

#####################

#' Calculate Multivariate Plasticity Index (MVPi)
#'
#' This function calculates the Multivariate Plasticity Index (MVPi) for a set of phenotypic traits,
#' quantifying plasticity across environments by computing the average Euclidean distance between
#' the centroids of environments in a reduced-dimensional space. Dimensionality is reduced using PCA
#' via \code{prcomp}, retaining the first \code{n_axes} axes.
#'
#' The formula is based on the average Euclidean distance between environment centroids:
#'
#' \deqn{MVPi = \text{mean}\Big\{ \sqrt{\sum_{k=1}^{n\_axes} \big(c_{ik} - c_{jk}\big)^2} \Big\}_{i < j}}
#'
#' @param trait_data A data structure (e.g., data frame, matrix, or numeric vector) containing phenotypic trait values.
#'                   Rows correspond to samples.
#' @param env Optional vector of environments (numeric, character, or factor) corresponding to each sample.
#'            If not provided, the function assumes equidistant environments using \code{seq_len(nrow)}.
#' @param n_axes Integer specifying the number of PCA axes to retain.
#'
#' @return A numeric value representing the Multivariate Plasticity Index (MVPi).
#'
#' @examples
#' # Example with explicit environments:
#' trait_data <- data.frame(
#'   Trait1 = c(1.2, 2.3, 1.9, 3.1, 2.0, 2.8),
#'   Trait2 = c(2.5, 3.0, 2.8, 3.6, 2.9, 3.3)
#' )
#' env <- factor(c(1, 1, 2, 2, 3, 3))
#' mvpi_value <- calculate_MVPi(trait_data, env, n_axes = 2)
#' print(mvpi_value)
#'
#' # Example without providing an environment vector (assumes equidistant environments):
#' mvpi_value <- calculate_MVPi(trait_data, n_axes = 2)
#' print(mvpi_value)
#'
#' @export
calculate_MVPi <- function(trait_data, env = NULL, n_axes = 1) {
  # Convert trait_data to a numeric matrix if not already one.
  if (!is.matrix(trait_data)) {
    trait_data <- as.matrix(trait_data)
  }
  if (!is.numeric(trait_data)) {
    stop("trait_data must be numeric.")
  }
  
  # Determine number of samples.
  n_samples <- nrow(trait_data)
  
  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_len(n_samples)
  } else {
    if (length(env) != n_samples) {
      stop("The length of env must equal the number of rows in trait_data.")
    }
  }
  
  # Convert env to a factor.
  env <- factor(env)
  
  # Perform PCA on the trait data.
  pca_result <- prcomp(trait_data, center = TRUE, scale. = FALSE)
  
  # Ensure that n_axes does not exceed the available dimensions.
  if (n_axes > ncol(pca_result$x)) {
    stop("n_axes exceeds the number of available principal components.")
  }
  
  # Retain the first n_axes of the PCA scores.
  scores <- pca_result$x[, 1:n_axes, drop = FALSE]
  
  # Compute centroids for each environment.
  env_levels <- levels(env)
  if (length(env_levels) < 2) {
    warning("Less than two environments found; MVPi is undefined. Returning NA.")
    return(NA)
  }
  centroids <- sapply(env_levels, function(x) {
    colMeans(scores[env == x, , drop = FALSE], na.rm = TRUE)
  })
  centroids <- t(centroids)  # rows = environment centroids
  
  # Compute pairwise Euclidean distances between centroids.
  distances <- as.vector(dist(centroids))
  
  # Calculate MVPi as the mean Euclidean distance.
  mvpi <- mean(distances, na.rm = TRUE)
  
  return(mvpi)
}



## test - to be tested with dataset from DOI: 10.1093/jxb/eraa545


#calculate_MVPi(synthetic_data1,1,c(3,4,5),3)


#######################

#' Calculate Standardized Plasticity Metric (SPM)
#'
#' This function calculates the Standardized Plasticity Metric (SPM), which quantifies phenotypic plasticity
#' by comparing trait values between resident and nonresident environments.
#'
#' The formula used is:
#' \deqn{SPM = \frac{|AinA - AinB|}{AinA} \text{ or } \frac{|BinB - BinA|}{BinB}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_col A string or numeric value indicating the column in the dataset containing the trait values.
#' @param resident_env A numeric or string value representing the identifier for the resident environment in `env_col`.
#' @param nonresident_env A numeric or string value representing the identifier for the nonresident environment in `env_col`.
#' @return A numeric value representing the Standardized Plasticity Metric (SPM).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_Values = c(10, 12, 15, 18),
#'   Environment = factor(c(1, 1, 2, 2))  # Resident = 1, Nonresident = 2
#' )
#' spm_value = calculate_SPM(synthetic_data, env_col = "Environment", trait_col = "Trait_Values", resident_env = 1, nonresident_env = 2)
#' print(spm_value)
#'
#' @export
calculate_SPM = function(data, env_col, trait_col, resident_env, nonresident_env) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Extract trait data
  trait_values = data[[trait_col]]
  
  # Calculate mean trait values for resident and nonresident environments
  trait_resident = mean(trait_values[env_col == resident_env], na.rm = TRUE)
  trait_nonresident = mean(trait_values[env_col == nonresident_env], na.rm = TRUE)
  
  # Calculate the Standardized Plasticity Metric (SPM)
  spm = abs(trait_resident - trait_nonresident) / trait_resident
  
  return(spm)
}

#
### test - passed on synthetic dataset 
#(mean(df_test2[df_test2[, 1] == 1, 2])-mean(df_test2[df_test2[, 1] == 2, 2]))/mean(df_test2[df_test2[, 1] == 1, 2])
#calculate_SPM(df_test2,1,2,1,2)
#
###test2
#(mean(synthetic_data1[synthetic_data1[, 1] == 1, 3])-mean(synthetic_data1[synthetic_data1[, 1] == 3, 3]))/mean(synthetic_data1[synthetic_data1[, 1] == 1, 3])
#calculate_SPM(synthetic_data1,1,3,1,3)




############################


#' Calculate SSpop/SStotal Plasticity Ratio for Multiple Traits
#'
#' This function calculates the Plasticity Ratio, which quantifies the relative contribution 
#' of population differences to the overall phenotypic variability across environments, 
#' for multiple trait columns. It is calculated as the ratio of the sum of squares between populations (SSpop) 
#' to the total sum of squares (SStotal) from a one-way ANOVA.
#'
#' The formula used is:
#' \deqn{Plasticity Ratio = \frac{SSpop}{SStotal}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string or numeric value indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A vector of strings or numeric values indicating the columns in the dataset containing the trait values.
#' @param pop_col A string or numeric value indicating the column in the dataset containing the population labels, or an external vector specifying population belonging.
#' @return A named numeric vector representing the SSpop/SStotal Plasticity Ratio for each trait.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_1 = c(10, 12, 15, 18, 11, 17),
#'   Trait_2 = c(20, 22, 25, 28, 21, 27),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3)),
#'   Population = factor(c("A", "A", "B", "B", "C", "C"))
#' )
#' plasticity_ratios = calculate_Plasticity_Ratio(synthetic_data, env_col = "Environment", trait_cols = c("Trait_1", "Trait_2"), pop_col = "Population")
#' print(plasticity_ratios)
#'
#' @export
calculate_Plasticity_Ratio = function(data, env_col, trait_cols, pop_col) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Handle pop_col
  if (is.numeric(pop_col) && length(pop_col) == 1) {
    pop_col = data[[pop_col]]
  } else if (is.vector(pop_col) && length(pop_col) == nrow(data)) {
    pop_col = pop_col
  } else {
    pop_col = data[[pop_col]]
  }
  
  # Initialize a vector to store Plasticity Ratios for each trait
  plasticity_ratios = numeric(length(trait_cols))
  names(plasticity_ratios) = trait_cols
  
  # Loop over each trait column
  for (i in seq_along(trait_cols)) {
    # Extract the trait values for the current trait
    trait_values = if (is.numeric(trait_cols[i])) data[[trait_cols[i]]] else data[[trait_cols[i]]]
    
    # Create a temporary data frame for the ANOVA
    temp_data = data.frame(trait_values, pop_col, env_col)
    
    # Perform a one-way ANOVA to obtain SSpop and SStotal
    anova_result = aov(trait_values ~ pop_col + env_col, data = temp_data)
    
    # Extract the sum of squares between populations (SSpop) and total sum of squares (SStotal)
    ss_total = sum(anova(anova_result)[, "Sum Sq"])
    ss_pop = anova(anova_result)["pop_col", "Sum Sq"]
    
    # Calculate and store the Plasticity Ratio
    plasticity_ratios[i] = ss_pop / ss_total
  }
  
  return(plasticity_ratios)
}


#pop=c(rep(1,150),rep(2,150))
#calculate_Plasticity_Ratio(synthetic_data1,1,3,pop)
#
### test - passed on synthetic dataset 
#
#overall_mean = mean(synthetic_data1[, 3])
#
#mean_pop1 = mean(synthetic_data1[1:150, 3])   # Mean for Population 1 (first 150 rows)
#mean_pop2 = mean(synthetic_data1[151:300, 3])  # Mean for Population 2 (next 150 rows)
#
#n_pop1 = 150  # Number of observations in Population 1
#n_pop2 = 150  # Number of observations in Population 2
#
#ss_pop1 = n_pop1 * (mean_pop1 - overall_mean)^2
#ss_pop2 = n_pop2 * (mean_pop2 - overall_mean)^2
#SS_pop = ss_pop1 + ss_pop2  # Total SS_pop
#
#SS_total = sum((synthetic_data1[, 3] - overall_mean)^2)
#
#plasticity_ratio = SS_pop / SS_total
#print(paste("Plasticity Ratio:", plasticity_ratio))



