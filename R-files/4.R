#Calculate Developmental Plasticity Index (DPI)(calculate_DPI),
#Calculate Coefficient of Environmental Variation (CEV)(calculate_CEV),
#Calculate Plasticity Response Index (PRI)(calculate_PRI),
#Calculate Phenotypic Flexibility Index (PFI)(calculate_PFI),
#Calculate Standardized Plasticity Index (SPI)(calculate_SPI),
#Calculate Absolute Plasticity Coefficient (APC)(calculate_APC),
#Calculate Stability Index (SI)(calculate_SI),
#Calculate Relative Stability Index (RSI)(calculate_RSI),
#Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS),
#Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS),
#Calculate Multivariate Plasticity Index (MVPi)(calculate_MVPi),
#Calculate Standardized Plasticity Metric (SPM)(calculate_SPM),
#Calculate SSpop/SStotal Plasticity Ratio(calculate_Plasticity_Ratio)


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

##################### datasets for testing

df_test1 = data.frame(Column1 = c(rep(4, 10), rep(2, 10)), Column2 = c(rep(10, 10), rep(1, 10)))
df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 10), rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
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

#' Calculate Developmental Plasticity Index (DPI)
#'
#' This function calculates the Developmental Plasticity Index (DPI) by quantifying how a trait value 
#' changes between two time points during development.
#'
#' The formula used is:
#' \deqn{DPI = \frac{Trait Value at Time2 - Trait Value at Time1}{Time Interval}}
#'
#' @param data A data frame containing the input data, with columns for trait values at different time points.
#' @param trait_col_time1 A string or numeric value indicating the column for trait values at Time 1.
#' @param trait_col_time2 A string or numeric value indicating the column for trait values at Time 2.
#' @param time_interval Numeric, the time interval between the two time points.
#' @return A numeric value representing the Developmental Plasticity Index (DPI).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_Time1 = c(10, 12, 14),
#'   Trait_Time2 = c(15, 18, 22)
#' )
#' 
#' # Calculate DPI
#' dpi_value = calculate_DPI(synthetic_data, "Trait_Time1", "Trait_Time2", 5)
#' print(dpi_value)
#'
#' @export
calculate_DPI = function(data, trait_col_time1, trait_col_time2, time_interval) {
  
  # Handle the input columns for time1 and time2
  trait_time1 = if (is.numeric(trait_col_time1)) data[[trait_col_time1]] else data[[trait_col_time1]]
  trait_time2 = if (is.numeric(trait_col_time2)) data[[trait_col_time2]] else data[[trait_col_time2]]
  
  # Calculate DPI
  dpi = (trait_time2 - trait_time1) / time_interval
  
  return(dpi)
}

###############################



#' Calculate Coefficient of Environmental Variation (CEV)
#'
#' This function calculates the Coefficient of Environmental Variation (CEV), 
#' which measures the variability of trait values across different environments.
#' 
#' The formula used is:
#' \deqn{CEV = \frac{Standard Deviation of Trait Values Across Environments}{Mean Trait Value Across Environments} \times 100}
#'
#' @param data A data frame containing the input data.
#' @param trait_col A string or numeric value indicating the column for trait values across environments.
#' @return A numeric value representing the Coefficient of Environmental Variation (CEV).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_Values = c(10, 15, 20, 25, 30)
#' )
#' 
#' # Calculate CEV
#' cev_value = calculate_CEV(synthetic_data, "Trait_Values")
#' print(cev_value)
#'
#' @export
calculate_CEV = function(data, trait_col) {
  
  # Handle the input column for trait values
  trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Calculate the mean and standard deviation of the trait values
  mean_trait_value = mean(trait_values, na.rm = TRUE)
  sd_trait_value = sd(trait_values, na.rm = TRUE)
  
  # Calculate CEV
  cev = (sd_trait_value / mean_trait_value) * 100
  
  return(cev)
}

################################


#' Calculate Plasticity Response Index (PRI)
#'
#' This function calculates the Plasticity Response Index (PRI), which quantifies the plasticity
#' by measuring the proportional difference in trait values between an extreme environment 
#' (e.g., the best or worst condition) and a control environment.
#' 
#' The formula used is:
#' \deqn{PRI = \frac{Trait Value in Extreme Environment - Trait Value in Control Environment}{Mean Trait Value Across All Environments}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_col A string or numeric value indicating the column in the dataset containing the trait values.
#' @param trait_extreme A numeric or string value representing the identifier for the extreme environment in `env_col`.
#' @param trait_control A numeric or string value representing the identifier for the control environment in `env_col`.
#' @return A numeric value representing the Plasticity Response Index (PRI).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_Values = c(10, 15, 20, 25, 30),
#'   Environment = factor(c(1, 1, 2, 2, 2))  # Control = 1, Extreme = 2
#' )
#' pri_value = calculate_PRI(synthetic_data, env_col = "Environment", trait_col = "Trait_Values", trait_extreme = 2, trait_control = 1)
#' print(pri_value)
#'
#' @export
calculate_PRI = function(data, env_col, trait_col, trait_extreme, trait_control) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Extract trait data
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Calculate the mean trait value across all environments
  mean_trait_value = mean(trait_col, na.rm = TRUE)
  
  # Calculate trait values for extreme and control environments
  trait_ex = mean(trait_col[env_col == trait_extreme], na.rm = TRUE)
  trait_con = mean(trait_col[env_col == trait_control], na.rm = TRUE)
  
  # Calculate the Plasticity Response Index (PRI)
  pri = (trait_ex - trait_con) / mean_trait_value
  
  return(pri)
}


############################

#' Calculate Phenotypic Flexibility Index (PFI)
#'
#' This function calculates the Phenotypic Flexibility Index (PFI), which measures the degree of 
#' flexibility in trait values over time or environments.
#'
#' The formula used is:
#' \deqn{PFI = \frac{|Maximum Deviation from Baseline|}{Baseline Trait Value}}
#'
#' @param data A data frame containing the input data.
#' @param baseline_col A numeric value indicating the baseline for trait values.
#' @param trait_values_col A string or numeric value indicating the column for trait values across environments.
#' @return A numeric value representing the Phenotypic Flexibility Index (PFI).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Baseline = 20,
#'   Trait_Values = c(18, 22, 25, 20, 19)
#' )
#' pfi_value = calculate_PFI(synthetic_data, "Baseline", "Trait_Values")
#' print(pfi_value)
#'
#' @export
calculate_PFI = function(data, baseline_col, trait_values_col) {
  
  # Handle input columns
  trait_values = if (is.numeric(trait_values_col)) data[[trait_values_col]] else data[[trait_values_col]]
  
  # Calculate the maximum deviation from the baseline
  max_deviation = max(abs(trait_values - baseline), na.rm = TRUE)
  
  # Calculate PFI
  pfi = max_deviation / baseline
  
  return(pfi)
}


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
calculate_SPI = function(data, env_col, trait_cols, trait_env1, trait_env2, reference_env) {
  
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
    trait_env1_value = mean(trait_values[env_col == trait_env1], na.rm = TRUE)
    trait_env2_value = mean(trait_values[env_col == trait_env2], na.rm = TRUE)
    
    # Calculate the standard deviation of the reference environment's trait values
    sd_reference = sd(trait_values[env_col == reference_env], na.rm = TRUE)
    
    # Calculate the Standardized Plasticity Index (SPI) for this trait
    spi = (trait_env2_value - trait_env1_value) / sd_reference
    
    # Store the result
    spi_results[i] = spi
  }
  
  return(spi_results)
}

calculate_SPI(synthetic_data1,1,c(3,4),1,2,3)

########################

#' Calculate Absolute Plasticity Coefficient (APC)
#'
#' This function calculates the Absolute Plasticity Coefficient (APC), 
#' which quantifies absolute plasticity by calculating the average absolute difference 
#' in trait values between environments. If the environments are not sequential, 
#' the function will use the order provided in the dataset.
#' 
#' The formula used is:
#' \deqn{APC = \frac{1}{n-1} \sum_{i=1}^{n-1} |Trait Value_{Envi+1} - Trait Value_{Envi}|}
#' 
#' @param data A data frame containing the input data.
#' @param env_col A string or numeric value indicating the column in `data` that contains the environment labels.
#' @param trait_cols A character vector or numeric vector indicating the columns in the dataset containing the trait values.
#' @param sequential_env A logical indicating whether the environments are in sequential order (e.g., increasing temperature).
#' @return A numeric vector of APC values, one for each trait column specified.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(10, 12, 15, 20),
#'   Trait2 = c(8, 11, 13, 17),
#'   Environment = c(1, 2, 3, 4)  # Sequential environments (e.g., increasing temperature)
#' )
#' apc_values = calculate_APC(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"))
#' print(apc_values)
#' 
#' @export
calculate_APC = function(data, env_col, trait_cols, sequential_env = TRUE) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  
  # Initialize a vector to store APC values
  apc_results = c()
  
  # Handle environment ordering if sequential_env is FALSE
  if (!sequential_env) {
    message("Environments are not sequential. Using the provided order in the dataset.")
  } else {
    # Ensure the environments are ordered if sequential
    ordered_indices = order(env_col)
    data = data[ordered_indices, ]
  }
  
  # Loop through each trait column and calculate the APC
  for (trait_col in trait_cols) {
    
    # Extract the trait data
    trait_values = data[[trait_col]]
    
    # Calculate the absolute differences between consecutive environments
    abs_diff = abs(diff(trait_values))
    
    # Calculate the APC
    apc = mean(abs_diff, na.rm = TRUE)
    
    # Store the result
    apc_results = c(apc_results, apc)
  }
  
  # Return the APC values for each trait
  names(apc_results) = trait_cols
  return(apc_results)
}


calculate_APC(synthetic_data1,1,c(3,4))

#############################


#' Calculate Stability Index (SI)
#'
#' This function calculates the Stability Index (SI), which quantifies the stability 
#' of a trait across environments. It measures the consistency of a genotype’s trait expression 
#' under varying conditions.
#'
#' The formula used is:
#' \deqn{SI = \frac{Variance across environments}{Mean trait value across environments}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A character vector or numeric vector indicating the columns in the dataset containing the trait values.
#' @return A numeric vector of Stability Index (SI) values, one for each trait column specified.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(10, 12, 15, 20),
#'   Trait2 = c(8, 11, 13, 17),
#'   Environment = c(1, 2, 3, 4)  # Different environments
#' )
#' si_values = calculate_SI(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"))
#' print(si_values)
#'
#' @export
calculate_SI = function(data, env_col, trait_cols) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else {
    env_col = env_col
  }
  
  # Initialize a vector to store SI values
  si_results = c()
  
  # Loop through each trait column and calculate the SI
  for (trait_col in trait_cols) {
    
    # Extract the trait data
    trait_values = data[[trait_col]]
    
    # Calculate the variance across environments
    variance_trait = var(trait_values, na.rm = TRUE)
    
    # Calculate the mean trait value across environments
    mean_trait_value = mean(trait_values, na.rm = TRUE)
    
    # Calculate the Stability Index (SI)
    si = variance_trait / mean_trait_value
    
    # Store the result
    si_results = c(si_results, si)
  }
  
  # Return the SI values for each trait
  names(si_results) = trait_cols
  return(si_results)
}


#######################


#' Calculate Relative Stability Index (RSI)
#'
#' This function calculates the Relative Stability Index (RSI), which measures the relative 
#' stability of a trait across environments, focusing on stability despite environmental fluctuations.
#'
#' The formula used is:
#' \deqn{RSI = 1 - \frac{Standard Deviation of Trait Values Across Environments}{Mean Trait Value Across Environments}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A character vector or numeric vector indicating the columns in the dataset containing the trait values.
#' @return A numeric vector of Relative Stability Index (RSI) values, one for each trait column specified.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(10, 12, 15, 20),
#'   Trait2 = c(8, 11, 13, 17),
#'   Environment = c(1, 2, 3, 4)  # Different environments
#' )
#' rsi_values = calculate_RSI(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"))
#' print(rsi_values)
#'
#' @export
calculate_RSI = function(data, env_col, trait_cols) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else {
    env_col = env_col
  }
  
  # Initialize a vector to store RSI values
  rsi_results = c()
  
  # Loop through each trait column and calculate the RSI
  for (trait_col in trait_cols) {
    
    # Extract the trait data
    trait_values = data[[trait_col]]
    
    # Calculate the standard deviation and mean trait value across environments
    sd_trait = sd(trait_values, na.rm = TRUE)
    mean_trait_value = mean(trait_values, na.rm = TRUE)
    
    # Calculate the Relative Stability Index (RSI)
    rsi = 1 - (sd_trait / mean_trait_value)
    
    # Store the result
    rsi_results = c(rsi_results, rsi)
  }
  
  # Return the RSI values for each trait
  names(rsi_results) = trait_cols
  return(rsi_results)
}

##########################

#' Calculate Environmental Variance Sensitivity (EVS)
#'
#' This function calculates the Environmental Variance Sensitivity (EVS), 
#' which quantifies how sensitive a trait is to changes in environmental variance.
#'
#' The formula used is:
#' \deqn{EVS = \frac{\Delta Trait Variance}{\Delta Environmental Variance}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A character or numeric vector indicating the columns in the dataset containing the trait values.
#' @return A numeric vector of Environmental Variance Sensitivity (EVS) values, one for each trait column specified.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(10, 12, 15, 18),
#'   Trait2 = c(8, 11, 13, 17),
#'   Environment = c(1, 2, 3, 4),  # Different environments
#'   Environmental_Variance = c(1.1, 2.2, 1.5, 2.5)  # Environmental variance
#' )
#' evs_values = calculate_EVS(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"), env_variance_col = "Environmental_Variance")
#' print(evs_values)
#'
#' @export
calculate_EVS = function(data, env_col, trait_cols, env_variance_col) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else {
    env_col = env_col
  }
  
  # Extract environmental variance data
  env_variance = data[[env_variance_col]]
  
  # Initialize a vector to store EVS values
  evs_results = c()
  
  # Loop through each trait column and calculate the EVS
  for (trait_col in trait_cols) {
    
    # Extract the trait data
    trait_values = data[[trait_col]]
    
    # Calculate the change in trait variance across environments
    delta_trait_variance = var(trait_values, na.rm = TRUE)
    
    # Calculate the change in environmental variance
    delta_env_variance = var(env_variance, na.rm = TRUE)
    
    # Calculate the Environmental Variance Sensitivity (EVS)
    evs = delta_trait_variance / delta_env_variance
    
    # Store the result
    evs_results = c(evs_results, evs)
  }
  
  # Return the EVS values for each trait
  names(evs_results) = trait_cols
  return(evs_results)
}
#################

#' Calculate Environmental Variance Sensitivity (EVS)
#'
#' This function calculates the Environmental Variance Sensitivity (EVS), 
#' which quantifies how sensitive a trait is to changes in environmental variance.
#'
#' The formula used is:
#' \deqn{EVS = \frac{\Delta Trait Variance}{\Delta Environmental Variance}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A character or numeric vector indicating the columns in the dataset containing the trait values.
#' @return A numeric vector of Environmental Variance Sensitivity (EVS) values, one for each trait column specified.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(10, 12, 15, 18),
#'   Trait2 = c(8, 11, 13, 17),
#'   Environment = c(1, 2, 3, 4)  # Different environments
#' )
#' evs_values = calculate_EVS(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"))
#' print(evs_values)
#'
#' @export
calculate_EVS = function(data, env_col, trait_cols) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else {
    env_col = env_col
  }
  
  # Initialize a vector to store EVS values
  evs_results = c()
  
  # Loop through each trait column and calculate the EVS
  for (trait_col in trait_cols) {
    
    # Extract the trait data
    trait_values = data[[trait_col]]
    
    # Calculate the variance of trait values across environments
    trait_variance = var(trait_values, na.rm = TRUE)
    
    # Calculate the variance of environments
    env_variance = var(env_col, na.rm = TRUE)
    
    # Calculate the Environmental Variance Sensitivity (EVS)
    evs = trait_variance / env_variance
    
    # Store the result
    evs_results = c(evs_results, evs)
  }
  
  # Return the EVS values for each trait
  names(evs_results) = trait_cols
  return(evs_results)
}


#####################

#' Calculate Multivariate Plasticity Index (MVPi)
#'
#' This function calculates the Multivariate Plasticity Index (MVPi) by using Canonical Variate Analysis (CVA)
#' to compute the Euclidean distances between phenotypic states across different environments.
#' The formula used is based on the Euclidean distance between phenotypic states.
#'
#' @param data A data frame containing phenotypic trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A numeric vector specifying the columns in the dataset containing phenotypic trait values.
#' @param n_axes Integer, the number of canonical variate axes (CVs) to consider.
#' @return A numeric value representing the Multivariate Plasticity Index (MVPi).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(1.2, 2.3, 1.9, 3.1, 2.0, 2.8),
#'   Trait2 = c(2.5, 3.0, 2.8, 3.6, 2.9, 3.3),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3))
#' )
#' mvpi_value = calculate_MVPi(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"), n_axes = 2)
#' print(mvpi_value)
#'
#' @export
calculate_MVPi = function(data, env_col, trait_cols, n_axes) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Extract phenotypic trait values
  trait_data = as.matrix(data[, trait_cols])
  
  # Ensure env_col is factor and create a design matrix for environments
  env_col = factor(env_col)
  env_design = model.matrix(~ env_col - 1)
  
  # Perform SVD for ordination
  svd_result = svd(trait_data)
  
  # Calculate scores for the samples in the reduced space
  scores = svd_result$u %*% diag(svd_result$d)
  
  # Select the first n_axes
  selected_scores = scores[, 1:n_axes]
  
  # Calculate Euclidean distances between environment centroids
  env_levels = unique(env_col)
  distances = matrix(0, nrow = length(env_levels), ncol = length(env_levels))
  
  for (i in 1:(length(env_levels) - 1)) {
    for (j in (i + 1):length(env_levels)) {
      # Calculate centroids in the reduced space for each environment
      centroid_i = colMeans(selected_scores[env_col == env_levels[i], ])
      centroid_j = colMeans(selected_scores[env_col == env_levels[j], ])
      
      # Calculate Euclidean distance between centroids
      distances[i, j] = sqrt(sum((centroid_i - centroid_j)^2))
    }
  }
  
  # Calculate the average Euclidean distance as MVPi
  mvpi = mean(distances[distances > 0])
  
  return(mvpi)
}

calculate_MVPi(synthetic_data1,1,c(3,4,5),3)


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


calculate_SPM(synthetic_data1,1,3,1,2)


############################


#' Calculate SSpop/SStotal Plasticity Ratio
#'
#' This function calculates the Plasticity Ratio, which quantifies the relative contribution 
#' of population differences to the overall phenotypic variability across environments. 
#' It is calculated as the ratio of the sum of squares between populations (SSpop) 
#' to the total sum of squares (SStotal) from a one-way ANOVA.
#'
#' The formula used is:
#' \deqn{Plasticity Ratio = \frac{SSpop}{SStotal}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string or  numeric indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_col A string or numeric value indicating the column in the dataset containing the trait values.
#' @param pop_col A string or numeric value indicating the column in the dataset containing the population labels, or an external vector specifying population belonging.
#' @return A numeric value representing the SSpop/SStotal Plasticity Ratio.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_Values = c(10, 12, 15, 18, 11, 17),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3)),
#'   Population = factor(c("A", "A", "B", "B", "C", "C"))
#' )
#' plasticity_ratio = calculate_Plasticity_Ratio(synthetic_data, env_col = "Environment", trait_col = "Trait_Values", pop_col = "Population")
#' print(plasticity_ratio)
#'
#' @export
calculate_Plasticity_Ratio = function(data, env_col, trait_col, pop_col) {
  
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
  
  trait_col = data[[trait_col]]
  
  
  # Perform a one-way ANOVA to obtain SSpop and SStotal
  anova_result = aov(trait_col ~ pop_col + env_col, data = data)
  print(anova_result)
  # Extract the sum of squares between populations (SSpop) and total sum of squares (SStotal)
  ss_total = sum(anova(anova_result)[, "Sum Sq"])
  ss_pop = anova(anova_result)["pop_col", "Sum Sq"]
  
  
  # Calculate the Plasticity Ratio
  plasticity_ratio = ss_pop / ss_total
  
  return(plasticity_ratio)
}

pop=c(rep(1,150),rep(2,150))
calculate_Plasticity_Ratio(synthetic_data1,1,3,pop)
