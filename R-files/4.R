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


source("~/CRC 1622 - Z2/R-files/2.R")

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
#' @param time_interval Numeric, the time interval between the two time points.
#'
#' @return A numeric vector of DPI values.
#' 
#' @examples
#' # Example usage with synthetic data
#' trait_time1 = c(10, 12, 14)
#' trait_time2 = c(15, 18, 22)
#'
#' # Calculate DPI
#' dpi_values = calculate_DPI(trait_time1, trait_time2, 5)
#' print(dpi_values)
#'
#' @export
calculate_DPI = function(trait_values_time1, trait_values_time2, time_interval) {
  # Ensure inputs are numeric vectors
  if (!is.numeric(trait_values_time1) || !is.numeric(trait_values_time2)) {
    stop("trait_values_time1 and trait_values_time2 must be numeric vectors.")
  }
  
  # Ensure the vectors are of equal length
  if (length(trait_values_time1) != length(trait_values_time2)) {
    stop("trait_values_time1 and trait_values_time2 must have the same length.")
  }
  
  # Ensure the time interval is a positive numeric value
  if (!is.numeric(time_interval) || time_interval <= 0) {
    stop("time_interval must be a positive numeric value.")
  }
  
  # Calculate DPI
  dpi = (trait_values_time2 - trait_values_time1) / time_interval
  
  return(dpi)
}


## test - passed on synthetic dataset - however it is not clear how one would structure a dataset with time resolved data (will the measurements at time x for one sample be in different columns all with equally separated time intervals or in rows?)
# mathematically the function does the wanted if one knows how to structure the data

#calculate_DPI(df_test2,2,3,3)



###############################



#' Calculate Coefficient of Environmental Variation (CEV) for Multiple Traits
#'
#' This function calculates the Coefficient of Environmental Variation (CEV) for multiple traits,
#' which measures the variability of each trait across different environments.
#'
#' The formula used is:
#' \deqn{CEV = \frac{Standard Deviation of Trait Values Across Environments}{Mean Trait Value Across Environments} \times 100}
#'
#' @param data A data frame containing the input data.
#' @param trait_cols A vector of strings or numeric values indicating the columns for trait values across environments.
#' @return A data frame with each trait and its corresponding Coefficient of Environmental Variation (CEV).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_1 = c(10, 15, 20, 25, 30),
#'   Trait_2 = c(5, 7, 9, 11, 13)
#' )
#'
#' # Calculate CEV for multiple traits
#' cev_values = calculate_CEV(synthetic_data, trait_cols = c("Trait_1", "Trait_2"))
#' print(cev_values)
#'
#' @export
calculate_CEV = function(data, trait_cols) {
  
  # Initialize a data frame to store CEV results for each trait
  CEV_results = data.frame(Trait = character(), CEV = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait column
  for (trait_col in trait_cols) {
    
    # Handle the input column for trait values
    trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Calculate the mean and standard deviation of the trait values
    mean_trait_value = mean(trait_values, na.rm = TRUE)
    sd_trait_value = sd(trait_values, na.rm = TRUE)
    
    # Calculate CEV
    cev = (sd_trait_value / mean_trait_value) * 100
    
    # Append the result to the data frame
    CEV_results = rbind(CEV_results, data.frame(Trait = colnames(data)[trait_col], CEV = cev))
  }
  
  return(CEV_results)
}



## tested - passed on synthetic dataset 

#calculate_CEV(df_test2,c(2,3))
#sd(df_test2[,3])/mean(df_test2[,3])*100

################################


#' Calculate Plasticity Response Index (PRI) for Multiple Traits
#'
#' This function calculates the Plasticity Response Index (PRI) for multiple traits, quantifying plasticity
#' by measuring the proportional difference in trait values between an extreme environment (e.g., the best or worst condition) and a control environment.
#'
#' The formula used is:
#' \deqn{PRI = \frac{Trait Value in Extreme Environment - Trait Value in Control Environment}{Mean Trait Value Across All Environments}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param trait_cols A vector of strings or numeric values indicating the columns in the dataset containing the trait values.
#' @param env_extreme A numeric or string value representing the identifier for the extreme environment in `env_col`.
#' @param env_control A numeric or string value representing the identifier for the control environment in `env_col`.
#' @return A data frame with each trait and its corresponding Plasticity Response Index (PRI).
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_1 = c(10, 15, 20, 25, 30),
#'   Trait_2 = c(5, 8, 12, 16, 20),
#'   Environment = factor(c(1, 1, 2, 2, 2))  # Control = 1, Extreme = 2
#' )
#' pri_values = calculate_PRI(synthetic_data, env_col = "Environment", trait_cols = c("Trait_1", "Trait_2"), env_extreme = 2, env_control = 1)
#' print(pri_values)
#'
#' @export
calculate_PRI = function(data, env_col, trait_cols, env_extreme, env_control) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a data frame to store PRI results for each trait
  PRI_results = data.frame(Trait = character(), PRI = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait column
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Calculate the mean trait value across all environments
    mean_trait_value = mean(trait_values, na.rm = TRUE)
    
    # Calculate trait values for extreme and control environments
    trait_ex = mean(trait_values[env_col == env_extreme], na.rm = TRUE)
    trait_con = mean(trait_values[env_col == env_control], na.rm = TRUE)
    
    # Calculate the Plasticity Response Index (PRI)
    pri = (trait_ex - trait_con) / mean_trait_value
    
    # Append the result to the data frame
    PRI_results = rbind(PRI_results, data.frame(Trait = colnames(data)[trait_col], PRI = pri))
  }
  
  return(PRI_results)
}

## test - passed on synthetic dataset

#calculate_PRI(df_test2,1,2,1,2)

############################

#' Calculate Phenotypic Flexibility Index (PFI) for Multiple Traits
#'
#' This function calculates the Phenotypic Flexibility Index (PFI) for each specified trait by finding the sample 
#' with the maximum absolute deviation from its respective baseline value, then dividing this deviation by that sample's baseline value.
#'
#' @param data A data frame containing the input data.
#' @param baseline_cols A vector of column names or indices indicating the baseline columns for each trait.
#' @param trait_cols A vector of column names or indices indicating the columns for trait values.
#' @return A named numeric vector containing the Phenotypic Flexibility Index (PFI) for each trait.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Baseline_Trait1 = c(20, 20, 20, 20, 20),
#'   Baseline_Trait2 = c(15, 15, 15, 15, 15),
#'   Trait_1 = c(18, 22, 25, 20, 19),
#'   Trait_2 = c(10, 15, 13, 18, 17)
#' )
#' pfi_values = calculate_PFI(synthetic_data, baseline_cols = c("Baseline_Trait1", "Baseline_Trait2"), trait_cols = c("Trait_1", "Trait_2"))
#' print(pfi_values)
#'
#' @export
calculate_PFI = function(data, baseline_cols, trait_cols) {
  
  # Initialize vector to store PFI values for each trait
  pfi_values = numeric(length(trait_cols))
  names(pfi_values) = trait_cols
  
  # Loop through each trait to calculate PFI independently
  for (i in seq_along(trait_cols)) {
    # Extract trait values and corresponding baseline
    trait_values = if (is.numeric(trait_cols[i])) data[[trait_cols[i]]] else data[[trait_cols[i]]]
    baseline_values = if (is.numeric(baseline_cols[i])) data[[baseline_cols[i]]] else data[[baseline_cols[i]]]
    
    # Calculate the absolute deviation from baseline for each sample
    deviations = abs(trait_values - baseline_values)
    
    # Find the sample with the maximum deviation
    max_deviation_index = which.max(deviations)
    max_deviation = deviations[max_deviation_index]

    # Calculate PFI as the maximum deviation divided by the corresponding baseline value
    pfi_values[i] = max_deviation / baseline_values[max_deviation_index]
  }
  
  return(pfi_values)
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

#' Calculate Absolute Plasticity Coefficient (APC)
#'
#' This function calculates the Absolute Plasticity Coefficient (APC), 
#' which quantifies absolute plasticity by calculating the average absolute difference 
#' in the mean trait values between environments. If the environments are not sequential, 
#' the function will use the order provided in the dataset.
#'
#' The formula used is:
#' \deqn{APC = \frac{1}{n-1} \sum_{i=1}^{n-1} |Mean(Trait Value)_{Envi+1} - Mean(Trait Value)_{Envi}|}
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
  
  # Calculate the APC for each trait column
  for (trait_col in trait_cols) {
    
    # Group data by environment and calculate the mean trait value for each environment
    env_means = aggregate(data[[trait_col]], by = list(env_col), FUN = mean)
    
    # Calculate the absolute differences between consecutive environments
    abs_diff = abs(diff(env_means$x))
  
    # Calculate the APC
    apc = mean(abs_diff, na.rm = TRUE)
    
    # Store the result
    apc_results = c(apc_results, apc)
  }
  
  # Return the APC values for each trait
  names(apc_results) = trait_cols
  return(apc_results)
}
## test - passed on synthetic dataset 
#calculate_APC(df_test2,1,2)

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

### test - passed on synthetic dataset 
#var(df_test2[,2])/mean(df_test2[,2])
#calculate_SI(df_test2,1,2)
#######################


#' Calculate Relative Stability Index (RSI)
#'
#' This function calculates the Relative Stability Index (RSI), which measures the relative 
#' stability of a trait across environments, focusing on stability despite environmental fluctuations. 
#' Users can optionally specify one or more environments for which to calculate the RSI. 
#' If no specific environments are provided, the RSI will be calculated over all available environments.
#'
#' The formula used is:
#' \deqn{RSI = 1 - \frac{Standard Deviation of Trait Values Across Environments}{Mean Trait Value Across Environments}}
#'
#' @param data A data frame containing the input data.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels.
#' @param trait_cols A character vector or numeric vector indicating the columns in the dataset containing the trait values.
#' @param specific_envs An optional vector of specific environments to include in the RSI calculation. 
#' If `NULL`, the function calculates RSI across all environments (default is `NULL`).
#' @return A numeric vector of Relative Stability Index (RSI) values, one for each trait column specified.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait1 = c(10, 12, 15, 20),
#'   Trait2 = c(8, 11, 13, 17),
#'   Environment = c(1, 2, 3, 4)  # Different environments
#' )
#'
#' # Calculate RSI for all environments
#' rsi_values = calculate_RSI(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"))
#' print(rsi_values)
#'
#' # Calculate RSI for specific environments (1 and 2)
#' rsi_values_specific = calculate_RSI(synthetic_data, env_col = "Environment", trait_cols = c("Trait1", "Trait2"), specific_envs = c(1, 2))
#' print(rsi_values_specific)
#'
#' @export
calculate_RSI = function(data, env_col, trait_cols, specific_envs = NULL) {
  
  # If specific environments are provided, subset the data accordingly
  if (!is.null(specific_envs)) {
    data = data[data[[env_col]] %in% specific_envs, ]
  }
  
  # Initialize a vector to store RSI values
  rsi_results = c()
  
  # Loop through each trait column and calculate the RSI
  for (trait_col in trait_cols) {
    
    # Extract the trait data for the specified environments (or all environments if none are specified)
    trait_values = data[[trait_col]]
    
    # Calculate the standard deviation and mean trait value across the selected environments
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
calculate_EVS = function(data, env_col, trait_cols ) {
 
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_variance = data[[env_col]]
  } else {
    env_variance = env_col
  }
  
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



