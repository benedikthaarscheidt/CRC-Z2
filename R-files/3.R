#this files contains the following indices/functions: 
#generate_synthetic_data,
#Phenotypic Stability Index (calculate_PSI),
#Relative Plasticity Index (calculate_RPI) - passed,
#Plasticity Quotient (calculate_PQ) - passed,
#Phenotypic Range (PR) (calculate_PR) - passed,
#Norm of reaction width (calculate_NRW) - passed,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - passed,
#Calculate Plasticity Differential (PD) (calculate_PD) - passed,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - passed,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - passed,

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

df_test6=data.frame(
  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   
  Column1 = c(rep(6, 10), rep(4, 10), rep(2, 10)),   
  Column2 = c(rep(1, 5),rep(2,5), rep(2, 10), rep(3, 10)),
  Column3 = c(rep(1, 10), rep(1, 10), rep(1, 10))
)



#########################################

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


#this files contains the following indices: 
#' Calculate Relative Plasticity Index (RPI)
#'
#' This function calculates the Relative Plasticity Index (RPI) for a given trait across two environments.
#' The RPI quantifies the relative change in a trait between two specified environments using the formula:
#' 
#' \deqn{RPI = \frac{|TraitValueEnv1 - TraitValueEnv2|}{TraitValueEnv1 + TraitValueEnv2}}
#' 
#' The user must specify the two environments either by providing the column index for the environment in 
#' the dataset or by passing an external vector indicating environmental belonging.
#'
#' @param data A data frame containing the input data, with one column for the trait and another column or external vector for the environment.
#' @param trait_col A string or numeric value indicating the column in `data` that contains the trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param env1 A value indicating the first environment to compare. Must be present in `env_col`.
#' @param env2 A value indicating the second environment to compare. Must be present in `env_col`.
#'
#' @return A numeric value representing the calculated RPI for the specified trait between the two environments.
#'
#' @details The RPI is calculated by taking the absolute difference between the trait values in the two environments
#' and dividing it by the sum of the trait values in the two environments. This index provides a relative measure of 
#' the trait's plasticity across the two environments.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait = c(100, 110, 120, 130, 140, 150),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3))
#' )
#' 
#' # Calculate RPI between environment 1 and 2
#' rpi_value = calculate_RPI(synthetic_data, trait_col = "Trait", env_col = "Environment", env1 = 1, env2 = 2)
#' print(rpi_value)
#'
#' @export
calculate_RPI = function(data, trait_col, env_col, env1, env2) {
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
  
  # Filter the data for the two specified environments
  env1_data = trait_col[env_col == env1]
  env2_data = trait_col[env_col == env2]
  
  # Ensure there are equal numbers of values for env1 and env2
  if (length(env1_data) != length(env2_data)) {
    stop("The number of trait values for the two environments must be equal.")
  }
  
  # Calculate RPI for each pair of values between the two environments
  RPI_values = abs(env1_data - env2_data) / (env1_data + env2_data)
  
  # Return the average RPI across all data points
  RPI = mean(RPI_values, na.rm = TRUE)
  
  return(RPI)
}


## test - passed on synthetic dataset

calculate_RPI(df_test2,2,1,1,2)

######################################

#' @title Calculate Plasticity Quotient (PQ)
#' @description This function calculates the Plasticity Quotient (PQ) based on the range of trait values across different environmental conditions and the corresponding environmental factor range.
#'
#' @param data A data frame containing the data for the analysis.
#' @param trait_col A column name (string) or numeric index representing the trait values in the data. This can also be a vector of trait values with the same length as the number of rows in the data.
#' @param env_col A column name (string) or numeric index representing the environmental variable in the data, or a vector of environmental values. This specifies the environmental conditions for calculating plasticity.
#' @param factor_col A column name (string) or numeric index representing a relevant environmental factor (e.g., temperature, moisture) in the data, or a vector of environmental factor values. This will be used for standardizing the trait range.
#'
#' @details 
#' The function first handles the provided column names or indices for `trait_col`, `env_col`, and `factor_col`, ensuring they are correctly extracted from the data. 
#' 
#' It then calculates the mean trait values and factor values across the different environmental conditions specified in `env_col`. The range of the mean trait values is computed, and this is divided by the range of the mean factor values to compute the Plasticity Quotient (PQ).
#'
#' @return A numeric value representing the Plasticity Quotient (PQ), which is the ratio of the range of trait values to the range of environmental factor values.
#'
#' @examples
#' # Example usage:
#' # Assuming 'data' is a data frame with columns 'trait', 'env', and 'factor':
#' PQ = calculate_PQ(data, 'trait', 'env', 'factor')
#' 
#' @export
calculate_PQ = function(data, env_col, trait_col, factor_col){
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Handle factor_col
  if (is.numeric(factor_col) && length(factor_col) == 1) {
    factor_col = data[[factor_col]]
  } else if (is.vector(factor_col) && length(factor_col) == nrow(data)) {
    factor_col = factor_col
  } else {
    factor_col = data[[factor_col]]
  }
  
  # Extract trait data
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Calculation of PQ
  mean_values = aggregate(cbind(trait_col, factor_col) ~ env_col, data = data, FUN = mean)
  range_trait = abs(max(mean_values$trait_col) - min(mean_values$trait_col))
  env_range = abs(max(mean_values$factor_col) - min(mean_values$factor_col))
  
  # Return PQ
  return(range_trait / env_range)
}

## test - passed on synthetic dataframe

calculate_PQ(df_test_6,1,3,2)

###########################################

#' @description This function calculates the Phenotypic Range (PR) for a trait across specified environments. PR is defined as the difference between the maximum and minimum trait values within each environment.
#'
#' @param data A data frame containing the data for the analysis.
#' @param env_col A column name (string), numeric index, or vector representing the environmental variable in the data. This specifies the environment for calculating the PR.
#' @param trait_col A column name (string) or numeric index representing the trait values in the data. This can also be a vector of trait values with the same length as the number of rows in the data.
#' @param env Optional. A vector specifying which environments to include in the calculation. If not provided, all unique environments from `env_col` are used.
#'
#' @details 
#' The function calculates the Phenotypic Range (PR) for each specified environment, defined as the range (maximum minus minimum) of trait values within each environment. If `env` is provided, only the specified environments are used; otherwise, it uses all unique environments in `env_col`. 
#'
#' @return A numeric vector containing the PR for each environment. The length of the vector is equal to the number of unique environments or the length of the `env` vector (if provided).
#'
#' @examples
#' # Example usage:
#' # Assuming 'data' is a data frame with columns 'trait' and 'env':
#' PR = calculate_PR(data, 'env', 'trait')
#' # Example with specific environments:
#' PR = calculate_PR(data, 'env', 'trait', env = c("Environment1", "Environment2"))
#'
#' @export
calculate_PR=function(data,env_col,trait_col,env=NULL, across=F){
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Extract trait data
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  if(across==T){
    # Calculate mean trait values for each environment
    means = tapply(trait_col, env_col, mean, na.rm = TRUE)
    # Identify maximum and minimum means
    max_mean = max(means, na.rm = TRUE)
    min_mean = min(means, na.rm = TRUE)
    return(abs(max_mean-min_mean))
  }else{
  
    if(is.null(env)){
      env=unique(env_col)
    }else{
      env=env
    }
    
    PR=rep(0,length(env))
    for(i in 1:length(env)){
      PR[i]=max(trait_col[which(env_col==env[i])])-min(trait_col[which(env_col==env[i])])
    }
    
    return(PR)
    }
  }



## test - passed on synthetic dataset 


calculate_PR(df_test6,1,3)
calculate_PR(df_test6,1,2,across=T)






##########################

#' @description This function calculates the Norm of Reaction Width (NRW) for a specific genotype or group, measuring the range of trait values across different environments. If no grouping factor is specified, it calculates the NRW for the entire dataset.
#'
#' @param data A data frame containing the data for the analysis.
#' @param trait_col A column name (string) or numeric index representing the trait values in the data. This can also be a vector of trait values with the same length as the number of rows in the data.
#' @param env_col A column name (string) or numeric index representing the environmental variable in the data, or a vector of environmental values. This specifies the environmental conditions for calculating NRW.
#' @param group_col Optional. A column name (string) or numeric index representing the grouping factor (e.g., genotype) in the data. If not provided, all data will be treated as one group.
#'
#' @details 
#' The function calculates the Norm of Reaction Width (NRW) by determining the range (maximum minus minimum) of trait values for each group (or the entire dataset if no group is specified) across all environments. It returns a data frame with the NRW for each group.
#'
#' @return A data frame with two columns: one for the group identifier and the other for the corresponding Norm of Reaction Width (NRW) for each group.
#'
#' @examples
#' # Example usage:
#' # Assuming 'data' is a data frame with columns 'trait', 'env', and 'genotype':
#' NRW = calculate_NRW(data, 'trait', 'env', 'genotype')
#' # If no grouping factor is provided:
#' NRW = calculate_NRW(data, 'trait', 'env')
#'
#' @export
calculate_NRW = function(data,env_col=NULL, trait_col) {
  
  if (is.numeric(trait_col)) trait_col = data[[trait_col]]
  if (is.null(env_col)){env_col=rep(1,nrow(data))}
  else{env_col = data[[env_col]]}
  
  # Calculate the range (NRW) for each group across environments
  NRW_df = aggregate(trait_col ~  env_col, data = data, FUN = function(x) max(x) - min(x))
  
  # Rename columns to make the result clearer
  colnames(NRW_df) = c("Group", "NRW")
  
  return(NRW_df)
}


## test - passed on synthetic dataset


calculate_NRW(df_test6,1,3)

##############################



calculate_ESP = function(data, env_col,trait_col,  env=NULL) {
  # Handle env_col properly (either column index, column name, or vector)
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # If env is NULL, calculate ESP for all unique environments
  if (is.null(env)) {
    env = unique(env_col)
  }
  
  # Extract trait data properly (either column index or column name)
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Calculate the mean trait value across all environments
  mean_trait_all_envs = mean(trait_col, na.rm = TRUE)
  
  # Initialize ESP vector
  ESP = rep(1, length(env))
  
  # Loop through each environment and calculate ESP
  for (i in 1:length(env)) {
    trait_values_env = trait_col[env_col == env[i]]
    ESP[i] = (mean(trait_values_env, na.rm = TRUE) - mean_trait_all_envs) / mean_trait_all_envs
  }
  
  return(ESP)
}

## test - passed on synthetic dataset 

calculate_ESP(df_test6,1,3,env=NULL)

######################### 01.10.2024

#' Calculate Plasticity Differential (PD)
#'
#' This function calculates the Plasticity Differential (PD) for a given trait across control and stress conditions.
#' The PD quantifies the absolute difference in trait values between the control and stress environments using the formula:
#'
#' \deqn{PD = |TraitValueStress - TraitValueControl|}
#'
#' The user must specify the control and stress environments either by providing the column index for the environment in
#' the dataset or by passing an external vector indicating environmental belonging.
#'
#' @param data A data frame containing the input data, with one column for the trait and another column or external vector for the environment.
#' @param trait_col A string or numeric value indicating the column in `data` that contains the trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param control_env A value indicating the control environment to compare. Must be present in `env_col`.
#' @param stress_env A value indicating the stress environment to compare. Must be present in `env_col`.
#'
#' @return A numeric value representing the calculated PD for the specified trait between the control and stress environments.
#'
#' @details The PD is calculated by taking the absolute difference between the trait values in the stress and control environments.
#' This index provides a measure of the change in trait value under stress relative to the control.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait = c(100, 110, 120, 130, 140, 150),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3)) # Assume 1=Control, 2=Stress
#' )
#' 
#' # Calculate PD between control (1) and stress (2)
#' pd_value = calculate_PD(synthetic_data, trait_col = "Trait", env_col = "Environment", control_env = 1, stress_env = 2)
#' print(pd_value)
#'
#' @export
calculate_PD = function(data, env_col, trait_col, control_env, stress_env) {
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
  
  # Filter the data for the control and stress environments
  control_data = trait_col[env_col == control_env]
  stress_data = trait_col[env_col == stress_env]
  
  # Ensure there are equal numbers of values for control and stress
  if (length(control_data) != length(stress_data)) {
    stop("The number of trait values for the control and stress environments must be equal.")
  }
  
  # Calculate PD for each pair of values between the control and stress environments
  PD_values = abs(stress_data - control_data)
  print(PD_values)
  
  # Return the average PD across all data points
  PD = mean(PD_values, na.rm = TRUE)
  
  return(PD)
}

## test - passed on synthetic dataset 

calculate_PD(df_test6,1,3,3,1)
##########################


#' Calculate Fitness Plasticity Index (FPI)
#'
#' This function calculates the Fitness Plasticity Index (FPI) for a fitness-related trait (e.g., survival, biomass)
#' across control and stress conditions. The FPI quantifies the relative change in fitness between the control
#' and stress environments using the formula:
#'
#' \deqn{FPI = \frac{FitnessStress - FitnessControl}{FitnessControl}}
#'
#' The user must specify the control and stress environments either by providing the column index for the environment
#' in the dataset or by passing an external vector indicating environmental belonging.
#'
#' @param data A data frame containing the input data, with one column for the fitness trait and another column or external vector for the environment.
#' @param trait_col A string or numeric value indicating the column in `data` that contains the fitness trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param control_env A value indicating the control environment to compare. Must be present in `env_col`.
#' @param stress_env A value indicating the stress environment to compare. Must be present in `env_col`.
#'
#' @return A numeric value representing the calculated FPI for the specified trait between the control and stress environments.
#'
#' @details The FPI is calculated by taking the difference between the fitness values in the stress and control environments,
#' divided by the fitness values in the control environment. This index provides a relative measure of the fitness trait's
#' plasticity across the two environments.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Fitness = c(10, 12, 14, 16),   # Fitness values
#'   Environment = factor(c(1, 1, 2, 2))  # Control = 1, Stress = 2
#' )
#'
#' # Calculate FPI between control (1) and stress (2)
#' fpi_value = calculate_FPI(synthetic_data, trait_col = "Fitness", env_col = "Environment", control_env = 1, stress_env = 2)
#' print(fpi_value)
#'
#' @export
calculate_FPI = function(data, env_col, trait_col, control_env, stress_env) {
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
  
  # Filter the data for the control and stress environments
  control_data = trait_col[env_col == control_env]
  stress_data = trait_col[env_col == stress_env]
  
  # Ensure there are equal numbers of values for control and stress
  if (length(control_data) != length(stress_data)) {
    stop("The number of trait values for the control and stress environments must be equal.")
  }
  
  # Calculate FPI for each pair of values between the control and stress environments
  FPI_values = (stress_data - control_data) / control_data
  
  # Return the average FPI across all data points
  FPI = mean(FPI_values, na.rm = TRUE)
  
  return(FPI)
}

## test - passed on synthetic dataset

calculate_FPI(df_test6,1,3,1,2)


#################################


#' Calculate Transplant Plasticity Score (TPS)
#'
#' This function calculates the Transplant Plasticity Score (TPS) for a given trait when an organism is transplanted to a different environment.
#' The TPS quantifies the relative change in a trait between the native environment and the transplanted environment using the formula:
#'
#' \deqn{TPS = \frac{TraitValueTransplantedEnv - TraitValueNativeEnv}{TraitValueNativeEnv}}
#'
#' The user must specify the native and transplanted environments either by providing the column index for the environment
#' in the dataset or by passing an external vector indicating environmental belonging.
#'
#' @param data A data frame containing the input data, with one column for the trait and another column or external vector for the environment.
#' @param trait_col A string or numeric value indicating the column in `data` that contains the trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param native_env A value indicating the native environment. Must be present in `env_col`.
#' @param transplanted_env A value indicating the transplanted environment. Must be present in `env_col`.
#'
#' @return A numeric value representing the calculated TPS for the specified trait between the native and transplanted environments.
#'
#' @details The TPS is calculated by taking the difference between the trait values in the transplanted and native environments,
#' divided by the trait values in the native environment. This index provides a relative measure of the trait's plasticity after transplantation.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait = c(50, 60, 70, 80),   # Trait values
#'   Environment = factor(c(1, 1, 2, 2))  # Native = 1, Transplanted = 2
#' )
#'
#' # Calculate TPS between native (1) and transplanted (2)
#' tps_value = calculate_TPS(synthetic_data, trait_col = "Trait", env_col = "Environment", native_env = 1, transplanted_env = 2)
#' print(tps_value)
#'
#' @export
calculate_TPS = function(data, env_col, trait_col, native_env, transplanted_env) {
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
  
  # Filter the data for the native and transplanted environments
  native_data = trait_col[env_col == native_env]
  transplanted_data = trait_col[env_col == transplanted_env]
  
  # Ensure there are equal numbers of values for native and transplanted environments
  if (length(native_data) != length(transplanted_data)) {
    stop("The number of trait values for the native and transplanted environments must be equal.")
  }
  
  # Calculate TPS for each pair of values between the native and transplanted environments
  TPS_values = (transplanted_data - native_data) / native_data
  
  # Return the average TPS across all data points
  TPS = mean(TPS_values, na.rm = TRUE)
  
  return(TPS)
}

## test - passed on synthetic dataset


calculate_TPS(df_test6,1,3,1,2)

###########################
