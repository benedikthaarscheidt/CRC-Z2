#this files contains the following indices/functions: 
#generate_synthetic_data,
#Phenotypic Stability Index (calculate_PSI),
#Relative Plasticity Index (calculate_RPI) - tested,
#Plasticity Quotient (calculate_PQ) - tested,
#Phenotypic Range (PR) (calculate_PR) - tested,
#Norm of reaction width (calculate_NRW) - tested,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
#Calculate Plasticity Differential (PD) (calculate_PD) - tested,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,

source("~/CRC 1622 - Z2/R-files/2.R")



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

#' Calculate the Phenotypic Stability Index (PSI) for Multiple Traits Using Linear Regression
#'
#' This function calculates the Phenotypic Stability Index (PSI) for multiple traits across different environments.
#' The PSI is derived from the regression coefficient obtained by regressing the trait values against the environmental factor.
#' A lower regression coefficient indicates higher phenotypic stability.
#'
#' The environmental factors can be provided either as a column in the data frame or as an external vector.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @param env_col The column number or name representing the environmental conditions within the data frame, or NULL if using an external vector.
#' @param external_env_factor (Optional) A numeric vector representing external environmental factors. If provided, this will override `env_col`.
#' @return A data frame with the PSI values for each trait.
#' @examples
#' df = data.frame(
#'   Environment = rep(1:3, each = 10),
#'   Trait_1 = c(rnorm(10, 100, 10), rnorm(10, 110, 15), rnorm(10, 120, 20)),
#'   Trait_2 = c(rnorm(10, 200, 8), rnorm(10, 195, 12), rnorm(10, 205, 10))
#' )
#' 
#' # Calculate PSI for multiple traits using the internal environment column
#' PSI_results = calculate_PSI(df, trait_cols = c("Trait_1", "Trait_2"), env_col = "Environment")
#' print(PSI_results)
#' 
#' # Calculate PSI for multiple traits using an external environmental factor vector
#' external_env = rep(1:3, each = 10)
#' PSI_results_external = calculate_PSI(df, trait_cols = c("Trait_1", "Trait_2"), env_col = NULL, external_env_factor = external_env)
#' print(PSI_results_external)
#' @export
calculate_PSI = function(data, trait_cols, env_col = NULL, external_env_factor = NULL) {
  
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
  
  # Initialize a list to store PSI values for each trait
  PSI_values_list = list()
  
  # Loop over each trait and calculate PSI
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Extract the trait values
    trait_values = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Perform a linear regression of trait values on environment values
    model = lm(trait_values ~ environment_values)
    
    # Extract the regression coefficient (slope)
    stability_coefficient = coef(model)[2]
    
    # Calculate the stability score
    stability_score = 1 / (1 + abs(stability_coefficient))
    
    # Store the PSI for the current trait
    PSI_values_list[[i]] = stability_score
    names(PSI_values_list)[i] = trait_column
  }
  
  return(PSI_values_list)
}



#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#PSI = calculate_PSI(synthetic_data1, trait_cols = 5, external_env_factor = external_light)
#print(PSI)

##################################



#' Calculate Relative Plasticity Index (RPI) for all environment combinations across multiple traits
#' 
#' This function calculates the Relative Plasticity Index (RPI) for multiple traits across all pairs of environments
#' if no specific environments are provided. Otherwise, RPI is calculated between two specified environments.
#' 
#' @param data A data frame containing the input data, with columns for the traits and another column or external vector for the environment.
#' @param trait_cols A vector of strings or numeric values indicating the columns in `data` that contain the trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param env1 A value indicating the first environment to compare. Must be present in `env_col`. Defaults to `NULL` for calculating all combinations.
#' @param env2 A value indicating the second environment to compare. Must be present in `env_col`. Defaults to `NULL` for calculating all combinations.
#' 
#' @return If env1 and env2 are provided, a data frame with the RPI values for each trait between the two environments.
#' If env1 and env2 are not provided, a data frame of RPI values for all possible environment pairs for each trait.
#' 
#' @details The RPI is calculated by taking the absolute difference between the trait values in the two environments
#' and dividing it by the sum of the trait values in the two environments. This is done for all specified traits.
#' 
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_1 = c(100, 110, 120, 130, 140, 150),
#'   Trait_2 = c(200, 210, 220, 230, 240, 250),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3))
#' )
#' 
#' # Calculate RPI for all environment combinations for multiple traits
#' rpi_values = calculate_RPI(synthetic_data, trait_cols = c("Trait_1", "Trait_2"), env_col = "Environment")
#' print(rpi_values)
#' 
#' @export
calculate_RPI = function(data, trait_cols, env_col, env1 = NULL, env2 = NULL) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Ensure trait_cols are either numeric or character and fetch the respective columns
  trait_cols = if (is.numeric(trait_cols)) names(data)[trait_cols] else trait_cols
  
  # Initialize a data frame to store results
  RPI_results_all = data.frame()
  
  # Loop through each trait
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_data = data[[trait_col]]
    
    # If no specific environments are provided, calculate RPI for all pairs of environments
    if (is.null(env1) && is.null(env2)) {
      
      # Get unique environments
      unique_envs = unique(env_col)
      env_combinations = combn(unique_envs, 2, simplify = TRUE)  # All pairs of environments
      
      # Initialize a data frame to store RPI values for all environment pairs
      RPI_results = data.frame(env1 = character(), env2 = character(), Trait = character(), RPI = numeric())
      
      # Loop through all environment pairs
      for (i in 1:ncol(env_combinations)) {
        env1_data = trait_data[env_col == env_combinations[1, i]]
        env2_data = trait_data[env_col == env_combinations[2, i]]
        
        # Ensure there are equal numbers of values for env1 and env2
        if (length(env1_data) != length(env2_data)) {
          stop("The number of trait values for the two environments must be equal.")
        }
        
        # Calculate RPI for each pair of values between the two environments
        RPI_values = abs(env1_data - env2_data) / (env1_data + env2_data)
        
        # Store the average RPI for this environment pair
        RPI_results = rbind(RPI_results, data.frame(
          env1 = env_combinations[1, i],
          env2 = env_combinations[2, i],
          Trait = trait_col,
          RPI = mean(RPI_values, na.rm = TRUE)
        ))
      }
      
      # Append to overall results
      RPI_results_all = rbind(RPI_results_all, RPI_results)
      
    } else {
      
      # If specific environments are provided, filter the data for the two specified environments
      env1_data = trait_data[env_col == env1]
      env2_data = trait_data[env_col == env2]
      
      # Ensure there are equal numbers of values for env1 and env2
      if (length(env1_data) != length(env2_data)) {
        stop("The number of trait values for the two environments must be equal.")
      }
      
      # Calculate RPI for each pair of values between the two environments
      RPI_values = abs(env1_data - env2_data) / (env1_data + env2_data)
      
      # Store the average RPI for this trait between the two environments
      RPI_results_all = rbind(RPI_results_all, data.frame(
        env1 = env1,
        env2 = env2,
        Trait = trait_col,
        RPI = mean(RPI_values, na.rm = TRUE)
      ))
    }
  }
  
  # Return the combined RPI results for all traits
  return(RPI_results_all)
}


## test - passed on synthetic dataset

#calculate_RPI(df_test2,trait_col=c(2,3),env_col=1)

######################################

#' @title Calculate Plasticity Quotient (PQ) for Multiple Traits
#' @description This function calculates the Plasticity Quotient (PQ) based on the range of trait values across different environmental conditions and the corresponding environmental factor range for multiple traits.
#' 
#' @param data A data frame containing the data for the analysis.
#' @param trait_cols A vector of column names (strings) or numeric indices representing the trait values in the data. Each will be analyzed separately.
#' @param env_col A column name (string) or numeric index representing the environmental variable in the data, or a vector of environmental values. This specifies the environmental conditions for calculating plasticity.
#' @param factor_col A column name (string) or numeric index representing a relevant environmental factor (e.g., temperature, moisture) in the data, or a vector of environmental factor values. This will be used for standardizing the trait range.
#' 
#' @details 
#' The function first handles the provided column names or indices for `trait_cols`, `env_col`, and `factor_col`, ensuring they are correctly extracted from the data. 
#' 
#' It then calculates the mean trait values and factor values across the different environmental conditions specified in `env_col`. The range of the mean trait values is computed, and this is divided by the range of the mean factor values to compute the Plasticity Quotient (PQ) for each trait.
#' 
#' @return A data frame where each row corresponds to a trait and the calculated Plasticity Quotient (PQ) for that trait.
#' 
#' @examples
#' # Example usage:
#' # Assuming 'data' is a data frame with columns 'Trait_1', 'Trait_2', 'env', and 'factor':
#' PQ = calculate_PQ(data, c('Trait_1', 'Trait_2'), 'env', 'factor')
#' print(PQ)
#' 
#' @export
calculate_PQ = function(data, env_col, trait_cols, factor_col) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col_data = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col_data = env_col
  } else {
    env_col_data = data[[env_col]]
  }
  
  # Handle factor_col
  if (is.numeric(factor_col) && length(factor_col) == 1) {
    factor_data = data[[factor_col]]
  } else if (is.vector(factor_col) && length(factor_col) == nrow(data)) {
    factor_data = factor_col
  } else {
    factor_data = data[[factor_col]]
  }
  
  # Initialize a list to store PQ values for each trait
  PQ_values_list = list()
  
  # Loop over each trait in trait_cols
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Extract trait data
    trait_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Calculate mean values for the current trait by environment and factor levels
    mean_values = aggregate(trait_data ~ env_col_data, data = data, FUN = mean)
    range_trait = abs(max(mean_values$trait_data) - min(mean_values$trait_data))
    env_range = abs(max(factor_data) - min(factor_data))
    
    # Calculate PQ for the current trait
    PQ_value = range_trait / env_range
    
    # Store the PQ value in the list
    PQ_values_list[[i]] = PQ_value
    names(PQ_values_list)[i] = trait_column
  }
  
  # Return the PQ values list
  return(PQ_values_list)
}

## test - passed on synthetic dataframe

#calculate_PQ(df_test6,1,3,2)

###########################################

#' @description This function calculates the Phenotypic Range (PR) for multiple traits across specified environments. PR is defined as the difference between the maximum and minimum trait values within each environment.
#' 
#' @param data A data frame containing the data for the analysis.
#' @param env_col A column name (string), numeric index, or vector representing the environmental variable in the data. This specifies the environment for calculating the PR.
#' @param trait_cols A vector of column names (strings) or numeric indices representing the trait values in the data. Each will be analyzed separately.
#' @param env Optional. A vector specifying which environments to include in the calculation. If not provided, all unique environments from `env_col` are used.
#' @param across Optional. A logical value indicating whether to calculate the PR across all environments combined. If `TRUE`, calculates the PR across all environments for each trait.
#' 
#' @details 
#' The function calculates the Phenotypic Range (PR) for each specified environment, defined as the range (maximum minus minimum) of trait values within each environment. If `env` is provided, only the specified environments are used; otherwise, it uses all unique environments in `env_col`. 
#' 
#' @return A data frame where each row corresponds to a trait and contains the PR for each environment.
#' 
#' @examples
#' # Example usage:
#' # Assuming 'data' is a data frame with columns 'Trait_1', 'Trait_2', and 'env':
#' PR = calculate_PR(data, 'env', c('Trait_1', 'Trait_2'))
#' 
#' # Example with specific environments:
#' PR = calculate_PR(data, 'env', c('Trait_1', 'Trait_2'), env = c("Environment1", "Environment2"))
#' 
#' @export
calculate_PR = function(data, env_col, trait_cols, env = NULL, across = FALSE) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a data frame to store the PR for each trait and environment
  PR_results = data.frame(Trait = character(), Environment = character(), PR = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait in trait_cols
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # If across = TRUE, calculate PR across all environments combined
    if (across == TRUE) {
      # Calculate mean trait values for each environment
      means = tapply(trait_data, env_col, mean, na.rm = TRUE)
      # Identify maximum and minimum means
      max_mean = max(means, na.rm = TRUE)
      min_mean = min(means, na.rm = TRUE)
      PR_across = abs(max_mean - min_mean)
      PR_results = rbind(PR_results, data.frame(Trait = colnames(data)[trait_col], Environment = "Across", PR = PR_across))
      
    } else {
      # If across = FALSE, calculate PR for each environment
      
      if (is.null(env)) {
        env = unique(env_col)
      }
      
      # Loop through each environment and calculate the PR
      for (e in env) {
        env_data = trait_data[env_col == e]
        PR_value = max(env_data, na.rm = TRUE) - min(env_data, na.rm = TRUE)
        PR_results = rbind(PR_results, data.frame(Trait = colnames(data)[trait_col], Environment = e, PR = PR_value))
      }
    }
  }
  
  return(PR_results)
}


## test - passed on synthetic dataset 


#calculate_PR(df_test6,1,3)
#calculate_PR(df_test6,1,trait_cols=c(2,3),across=T)






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
calculate_NRW = function(data, env_col = NULL, trait_cols) {
  
  # Handle env_col
  if (is.null(env_col)) {
    env_col_data = rep(1, nrow(data))  # Treat as one environment if not provided
  } else if (is.numeric(env_col)) {
    env_col_data = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col_data = env_col
  } else {
    env_col_data = data[[env_col]]
  }
  
  # Initialize a vector to store NRW values for each trait
  NRW_values = numeric(length(trait_cols))
  names(NRW_values) = trait_cols  # Name the vector by trait columns
  
  # Loop through each trait in trait_cols
  for (i in seq_along(trait_cols)) {
    trait_col = trait_cols[i]
    
    # Extract trait data
    trait_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Calculate the range (NRW) for each environment or group
    NRW_df = aggregate(trait_data ~ env_col_data, data = data, FUN = function(x) max(x) - min(x))
    
    # Calculate the mean NRW across all groups for the current trait
    NRW_values[i] = mean(NRW_df$trait_data, na.rm = TRUE)
  }
  
  return(NRW_values)
}

## test - passed on synthetic dataset


#calculate_NRW(df_test6,trait_cols=3)

##############################


#' @title Calculate Environmental Sensitivity Performance (ESP) for Multiple Traits
#' @description This function calculates the Environmental Sensitivity Performance (ESP) for multiple traits across specified environments. 
#' ESP is calculated as the relative change in trait values within each environment compared to the overall mean trait value across all environments.
#' 
#' @param data A data frame containing the data for the analysis.
#' @param env_col A column name (string), numeric index, or vector representing the environmental variable in the data. This specifies the environment for calculating ESP.
#' @param trait_cols A vector of column names (strings) or numeric indices representing the trait values in the data. This can also be a vector of trait values with the same length as the number of rows in the data.
#' @param env Optional. A vector specifying which environments to include in the calculation. If not provided, all unique environments from `env_col` are used.
#' 
#' @details 
#' The function calculates the Environmental Sensitivity Performance (ESP) for each trait specified in `trait_cols` by determining the difference between the mean trait value for each environment and the overall mean trait value across all environments, divided by the overall mean. 
#' This results in a measure of how sensitive the trait is to environmental conditions, with values closer to zero indicating low sensitivity and larger values indicating higher sensitivity.
#' 
#' @return A data frame containing the ESP values for each trait across each environment. The data frame includes the following columns:
#' \item{Trait}{The name of the trait for which ESP was calculated.}
#' \item{Environment}{The environment in which the trait's ESP was calculated.}
#' \item{ESP}{The calculated Environmental Sensitivity Performance value.}
#' 
#' @examples 
#' # Example usage:
#' df_test <- data.frame(
#'   Environment = rep(1:3, each = 10),
#'   Trait_1 = c(rnorm(10, 100, 10), rnorm(10, 110, 15), rnorm(10, 120, 20)),
#'   Trait_2 = c(rnorm(10, 200, 8), rnorm(10, 195, 12), rnorm(10, 205, 10))
#' )
#' 
#' # Calculate ESP for multiple traits
#' ESP_values <- calculate_ESP(df_test, env_col = "Environment", trait_cols = c(2, 3))
#' print(ESP_values)
#' 
#' @export
calculate_ESP = function(data, env_col, trait_cols, env = NULL) {
  
  # Handle env_col
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
  
  # Initialize a data frame to store the ESP results for each trait
  ESP_results = data.frame(Trait = character(), Environment = character(), ESP = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait in trait_cols
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Calculate the mean trait value across all environments
    mean_trait_all_envs = mean(trait_data, na.rm = TRUE)
    
    # Loop through each environment and calculate ESP
    for (i in 1:length(env)) {
      trait_values_env = trait_data[env_col == env[i]]
      ESP_value = (mean(trait_values_env, na.rm = TRUE) - mean_trait_all_envs) / mean_trait_all_envs
      
      # Store the results in the data frame
      ESP_results = rbind(ESP_results, data.frame(
        Trait = colnames(data)[trait_col],
        Environment = env[i],
        ESP = ESP_value
      ))
    }
  }
  
  return(ESP_results)
}


## test - passed on synthetic dataset 

#calculate_ESP(df_test6,1,trait_cols=c(2,3),env=NULL)

######################### 01.10.2024

#' @title Calculate Plasticity Differential (PD) for Multiple Traits
#' @description This function calculates the Plasticity Differential (PD) for multiple traits across control and stress conditions.
#' PD quantifies the absolute difference in trait values between control and stress environments using the formula:
#' \deqn{PD = |TraitValueStress - TraitValueControl|}
#'
#' @param data A data frame containing the input data, with one or more columns for traits and another column or external vector for the environment.
#' @param trait_cols A vector of strings or numeric values indicating the columns in `data` that contain the trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param control_env A value indicating the control environment to compare. Must be present in `env_col`.
#' @param stress_env A value indicating the stress environment to compare. Must be present in `env_col`.
#'
#' @return A data frame with columns for each trait and the calculated PD between the control and stress environments.
#'
#' @details The PD is calculated for each trait by taking the absolute difference between the trait values in the stress and control environments.
#' This index provides a measure of the change in trait value under stress relative to the control for each trait.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_1 = c(100, 110, 120, 130, 140, 150),
#'   Trait_2 = c(200, 210, 220, 230, 240, 250),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3)) # Assume 1=Control, 2=Stress
#' )
#' 
#' # Calculate PD between control (1) and stress (2) for multiple traits
#' pd_values = calculate_PD(synthetic_data, trait_cols = c("Trait_1", "Trait_2"), env_col = "Environment", control_env = 1, stress_env = 2)
#' print(pd_values)
#'
#' @export
calculate_PD = function(data, env_col, trait_cols, control_env, stress_env) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a data frame to store the PD results for each trait
  PD_results = data.frame(Trait = character(), PD = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait column
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Filter the data for the control and stress environments
    control_data = trait_data[env_col == control_env]
    stress_data = trait_data[env_col == stress_env]
    
    # Ensure there are equal numbers of values for control and stress
    if (length(control_data) != length(stress_data)) {
      stop("The number of trait values for the control and stress environments must be equal.")
    }
    
    # Calculate PD for each pair of values between the control and stress environments
    PD_values = abs(stress_data - control_data)
    
    # Calculate and store the average PD for this trait
    PD = mean(PD_values, na.rm = TRUE)
    
    # Append the result to the data frame
    PD_results = rbind(PD_results, data.frame(Trait = colnames(data)[trait_col], PD = PD))
  }
  
  return(PD_results)
}

## test - passed on synthetic dataset 

#calculate_PD(df_test6,1,c(2,3),3,1)
##########################


#' @title Calculate Fitness Plasticity Index (FPI) for Multiple Traits
#' @description This function calculates the Fitness Plasticity Index (FPI) for multiple fitness-related traits
#' across control and stress conditions. The FPI quantifies the relative change in fitness between control
#' and stress environments for each trait.
#'
#' @param data A data frame containing the input data, with multiple columns for fitness traits and another column or external vector for the environment.
#' @param trait_cols A vector of strings or numeric values indicating the columns in `data` that contain the fitness trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param control_env A value indicating the control environment to compare. Must be present in `env_col`.
#' @param stress_env A value indicating the stress environment to compare. Must be present in `env_col`.
#'
#' @return A data frame with each trait and the calculated FPI for each trait between the control and stress environments.
#'
#' @details The FPI is calculated for each trait by taking the difference between the fitness values in the stress and control environments,
#' divided by the fitness values in the control environment. This index provides a relative measure of each fitness trait's
#' plasticity across the two environments.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Fitness_1 = c(10, 12, 14, 16),
#'   Fitness_2 = c(20, 22, 24, 26),
#'   Environment = factor(c(1, 1, 2, 2))  # Control = 1, Stress = 2
#' )
#' 
#' # Calculate FPI between control (1) and stress (2) for multiple traits
#' fpi_values = calculate_FPI(synthetic_data, trait_cols = c("Fitness_1", "Fitness_2"), env_col = "Environment", control_env = 1, stress_env = 2)
#' print(fpi_values)
#'
#' @export
calculate_FPI = function(data, env_col, trait_cols, control_env, stress_env) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a data frame to store FPI results for each trait
  FPI_results = data.frame(Trait = character(), FPI = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait column
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Filter the data for the control and stress environments
    control_data = trait_data[env_col == control_env]
    stress_data = trait_data[env_col == stress_env]
    
    # Ensure there are equal numbers of values for control and stress
    if (length(control_data) != length(stress_data)) {
      stop("The number of trait values for the control and stress environments must be equal.")
    }
    
    # Calculate FPI for each pair of values between the control and stress environments
    FPI_values = (stress_data - control_data) / control_data
    
    # Calculate and store the average FPI for this trait
    FPI = mean(FPI_values, na.rm = TRUE)
    
    # Append the result to the data frame
    FPI_results = rbind(FPI_results, data.frame(Trait = colnames(data)[trait_col], FPI = FPI))
  }
  
  return(FPI_results)
}


## test - passed on synthetic dataset

#calculate_FPI(df_test6,1,c(2,3),1,2)


#################################


#' @title Calculate Transplant Plasticity Score (TPS) for Multiple Traits
#' @description This function calculates the Transplant Plasticity Score (TPS) for multiple traits when an organism is transplanted to a different environment.
#' The TPS quantifies the relative change in a trait between the native environment and the transplanted environment.
#'
#' @param data A data frame containing the input data, with multiple columns for trait values and another column or external vector for the environment.
#' @param trait_cols A vector of strings or numeric values indicating the columns in `data` that contain the trait values.
#' @param env_col A string, numeric, or vector indicating the column in `data` that contains the environment labels, or an external vector specifying environment belonging.
#' @param native_env A value indicating the native environment. Must be present in `env_col`.
#' @param transplanted_env A value indicating the transplanted environment. Must be present in `env_col`.
#'
#' @return A data frame with each trait and the calculated TPS for each trait between the native and transplanted environments.
#'
#' @details The TPS is calculated for each trait by taking the difference between the trait values in the transplanted and native environments,
#' divided by the trait values in the native environment. This index provides a relative measure of each trait's plasticity after transplantation.
#'
#' @examples
#' # Example usage with synthetic data
#' synthetic_data = data.frame(
#'   Trait_1 = c(50, 60, 70, 80),
#'   Trait_2 = c(30, 40, 50, 60),
#'   Environment = factor(c(1, 1, 2, 2))  # Native = 1, Transplanted = 2
#' )
#'
#' # Calculate TPS between native (1) and transplanted (2) for multiple traits
#' tps_values = calculate_TPS(synthetic_data, trait_cols = c("Trait_1", "Trait_2"), env_col = "Environment", native_env = 1, transplanted_env = 2)
#' print(tps_values)
#'
#' @export
calculate_TPS = function(data, env_col, trait_cols, native_env, transplanted_env) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a data frame to store TPS results for each trait
  TPS_results = data.frame(Trait = character(), TPS = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each trait column
  for (trait_col in trait_cols) {
    
    # Extract trait data
    trait_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    
    # Filter the data for the native and transplanted environments
    native_data = trait_data[env_col == native_env]
    transplanted_data = trait_data[env_col == transplanted_env]
    
    # Ensure there are equal numbers of values for native and transplanted environments
    if (length(native_data) != length(transplanted_data)) {
      stop("The number of trait values for the native and transplanted environments must be equal.")
    }
    
    # Calculate TPS for each pair of values between the native and transplanted environments
    TPS_values = (transplanted_data - native_data) / native_data
    
    # Calculate and store the average TPS for this trait
    TPS = mean(TPS_values, na.rm = TRUE)
    
    # Append the result to the data frame
    TPS_results = rbind(TPS_results, data.frame(Trait = colnames(data)[trait_col], TPS = TPS))
  }
  
  return(TPS_results)
}

## test - passed on synthetic dataset


#calculate_TPS(df_test6,1,c(2,3),1,2)

###########################
