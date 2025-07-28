#this files contains the following indices/functions: 
#generate_synthetic_data,
#check_and_install_packages,
#RDPI	(rdpi_calculation) - tested,
#RDPIs (rdpi_mean_calculation) - tested,
#ESPI (calculate_ESPI) - tested,
#ESPIid (espiid_calculation) - tested,
#evwpi_calculation (idea from Benedikt)


##### datasets for testing

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



if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2")
  library(roxygen2)
}
############################


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



#' @title Calculate Relative Distance Plasticity Index (RDPI)
#' 
#' @description
#' Calculates the Relative Distance Plasticity Index (RDPI) for trait values across different 
#' environmental conditions. The calculation follows Valladares et al. (2006) methodology.
#' 
#' @param trait_values Numeric vector containing trait measurements (already averaged across replicates)
#' @param env_values Optional vector of environmental conditions. If NULL, equidistant steps are assumed
#' 
#' @return Numeric value representing the RDPI score
#' 
#' @examples
#' # With equidistant environments
#' traits = c(10, 12, 15, 18, 20)
#' rdpi = calculate_rdpi(traits)
#' 
#' # With specified environments
#' traits = c(10, 12, 15, 18, 20)
#' envs = c(1, 2, 4, 6, 8)  # Non-equidistant environments
#' rdpi = calculate_rdpi(traits, envs)
#'
#' @export
calculate_rdpi = function(trait_values, env_values = NULL) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be numeric")
  }
  
  n_envs = length(trait_values)
  
  # If no environment values provided, create sequential environments
  if (is.null(env_values)) {
    env_values = seq_len(n_envs)
  }
  
  # Ensure same length
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }
  
  # Get all pairs of environment indices
  env_pairs = combn(n_envs, 2)
  n_pairs = ncol(env_pairs)
  
  # Calculate RDPIs for each pair
  rdpis = numeric(n_pairs)
  
  for (i in seq_len(n_pairs)) {
    idx1 = env_pairs[1, i]
    idx2 = env_pairs[2, i]
    
    # Calculate relative distance for this pair
    abs_diff = abs(trait_values[idx2] - trait_values[idx1])
    sum_vals = trait_values[idx1] + trait_values[idx2]
    
    # Handle potential division by zero
    if (sum_vals == 0) {
      rdpis[i] = 0
    } else {
      rdpis[i] = abs_diff / sum_vals
    }
  }
  
  # Calculate final RDPI as mean of all RDPIs
  rdpi = mean(rdpis, na.rm = TRUE)
  
  return(rdpi)
}


## Example usage with synthetic data
#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#external_water = sample(rep(c("Low", "High"), each = 150))
#
## Calculate RDPI with factors not in the dataframe
#RDPI= rdpi_calculation(synthetic_data1, sp = NULL, trait_cols = c(3, 4), factors = 2, factors_not_in_dataframe = list(external_light, external_water))
#
#
#### test - passed on synthetic dataset and crosschecked with data from package
#
#
#df_test20 = data.frame(
#  Column0 = as.factor(c(rep(1, 20), rep(2, 20))),  # Two groups/species, 20 observations each
#  Column1 = as.factor(c(rep(1, 10), rep(2, 10), rep(1, 10), rep(2, 10))),  # Two environments per group
#  Column2 = c(rep(2, 10), rep(1, 10), rep(2, 10), rep(1, 10))  # Numeric trait values
#)
##rdpi_calculation(df_test20,sp=1, trait_cols=3,factors = 2)


#####################################################





#' Calculate the Relative Distance Plasticity Index for Mean Phenotypic Values (RDPIs)
#'
#' This function calculates the Relative Distance Plasticity Index (RDPI) for mean phenotypic values. 
#' RDPI is defined as the absolute phenotypic distance between the means of the same genotype 
#' in different environments, divided by the smaller of the two mean phenotypic values.
#' 
#' The function is designed to handle multiple traits, genotypes, and environmental combinations, 
#' providing a comprehensive analysis of phenotypic plasticity across environmental gradients.
#'
#' @param dataframe A data frame containing the phenotypic data, including traits, genotype identifiers, 
#' and environmental factors.
#' @param trait_cols A vector specifying the column indices or names of the traits to be analyzed.
#' @param sp (Optional) A column index or name indicating the species or genotype identifier. If not provided, 
#' the data is treated as a single group.
#' @param factors (Optional) A vector of column indices or names specifying the internal environmental factors 
#' within the data frame that should be combined to create unique environmental combinations.
#' @param factors_not_in_dataframe (Optional) A list of vectors representing external environmental factors 
#' that are not included in the data frame. Each vector must have a length equal to the number of rows in the data frame.
#' 
#' @return A list containing:
#' \item{all_results}{A list where each entry corresponds to a genotype (or species) and contains the RDPI values 
#' and trait data for the specified traits.}
#' \item{trait_boxplots}{A ggplot object representing boxplots of trait values across environmental combinations.}
#' \item{anova_results}{A list of ANOVA summaries for each trait across the environmental combinations.}
#' \item{tukey_results}{A list of Tukey's HSD test results for each trait, assessing significant differences 
#' between environmental combinations.}
#'
#' @details 
#' The RDPI provides a measure of phenotypic plasticity by comparing the mean phenotypic values of the same genotype 
#' across different environments. This metric is particularly useful for datasets with multiple environments and 
#' genotypes, where it can reveal the degree of plasticity in response to environmental changes.
#'
#' The function first combines internal and external environmental factors to create unique environmental 
#' combinations. For each trait, it then calculates the absolute differences between the means of the trait across 
#' all pairs of environments for the same genotype. These differences are normalized by dividing by the smaller of 
#' the two means. The final RDPI value for each genotype is the average of these normalized differences across 
#' all environment pairs.
#'
#' ANOVA is performed on the trait data to assess the significance of the environmental factors, and Tukey's HSD test 
#' is used to determine which specific environmental combinations differ significantly. The function also generates 
#' boxplots to visualize the distribution of trait values across the different environmental combinations.
#' 
#' @examples
#' # Example usage
#' df = data.frame(
#'   Genotype = rep(c("G1", "G2"), each = 6),
#'   Environment = rep(c("Env1", "Env2", "Env3"), times = 4),
#'   Height = c(10, 12, 11, 14, 16, 15, 20, 22, 21, 23, 25, 24)
#' )
#' 
#' # Calculate RDPIs for the 'Height' trait across genotypes and environments
#' results = rdpi_mean_calculation(
#'   dataframe = df, 
#'   trait_cols = "Height", 
#'   sp = "Genotype", 
#'   factors = "Environment"
#' )
#' 
#' # Access the RDPI values and boxplots
#' print(results$all_results)
#' print(results$trait_boxplots)
#' 
#' @export
rdpi_mean_calculation = function(dataframe, trait_cols, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL, stat_analysis = NULL) {
  
  # List of required packages
  required_packages = c("ggplot2", "agricolae", "dplyr", "reshape2")
  
  # Check and install required packages
  check_and_install_packages(required_packages)
  
  # Convert column indices to names if necessary
  if (!is.null(sp)) {
    sp = if (is.numeric(sp)) names(dataframe)[sp] else sp
  }
  traits = if (is.numeric(trait_cols)) names(dataframe)[trait_cols] else trait_cols
  
  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    external_factors_df = as.data.frame(factors_not_in_dataframe)
    if (!is.null(factors)) {
      factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
      combined_factors_df = cbind(dataframe[factors], external_factors_df)
    } else {
      combined_factors_df = external_factors_df
    }
    dataframe$Combined_Factors = interaction(combined_factors_df, drop = TRUE)
  } else if (!is.null(factors)) {
    factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
    dataframe$Combined_Factors = interaction(dataframe[factors], drop = TRUE)
  } else {
    stop("You must provide either internal factors, external factors, or both.")
  }
  
  dataframe$Combined_Factors = as.factor(dataframe$Combined_Factors)
  
  all_rdpi_results = data.frame()  # Initialize a dataframe to store all RDPI values
  
  if (is.null(sp)) {
    unique_species = list("Single_Group" = dataframe)
  } else {
    unique_species = split(dataframe, dataframe[[sp]])
  }
  
  for (species_name in names(unique_species)) {
    species_data = unique_species[[species_name]]
    
    for (trait in traits) {
      RDPI_values = data.frame(sp = character(), env1 = character(), env2 = character(), rdpi = numeric())
      
      # Get unique environment combinations
      env_levels = levels(species_data$Combined_Factors)
      n_env_levels = length(env_levels)
      
      mean_values = aggregate(species_data[[trait]], list(species_data$Combined_Factors), mean)
      colnames(mean_values) = c("Env_Combination", "Mean_Trait")
      
      for (i in 1:(n_env_levels - 1)) {
        for (j in (i + 1):n_env_levels) {
          mean_i = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[i]]
          mean_j = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[j]]
          
          # Calculate RDPI for this trait between the two environments
          rdpi_value = abs(mean_i - mean_j) / min(mean_i, mean_j)
          
          # Append the RDPI value for this species, trait, and environment pair
          RDPI_values = rbind(RDPI_values, data.frame(sp = species_name, env1 = env_levels[i], env2 = env_levels[j], rdpi = rdpi_value))
        }
      }
      
      all_rdpi_results = rbind(all_rdpi_results, RDPI_values)
    }
  }
  
  # If statistical analysis is requested, perform ANOVA and Tukey's HSD
  if (!is.null(stat_analysis)) {
    all_trait_data = data.frame()
    
    for (species_name in names(unique_species)) {
      species_data = unique_species[[species_name]]
      
      for (trait in traits) {
        trait_data = data.frame(
          sp = species_name,
          trait = trait,
          Combined_Factors = species_data$Combined_Factors,
          Trait_Value = species_data[[trait]]
        )
        
        all_trait_data = rbind(all_trait_data, trait_data)
      }
    }
    
    # Perform ANOVA and Tukey's HSD test
    anova_results = list()
    tukey_results = list()
    
    for (trait in traits) {
      fit = aov(Trait_Value ~ Combined_Factors, data = subset(all_trait_data, trait == trait))
      anova_results[[trait]] = summary(fit)
      
      # Perform Tukey's HSD test
      Tukey = agricolae::HSD.test(fit, trt = "Combined_Factors")
      tukey_results[[trait]] = Tukey
    }
    
    # Create boxplots for the traits across environmental combinations
    boxplot_traits = ggplot2::ggplot(all_trait_data, ggplot2::aes(x = Combined_Factors, y = Trait_Value, fill = trait)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~trait, scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::xlab("Environmental Combinations") +
      ggplot2::ylab("Trait Values") +
      ggplot2::ggtitle("Boxplots of Trait Values Across Environmental Combinations") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    return(list(
      rdpi_results = all_rdpi_results,
      trait_boxplots = boxplot_traits,
      anova_results = anova_results,
      tukey_results = tukey_results
    ))
  }
  
  # If stat_analysis is NULL, return only the RDPI values
  return(all_rdpi_results)
}

  


## Example usage with synthetic data
#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#external_water =rep(c("Low", "High"), each = 150)
#
### test - passed on synthetic dataset
#
#rdpi_mean_calculation(df_test20,sp=1, trait_cols=3,factors = 2)


#########################################


#' Calculate Environmental Sensitivity Performance Index (ESPI)
#'
#' This function calculates the Environmental Sensitivity Performance Index (ESPI)
#' for one or more traits across different environmental conditions. ESPI is calculated
#' as:
#' \deqn{ESPI = \frac{(\text{Maximum mean} - \text{Minimum mean})}{\text{Absolute distance between environmental values}}}
#'
#' @param trait_values A numeric vector for a single trait or a data frame with one column per trait.
#' @param env_values (Optional) A numeric vector of environmental values corresponding to `trait_values`. 
#' If `NULL`, equidistant environments are generated automatically.
#'
#' @return A numeric value for a single trait or a named numeric vector for multiple traits, with the ESPI
#' value for each trait.
#'
#' @details
#' The ESPI measures how sensitive a trait is to changes in environmental conditions. It is based on the
#' maximum and minimum mean trait values across different environments, normalized by the absolute range
#' of the environmental values.
#'
#' @examples
#' # Example 1: Single Trait
#' trait_values = c(10, 15, 20, 25, 30)
#' env_values = c(1, 2, 3, 4, 5)
#' calculate_ESPI(trait_values, env_values)
#'
#' # Example 2: Multiple Traits
#' trait_values = data.frame(
#'   Height = c(10, 15, 20, 25, 30),
#'   Weight = c(5, 7, 9, 11, 13)
#' )
#' env_values = c(1, 2, 3, 4, 5)
#' calculate_ESPI(trait_values, env_values)
#'
#' # Example 3: Without env_values
#' trait_values = c(10, 15, 20, 25, 30)
#' calculate_ESPI(trait_values)
#'
#' @export
calculate_ESPI = function(trait_values, env_values = NULL) {
  
  if (is.null(env_values)) {
    env_values = seq_along(trait_values)
  
  }
  
  if (!is.vector(trait_values) && !is.data.frame(trait_values)) {
    stop("trait_values must be a numeric vector or a data frame where each column is a trait.")
  }
  
  if (length(env_values) != length(trait_values)) {
    stop("env_values must be the same length as trait_values.")
  }
  
  # Function to calculate ESPI for a single trait
  calculate_single_espi = function(single_trait) {
    means = tapply(single_trait, env_values, mean, na.rm = TRUE)
    max_mean = max(means, na.rm = TRUE)
    min_mean = min(means, na.rm = TRUE)
    
    abs_env_distance = abs(max(as.numeric(env_values), na.rm = TRUE) - min(as.numeric(env_values), na.rm = TRUE))
    if (abs_env_distance > 0) {
      return((max_mean - min_mean) / abs_env_distance)
    } else {
      return(NA)
    }
  }
  
  # Handle single trait
  if (is.vector(trait_values)) {
    return(calculate_single_espi(trait_values))
  } else if (is.data.frame(trait_values)) {
    # Handle multiple traits
    espi_results = sapply(trait_values, calculate_single_espi)
    names(espi_results) = colnames(trait_values)
    return(espi_results)
  }
}
## test - passed on synthetic dataset 


#df_test2.2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 5),rep(4,5) ,rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
#calculate_ESPI(df_test2.2,trait_cols = c(2,3),env_col = 1)


#####################################


#' Calculate Environmental Sensitivity Performance Index for Individual Differences (ESPIID)
#'
#' This function calculates the Environmental Sensitivity Performance Index for Individual Differences (ESPIID)
#' for a single trait across different environmental conditions. ESPIID is calculated as:
#' \deqn{ESPIID = \frac{\text{Mean or Median Absolute Phenotypic Distance}}{\text{Absolute Environmental Distance}}}
#'
#' @param trait_values A numeric vector containing the trait values.
#' @param env_values (Optional) A vector representing the environmental conditions for the trait values.
#' If `NULL`, equidistant environments are assumed.
#' @param use_median (Optional) A logical value indicating whether to use the median instead of the mean
#' for calculating the absolute phenotypic distances. Defaults to `FALSE` (use mean).
#'
#' @return A numeric value representing the ESPIID for the given trait.
#'
#' @details
#' The ESPIID quantifies the sensitivity of a trait to environmental changes by comparing
#' the phenotypic distances between individuals across pairs of environments. The phenotypic
#' distances are normalized by the absolute distance between the environmental values.
#'
#' The function assumes that the input trait values and environmental values are of equal length.
#' If environmental values are not provided, equidistant environments are generated automatically.
#'
#' @examples
#' # Example 1: Single Trait with Specified Environmental Values
#' trait_values = c(10, 12, 15, 18, 20, 25)
#' env_values = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3")
#' espiid = calculate_espiid(trait_values, env_values)
#' print(espiid)
#'
#' # Example 2: Single Trait with Equidistant Environments
#' trait_values = c(10, 12, 15, 18, 20, 25)
#' espiid = calculate_espiid(trait_values)
#' print(espiid)
#'
#' # Example 3: Use Median for ESPIID Calculation
#' trait_values = c(10, 12, 15, 18, 20, 25)
#' env_values = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3")
#' espiid = calculate_espiid(trait_values, env_values, use_median = TRUE)
#' print(espiid)
#'
#' @export
calculate_espiid = function(trait_values, use_median = FALSE) {


  env_values = seq_along(trait_values) # Create equidistant environment indices
  
  
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length.")
  }
  
  env_values = as.factor(env_values) # Convert to factor
  
  # Get all unique pairs of environments
  env_levels = levels(env_values)
  n_envs = length(env_levels)
  
  if (n_envs < 2) {
    stop("At least two environments are required to calculate ESPIID.")
  }
  
  # Initialize storage for ESPIID values
  espiid_values = numeric()
  
  # Loop over all environment pairs
  for (i in 1:(n_envs - 1)) {
    for (j in (i + 1):n_envs) {
      # Extract trait values for the two environments
      trait_i = trait_values[env_values == env_levels[i]]
      trait_j = trait_values[env_values == env_levels[j]]
      
      # Calculate absolute differences between all pairs of individuals
      abs_diff = abs(outer(trait_i, trait_j, "-"))
      
      # Calculate mean or median absolute phenotypic distance
      mean_abs_diff = if (use_median) median(abs_diff, na.rm = TRUE) else mean(abs_diff, na.rm = TRUE)
      
      # Calculate absolute distance between environmental values
      abs_env_distance = abs(as.numeric(env_levels[i]) - as.numeric(env_levels[j]))
      
      # Handle potential division by zero
      espiid_value = if (abs_env_distance > 0) mean_abs_diff / abs_env_distance else NA
      
      # Store ESPIID value for this environment pair
      espiid_values = c(espiid_values, espiid_value)
    }
  }
  
  # Return the mean ESPIID across all environment pairs
  return(mean(espiid_values, na.rm = TRUE))
}

#l=list(external_light)
#
#external_factors=list(external_light,external_water)
#
### test - passed on synthetic dataset 
#
#espiid_calculation(df_test2.2, trait_cols=2,factors=1)


############################################## dont mind the following 

#' Calculate Environmental Variance Weighted Plasticity Index (EVWPI)
#'
#' This function calculates the Environmental Variance Weighted Plasticity Index (EVWPI) 
#' based on user-specified environments and weighted environmental factors.
#'
#' @param dataframe A data frame containing the phenotypic data, including traits, genotype identifiers, and environmental factors.
#' @param trait_cols A vector specifying the column indices or names of the traits to be analyzed.
#' @param sp (Optional) An integer or string specifying the column that identifies the species or group for which EVWPI will be calculated. 
#' This can be the name or index of the column in the `dataframe`. If omitted, the function assumes the data pertains to a single species.
#' @param env_col A vector specifying the column indices or names of the pre-defined environments or a list of vectors representing the environments.
#' @param factors A vector specifying the column indices or names of the environmental factors to be weighted, or a list of vectors representing the environmental factors.
#' @param use_median (Optional) A logical value indicating whether to use the median instead of the mean for the EVWPI calculation. Defaults to FALSE (use mean).
#'
#' @return A list containing:
#' \item{all_results}{A list where each entry corresponds to a genotype (or species) and contains the EVWPI values for the specified traits.}
#'
#' @examples
#' # Example usage with column indices
#' df = data.frame(
#'   Genotype = rep(c("G1", "G2"), each = 6),
#'   Temperature = c(10, 20, 30, 10, 20, 30),
#'   pH = c(5, 6.5, 8, 5, 6.5, 8),
#'   Height = c(10, 12, 11, 14, 16, 15, 20, 22, 21, 23, 25, 24),
#'   Environment1 = c("EnvA", "EnvA", "EnvA", "EnvB", "EnvB", "EnvB")
#' )
#' 
#' # Specify the environments directly
#' environments = df$Environment1
#' factors = list(df$Temperature, df$pH)
#' 
#' # Calculate EVWPI for the 'Height' trait across genotypes and environments
#' results = evwpi_calculation(
#'   dataframe = df, 
#'   traits = "Height", 
#'   sp = "Genotype", 
#'   environments = environments,
#'   factors = factors
#' )
#' 
#' # Access the EVWPI values
#' print(results$all_results)
#' 
#' @export
evwpi_calculation = function(dataframe, trait_cols, sp = NULL, env_col, factors, use_median = FALSE) {
  #NOTE:this was an idea but needs to be worked on
  # Convert traits and sp from column indices to names if necessary
  if (!is.null(sp) && is.numeric(sp)) {
    sp = names(dataframe)[sp]
  }
  traits = if (is.numeric(trait_cols)) names(dataframe)[trait_cols] else trait_cols
  
  # Handle environments: if it's a vector of indices, extract the columns; if it's a list of vectors, use directly
  if (is.numeric(env_col)) {
    environments = dataframe[[names(dataframe)[env_col]]]
  }
  
  # Handle factors: if factors is a vector of indices/names, extract the columns; if it's a list of vectors, use directly
  if (is.numeric(factors) || is.character(factors)) {
    factors = lapply(factors, function(f) dataframe[[f]])
  }
  
  # Ensure the factors are correctly handled and match the length of the dataframe
  if (any(sapply(factors, length) != nrow(dataframe))) {
    stop("Length of factors must match the number of rows in the dataframe.")
  }
  
  # Step 1: Calculate the variance induced by each factor
  variance_scores = sapply(factors, function(factor) {
    var(factor, na.rm = TRUE)
  })
  
  # Normalize the variance scores to get the weights
  total_variance = sum(variance_scores, na.rm = TRUE)
  weights = variance_scores / total_variance
  
  all_results = list()
  
  if (is.null(sp)) {
    unique_species = list("Single_Group" = dataframe)
  } else {
    unique_species = split(dataframe, dataframe[[sp]])
  }
  
  # Loop over each species or group
  for (species_name in names(unique_species)) {
    species_data = unique_species[[species_name]]
    
    EVWPI_results = list()
    
    # Initialize a data frame to store EVWPI values for each pair of environments
    EVWPI_values = data.frame(sp = character(), trait = character(), evwpi = numeric())
    unique_envs=unique(environments)
    
    # Loop over each trait
    for (trait in traits) {
      for (i in 1:(length(unique_envs) - 1)) {
        for (j in (i + 1):length(unique_envs)) {
          
          env_i = unique_envs[i]
          env_j = unique_envs[j]
          
          trait_i = species_data[[trait]][environments == env_i]
          trait_j = species_data[[trait]][environments == env_j]
          abs_diff = abs(outer(trait_i, trait_j, "-"))
          
          # Calculate mean or median absolute phenotypic distance
          mean_abs_diff = if (use_median) median(abs_diff, na.rm = TRUE) else mean(abs_diff, na.rm = TRUE)
          
          # Calculate the weighted environmental distance using the variance-based weights
          weighted_env_distance = sum(weights)
          
          
          # Calculate EVWPI for this pair of environments
          evwpi_value = if (weighted_env_distance != 0) mean_abs_diff / weighted_env_distance else NaN
        
          
          # Store the results
          EVWPI_values = rbind(EVWPI_values, data.frame(sp = species_name, trait = trait, evwpi = evwpi_value))
        }
      }
    }
    
    all_results[[species_name]] = EVWPI_values
  }
  
  return(all_results)
}

#external_ph = sample(rep(c(1, 6, 18), each = 100))
#evwpi=evwpi_calculation(synthetic_data1,trait_cols = c(2,3),env_col=1, factors=1)
#print(evwpi)