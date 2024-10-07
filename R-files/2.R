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
  c(10, 11, 12,  
    51, 52, 53,  
    100, 101, 102), 
  nrow = 3, ncol = 3, byrow = TRUE
)

# Define within-environment variance for each trait in each environment
within_variance = matrix(
  c(1, 5, 5,   
    8, 12, 10,    
    5, 7, 6),     
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



#' @title Calculate Relative Distance Plasticity Index (RDPI) Across Multiple Traits and Environmental Conditions
#'
#' @description
#' This function computes the RDPI (Relative Distance Plasticity Index) for one or more traits within a single species or across multiple species/groups. 
#' RDPI quantifies the phenotypic plasticity of traits across different environmental conditions. The function allows users to specify 
#' which columns should be used for species (optional), traits, and environmental factors, either by name or by index.
#' The function generates boxplots to visualize RDPI distributions and performs ANOVA to determine significant differences in RDPI values 
#' across environmental groups and traits, followed by Tukey's HSD for pairwise comparisons.
#' 
#' RDPI Calculation:
#' For each pair of environmental combinations, it calculates the absolute difference between their mean trait values.
#' This difference is normalized by dividing it by the overall mean trait value across all environmental combinations.
#' The RDPI for a given trait is the average of these normalized differences across all pairs of environmental combinations.

#'
#' @param dataframe A data frame containing the data to be analyzed. This should include the traits of interest and the environmental factors. 
#' The species identifier column is optional and only necessary if multiple species/groups are being analyzed.
#' @param sp (Optional) An integer or string specifying the column that identifies the species or group for which RDPI will be calculated. 
#' This can be the name or index of the column in the `dataframe`. If omitted, the function assumes the data pertains to a single species.
#' @param traits A vector of integers or strings specifying the columns that contain the trait values to be analyzed. 
#' These can be the names or indices of the columns in the `dataframe`.
#' @param factors (Optional) A vector of integers or strings specifying the columns that contain the environmental factors. 
#' RDPI is calculated by comparing trait values across these environmental factors. This can be the names or indices of the columns in the `dataframe`.
#' If `factors_not_in_dataframe` is provided, this argument can be omitted or set to `NULL`.
#' @param factors_not_in_dataframe (Optional) A list of vectors representing environmental factors that are not included in the `dataframe`. 
#' These vectors should have the same length as the `dataframe` and correspond to the environmental conditions for each observation. 
#' If `factors` is provided, this argument can be omitted or set to `NULL`.
#' 
#' @return A nested list containing:
#'   \item{all_results}{A list where each element corresponds to a species or group. Each species/group contains a list of RDPI results for each trait. 
#'   Each trait includes RDPI values, summary statistics, and Tukey's HSD results if applicable.}
#'   \item{test_result}{The result of an ANOVA test comparing RDPI values across traits and environmental groups.}
#'   \item{boxplot}{A ggplot object representing the boxplot of RDPI values across all specified traits.}
#' 
#' @details
#' This function computes RDPI for each species or group based on the trait values across different levels of environmental factors. 
#' If the `sp` parameter is omitted, the function treats the data as pertaining to a single species. The results include a summary of RDPI values, 
#' a boxplot for visualization, and the result of a statistical test (ANOVA) comparing RDPI values across traits and environmental groups. 
#' Tukey's HSD test is performed for pairwise comparisons if applicable.
#' The RDPI is calculated by taking the absolute differences in trait values between pairs of environmental levels, standardizing these differences by dividing by 
#' the overall mean trait value across all levels, and averaging these standardized differences.
#'
#' @examples
#' # Example with a custom dataset using column names, without specifying a species identifier
#' df = data.frame(
#'   Height = rnorm(36, mean = 50, sd = 10),
#'   LeafArea = rnorm(36, mean = 25, sd = 5),
#'   RootDepth = rnorm(36, mean = 15, sd = 3),
#'   Light = rep(c("Low", "High"), times = 18),
#'   Water = rep(c("Low", "High"), each = 6, times = 3)
#' )
#' 
#' # Calculate RDPI with factors from the dataframe
#' result = rdpi_calculation(df, traits = c("Height", "LeafArea", "RootDepth"), factors = c("Light", "Water"))
#'
#' # View the RDPI values, summary statistics, boxplot, and test results for a specific trait
#' print(result$all_results$Single_Group$Height$RDPI_values)
#' print(result$all_results$Single_Group$Height$summary)
#' print(result$boxplot)
#' print(result$test_result)
#'
#' # Example using external factors not in the dataframe, without specifying a species identifier
#' external_light = rep(c("Low", "High"), times = 18)
#' external_water = rep(c("Low", "High"), each = 6, times = 3)
#' 
#' result_external = rdpi_calculation(df, traits = c("Height", "LeafArea", "RootDepth"), factors_not_in_dataframe = list(external_light, external_water))
#' 
#' print(result_external$all_results$Single_Group$Height$RDPI_values)
#' print(result_external$all_results$Single_Group$Height$summary)
#' print(result_external$boxplot)
#' print(result_external$test_result)
#'
#' # Example using the species identifier
#' df_with_species = data.frame(
#'   Species = rep(c("A", "B"), each = 18),
#'   Height = rnorm(36, mean = 50, sd = 10),
#'   LeafArea = rnorm(36, mean = 25, sd = 5),
#'   RootDepth = rnorm(36, mean = 15, sd = 3),
#'   Light = rep(c("Low", "High"), times = 18),
#'   Water = rep(c("Low", "High"), each = 6, times = 3)
#' )
#' 
#' result_with_species = rdpi_calculation(df_with_species, sp = "Species", traits = c("Height", "LeafArea", "RootDepth"), factors = c("Light", "Water"))
#' 
#' print(result_with_species$all_results$A$Height$RDPI_values)
#' print(result_with_species$all_results$A$summary)
#' print(result_with_species$boxplot)
#' print(result_with_species$test_result)
#' 
#' @export
rdpi_calculation <- function(dataframe, traits, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL) {
  
  # List of required packages
  required_packages <- c("ggplot2", "agricolae", "dplyr", "reshape2")
  
  # Check and install required packages
  check_and_install_packages(required_packages)
  
  # Convert column indices to names if necessary
  if (!is.null(sp)) {
    sp <- if (is.numeric(sp)) names(dataframe)[sp] else sp
  }
  traits <- if (is.numeric(traits)) names(dataframe)[traits] else traits
  
  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    external_factors_df <- as.data.frame(factors_not_in_dataframe)
    if (!is.null(factors)) {
      factors <- if (is.numeric(factors)) names(dataframe)[factors] else factors
      combined_factors_df <- cbind(dataframe[factors], external_factors_df)
    } else {
      combined_factors_df <- external_factors_df
    }
    dataframe$Combined_Factors <- interaction(combined_factors_df, drop = TRUE)
  } else if (!is.null(factors)) {
    factors <- if (is.numeric(factors)) names(dataframe)[factors] else factors
    dataframe$Combined_Factors <- interaction(dataframe[factors], drop = TRUE)
  } else {
    stop("You must provide either internal factors, external factors, or both.")
  }
  
  dataframe$Combined_Factors <- as.factor(dataframe$Combined_Factors)
  
  all_results <- list()
  
  if (is.null(sp)) {
    unique_species <- list("Single_Group" = dataframe)
  } else {
    unique_species <- split(dataframe, dataframe[[sp]])
  }
  
  # Initialize a dataframe to collect all trait data for statistical analysis
  all_trait_data <- data.frame()
  
  for (species_name in names(unique_species)) {
    species_data <- unique_species[[species_name]]
    
    RDPI_results <- list()
    
    for (trait in traits) {
      RDPI <- data.frame(sp = character(0), rdpi = numeric(0))
      
      # Compute pairwise Canberra distance for all individuals in the dataset
      RDPI_temp <- as.matrix(dist(x = species_data[[trait]], method = "canberra"))
      
      # Generate a matrix where value is "TRUE" only if observation i and j belong to different levels of the factor
      filter_frame <- matrix(NA, nrow(species_data), nrow(species_data))
      
      for (i in 1:nrow(filter_frame)) {
        for (j in 1:ncol(filter_frame)) {
          if (species_data$Combined_Factors[i] == species_data$Combined_Factors[j]) {
            filter_frame[i, j] <- FALSE
          } else {
            filter_frame[i, j] <- TRUE
          }
        }
      }
      
      # Keep only the lower triangle (no need to compare twice)
      filter_frame[upper.tri(filter_frame, diag = TRUE)] <- FALSE
      
      # Subset RDPI so that it only includes comparisons between individuals from different environments
      rdpi_values <- RDPI_temp[filter_frame == TRUE]
      
      # Calculate the mean RDPI value for this trait
      RDPI_sp_value <- mean(rdpi_values, na.rm = TRUE)
      RDPI_sp <- data.frame(sp = species_name, rdpi = RDPI_sp_value)
      RDPI <- rbind(RDPI, RDPI_sp)
      
      # Collect all trait data for statistical analysis
      trait_data <- data.frame(
        sp = species_name,
        trait = trait,
        Combined_Factors = species_data$Combined_Factors,
        Trait_Value = species_data[[trait]]
      )
      
      all_trait_data <- rbind(all_trait_data, trait_data)
      
      RDPI_results[[trait]] <- list(
        RDPI_values = RDPI,
        trait_data = trait_data
      )
    }
    
    all_results[[species_name]] <- RDPI_results
  }
  
  # Perform ANOVA on the trait data
  anova_results <- list()
  tukey_results <- list()
  
  for (trait in traits) {
    fit <- aov(Trait_Value ~ Combined_Factors, data = subset(all_trait_data, trait == trait))
    anova_results[[trait]] <- summary(fit)
    
    # Perform Tukey's HSD test
    Tukey <- agricolae::HSD.test(fit, trt = "Combined_Factors")
    tukey_results[[trait]] <- Tukey
  }
  
  # Create boxplots for the traits across environmental combinations
  boxplot_traits <- ggplot2::ggplot(all_trait_data, ggplot2::aes(x = Combined_Factors, y = Trait_Value, fill = trait)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~trait, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Environmental Combinations") +
    ggplot2::ylab("Trait Values") +
    ggplot2::ggtitle("Boxplots of Trait Values Across Environmental Combinations") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Return the RDPI results, trait boxplots, ANOVA, and Tukey's HSD test results
  return(list(
    all_results = all_results,
    trait_boxplots = boxplot_traits,
    anova_results = anova_results,
    tukey_results = tukey_results
  ))
}


# Example usage with synthetic data
external_light = rep(c(0.4, 0.6, 0.8), each = 100)
external_water = sample(rep(c("Low", "High"), each = 150))

# Calculate RDPI with factors not in the dataframe
RDPI= rdpi_calculation(synthetic_data1, sp = NULL, traits = c(3, 4), factors = 2, factors_not_in_dataframe = list(external_light, external_water))


### test - passed on synthetic dataset and crosschecked with data from package


df_test20 <- data.frame(
  Column0 = as.factor(c(rep(1, 20), rep(2, 20))),  # Two groups/species, 20 observations each
  Column1 = as.factor(c(rep(1, 10), rep(2, 10), rep(1, 10), rep(2, 10))),  # Two environments per group
  Column2 = c(rep(2, 10), rep(1, 10), rep(2, 10), rep(1, 10))  # Numeric trait values
)
rdpi_calculation(df_test20,sp=1, traits=3,factors = 2)


#####################################################





#' Calculate the Relative Distance Plasticity Index for Mean Phenotypic Values (RDPIs)
#'
#' This function calculates the Relative Distance Plasticity Index (RDPI) for mean phenotypic values. 
#' RDPI is defined as the absolute phenotypic distance between the means of the same genotype 
#' in different environments, divided by one of the two mean phenotypic values (typically the smaller one).
#' 
#' The function is designed to handle multiple traits, genotypes, and environmental combinations, 
#' providing a comprehensive analysis of phenotypic plasticity across environmental gradients.
#'
#' @param dataframe A data frame containing the phenotypic data, including traits, genotype identifiers, 
#' and environmental factors.
#' @param traits A vector specifying the column indices or names of the traits to be analyzed.
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
#'   traits = "Height", 
#'   sp = "Genotype", 
#'   factors = "Environment"
#' )
#' 
#' # Access the RDPI values and boxplots
#' print(results$all_results)
#' print(results$trait_boxplots)
#' 
#' @export
rdpi_mean_calculation = function(dataframe, traits, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL) {
#NOTE: would it make sense to calculate this score not only pairwise for 2 different env and 1 trait, but as a sum over all traits those 2 envs share?
    
    # List of required packages
    required_packages = c("ggplot2", "agricolae", "dplyr", "reshape2")
    
    # Check and install required packages
    check_and_install_packages(required_packages)
    
    # Convert column indices to names if necessary
    if (!is.null(sp)) {
      sp = if (is.numeric(sp)) names(dataframe)[sp] else sp
    }
    traits = if (is.numeric(traits)) names(dataframe)[traits] else traits
    
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
    
    all_results = list()
    
    if (is.null(sp)) {
      unique_species = list("Single_Group" = dataframe)
    } else {
      unique_species = split(dataframe, dataframe[[sp]])
    }
    
    # Loop over each species or group
    for (species_name in names(unique_species)) {
      species_data = unique_species[[species_name]]
      
      RDPI_results = list()
      
      # Initialize a data frame to store RDPI values for each pair of environments
      RDPI_values = data.frame(sp = character(), env1 = character(), env2 = character(), rdpi = numeric())
      
      # Get unique environment combinations
      env_levels = levels(species_data$Combined_Factors)
      n_env_levels = length(env_levels)
      
      # Loop over each trait
      for (trait in traits) {
        mean_values = aggregate(species_data[[trait]], list(species_data$Combined_Factors), mean)
        colnames(mean_values) = c("Env_Combination", "Mean_Trait")
        
        for (i in 1:(n_env_levels - 1)) {
          for (j in (i + 1):n_env_levels) {
            mean_i = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[i]]
            mean_j = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[j]]
            
            # Calculate RDPI for this trait between the two environments
            rdpi_value = abs(mean_i - mean_j) / min(mean_i, mean_j)
            
            # Sum the RDPI values across all traits
            RDPI_values = rbind(RDPI_values, data.frame(sp = species_name, env1 = env_levels[i], env2 = env_levels[j], rdpi = rdpi_value))
          }
        }
      }
      
      all_results[[species_name]] = RDPI_values
    }
    
    return(all_results)
  }
  


# Example usage with synthetic data
external_light = rep(c(0.4, 0.6, 0.8), each = 100)
external_water =rep(c("Low", "High"), each = 150)

## test - passed on synthetic dataset

rdpi_mean_calculation(df_test20,sp=1, traits=3,factors = 2)


#########################################


#' Calculate Environmental Sensitivity Performance Index (ESPI)
#'
#' This function calculates the Environmental Sensitivity Performance Index (ESPI) 
#' for a given trait across different environmental conditions. ESPI is calculated 
#' using the formula: 
#' \deqn{ESPI = \frac{(\text{Maximum mean} - \text{Minimum mean})}{\text{Absolute distance between environmental values}}}
#'
#' @param data A data frame containing the trait and environmental variables.
#' @param trait_col A numeric or character vector specifying the column index or 
#' name of the trait in the data frame.
#' @param env_col A numeric or character vector specifying the column index or 
#' name of the environmental variable in the data frame.
#' 
#' @return A numeric value representing the Environmental Sensitivity Performance Index (ESPI).
#' 
#' @details 
#' The ESPI is a measure of how sensitive a particular trait is to changes 
#' in environmental conditions. It is calculated by taking the difference between the 
#' maximum and minimum mean values of the trait across different environments and 
#' dividing this by the absolute distance between the maximum and minimum environmental 
#' values. This index assumes normality in the data, and the choice of an appropriate 
#' environmental range is crucial for accurate calculations.
#' 
#' @examples 
#' \dontrun{
#' # Example usage:
#' df = data.frame(
#'   trait = rnorm(100),
#'   environment = seq(1, 100)
#' )
#' 
#' ESPI_value = calculate_ESPI(df, trait_col = "trait", env_col = "environment")
#' print(ESPI_value)
#' }
#' 
#' @seealso 
#' \code{\link{tapply}}, \code{\link{mean}}
#' 
#' @export
calculate_ESPI = function(data, trait_col, env_col) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Handle trait_col
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Calculate mean trait values for each environment
  means = tapply(trait_col, env_col, mean, na.rm = TRUE)
  
  # Identify maximum and minimum means
  max_mean = max(means, na.rm = TRUE)
  min_mean = min(means, na.rm = TRUE)
  
  # Calculate the absolute distance between the maximum and minimum environmental values
  abs_env_distance = abs(max(env_col, na.rm = TRUE) - min(env_col, na.rm = TRUE))
  
  # Calculate ESPI
  ESPI = (max_mean - min_mean) / abs_env_distance
  
  return(ESPI)
}

## test - passed on synthetic dataset 


df_test2.2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 5),rep(4,5) ,rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
calculate_ESPI(df_test2.2,trait_col = 2,env_col = 1)


#####################################


#' Calculate Environmental Sensitivity Performance Index for Individual Differences (ESPIID)
#'
#' This function calculates the Environmental Sensitivity Performance Index for Individual Differences (ESPIID). 
#' ESPIID is defined as the absolute phenotypic distance between individuals of the same genotype in different environments, 
#' divided by the absolute distance between environmental values.
#' 
#' The function is designed to handle multiple traits, genotypes, and environmental combinations, providing a comprehensive analysis 
#' of phenotypic sensitivity to environmental changes across individual plants.
#'
#' @param dataframe A data frame containing the phenotypic data, including traits, genotype identifiers, and environmental factors.
#' @param traits A vector specifying the column indices or names of the traits to be analyzed.
#' @param sp (Optional) A column index or name indicating the species or genotype identifier. If not provided, the data is treated as a single group.
#' @param factors (Optional) A vector of column indices or names specifying the internal environmental factors within the data frame that should be combined to create unique environmental combinations.
#' @param factors_not_in_dataframe (Optional) A list of vectors representing external environmental factors that are not included in the data frame. Each vector must have a length equal to the number of rows in the data frame.
#' @param use_median (Optional) A logical value indicating whether to use the median instead of the mean for the ESPIID calculation. Defaults to FALSE (use mean).
#' 
#' @return A list containing:
#' \item{all_results}{A list where each entry corresponds to a genotype (or species) and contains the ESPIID values for the specified traits.}
#'
#' @details 
#' The ESPIID provides a measure of phenotypic sensitivity by comparing the individual phenotypic values of the same genotype across different environments. 
#' For each pair of environmental combinations, it calculates the absolute differences between all pairs of individuals in those environments using the `outer()` function. 
#' The ESPIID is the mean or median of these absolute phenotypic distances, normalized by dividing by the absolute distance between the environmental values.
#' 
#' When to use ESPIID vs. ESPI:
#' - **ESPIID** focuses on the phenotypic sensitivity of individual plants within the same genotype across environments, making it useful for datasets where individual variation is important.
#' - **ESPI** provides a broader measure of phenotypic plasticity by comparing the mean trait values of genotypes or groups across environments. ESPI is more suitable when you want to assess general plasticity across groups rather than individual sensitivity.
#'
#' Use ESPIID when you are interested in within-genotype variation and want to evaluate the spread of individual phenotypic responses. ESPI, on the other hand, is more suitable when you are comparing the overall average plasticity across environmental conditions.
#'
#' @examples
#' # Example usage
#' df = data.frame(
#'   Genotype = rep(c("G1", "G2"), each = 6),
#'   Environment = rep(c("Env1", "Env2", "Env3"), times = 4),
#'   Height = c(10, 12, 11, 14, 16, 15, 20, 22, 21, 23, 25, 24)
#' )
#' 
#' # Calculate ESPIID for the 'Height' trait across genotypes and environments
#' results = espiid_calculation(
#'   dataframe = df, 
#'   traits = "Height", 
#'   sp = "Genotype", 
#'   factors = "Environment"
#' )
#' 
#' # Access the ESPIID values
#' print(results$all_results)
#' 
#' @export
espiid_calculation = function(dataframe, traits, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL, use_median = FALSE) {
  
  # List of required packages
  required_packages = c("dplyr", "reshape2")
  
  # Check and install required packages
  check_and_install_packages(required_packages)
  
  # Convert column indices to names if necessary
  if (!is.null(sp)) {
    sp = if (is.numeric(sp)) names(dataframe)[sp] else sp
  }
  traits = if (is.numeric(traits)) names(dataframe)[traits] else traits
  
  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    # Check if Combined_Factors is continuous
    for(i in 1:length(factors_not_in_dataframe)){
      if(!is.numeric(factors_not_in_dataframe[[i]])){
        stop("The combined environmental factors are not continuous. ESPIID assumes continuous environmental variables. Results may not be accurate.")
      }
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
  
  all_results = list()
  
  if (is.null(sp)) {
    unique_species = list("Single_Group" = dataframe)
  } else {
    unique_species = split(dataframe, dataframe[[sp]])
  }
  
  # Loop over each species or group
  for (species_name in names(unique_species)) {
    species_data = unique_species[[species_name]]
    
    ESPIID_results = list()
    
    # Initialize a data frame to store ESPIID values for each pair of environments
    ESPIID_values = data.frame(sp = character(), env1 = character(), env2 = character(), espiid = numeric())
    
    # Get unique environment combinations
    env_levels = levels(species_data$Combined_Factors)
    n_env_levels = length(env_levels)
    
    # Loop over each trait
    for (trait in traits) {
      for (i in 1:(n_env_levels - 1)) {
        for (j in (i + 1):n_env_levels) {
          trait_i = species_data[[trait]][species_data$Combined_Factors == env_levels[i]]
          trait_j = species_data[[trait]][species_data$Combined_Factors == env_levels[j]]
          
          # Calculate absolute phenotypic distances between all pairs of individuals in env_i and env_j
          abs_diff = abs(outer(trait_i, trait_j, "-"))
          
          # Calculate mean or median absolute phenotypic distance
          mean_abs_diff = if (use_median) median(abs_diff, na.rm = TRUE) else mean(abs_diff, na.rm = TRUE)
          
          # Calculate absolute distance between environmental values
          abs_env_distance = abs(as.numeric(env_levels[i]) - as.numeric(env_levels[j]))
          
          # Calculate ESPIID for this pair of environments
          espiid_value = mean_abs_diff / abs_env_distance
          
          # Store the results
          ESPIID_values = rbind(ESPIID_values, data.frame(sp = species_name, env1 = env_levels[i], env2 = env_levels[j], espiid = espiid_value))
        }
      }
    }
    
    all_results[[species_name]] = ESPIID_values
  }
  
  return(all_results)
}

l=list(external_light)

external_factors=list(external_light,external_water)

## test - passed on synthetic dataset 

espiid_calculation(df_test2.2, traits=2,factors=1)


############################################## dont mind the following 

#' Calculate Environmental Variance Weighted Plasticity Index (EVWPI)
#'
#' This function calculates the Environmental Variance Weighted Plasticity Index (EVWPI) 
#' based on user-specified environments and weighted environmental factors.
#'
#' @param dataframe A data frame containing the phenotypic data, including traits, genotype identifiers, and environmental factors.
#' @param traits A vector specifying the column indices or names of the traits to be analyzed.
#' @param sp (Optional) An integer or string specifying the column that identifies the species or group for which EVWPI will be calculated. 
#' This can be the name or index of the column in the `dataframe`. If omitted, the function assumes the data pertains to a single species.
#' @param environments A vector specifying the column indices or names of the pre-defined environments or a list of vectors representing the environments.
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
evwpi_calculation = function(dataframe, traits, sp = NULL, environments, factors, use_median = FALSE) {
  #NOTE:this was an idea but needs to be worked on
  # Convert traits and sp from column indices to names if necessary
  if (!is.null(sp) && is.numeric(sp)) {
    sp = names(dataframe)[sp]
  }
  traits = if (is.numeric(traits)) names(dataframe)[traits] else traits
  
  # Handle environments: if it's a vector of indices, extract the columns; if it's a list of vectors, use directly
  if (is.numeric(environments)) {
    environments = dataframe[[names(dataframe)[environments]]]
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

external_ph = sample(rep(c(1, 6, 18), each = 100))
evwpi=evwpi_calculation(synthetic_data1,traits = 3,environments=1, factors=list(external_ph))
print(evwpi)