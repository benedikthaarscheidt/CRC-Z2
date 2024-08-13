
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







#' @title Calculate Relative Distance Plasticity Index (RDPI) Across Multiple Traits and Environmental Conditions
#'
#' @description
#' This function computes the RDPI (Relative Distance Plasticity Index) for one or more traits within a single species or across multiple species/groups. 
#' RDPI quantifies the phenotypic plasticity of traits across different environmental conditions. The function allows users to specify 
#' which columns should be used for species (optional), traits, and environmental factors, either by name or by index.
#' The function generates boxplots to visualize RDPI distributions and performs ANOVA to determine significant differences in RDPI values 
#' across environmental groups and traits, followed by Tukey's HSD for pairwise comparisons.
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
  
  # Function to check and install missing packages
  check_and_install_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
      } else {
        library(pkg, character.only = TRUE)
      }
    }
  }
  
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
      
      mean_values <- aggregate(species_data[[trait]], list(species_data$Combined_Factors), mean)
      colnames(mean_values) <- c("Env_Combination", "Mean_Trait")
      overall_mean_trait_value <- mean(mean_values$Mean_Trait, na.rm = TRUE)
      env_levels <- levels(species_data$Combined_Factors)
      n_env_levels <- length(env_levels)
      
      rdpi_values <- c()
      
      for (i in 1:(n_env_levels - 1)) {
        for (j in (i + 1):n_env_levels) {
          mean_i <- mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[i]]
          mean_j <- mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[j]]
          
          diff <- abs(mean_i - mean_j)
          
          if (overall_mean_trait_value != 0) {
            rdpi_values <- c(rdpi_values, diff / overall_mean_trait_value)
          }
        }
      }
      
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
external_light <- rep(c(0.4, 0.6, 0.8), each = 100)
external_water <- sample(rep(c("Low", "High"), each = 150))

# Calculate RDPI with factors not in the dataframe
result <- rdpi_calculation(synthetic_data1, sp = NULL, traits = c(3, 4), factors = 2, factors_not_in_dataframe = list(external_light, external_water))

# Print the RDPI results and plot the boxplots
print(result$all_results)
print(result$trait_boxplots)

# View ANOVA and Tukey's HSD results
print(result$anova_results)
print(result$tukey_results)
