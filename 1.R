#this files contains the following indices: 
#coefficient-of variation total 
#slope of norm reaction,
# slope of plastic response (D),
# response coefficient (RC),

################################

if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2")
  library(roxygen2)
}

################################

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
generate_synthetic_data = function(n_plants, baseline_values, within_variance) {
  # Number of traits is determined by the number of columns in baseline_values
  n_traits = ncol(baseline_values)
  n_environments = nrow(baseline_values)
  
  # Initialize empty data frame to store the results
  synthetic_data = data.frame(matrix(nrow = n_plants * n_environments, ncol = n_traits+1))
  
  # Generate row names where the the integer part is the plant id and the fractional part is the environment indicator 
  # this can potentially be left out but I thought is nice to have
  rownames(synthetic_data) = rep(1:n_plants, times = n_environments) + rep(seq(0.1, (n_environments - 1) / 10 + 0.1, by = 0.1), each = n_plants)
  
  # Set column names for the traits
  colnames(synthetic_data) = c("Env_indicator",paste("Trait", 1:n_traits, sep = "_"))
  
  for (env in 1:n_environments) {
    # Calculate the row indices for the current environment. Like this it is easier to batch generate the normally distributed values.
    start_idx = (env - 1) * n_plants + 1
    end_idx = env * n_plants
    for (i in 1:n_traits) {
      # Generate random data for the current trait and environment with specific variance
      synthetic_data[start_idx:end_idx, 1] = env
      synthetic_data[start_idx:end_idx, i+1] = rnorm(n_plants, mean = baseline_values[env, i], sd = within_variance[env, i])
    }
  }
  
  
  return(synthetic_data)
}



n_plants = 100

# Define baseline values for each trait in each environment
baseline_values = matrix(
  c(100, 110, 120,  
    50, 195, 205,  
    70, 310, 290), 
  nrow = 3, ncol = 3, byrow = TRUE
)

# Define within-environment variance for each trait in each environment
within_variance = matrix(
  c(10, 15, 20,   
    8, 12, 10,    
    5, 7, 6),     
  nrow = 3, ncol = 3, byrow = TRUE
)

set.seed(12345)
synthetic_data1 = generate_synthetic_data(n_plants, baseline_values, within_variance)
print(class(synthetic_data1))

################################


#' Impute Missing Values in a Data Frame or Matrix
#'
#' This function imputes missing values in a data frame or matrix using various methods: median, mean, rpart, or k-NN.
#'
#' @param mt A data frame or matrix containing the data to be imputed.
#' @param mode The method of imputation. Options are "median", "mean", "rpart", or "knn".
#' @return The imputed data frame or matrix.
#' @details If the "rpart" method is selected, and the data frame is large, a warning is issued that the process might take a while.
#' If the "knn" method is selected and the VIM package is not installed, the user is prompted to install it.
#' @examples
#' df = data.frame(
#'   Trait_1 = c(1, 2, NA, 4, 5),
#'   Trait_2 = c(NA, 3, 3, 4, 5),
#'   Trait_3 = c(10, NA, 12, 13, 14)
#' )
#' imputed_df = impute(df, mode = "median")
#' @export
impute = function(mt, mode = "median") {
  if (!is.matrix(mt) && !is.data.frame(mt)) {
    stop("Input must be a matrix or data frame")
  }
  
  ismt = is.matrix(mt)
  if (ismt) {
    mt = as.data.frame(mt)
  }
  
  if (mode == "median") {
    for (i in 1:ncol(mt)) {
      if (!is.factor(mt[, i])) {
        mt[is.na(mt[, i]), i] = median(mt[, i], na.rm = TRUE)
      }
    }
  } else if (mode == "mean") {
    for (i in 1:ncol(mt)) {
      if (!is.factor(mt[, i])) {
        mt[is.na(mt[, i]), i] = mean(mt[, i], na.rm = TRUE)
      }
    }
  } else if (mode == "rpart") {
    if (nrow(mt) > 1000) {
      warning("Imputation using rpart may take a while for large datasets.")
    }
    for (i in 1:ncol(mt)) {
      midx = which(is.na(mt[, i]))
      if (length(midx) == 0) next
      idx = which(!is.na(mt[, i]))
      colname = colnames(mt)[i]
      frm = as.formula(paste(colname, "~ ."))
      mod = rpart(frm, data = mt[idx, ], method = ifelse(is.factor(mt[, i]), "class", "anova"))
      vals = predict(mod, newdata = mt[midx, ], type = ifelse(is.factor(mt[, i]), "class", "vector"))
      mt[midx, i] = vals
    }
  } else if (mode == "knn") {
    if (!requireNamespace("VIM", quietly = TRUE)) {
      choice = utils::menu(c("Yes", "No"), title = "Package 'VIM' is required for k-NN imputation. Would you like to install it?")
      if (choice == 1) {
        install.packages("VIM")
        library(VIM)
      } else {
        stop("k-NN imputation requires the 'VIM' package. Please install it to proceed.")
      }
    }
    mt = VIM::kNN(mt, k = 5, imp_var = FALSE)
  } else {
    stop(paste("mode", mode, "not yet implemented"))
  }
  
  if (ismt) {
    mt = as.matrix(mt)
  }
  
  return(mt)
}
################################

#' Calculate the Coefficient of Variation for Plasticity (CVt)
#'
#' This function calculates the Coefficient of Variation (CVt) for a given dataset. 
#' The CVt can be calculated either for each trait separately or as an overall value across all traits, 
#' depending on the `traitwise` parameter.
#'
#' @param data A numeric data frame, matrix or vector containing the trait data. Each column represents a trait, and each row represents an observation. 
#' @param traitwise A logical parameter. If `TRUE`, the CVt is calculated separately for each trait. If `FALSE`, the CVt is calculated across all traits combined.
#' 
#' @return A numeric vector if `traitwise = TRUE`, where each element represents the CVt for a single trait. 
#' If `traitwise = FALSE`, returns a single numeric value representing the overall CVt across all traits.
#' 
#' @details The Coefficient of Variation (CVt) is calculated as the ratio of the standard deviation to the mean. 
#' This function provides flexibility in calculating the CVt for each trait individually or for the entire dataset as a whole.
#' 
#' @examples
#' # Example dataset
#' synthetic_data = data.frame(
#'   Trait_1 = c(100.2, 92, 83.5, 96.5, 103.1),
#'   Trait_2 = c(101, 127, 139, 103, 113),
#'   Trait_3 = c(127.8, 103.5, 83.6, 92.5, 101.9)
#' )
#' 
#' # Calculate CVt for each trait separately
#' CVt_per_trait = calculate_CVt(synthetic_data1, traitwise = TRUE)
#' print(CVt_per_trait)
#' 
#' # Calculate overall CVt across all traits
#' overall_CVt = calculate_CVt(synthetic_data1, traitwise = FALSE)
#' print(overall_CVt)
#' 
#' @export
calculate_CVt = function(data, traitwise = TRUE) {
  if (traitwise) {
    # Calculate CVt for each trait separately
    CVt_values = apply(data, 2, function(x) sd(x) / mean(x))
    names(CVt_values) = colnames(data)
    return(CVt_values)
  } else {
    # Calculate overall CVt across all traits
    flattened_data = as.numeric(unlist(data))  # Flatten the data frame to a numeric vector
    overall_CVt = sd(flattened_data) / mean(flattened_data)
    return(overall_CVt)
  }
}

CVt_per_trait = calculate_CVt(synthetic_data1[,-1], traitwise = TRUE)

overall_CVt = calculate_CVt(synthetic_data1[,-1], traitwise = FALSE)


################################

#' Calculate the Reaction Norm Slope for a Specific Trait
#'
#' This function calculates the reaction norm slope for a specified trait in a data frame, 
#' using the first column as the environment indicator. It also provides an option to plot the reaction norm.
#'
#' @param data A data frame containing the environment indicators in the first column and the trait data in the remaining columns.
#' @param trait_col The column number of the trait for which the reaction norm slope is to be calculated. Defaults to 2 if not specified (since 1 is the environment indicator).
#' @param plot A logical flag indicating whether to plot the reaction norm. Defaults to FALSE.
#' @return The slope of the reaction norm for the specified trait.
#' @examples
#' df = data.frame(
#'   Environment = c(1, 1, 2, 2, 3, 3),
#'   Trait_1 = c(1, 2, 3, 4, 5, 6),
#'   Trait_2 = c(2, 4, 6, 8, 10, 12)
#' )
#' slope = calculate_reaction_norm_slope(df, trait_col = 2, plot = TRUE)
#' @export
calculate_reaction_norm_slope = function(data, env_col=1 ,trait_col = 2, plot = FALSE) {
  if (!is.numeric(data[[trait_col]])) {
    stop("The specified column is not numeric.")
  }
  
  if (any(is.na(data[[trait_col]]))) {
    stop("There are missing values in your data. Consider using the function impute().")
  }
  
  # Extract the environment indicators from the first column
  env_indicators = data[[env_col]]
  
  if (length(env_indicators) != nrow(data)) {
    stop("Length of env_indicators must match the number of rows in the data.")
  }
  
  # Fit a linear model to estimate the slope
  model = lm(data[[trait_col]] ~ env_indicators)
  
  # Extract the slope of the reaction norm
  slope = coef(model)[["env_indicators"]]
  
  # Plot the reaction norm if requested
  if (plot) {
    plot(data[[trait_col]] ~ env_indicators, 
         xlab = "Environment Indicator", 
         ylab = colnames(data)[trait_col],
         main = paste("Reaction Norm for", colnames(data)[trait_col]),
         pch = 19, col = "blue")
    abline(model, col = "red", lwd = 2)
    legend("topleft", legend = c("Observed", "Fitted Line"), col = c("blue", "red"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
  }
  
  return(slope)
}

slope=calculate_reaction_norm_slope(synthetic_data1,plot=F)
print(slope)

################################


#NOTE: if the resource availability is an actual measurement of a metabolite which is being used by the plant then the grouping of the plants into high vs low resource availability should be done by clustering.
# Sadly I am missing the isight into the common practices in the field.


#' Calculate the D Slope for a Specific Trait Between High and Low Resource Availability
#'
#' This function calculates the D slope, which is the difference in the mean value of a specified trait 
#' between high and low resource availability conditions. The D slope quantifies the scope of the plastic response.
#' 
#' By default, the function assumes that the environmental indicator is in the first column of the data frame. 
#' The first column should be numeric and sorted from low to high resource availability. 
#' The function will then automatically assume that the lowest values represent the low resource environment and 
#' the highest values represent the high resource environment. 
#' 
#' Alternatively, if a separate vector containing environmental information is provided, the function will use this vector 
#' instead of the first column to determine the environment categories.
#'
#' @param data A data frame containing the environmental indicators and trait data. By default, the first column is assumed to be the environmental indicator.
#' @param env_col The column number or name of the environmental indicator. Defaults to 1 if not specified.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_vector An optional vector specifying the environmental resource availability ("Low" or "High"). If provided, this vector is used instead of the environmental indicator column in the data frame.
#' @return The D slope, representing the difference in mean trait values between high and low resource conditions.
#' @examples
#' df = data.frame(
#'   Environment = c(1, 2, 3, 4),
#'   Height = c(10, 12, 20, 22)
#' )
#' # Assuming Environment is sorted from low to high resource availability
#' D_slope = calculate_D_slope(df, trait_col = "Height")
#' print(D_slope)
#' 
#' # With an explicit environment vector
#' env_vector = c("Low", "Low", "High", "High")
#' D_slope = calculate_D_slope(df, trait_col = "Height", env_vector = env_vector)
#' print(D_slope)
#' @export
calculate_D_slope = function(data, env_col = 1, trait_col, env_vector = NULL) {
  if (is.null(env_vector)) {
    # Use the specified or default environmental column
    env_data = data[[env_col]]
    # Sort the data based on the environmental indicator
    sorted_data = data[order(env_data), ]
    
    # Assume the lowest values represent "Low" and the highest values represent "High"
    num_rows = nrow(sorted_data)
    mid_point = ceiling(num_rows / 2)
    env_vector = c(rep("Low", mid_point), rep("High", num_rows - mid_point))
  } else {
    # Ensure the length of env_vector matches the number of rows in data
    if (length(env_vector) != nrow(data)) {
      stop("Length of env_vector must match the number of rows in the data.")
    }
  }
  
  # Extract the relevant data for the "High" and "Low" conditions
  trait_high = data[[trait_col]][env_vector == "High"]
  trait_low = data[[trait_col]][env_vector == "Low"]
  
  # Calculate the mean trait value for each condition
  mean_high = mean(trait_high, na.rm = TRUE)
  mean_low = mean(trait_low, na.rm = TRUE)
  
  # Calculate the D slope
  D_slope = mean_high - mean_low
  
  return(D_slope)
}

D=calculate_D_slope(synthetic_data1,env_col=1, trait_col = 2)
print(D)

################################

#' Calculate the Response Coefficient (RC) for a Specific Trait
#'
#' This function calculates the Response Coefficient (RC), which is the ratio of the mean value 
#' of a specified trait between high and low resource availability conditions.
#' 
#' By default, the function assumes that the environmental indicator is in the first column of the data frame. 
#' The first column should be numeric and sorted from low to high resource availability. 
#' The function will then automatically assume that the lowest values represent the low resource environment and 
#' the highest values represent the high resource environment. 
#' 
#' Alternatively, if a separate vector containing environmental information is provided, the function will use this vector 
#' instead of the first column to determine the environment categories.
#'
#' @param data A data frame containing the environmental indicators and trait data. By default, the first column is assumed to be the environmental indicator.
#' @param env_col The column number or name of the environmental indicator. Defaults to 1 if not specified.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_vector An optional vector specifying the environmental resource availability ("Low" or "High"). If provided, this vector is used instead of the environmental indicator column in the data frame.
#' @return The Response Coefficient (RC), representing the ratio of mean trait values between high and low resource conditions.
#' @examples
#' df = data.frame(
#'   Environment = c(1, 2, 3, 4),
#'   Height = c(10, 12, 20, 22)
#' )
#' # Assuming Environment is sorted from low to high resource availability
#' RC = calculate_RC(df, trait_col = "Height")
#' print(RC)
#' 
#' # With an explicit environment vector
#' env_vector = c("Low", "Low", "High", "High")
#' RC = calculate_RC(df, trait_col = "Height", env_vector = env_vector)
#' print(RC)
#' @export
calculate_RC = function(data, env_col = 1, trait_col, env_vector = NULL) {
  if (is.null(env_vector)) {
    # Use the specified or default environmental column
    env_data = data[[env_col]]
    # Sort the data based on the environmental indicator
    sorted_data = data[order(env_data), ]
    # Assume the lowest values represent "Low" and the highest values represent "High"
    num_rows = nrow(sorted_data)
    mid_point = ceiling(num_rows / 2)
    env_vector = c(rep("Low", mid_point), rep("High", num_rows - mid_point))
  } else {
    # Ensure the length of env_vector matches the number of rows in data
    if (length(env_vector) != nrow(data)) {
      stop("Length of env_vector must match the number of rows in the data.")
    }
  }
  
  # Extract the relevant data for the "High" and "Low" conditions
  trait_high = data[[trait_col]][env_vector == "High"]
  trait_low = data[[trait_col]][env_vector == "Low"]
  
  # Calculate the mean trait value for each condition
  mean_high = mean(trait_high, na.rm = TRUE)
  mean_low = mean(trait_low, na.rm = TRUE)
  
  # Calculate the Response Coefficient (RC)
  RC = mean_high / mean_low
  
  return(RC)
}

RC= calculate_RC(synthetic_data1,trait_col=2)
print(RC)

################################

#' Calculate the Coefficient of Variation of Means (CVm)
#'
#' This function calculates the Coefficient of Variation of Means (CVm), 
#' which is the standard deviation of the means divided by the mean of the means across different environments.
#'
#' The function allows for flexible grouping of the data. By default, it groups data using a specified column from the data frame 
#' (with the first column as the default). Alternatively, an external grouping vector can be provided. If an external vector is provided, 
#' it will be used for grouping instead of any column from the data frame.
#'
#' @param data A data frame containing the trait data.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number or name of the environmental indicator, used to group the data. Defaults to the first column (`data[,1]`).
#' @param env_vector An optional vector specifying the grouping (e.g., environmental conditions). If provided, this vector is used instead of `env_col` for grouping.
#' @return The CVmd, representing the ratio of the standard deviation of medians to the mean of medians.
#' @examples
#' df = data.frame(
#'   Environment = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3"),
#'   Height = c(10, 12, 20, 22, 15, 18)
#' )
#' CVm = calculate_CVm(df, env_col = "Environment", trait_col = "Height")
#' print(CVm)
#' @export
calculate_CVm = function(data, trait_col, env_col = data[,1], env_vector = NULL) {
  # Use the external vector if provided; otherwise, use the specified column from the data
  if (!is.null(env_vector)) {
    group_var = env_vector
  } else {
    group_var = env_col
  }
  means = tapply(data[[trait_col]], group_var, mean, na.rm = TRUE)
  # Calculate the standard deviation of the means
  sd_of_means = sd(means, na.rm = TRUE)

  # Calculate the mean of the means
  mean_of_means = mean(means, na.rm = TRUE)
  
  # Calculate the CVm
  CVm = sd_of_means / mean_of_means
  
  return(CVm)
}

CVm=calculate_CVm(synthetic_data1,trait_col = 3)
print(CVm)


###############################


#' Calculate the Coefficient of Variation of Medians (CVmd)
#'
#' This function calculates the Coefficient of Variation of Medians (CVmd), 
#' which is the standard deviation of the medians divided by the mean of the medians across different environments or groups.
#'
#' The function allows for flexible grouping of the data. By default, it groups data using a specified column from the data frame 
#' (with the first column as the default). Alternatively, an external grouping vector can be provided. If an external vector is provided, 
#' it will be used for grouping instead of any column from the data frame.
#'
#' @param data A data frame containing the trait data.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number or name of the environmental indicator, used to group the data. Defaults to the first column (`data[,1]`).
#' @param env_vector An optional vector specifying the grouping (e.g., environmental conditions). If provided, this vector is used instead of `env_col` for grouping.
#' @return The CVmd, representing the ratio of the standard deviation of medians to the mean of medians.
#' @examples
#' df = data.frame(
#'   Environment = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3"),
#'   Height = c(10, 12, 20, 22, 15, 18)
#' )
#' 
#' # Calculate CVmd using the default environment column
#' CVmd_default = calculate_CVmd(df, trait_col = "Height")
#' print(CVmd_default)
#' 
#' # Calculate CVmd using an external grouping vector
#' env_vector = c("Low", "Low", "High", "High", "Medium", "Medium")
#' CVmd_vector = calculate_CVmd(df, trait_col = "Height", env_vector = env_vector)
#' print(CVmd_vector)
#' @export
calculate_CVmd = function(data, trait_col, env_col = data[,1], env_vector = NULL) {
  
  # Use the external vector if provided; otherwise, use the specified column from the data
  if (!is.null(env_vector)) {
    group_var = env_vector
  } else {
    group_var = env_col
  }
  
  # Calculate the medians for each group
  medians = tapply(data[[trait_col]], group_var, median, na.rm = TRUE)
  
  # Calculate the standard deviation of the medians
  sd_of_medians = sd(medians, na.rm = TRUE)
  
  # Calculate the mean of the medians
  mean_of_medians = mean(medians, na.rm = TRUE)
  
  # Calculate the CVmd
  CVmd = sd_of_medians / mean_of_medians
  
  return(CVmd)
}

CVmd=calculate_CVmd(synthetic_data1,trait_col = 3)
print(CVmd)

#note that if the distribution of the trait values within each environment is roughly symmetric, then the mean and median of the data will be close to each other.
#As a result, both the CVm and CVmd will likely be close because the standard deviation and mean calculated from the means will be similar to those calculated 
#from the medians. If the data within each group is skewed, the mean and median will differ more significantly.



###########################
#' @title Calculate Relative Distance Plasticity Index (RDPI)
#'
#' @description
#' This function computes the RDPI (Relative Distance Plasticity Index) for a given dataset. 
#' RDPI quantifies the phenotypic plasticity of a trait across different environmental factors. 
#' The function allows flexibility in the structure of the input data, enabling the user to specify 
#' which columns should be used for species, traits, and environmental factors, either by name or by index.
#'
#' @param dataframe A data frame containing the data.
#' @param sp The name or index of the column that identifies the group or species for which RDPI will be calculated.
#' @param trait The name or index of the column that contains the trait values. This column must be numeric.
#' @param factor The name or index of the column that contains the environmental factor levels.
#' RDPI is calculated by comparing trait values across these levels.
#' 
#' @return A list containing:
#'   \item{RDPI_values}{A data frame with RDPI values for each species or group.}
#'   \item{summary}{A summary statistics table for the RDPI values, including mean, standard deviation, and standard error per species.}
#'   \item{plot}{A ggplot object representing the boxplot of RDPI values across species or groups.}
#'   \item{test_result}{The result of a statistical test (t-test or ANOVA) comparing RDPI across species.}
#' 
#' @details
#' This function computes RDPI for each species or group based on the trait values across different levels 
#' of an environmental factor. The results include a summary of RDPI values, a boxplot for visualization, 
#' and the result of a statistical test (t-test if there are two groups, ANOVA if there are more than two).
#'
#' @examples
#' # Example with a custom dataset
#' df = data.frame(
#'   Species = rep(c("A", "B", "C"), each = 10),
#'   Trait = rnorm(30, mean = 50, sd = 10),
#'   Env = rep(c("Low", "High"), times = 15)
#' )
#' result = rdpi(df, sp = "Species", trait = "Trait", factor = "Env")
#' 
#' # Example using column indices
#' result_indices = rdpi(df, sp = 1, trait = 2, factor = 3)
#'
#' # Accessing results
#' print(result$RDPI_values)
#' print(result$summary)
#' print(result$plot)
#' print(result$test_result)
#'
#' @export
rdpi = function(dataframe, sp, trait, factor) {
  library(dplyr)
  # Load only the required functions from necessary packages
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("agricolae", quietly = TRUE)
  requireNamespace("psych", quietly = TRUE)
  
  # Convert column indices to names if necessary
  sp = if (is.numeric(sp)) names(dataframe)[sp] else sp
  trait = if (is.numeric(trait)) names(dataframe)[trait] else trait
  factor = if (is.numeric(factor)) names(dataframe)[factor] else factor
  
  # Create the object (an empty data frame) that will store the results of RDPI
  RDPI = data.frame(sp = character(0), rdpi = numeric(0))
  
  # Get the unique levels of the species (or group) column
  levels_dataframe = levels(droplevels(dataframe %>% pull({{sp}})))
  
  for (a in levels_dataframe) {
    
    # Subset the data for the current species or group
    data_sp = dataframe |> dplyr::filter(.data[[sp]] == a)
    
    # Compute RDPI using a custom function (assuming rdpi_matrix is available)
    RDPI_temp = tryCatch({
      rdpi_matrix(data_sp, .data[[trait]], .data[[factor]])
    }, error = function(e) {
      warning(sprintf("Error in RDPI calculation for %s: %s", a, e$message))
      return(NA)
    })
    
    RDPI_sp = data.frame(sp = as.character(a),
                          rdpi = RDPI_temp)
    
    RDPI = dplyr::bind_rows(RDPI, RDPI_sp) %>%
      dplyr::mutate(sp = as.factor(sp))
  }
  
  # Summary statistics for each species or group
  summary = RDPI %>%
    dplyr::group_by(sp) %>%
    dplyr::summarise(mean = mean(rdpi, na.rm = TRUE),
                     sd = sd(rdpi, na.rm = TRUE),
                     se = psych::se(rdpi, na.rm = TRUE))
  
  # Boxplot for RDPI values across species or groups
  boxplot_rdpi = ggplot2::ggplot(RDPI) +
    ggplot2::geom_boxplot(ggplot2::aes(x = sp, y = rdpi)) +
    ggplot2::ylab("RDPI") +
    ggplot2::xlab(sp)
  
  print(boxplot_rdpi)
  
  # Perform a statistical test (t-test or ANOVA) depending on the number of levels in sp
  if (nlevels(RDPI$sp) < 3) {
    fit = try(t.test(RDPI$rdpi ~ RDPI$sp), silent = TRUE)
    test_result = if (class(fit) != "try-error") fit else "To compute a t-test, the grouping factor must have exactly 2 levels"
  } else {
    fit = aov(RDPI$rdpi ~ RDPI$sp)
    test_result = summary(fit)
    Tuk = agricolae::HSD.test(fit, trt = "RDPI$sp")
    Tuk$groups$sp = as.factor(row.names(Tuk$groups))
    summary = dplyr::left_join(summary, Tuk$groups) %>%
      dplyr::select(-`RDPI$rdpi`)
  }
  
  # Return the RDPI data, summary, plot, and test result
  return(list(RDPI_values = RDPI, summary = summary, plot = boxplot_rdpi, test_result = test_result))
}

result_indices = rdpi(synthetic_data1, sp = 1, trait = 2, factor = 3)
print(result_indices)













