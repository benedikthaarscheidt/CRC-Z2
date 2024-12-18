#this files contains the following indices/functions: 
#generate_synthetic_data,
#impute,
#coefficient-of variation total (calculate_CVt) - tested,
#slope of norm reaction (calculate_reaction_norm_slope) - tested,
#slope of plastic response (D) (calculate_D_slope)- tested,
#response coefficient (RC) (calculate_RC)- tested,
#Standard deviation of means (CVm) (calculate_CVm)- tested,
#Standard deviation of medians (CVmd)(calculate_CVmd)- tested,
#Grand plasticity (calculate_grand_plasticity)- tested,
#combine_factors- tested,
#Phenotypic Plasticity Index (calculate_PPF)- tested,
#Phenotypic Plasticity Index (calculate_Phenotypic_Plasticity_Index)- tested,
#PImd (calculate_PImd)- tested,
#PILSM (calculate_PILSM)- tested,
#RTR (calculate_RTR)- tested,
#PIR (calculate_PIR) - tested


################################

if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2")
  library(roxygen2)
}

source("~/CRC 1622 - Z2/R-files/2.R")

################################

##' Generate Synthetic Plant Plasticity Data with Sequential Environments
##'
##' This function generates a synthetic dataset for plant plasticity measurements across multiple environments. 
##' It simulates traits for a specified number of plants in different environments, with trait values normally 
##' distributed around given baseline values. The generated dataset is structured such that the first set of rows 
##' corresponds to all plants in the first environment, followed by all plants in the second environment, and so on. 
##'
##' @param n_plants An integer specifying the number of plants (rows) to generate in each environment.
##' @param baseline_values A matrix or data frame where each row represents an environment and each column 
##' represents a trait. The values in this matrix provide the baseline (mean) values for each trait in each environment.
##' @param within_variance A matrix or data frame where each row represents an environment and each column 
##' represents a trait. The values in this matrix specify the standard deviation for each trait in each environment, 
##' determining the variability around the baseline values.
##' 
##' @return A data frame containing the synthetic dataset. The data frame will have `n_plants * n_environments` rows 
##' and `n_traits` columns. The row names indicate the plant number and environment (e.g., 1.1 for plant 1 in environment 1).
##' The column names represent the traits (e.g., `Trait_1`, `Trait_2`, etc.).
##' 
##' @details This function can be used for simulating data where different environments may have different 
##' levels of variability for the same trait. The data is structured sequentially by environment, making it easy to analyze 
##' the effects of different environmental conditions on plant traits.
##' 
##' @examples
##' # Define the parameters
##' n_plants = 100
##' 
##' # Define baseline values for each trait in each environment
##' # Rows are environments, columns are traits
##' baseline_values = matrix(
##'   c(100, 110, 120,  # Trait 1 baseline values in Env 1, 2, 3
##'     200, 195, 205,  # Trait 2 baseline values in Env 1, 2, 3
##'     300, 310, 290), # Trait 3 baseline values in Env 1, 2, 3
##'   nrow = 3, ncol = 3, byrow = TRUE
##' )
##' 
##' # Define within-environment variance for each trait in each environment
##' within_variance = matrix(
##'   c(10, 15, 20,   # Trait 1 variance in Env 1, 2, 3
##'     8, 12, 10,    # Trait 2 variance in Env 1, 2, 3
##'     5, 7, 6),     # Trait 3 variance in Env 1, 2, 3
##'   nrow = 3, ncol = 3, byrow = TRUE
##' )
##' 
##' # Generate the synthetic dataset
##' synthetic_data = generate_synthetic_data(n_plants, baseline_values, within_variance)
##' 
##' # View the first few rows of the generated dataset
##' head(synthetic_data)
##' 
##' @export
#generate_synthetic_data = function(n_plants, baseline_values, within_variance) {
#  # Number of traits is determined by the number of columns in baseline_values
#  n_traits = ncol(baseline_values)
#  n_environments = nrow(baseline_values)
#  
#  # Initialize empty data frame to store the results
#  synthetic_data = data.frame(matrix(nrow = n_plants * n_environments, ncol = n_traits+1))
#  
#  # Generate row names where the the integer part is the plant id and the fractional part is the environment indicator 
#  # this can potentially be left out but I thought is nice to have
#  rownames(synthetic_data) = rep(1:n_plants, times = n_environments) + rep(seq(0.1, (n_environments - 1) / 10 + 0.1, by = 0.1), each = n_plants)
#  
#  # Set column names for the traits
#  colnames(synthetic_data) = c("Env_indicator",paste("Trait", 1:n_traits, sep = "_"))
#  
#  for (env in 1:n_environments) {
#    # Calculate the row indices for the current environment. Like this it is easier to batch generate the normally distributed values.
#    start_idx = (env - 1) * n_plants + 1
#    end_idx = env * n_plants
#    for (i in 1:n_traits) {
#      # Generate random data for the current trait and environment with specific variance
#      synthetic_data[start_idx:end_idx, 1] = env
#      synthetic_data[start_idx:end_idx, i+1] = rnorm(n_plants, mean = baseline_values[env, i], sd = within_variance[env, i])
#    }
#  }
#  
#  
#  return(synthetic_data)
#}
#
#
#
#n_plants = 100
#
## Define baseline values for each trait in each environment
#baseline_values = matrix(
#  c(100, 110, 120,  
#    50, 195, 205,  
#    70, 310, 290), 
#  nrow = 3, ncol = 3, byrow = TRUE
#)
#
## Define within-environment variance for each trait in each environment
#within_variance = matrix(
#  c(10, 15, 10,   
#    8, 12, 10,    
#    5, 7, 10),     
#  nrow = 3, ncol = 3, byrow = TRUE
#)
#
#
#set.seed(12345)
#synthetic_data1 = generate_synthetic_data(n_plants, baseline_values, within_variance)
#print(class(synthetic_data1))
#
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
calculate_CVt = function(data, exclude_col=NULL, traitwise = TRUE) {
  if (any(is.na(data))) {
    stop("There are missing values in your data. Consider using the function impute().")
  }
  
  # Exclude specified columns if needed
  if (!is.null(exclude_col)) {
    data = data[ , -exclude_col]
  }
  
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


### test - passed in by hand calculation with this dataset

#df_test1 = data.frame(Column1 = c(rep(4, 10), rep(2, 10)), Column2 = c(rep(10, 10), rep(1, 10)))
#
##traitwise
#calculate_CVt(df_test1)
#print(sd(df_test1[,1])/mean(df_test1[,1]))
#print(sd(df_test1[,2])/mean(df_test1[,2]))
#
##total
#df_test1_flattened=as.numeric(unlist(df_test1))
#calculate_CVt(df_test1,traitwise = F)
#print(sd(df_test1_flattened)/mean(df_test1_flattened))





################################

#' Calculate the Reaction Norm Slope for a Specific Trait
#'
#' This function calculates the reaction norm slope for a specified trait in a data frame, 
#' using the first column as the environment indicator. It also provides an option to plot the reaction norm.
#'
#' @param data A data frame containing the environment indicators in the first column and the trait data in the remaining columns.
#' @param trait_col The column number of the trait for which the reaction norm slope is to be calculated. This can also be multiple indicator numbers in a vector which then will be used separately to calculate the slope.
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
calculate_reaction_norm_slope = function(data, env_col, trait_cols, plot = FALSE) {
  
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_indicators = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_indicators = env_col
  } else {
    env_indicators = data[[env_col]]
  }
  
  # Initialize an empty list to store slopes for multiple traits
  slopes_list = list()
  
  # Loop over each specified trait column
  for (trait_col in trait_cols) {
    
    if (!is.numeric(data[[trait_col]])) {
      stop(paste("The column", colnames(data)[trait_col], "is not numeric."))
    }
    
    if (any(is.na(data[[trait_col]]))) {
      stop(paste("There are missing values in column", colnames(data)[trait_col], ". Consider using the function impute()."))
    }
    
    if (length(env_indicators) != nrow(data)) {
      stop("Length of env_indicators must match the number of rows in the data.")
    }
    
    # Fit a linear model to estimate the slope for the trait column
    model = lm(data[[trait_col]] ~ env_indicators)
    
    # Extract the slope of the reaction norm
    slope = coef(model)[["env_indicators"]]
    
    # Store the slope in the list with the trait name
    slopes_list[[colnames(data)[trait_col]]] = slope
    
    # Plot the reaction norm if requested
    if (plot) {
      par(mai = c(1, 1, 0.5, 0.5))  # Adjust the margins for plotting
      plot(data[[trait_col]] ~ env_indicators, 
           xlab = "Environment Indicator", 
           ylab = colnames(data)[trait_col],
           main = paste("Reaction Norm for", colnames(data)[trait_col]),
           pch = 19, col = "blue")
      abline(model, col = "red", lwd = 2)
      legend("topleft", legend = c("Observed", "Fitted Line"), col = c("blue", "red"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
    }
  }
  
  return(slopes_list)  # Return the list of slopes for each trait
}

#test - passed with synthetic dataset

#df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 10), rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
#calculate_reaction_norm_slope(df_test2,env_col = 1,trait_cols = c(2,3),plot = T)

################################


#NOTE: if the resource availability is an actual measurement of a metabolite which is being used by the plant then the grouping of the plants into high vs low resource availability should be done by clustering.
# Sadly I am missing the insight into the common practices in the field.


#' Calculate the D Slope for a Specific Trait Between High and Low Resource Availability
#'
#' This function calculates the D slope, which is the difference in the mean value of a specified trait 
#' between high and low resource availability conditions. The D slope quantifies the scope of the plastic response.
#' 
#' The function assumes that the environmental indicator is in the `env_col` column of the data frame or passed as a vector. 
#' The function will automatically assume that the lowest values represent the low resource environment and 
#' the highest values represent the high resource environment.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param env_col The column number, name, or a vector representing the environmental conditions. Defaults to 1 if not specified.
#' @param trait_col The column number or name of the trait to analyze.
#' @return The D slope, representing the difference in mean trait values between high and low resource conditions.
#' @examples
#' df = data.frame(
#'   Environment = c(1, 2, 3, 4),
#'   Height = c(10, 12, 20, 22)
#' )
#' D_slope = calculate_D_slope(df, trait_col = "Height")
#' print(D_slope)
#' 
#' # With an explicit environment vector
#' env_vector = c("Low", "Low", "High", "High")
#' D_slope = calculate_D_slope(df, trait_col = "Height", env_col = env_vector)
#' print(D_slope)
#' @export
calculate_D_slope = function(data, env_col, trait_cols) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    # If env_col is a vector, keep it as is
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Combine env_col with the data to maintain alignment after sorting
  sorted_data = data[order(env_col), ]
  sorted_env_col = env_col[order(env_col)]
  
  # Determine the midpoint for splitting
  num_rows = nrow(sorted_data)
  mid_point = ceiling(num_rows / 2)
  
  # Assign "Low" to the first half and "High" to the second half
  env_vector = c(rep("Low", mid_point), rep("High", num_rows - mid_point))
  
  # Initialize a vector to store the D slope for each trait
  D_slope_values = numeric(length(trait_cols))
  names(D_slope_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column to calculate the D slope
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Extract the trait data aligned with "High" and "Low" conditions
    trait_high = sorted_data[env_vector == "High", trait_column]
    trait_low = sorted_data[env_vector == "Low", trait_column]
    
    # Calculate the mean trait value for each condition
    mean_high = mean(trait_high, na.rm = TRUE)
    mean_low = mean(trait_low, na.rm = TRUE)
    
    # Calculate the D slope for the current trait
    D_slope_values[i] = mean_high - mean_low
  }
  
  return(D_slope_values)
}

#test - passed with synthetic 

#calculate_D_slope(df_test2,env_col = 1,trait_cols = 2)
#df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
#calculate_D_slope(df_test3,env_col = 1,trait_cols = 3)
#mean(df_test3[1:15,3])-mean(df_test3[16:nrow(df_test3),3])




################################

#' Calculate the Response Coefficient (RC) for One or Multiple Traits
#'
#' This function calculates the Response Coefficient (RC), which is the ratio of the mean value
#' of specified traits between high and low resource availability conditions. It supports calculating
#' the RC for one or more traits at once.
#'
#' The function assumes that the environmental indicator is in the `env_col` column of the data frame or passed as a vector.
#' The function automatically assumes that the lowest values represent the low resource environment and
#' the highest values represent the high resource environment.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param env_col The column number, name, or a vector representing the environmental conditions. Defaults to 1 if not specified.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @return A named numeric vector of Response Coefficients (RC), representing the ratio of mean trait values between high and low resource conditions for each trait.
#' @examples
#' df = data.frame(
#'   Environment = c(1, 2, 1, 2),
#'   Height = c(10, 12, 20, 22),
#'   Weight = c(30, 40, 50, 60)
#' )
#' RC_values = calculate_RC(df, env_col = "Environment", trait_cols = c("Height", "Weight"))
#' print(RC_values)
#'
#' # With an explicit environment vector
#' env_vector = c("Low", "Low", "High", "High")
#' RC_values = calculate_RC(df, env_col = env_vector, trait_cols = c("Height", "Weight"))
#' print(RC_values)
#' @export
calculate_RC = function(data, env_col = 1, trait_cols) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Sort the data based on the environmental indicator
  sorted_indices = order(env_col)
  sorted_data = data[sorted_indices, ]
  
  # Determine the midpoint for splitting
  num_rows = nrow(sorted_data)
  mid_point = ceiling(num_rows / 2)
  
  # Assign "Low" to the first half and "High" to the second half
  env_vector = c(rep("Low", mid_point), rep("High", num_rows - mid_point))
  
  # Initialize a vector to store the RC for each trait
  RC_values = numeric(length(trait_cols))
  names(RC_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column to calculate the RC
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Extract the relevant data for the "High" and "Low" conditions
    trait_high = sorted_data[[trait_column]][env_vector == "High"]
    trait_low = sorted_data[[trait_column]][env_vector == "Low"]
    
    # Calculate the mean trait value for each condition
    mean_high = mean(trait_high, na.rm = TRUE)
    mean_low = mean(trait_low, na.rm = TRUE)
    
    # Calculate the Response Coefficient (RC) for the current trait
    RC_values[i] = mean_high / mean_low
  }
  
  return(RC_values)
}


# test - passed with synthetic dataset

#calculate_RC(df_test3,env_col = 1,trait_cols = 3)
#
#mean(df_test3[1:15,3])/mean(df_test3[16:nrow(df_test3),3])


################################

#' Calculate the Coefficient of Variation of Means (CVm) for One or Multiple Traits
#'
#' This function calculates the Coefficient of Variation of Means (CVm),
#' which is the standard deviation of the means divided by the mean of the means across different environments.
#'
#' The function allows for flexible grouping of the data. The grouping can be specified using `env_col`, which can be
#' either a column index/name from the data or an external grouping vector. It supports calculating the CVm for one or
#' multiple traits at once.
#'
#' @param data A data frame containing the trait data.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @param env_col The column number, name, or a vector representing the grouping (e.g., environmental conditions).
#' @return A named numeric vector where each element represents the CVm for a single trait.
#' @examples
#' df = data.frame(
#'   Environment = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3"),
#'   Height = c(10, 12, 20, 22, 15, 18),
#'   Weight = c(100, 105, 110, 120, 115, 125)
#' )
#' CVm_values = calculate_CVm(df, trait_cols = c("Height", "Weight"), env_col = "Environment")
#' print(CVm_values)
#' @export
calculate_CVm = function(data,  env_col = 1,trait_cols) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a vector to store the CVm for each trait
  CVm_values = numeric(length(trait_cols))
  names(CVm_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column to calculate the CVm
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Calculate the means for each group (environment)
    means = tapply(data[[trait_column]], env_col, mean, na.rm = TRUE)
    
    # Calculate the standard deviation of the means
    sd_of_means = sd(means, na.rm = TRUE)
    
    # Calculate the mean of the means
    mean_of_means = mean(means, na.rm = TRUE)
    
    # Calculate the CVm for the current trait
    CVm_values[i] = sd_of_means / mean_of_means
  }
  
  return(CVm_values)
}


## test - passed with synthetic dataset

#calculate_CVm(df_test3,env_col = 1,trait_col = 3)
#
#sd(c(mean(df_test3[1:10,3]),mean(df_test3[11:20,3]),mean(df_test3[21:30,3])))/mean(mean(df_test3[,3]))


###############################


#' Calculate the Coefficient of Variation of Medians (CVmd) for One or Multiple Traits
#'
#' This function calculates the Coefficient of Variation of Medians (CVmd),
#' which is the standard deviation of the medians divided by the mean of the medians across different environments or groups.
#'
#' The function allows for flexible grouping of the data. The grouping can be specified using `env_col`, which can be
#' either a column index/name from the data or an external grouping vector. It supports calculating the CVmd for one or
#' multiple traits at once.
#'
#' @param data A data frame containing the trait data.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @param env_col The column number, name, or a vector representing the grouping (e.g., environmental conditions).
#' @return A named numeric vector where each element represents the CVmd for a single trait.
#' @examples
#' df = data.frame(
#'   Environment = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3"),
#'   Height = c(10, 12, 20, 22, 15, 18),
#'   Weight = c(100, 105, 110, 120, 115, 125)
#' )
#'
#' # Calculate CVmd for both "Height" and "Weight" using the "Environment" column
#' CVmd_values = calculate_CVmd(df, trait_cols = c("Height", "Weight"), env_col = "Environment")
#' print(CVmd_values)
#'
#' @export
calculate_CVmd = function(data, trait_cols, env_col) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Initialize a vector to store the CVmd for each trait
  CVmd_values = numeric(length(trait_cols))
  names(CVmd_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column to calculate the CVmd
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Calculate the medians for each group (environment)
    medians = tapply(data[[trait_column]], env_col, median, na.rm = TRUE)
    
    # Calculate the standard deviation of the medians
    sd_of_medians = sd(medians, na.rm = TRUE)
    
    # Calculate the mean of the medians
    mean_of_medians = mean(medians, na.rm = TRUE)
    
    # Calculate the CVmd for the current trait
    CVmd_values[i] = sd_of_medians / mean_of_medians
  }
  
  return(CVmd_values)
}


# test - passed on synthetic dataset 

#calculate_CVmd(df_test3,env_col = 1,trait_col = 2)
#
#sd(c(median(df_test3[1:10,2]),median(df_test3[11:20,2]),median(df_test3[21:30,2])))/mean(c(median(df_test3[1:10,2]),median(df_test3[11:20,2]),median(df_test3[21:30,2])))

###############################

#' Calculate Grand Plasticity (Pi) using Adjusted Means for One or Multiple Traits
#'
#' This function calculates the grand plasticity (Pi) for specified traits across different treatments,
#' correcting for a covariate (e.g., size effect like biomass) using adjusted means. The Pi is calculated
#' as the coefficient of variation (CV) of the adjusted means (standard deviation of adjusted means
#' divided by the grand mean of adjusted means) across the treatments.
#'
#' @param data A data frame containing the trait data, environmental conditions, and covariate.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @param env_col The column number or name representing the environmental conditions (must include control).
#' @param covariate_col The column number, name, or a vector representing the covariate (e.g., biomass).
#' @param control_env The label or value in `env_col` that represents the control environment.
#' @return A named numeric vector where each element represents the grand plasticity (Pi) for a single trait.
#' @examples
#' df = data.frame(
#'   Environment = c("Control", "Control", "Treatment1", "Treatment2"),
#'   Height = c(10, 12, 20, 22),
#'   Weight = c(5, 6, 7, 8),
#'   Biomass = c(5, 6, 7, 8)
#' )
#' grand_plasticity = calculate_grand_plasticity(df, trait_cols = c("Height", "Weight"), env_col = "Environment",
#'   covariate_col = "Biomass", control_env = "Control")
#' print(grand_plasticity)
#' @export
calculate_grand_plasticity = function(data, trait_cols, env_col, covariate_col, control_env) {
  
  # Ensure required libraries are loaded
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("The 'emmeans' package is required but not installed.")
  }
  
  # Handle env_col (as column index, name, or vector)
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col_data = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col_data = env_col
  } else {
    env_col_data = data[[env_col]]
  }
  
  env_col_data = factor(env_col_data)
  
  # Handle covariate_col (as column index, name, or vector)
  if (is.numeric(covariate_col) && length(covariate_col) == 1) {
    covariate_data = data[[covariate_col]]
  } else if (is.vector(covariate_col) && length(covariate_col) == nrow(data)) {
    covariate_data = covariate_col
  } else {
    covariate_data = data[[covariate_col]]
  }
  
  
  # Initialize a vector to store grand plasticity for each trait
  grand_plasticity_values = numeric(length(trait_cols))
  names(grand_plasticity_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column to calculate grand plasticity
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Handle trait_col dynamically (whether it's a column number or name)
    trait_col_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Fit a linear model to account for the covariate (adjusted means)
    model = lm(trait_col_data ~ covariate_data + env_col_data, data = data)
    
    # Use emmeans to get the adjusted means for each environment
    adjusted_means = emmeans::emmeans(model, ~ env_col_data)
    
    # Extract the adjusted means for treatments (excluding control)
    adjusted_means_summary = summary(adjusted_means)
    treatment_means = adjusted_means_summary[,2]
    
    # Check if we successfully extracted the means
    if (length(treatment_means) == 0) {
      stop(paste("Could not extract treatment means for trait:", trait_column))
    }
    
    # Calculate the standard deviation and grand mean of the adjusted means
    sd_means = sd(treatment_means, na.rm = TRUE)
    
    grand_mean = mean(treatment_means, na.rm = TRUE)
    
    # Calculate grand plasticity (Pi) as the Coefficient of Variation (CV)
    grand_plasticity_values[i] = sd_means / grand_mean
  }
  
  return(grand_plasticity_values)
}

## test - passed on synthetic dataset

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
  Column1 = c(1, 1, 2, 2),      # Control (2) and Treatment (3)
  Column0 = c(10, 12, 20, 22),  # Response variable (trait)
  Column2 = c(1, 1, 2, 3)       # Covariate (e.g., biomass)
)

#calculate_grand_plasticity(df_test4,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)
#calculate_grand_plasticity(df_test5,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)
#calculate_grand_plasticity(df_test_simple,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)

##########################################
  
#' Combine Internal and External Factors into a Single Dataframe and Return Levels
#'
#' This function combines internal and external factors into a single dataframe by creating 
#' an interaction between them. It returns a modified dataframe with a new column `Combined_Factors` 
#' that represents the interaction of these factors and a list of levels for the combined factors.
#'
#' @param dataframe A data frame containing the data to be analyzed.
#' @param factors (Optional) A vector of integers or strings specifying the columns that contain the internal factors in the dataframe. 
#' These can be the names or indices of the columns in the `dataframe`.
#' @param factors_not_in_dataframe (Optional) A list of vectors representing external factors that are not included in the `dataframe`. 
#' These vectors should have the same length as the `dataframe` and correspond to the environmental conditions for each observation. 
#' @return A list containing:
#'   \item{dataframe}{The modified dataframe with a new column `Combined_Factors` representing the interaction of the provided factors.}
#'   \item{levels}{A list of the levels for the combined factors.}
#' @examples
#' # Example with internal factors
#' df = data.frame(
#'   Light = rep(c("Low", "High"), times = 6),
#'   Water = rep(c("Low", "High"), each = 6)
#' )
#' result = combine_factors(df, factors = c("Light", "Water"))
#' head(result$dataframe)
#' print(result$levels)
#' 
#' # Example with external factors
#' external_light = rep(c(0.4, 0.6, 0.8), each = 4)
#' result_external = combine_factors(df, factors_not_in_dataframe = list(external_light))
#' head(result_external$dataframe)
#' print(result_external$levels)
#' @export
combine_factors = function(dataframe, factors = NULL, factors_not_in_dataframe = NULL) {
  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    # Ensure the lengths match
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    # Create a data frame for external factors
    external_factors_df = as.data.frame(factors_not_in_dataframe)
    
    # If there are internal factors, combine them with external factors
    if (!is.null(factors)) {
      factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
      combined_factors_df = cbind(dataframe[factors], external_factors_df)
    } else {
      combined_factors_df = external_factors_df
    }
    
    # Create a combined factor interaction
    dataframe$Combined_Factors = interaction(combined_factors_df, drop = TRUE)
  } else if (!is.null(factors)) {
    # If only internal factors are provided
    factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
    dataframe$Combined_Factors = interaction(dataframe[factors], drop = TRUE)
  } else {
    stop("You must provide either internal factors, external factors, or both.")
  }
  
  # Ensure Combined_Factors is treated as a factor
  dataframe$Combined_Factors = as.factor(dataframe$Combined_Factors)
  levels=levels(dataframe$Combined_Factors)
  print(levels)
  return(dataframe)
}



## Example usage with synthetic data
#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#external_water = sample(rep(c("Low", "High"), each = 150))
#
#combined_factors= combine_factors(synthetic_data1, factors_not_in_dataframe = list(external_light,external_water), factors=1)

##########################################

#' Calculate the Phenotypic Plasticity Index (PPF) Based on Least Square Means
#'
#' This function calculates the Phenotypic Plasticity Index (PPF), which quantifies the phenotypic plasticity 
#' of one or multiple traits between two environments or environmental groupings. The PPF is calculated using 
#' the least square means (LSMs) of the traits in each environment, with optional adjustment for covariates 
#' that may influence the traits.
#'
#' Covariates can be included to adjust the model for additional environmental influences or biological factors 
#' that might affect the traits. If `env_pairs` is not provided, the function calculates the PPF for all possible 
#' combinations of environments.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @param env_col The column number, name, or a vector representing the environmental conditions.
#' @param covariate_col (Optional) The column number, name, or vector to include as a covariate in the model.
#' @param env_pairs (Optional) A list of environment pairs to calculate the PPF for. Each pair should be a list of two elements (e.g., `list(2, 3)`).
#' If `env_pairs` is NULL (default), the function calculates the PPF for all possible combinations of environments.
#' @return A named list of PPF values for each trait, with the environment pair as the names of each entry.
#' @examples
#' df = data.frame(
#'   Environment = rep(c("Env1", "Env2", "Env3"), each = 10),
#'   Height = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 20, 22, 21, 23, 24, 25, 23, 24, 22, 23, 18, 19, 21, 22, 20, 23, 24, 25, 23, 22),
#'   Biomass = c(3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
#' )
#' 
#' # Calculate PPF for all environment combinations
#' PPF_all = calculate_PPF(df, trait_cols = "Height", env_col = "Environment", covariate_col = "Biomass")
#' print(PPF_all)
#' 
#' # Calculate PPF for specific environment pairs
#' PPF_specific = calculate_PPF(df, trait_cols = "Height", env_col = "Environment", 
#'                              covariate_col = "Biomass", env_pairs = list(list("Env1", "Env2"), list("Env2", "Env3")))
#' print(PPF_specific)
#' @export
calculate_PPF = function(data, trait_cols, env_col, covariate_col = NULL, env_pairs = NULL) {
  
  # Handle env_col (as column index, name, or vector)
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col_data = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col_data = env_col
  } else {
    env_col_data = data[[env_col]]
  }
  env_col_data = as.factor(env_col_data)
  
  # Handle covariate_col (either as column or external vector)
  if (!is.null(covariate_col)) {
    if (is.numeric(covariate_col) && length(covariate_col) == 1) {
      covariate_data = data[[covariate_col]]
    } else if (is.vector(covariate_col) && length(covariate_col) == nrow(data)) {
      covariate_data = covariate_col
    } else {
      covariate_data = data[[covariate_col]]
    }
  }
  
  # Default: if no env_pairs are passed, calculate for all pairs of environments
  if (is.null(env_pairs)) {
    env_pairs = combn(levels(env_col_data), 2, simplify = FALSE)
  } else if (is.vector(env_pairs) && length(env_pairs) == 2) {
    env_pairs = list(env_pairs)  # Convert a single pair into a list of one pair
  }
  
  # Initialize a list to store PPF values for each environment pair
  PPF_values_list = list()
  # Loop over each trait
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    trait_col_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Store PPF values for this trait
    PPF_values_for_trait = numeric(length(env_pairs))
    names(PPF_values_for_trait) = sapply(env_pairs, function(pair) paste(pair, collapse = "-"))
    
    # Loop through each environment pair
    for (j in seq_along(env_pairs)) {
      env_pair = env_pairs[[j]]
      
      # Subset the data to only include the current environment pair
      subset_idx = which(env_col_data %in% env_pair)
      subset_data = data[subset_idx, ]
      subset_env_col = env_col_data[subset_idx]
      subset_trait_col_data = trait_col_data[subset_idx]
      
      # Fit the linear model
      if (is.null(covariate_col)) {
        model = lm(subset_trait_col_data ~ subset_env_col, data = subset_data)
      } else {
        subset_covariate_data = covariate_data[subset_idx]
        model = lm(subset_trait_col_data ~ subset_env_col + subset_covariate_data, data = subset_data)
      }
      
      # Calculate least square means (LSMs) for the current environment pair
      lsm = as.data.frame(emmeans::emmeans(model, ~ subset_env_col))
      
      # Calculate PPF for the current pair using the formula: 100 x ((LSM in one environment - LSM in the other) / LSM in the first environment)
      PPF_values_for_trait[j] = 100 * abs((lsm[1, "emmean"] - lsm[2, "emmean"]) / lsm[1, "emmean"])
    }
    
    # Add PPF values for the current trait to the list
    PPF_values_list[[i]] = PPF_values_for_trait
  }
  
  return(PPF_values_list)
}

#
#synthetic_data2=combine_factors(synthetic_data1,factors=NULL, factors_not_in_dataframe=list(external_water))
#
### test - passed on synthetic dataset
#
#df_test6 = data.frame(
#  Column1 = c(rep(3, 15), rep(2, 15)),   # Response variable
#  Column2 = c(rep(4, 15), rep(2, 15)),   # Control (2) and Treatment (3)
#  Column3 = c(rep(3, 10), rep(2, 10), rep(1, 10))    # Covariate (matches values of Column0)
#)
#df_test6$Column1 = as.factor(df_test6$Column1)
#calculate_PPF(df_test6,env_col = 1, trait_cols = 2, env_pairs = list(2,3))
#
#model = lm(Column2 ~ Column1 , data = df_test6)
#lsmeans_env = emmeans(model, ~ Column1)
#summary_lsmeans = summary(lsmeans_env)
## Extract only the LSMeans (adjusted means)
#lsmeans_values = summary_lsmeans$emmean
#100*abs((lsmeans_values[[1]]-lsmeans_values[[2]])/lsmeans_values[[1]])
#
#
#model = lm(Column3 ~ Column1 , data = df_test6)
#lsmeans_env = emmeans(model, ~ Column1)
#summary_lsmeans = summary(lsmeans_env)
## Extract only the LSMeans (adjusted means)
#lsmeans_values = summary_lsmeans$emmean
#100*abs((lsmeans_values[[1]]-lsmeans_values[[2]])/lsmeans_values[[1]])


###########################################


#' Calculate the Phenotypic Plasticity Index (Pi)
#'
#' This function calculates the Phenotypic Plasticity Index (Pi), defined as the 
#' difference between the maximum and minimum values of a trait divided by the maximum value.
#' The function can be applied to multiple traits at once if provided with multiple trait columns.
#'
#' @param data A numeric vector, data frame, or matrix containing the trait data. 
#' If a data frame or matrix is provided, each column represents a trait.
#' @param trait_cols (Optional) A vector of column numbers or names of the traits to analyze 
#' if `data` is a data frame or matrix. If not specified, it calculates Pi for the entire data vector.
#' @return A named numeric vector of Phenotypic Plasticity Index (Pi) values for each trait.
#' @examples
#' df = data.frame(
#'   Trait1 = c(10, 15, 12, 14, 18),
#'   Trait2 = c(20, 25, 22, 23, 27),
#'   Trait3 = c(30, 35, 32, 34, 37)
#' )
#' Pi_values = calculate_Phenotypic_Plasticity_Index(df, trait_cols = c("Trait1", "Trait2"))
#' print(Pi_values)
#' @export
calculate_Phenotypic_Plasticity_Index = function(data, trait_cols = NULL) {
  if (is.null(trait_cols)) {
    # If trait_col is not provided, assume data is a numeric vector
    max_value = max(data, na.rm = TRUE)
    min_value = min(data, na.rm = TRUE)
    
    # Calculate Pi for the single vector
    Pi = (max_value - min_value) / max_value
    return(Pi)
  } else {
    # If trait_cols are provided, assume data is a data frame or matrix
    Pi_values = numeric(length(trait_cols))  # Initialize vector to store Pi for each trait
    names(Pi_values) = trait_cols  # Name the vector by trait columns
    
    for (i in seq_along(trait_cols)) {
      trait_column = trait_cols[i]
      trait_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
      
      max_value = max(trait_data, na.rm = TRUE)
      min_value = min(trait_data, na.rm = TRUE)
      
      # Calculate Pi as (Max - Min) / Max
      Pi_values[i] = (max_value - min_value) / max_value
    }
    
    return(Pi_values)
  }
}


##test - passed on a synthetic dataset (look at the dataset for confirmation)

#calculate_Phenotypic_Plasticity_Index(df_test6,trait_col = 2)




####################################

#' Calculate the Proportional Inter-Median Difference (PImd) for One or Multiple Traits
#'
#' This function calculates the Proportional Inter-Median Difference (PImd), 
#' defined as the difference between the maximum and minimum medians of a trait 
#' across different environmental conditions divided by the maximum median. 
#' It provides a measure of the relative variability in the trait across the specified conditions.
#' The function can calculate PImd for multiple traits at once.
#'
#' @param data A data frame containing the trait and environmental condition data.
#' @param trait_cols A vector of column numbers or names of the traits to analyze. 
#' It can be a column index (integer) or a column name (string) within the data frame.
#' @param env_col The column number, name, or a vector representing the environmental conditions. 
#' It can be a column index (integer), a column name (string), or a vector of values.
#' @return A named numeric vector where each element represents the PImd value for a single trait.
#' @examples
#' # Example usage with a data frame
#' data = data.frame(
#'   Trait1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'   Trait2 = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
#'   Env = factor(c("A", "A", "B", "B", "C", "C", "A", "B", "C", "A"))
#' )
#' result = calculate_PImd(data, trait_cols = c("Trait1", "Trait2"), env_col = "Env")
#' print(result)
#' @export
calculate_PImd = function(data, trait_cols, env_col) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    # If env_col is a vector, keep it as is
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Ensure env_col is a factor to use levels
  if (!is.factor(env_col)) {
    env_col = factor(env_col)
  }
  
  # Initialize a vector to store PImd values for each trait
  PImd_values = numeric(length(trait_cols))
  names(PImd_values) = trait_cols
  
  # Loop over each trait column to calculate PImd
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    trait_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Get levels of env_col
    levels = levels(env_col)
    medians = c()
    
    # Calculate medians for each environment
    for (level in levels) {
      medians = c(medians, median(trait_data[env_col == level], na.rm = TRUE))
    }
    
    # Calculate PImd as (Max - Min) / Max
    PImd_values[i] = (max(medians) - min(medians)) / max(medians)
  }
  
  return(PImd_values)
}

##test - passed on synthetic dataset (check  dataset for confirmation)

#calculate_PImd(df_test6,trait_cols = c(2,3),env_col = 1)

###############################################



#' Calculate the Proportional Inter-Least Square Mean Difference (PILSM)
#' for One or Multiple Traits and Optionally Plot LSMs Against the Trait Data
#'
#' This function calculates the Proportional Inter-Least Square Mean Difference (PILSM), 
#' defined as the difference between the maximum and minimum least square means (LSMs) 
#' of a trait across different environments, divided by the maximum LSM.
#' Optionally, it can plot the LSMs against the original trait data.
#'
#' @param data A data frame containing the trait data and environmental indicators.
#' @param trait_cols A vector of column numbers or names of the traits to analyze.
#' @param env_col The column number, name, or a vector representing the environmental conditions.
#' @param covariates (Optional) A vector of column names or indices to include as covariates in the model.
#' @param plot (Optional) Logical indicating whether to plot the LSMs against the original data (default is FALSE).
#' @return A list with the PILSM values for each trait and the optional LSM plot.
#' @examples
#' df = data.frame(
#'   Environment = rep(c("Env1", "Env2", "Env3"), each = 10),
#'   Height = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'             20, 22, 21, 23, 24, 25, 23, 24, 22, 23,
#'             30, 32, 31, 33, 34, 35, 33, 34, 32, 33),
#'   SoilQuality = c(rep(3, 10), rep(4, 10), rep(5, 10))
#' )
#' 
#' # Calculate PILSM with a covariate and plot
#' result = calculate_PILSM(df, trait_cols = c("Height"), env_col = "Environment", covariates = "SoilQuality", plot = TRUE)
#' print(result)
#' @export
calculate_PILSM = function(data, trait_cols, env_col, covariates = NULL, plot = FALSE) {
  
  # Initialize a vector to store PILSM values for each trait
  PILSM_values = numeric(length(trait_cols))
  names(PILSM_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column
  for (i in seq_along(trait_cols)) {
    trait_col = trait_cols[i]
    
    # Handle env_col
    if (is.numeric(env_col) && length(env_col) == 1) {
      env_col_data = data[[env_col]]
    } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
      env_col_data = env_col
    } else {
      env_col_data = data[[env_col]]
    }
    
    trait_col_data = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
    env_col_data = as.factor(env_col_data)
    
    # Fit the linear model
    if (is.null(covariates)) {
      model = lm(trait_col_data ~ env_col_data, data = data)
    } else {
      covariates = if (is.numeric(covariates)) names(data)[covariates] else covariates
      formula = as.formula(paste("trait_col_data ~ env_col_data +", paste(covariates, collapse = " + ")))
      model = lm(formula, data = data)
    }
    
    # Calculate least square means (LSMs) for each environment
    lsm = as.data.frame(emmeans::emmeans(model, ~ env_col_data))
    
    # Calculate PILSM
    max_lsm = max(lsm$emmean, na.rm = TRUE)
    min_lsm = min(lsm$emmean, na.rm = TRUE)
    PILSM_values[i] = (max_lsm - min_lsm) / max_lsm
    
    # Optionally plot (stored in list if plot is requested)
    if (plot) {
      plot_data = data.frame(Environment = env_col_data, Trait = trait_col_data)
      lsm_plot = ggplot2::ggplot(plot_data, ggplot2::aes(x = Environment, y = Trait)) +
        ggplot2::geom_point(color = "blue", alpha = 0.5) +
        ggplot2::geom_point(data = lsm, ggplot2::aes(x = env_col_data, y = emmean), color = "red", size = 3) +
        ggplot2::geom_line(data = lsm, ggplot2::aes(x = env_col_data, y = emmean), color = "red", linetype = "dashed") +
        ggplot2::labs(title = paste("LSMs for", trait_col), x = "Environment", y = trait_col) +
        ggplot2::theme_minimal()
      
      # Print the plot for each trait
      print(lsm_plot)
    }
  }
  
  return(PILSM_values)
}



## test - passed on synthetic dataset 


#calculate_PILSM(df_test6,trait_col=c(2,3),env_col=1,plot = F)
#
#model = lm(Column2 ~ Column1 , data = df_test6)
#lsmeans_env = emmeans(model, ~ Column1)
#summary_lsmeans = summary(lsmeans_env)
## Extract only the LSMeans (adjusted means)
#lsmeans_values = summary_lsmeans$emmean
#(max(lsmeans_values)-min(lsmeans_values))/max(lsmeans_values)

################################################


#' Calculate the Relative Trait Response (RTR) for One or Multiple Traits
#'
#' This function calculates the Relative Trait Response (RTR) score, defined as the 
#' difference between the mean trait value at one end of an environmental gradient 
#' and the mean trait value at the opposite end, divided by the absolute maximum value of the trait.
#'
#' @param data A data frame containing the trait data and environmental conditions.
#' @param trait_cols A vector of column names or numbers for the traits to analyze.
#' @param env_col The column name or number for the environmental conditions.
#' @param env_low The value of the environmental condition representing one end of the gradient.
#' @param env_high The value of the environmental condition representing the opposite end of the gradient.
#' @return A named numeric vector where each element represents the RTR value for a single trait.
#' @examples
#' df = data.frame(
#'   Environment = rep(c("Low", "High"), each = 10),
#'   Height = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'             20, 22, 21, 23, 24, 25, 23, 24, 22, 23),
#'   Weight = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
#'             15, 16, 17, 18, 19, 20, 21, 22, 23, 24)
#' )
#' 
#' RTR_values = calculate_RTR(df, trait_cols = c("Height", "Weight"), env_col = "Environment", env_low = "Low", env_high = "High")
#' print(RTR_values)
#' @export
calculate_RTR = function(data, trait_cols, env_col, env_low, env_high) {
  
  # Initialize a vector to store RTR values for each trait
  RTR_values = numeric(length(trait_cols))
  names(RTR_values) = trait_cols  # Name the vector by trait columns
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Loop over each trait column to calculate RTR
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Handle trait_col dynamically (whether it's a column number or name)
    trait_col_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Subset data for the specified environmental conditions
    data_low = data[env_col == env_low, ]
    data_high = data[env_col == env_high, ]
    
    # Calculate mean trait values for each end of the gradient
    mean_low = mean(data_low[[trait_column]], na.rm = TRUE)
    mean_high = mean(data_high[[trait_column]], na.rm = TRUE)
    
    # Calculate the RTR value
    RTR_values[i] = (mean_high - mean_low) / max(abs(trait_col_data), na.rm = TRUE)
  }
  
  return(RTR_values)
}



## test - passed on synthetic dataset
#calculate_RTR(df_test2,trait_col=2,env_col=1,env_low=1,env_high=2)


######################################

#' @title Calculate Phenotypic Instability Ratio (PIR) for One or Multiple Traits
#'
#' @description This function calculates the Phenotypic Instability Ratio (PIR) for a given trait across different environments,
#' following the method described by Robinson (1989). PIR is calculated as the ratio of the difference between the maximum and minimum
#' mean trait values across environments to the mean trait value in the environment with the maximum relative growth rate (RGR).
#' The PIR metric provides an intermediate measure of plasticity and assumes normality in the data. It also requires prior knowledge of the
#' relative growth rates (RGRs) for each environment. Note that the method has statistical limitations.
#'
#' @param data A data frame containing the data for trait values and corresponding environment and RGR information.
#' @param trait_cols A vector of column names or numeric indices indicating the trait values. 
#' @param env_col A column name or numeric index indicating the environment identifiers. Can be a numeric index of the column in the data frame or a vector.
#' @param rgr_col A column name or numeric index indicating the relative growth rate (RGR) values. Can be a numeric index of the column in the data frame or a vector.
#'
#' @details The function first converts the environment column into a factor and computes the mean trait values for each environment.
#' It then identifies the environment with the maximum relative growth rate and uses the trait's mean value in that environment
#' to compute the PIR score. The calculation follows the formula:
#' \deqn{PIR = \frac{MaxMean - MinMean}{Mean_{MaxRGR}}}
#' where \code{MaxMean} is the maximum mean trait value, \code{MinMean} is the minimum mean trait value,
#' and \code{Mean_{MaxRGR}} is the mean trait value in the environment with the maximum RGR.
#'
#' @return A named numeric vector where each element represents the PIR for a specific trait.
#'
#' @references
#' Robinson, D. (1989). \emph{Plasticity in plant growth and resource use as a trait for crop breeding}. Field Crops Research, 11(2-3), 153-159.
#'
#' @examples
#' df = data.frame(
#'   Environment = rep(c("Env1", "Env2", "Env3"), each = 10),
#'   Height = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'             20, 22, 21, 23, 24, 25, 23, 24, 22, 23,
#'             30, 32, 31, 33, 34, 35, 33, 34, 32, 33),
#'   Weight = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
#'             15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
#'   RGR = c(1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
#'           1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 
#'           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
#' )
#' 
#' PIR_values = calculate_PIR(df, trait_cols = c("Height", "Weight"), env_col = "Environment", rgr_col = "RGR")
#' print(PIR_values)
#' @export
calculate_PIR = function(data, trait_cols, env_col, rgr_col) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Handle rgr_col
  if (is.numeric(rgr_col) && length(rgr_col) == 1) {
    rgr_col = data[[rgr_col]]
  } else if (is.vector(rgr_col) && length(rgr_col) == nrow(data)) {
    rgr_col = rgr_col
  } else {
    rgr_col = data[[rgr_col]]
  }
  
  # Convert env_col to a factor
  env_col = as.factor(env_col)
  
  # Initialize a vector to store PIR values for each trait
  PIR_values = numeric(length(trait_cols))
  names(PIR_values) = trait_cols  # Name the vector by trait columns
  
  # Loop over each trait column to calculate PIR
  for (i in seq_along(trait_cols)) {
    trait_column = trait_cols[i]
    
    # Handle trait_col dynamically (whether it's a column number or name)
    trait_col_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]
    
    # Calculate mean trait values for each environment
    means = tapply(trait_col_data, env_col, mean, na.rm = TRUE)
    
    # Identify maximum and minimum means
    max_mean = max(means, na.rm = TRUE)
    min_mean = min(means, na.rm = TRUE)
    
    # Identify the environment with the maximum growth rate
    max_rgr_env = levels(env_col)[which.max(tapply(rgr_col, env_col, mean, na.rm = TRUE))]
    
    # Find the mean of the trait at the environment where the growth rate is maximum
    mean_at_max_rgr = means[max_rgr_env]
    
    # Calculate PIR
    PIR_values[i] = (max_mean - min_mean) / mean_at_max_rgr
  }
  
  return(PIR_values)
}

## test - passed on synthetic dataset
#specific_growthrate=c(rep(10,10),rep(20,10))
#
#calculate_PIR(df_test2 , trait_col = 2 , env_col = 1, rgr_col = specific_growthrate)




