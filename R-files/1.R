#this files contains the following indices/functions: 
#generate_synthetic_data,
#impute,
#coefficient-of variation total (calculate_CVt) - tested,
#slope of norm reaction (calculate_reaction_norm_slope) - tested,
#slope of plastic response (D) (calculate_D_slope)- tested,
#response coefficient (RC) (calculate_RC)- tested,
#Standard deviation of means (CVm) (calculate_CVm)- tested,
#Standard deviation of medians (CVmd)(calculate_CVmd)- tested,
#Grand plasticity (calculate_GPi)- tested,
#combine_factors- tested,
#Phenotypic Plasticity Index (calculate_GPi)- tested,
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
  c(10, 15, 10,   
    8, 12, 10,    
    5, 7, 10),     
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
  if (any(is.na(data))) {
    stop("There are missing values in your data. Consider using the function impute().")
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

df_test1 = data.frame(Column1 = c(rep(4, 10), rep(2, 10)), Column2 = c(rep(10, 10), rep(1, 10)))

#traitwise
calculate_CVt(df_test1)
print(sd(df_test1[,1])/mean(df_test1[,1]))
print(sd(df_test1[,2])/mean(df_test1[,2]))

#total
df_test1_flattened=as.numeric(unlist(df_test1))
calculate_CVt(df_test1,traitwise = F)
print(sd(df_test1_flattened)/mean(df_test1_flattened))





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

df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 10), rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
calculate_reaction_norm_slope(df_test2,env_col = 1,trait_cols = c(2,3),plot = T)

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
calculate_D_slope = function(data, env_col, trait_col) {
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
  
  # Extract the trait data aligned with "High" and "Low" conditions
  trait_high = sorted_data[[trait_col]][env_vector == "High"]
  trait_low = sorted_data[[trait_col]][env_vector == "Low"]
  
  # Calculate the mean trait value for each condition
  mean_high = mean(trait_high, na.rm = TRUE)
  mean_low = mean(trait_low, na.rm = TRUE)
  
  # Calculate the D slope
  D_slope = mean_high - mean_low
  
  return(D_slope)
}

#test - passed with synthetic 

calculate_D_slope(df_test2,env_col = 1,trait_col = 2)
df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
calculate_D_slope(df_test3,env_col = 1,trait_col = 3)
mean(df_test3[1:15,3])-mean(df_test3[16:nrow(df_test3),3])




################################

#' Calculate the Response Coefficient (RC) for a Specific Trait
#'
#' This function calculates the Response Coefficient (RC), which is the ratio of the mean value 
#' of a specified trait between high and low resource availability conditions.
#' 
#' The function assumes that the environmental indicator is in the `env_col` column of the data frame or passed as a vector.
#' The function will automatically assume that the lowest values represent the low resource environment and 
#' the highest values represent the high resource environment.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param env_col The column number, name, or a vector representing the environmental conditions. Defaults to 1 if not specified.
#' @param trait_col The column number or name of the trait to analyze.
#' @return The Response Coefficient (RC), representing the ratio of mean trait values between high and low resource conditions.
#' @examples
#' df = data.frame(
#'   Environment = c(1, 2, 3, 4),
#'   Height = c(10, 12, 20, 22)
#' )
#' RC = calculate_RC(df, trait_col = "Height")
#' print(RC)
#' 
#' # With an explicit environment vector
#' env_vector = c("Low", "Low", "High", "High")
#' RC = calculate_RC(df, trait_col = "Height", env_col = env_vector)
#' print(RC)
#' @export
calculate_RC = function(data, env_col = 1, trait_col) {
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
  
  # Extract the relevant data for the "High" and "Low" conditions
  trait_high = sorted_data[[trait_col]][env_vector == "High"]
  trait_low = sorted_data[[trait_col]][env_vector == "Low"]
  
  # Calculate the mean trait value for each condition
  mean_high = mean(trait_high, na.rm = TRUE)
  mean_low = mean(trait_low, na.rm = TRUE)
  
  # Calculate the Response Coefficient (RC)
  RC = mean_high / mean_low
  
  return(RC)
}


# test - passed with synthetic dataset

calculate_RC(df_test3,env_col = 1,trait_col = 3)

mean(df_test3[1:15,3])/mean(df_test3[16:nrow(df_test3),3])


################################

#' Calculate the Coefficient of Variation of Means (CVm)
#'
#' This function calculates the Coefficient of Variation of Means (CVm), 
#' which is the standard deviation of the means divided by the mean of the means across different environments.
#'
#' The function allows for flexible grouping of the data. The grouping can be specified using `env_col`, which can be 
#' either a column index/name from the data or an external grouping vector.
#'
#' @param data A data frame containing the trait data.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number, name, or a vector representing the grouping (e.g., environmental conditions).
#' @return The CVm, representing the ratio of the standard deviation of means to the mean of means.
#' @examples
#' df = data.frame(
#'   Environment = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3"),
#'   Height = c(10, 12, 20, 22, 15, 18)
#' )
#' CVm = calculate_CVm(df, trait_col = "Height", env_col = "Environment")
#' print(CVm)
#' @export
calculate_CVm = function(data, trait_col, env_col = 1) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Calculate the means for each group
  means = tapply(data[[trait_col]], env_col, mean, na.rm = TRUE)
  
  # Calculate the standard deviation of the means
  sd_of_means = sd(means, na.rm = TRUE)
  
  # Calculate the mean of the means
  mean_of_means = mean(means, na.rm = TRUE)
  
  # Calculate the CVm
  CVm = sd_of_means / mean_of_means
  
  return(CVm)
}

## test - passed with synthetic dataset

calculate_CVm(df_test3,env_col = 1,trait_col = 3)

sd(c(mean(df_test3[1:10,3]),mean(df_test3[11:20,3]),mean(df_test3[21:30,3])))/mean(mean(df_test3[,3]))


###############################


#' Calculate the Coefficient of Variation of Medians (CVmd)
#'
#' This function calculates the Coefficient of Variation of Medians (CVmd), 
#' which is the standard deviation of the medians divided by the mean of the medians across different environments or groups.
#'
#' The function allows for flexible grouping of the data. The grouping can be specified using `env_col`, which can be 
#' either a column index/name from the data or an external grouping vector.
#'
#' @param data A data frame containing the trait data.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number, name, or a vector representing the grouping (e.g., environmental conditions).
#' @return The CVmd, representing the ratio of the standard deviation of medians to the mean of medians.
#' @examples
#' df = data.frame(
#'   Environment = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3"),
#'   Height = c(10, 12, 20, 22, 15, 18)
#' )
#' 
#' # Calculate CVmd using the environment column
#' CVmd_default = calculate_CVmd(df, trait_col = "Height", env_col = "Environment")
#' print(CVmd_default)
#' 
#' # Calculate CVmd using an external grouping vector
#' env_vector = c("Low", "Low", "High", "High", "Medium", "Medium")
#' CVmd_vector = calculate_CVmd(df, trait_col = "Height", env_col = env_vector)
#' print(CVmd_vector)
#' @export
calculate_CVmd = function(data, trait_col, env_col) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Calculate the medians for each group
  medians = tapply(data[[trait_col]], env_col, median, na.rm = TRUE)
  
  # Calculate the standard deviation of the medians
  sd_of_medians = sd(medians, na.rm = TRUE)
  
  # Calculate the mean of the medians
  mean_of_medians = mean(medians, na.rm = TRUE)
  
  # Calculate the CVmd
  CVmd = sd_of_medians / mean_of_medians
  
  return(CVmd)
}

# test - passed on synthetic dataset 

calculate_CVmd(df_test3,env_col = 1,trait_col = 2)

sd(c(median(df_test3[1:10,2]),median(df_test3[11:20,2]),median(df_test3[21:30,2])))/mean(median(df_test3[1:10,2]),median(df_test3[11:20,2]),median(df_test3[21:30,2]))

###############################

#' Calculate Plasticity (Pi) using Adjusted Means
#'
#' This function calculates the Plasticity (Pi), defined as the relative change
#' in a specified trait under different environmental conditions (treatment and control).
#' The plasticity is calculated as (Trait_treatment - Trait_control) / Trait_control,
#' using the adjusted means from a linear model that controls for a covariate (e.g., biomass).
#'
#' @param data A data frame containing the trait data, environmental conditions, and covariate.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number or name representing the environmental conditions (must include control).
#' @param covariate_col The column number or name representing the covariate (e.g., biomass).
#' @param control_env The label or value in `env_col` that represents the control environment.
#' @return The plasticity (Pi) value for the specified trait.
#' @examples
#' df = data.frame(
#'   Environment = c("Control", "Control", "Treatment", "Treatment"),
#'   Height = c(10, 12, 20, 22),
#'   Biomass = c(5, 6, 7, 8)
#' )
#' Pi = calculate_Pi_with_adjusted_means(df, trait_col = "Height", env_col = "Environment",
#'   covariate_col = "Biomass", control_env = "Control")
#' print(Pi)
#' @export
calculate_Pi_with_adjusted_means = function(data, env_col , trait_col, covariate_col, control_env) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  covariate_col = if (is.numeric(covariate_col)) data[[covariate_col]] else data[[covariate_col]]
  
  # Fit a linear model to account for the covariate (adjusted means)
  model = lm(trait_col ~ covariate_col + env_col, data = data)
  
  # Use emmeans to get the adjusted means for each environment
  library(emmeans)
  adjusted_means = emmeans(model, ~ env_col)
  
  # Extract the adjusted means for control and treatment
  trait_control_mean = summary(adjusted_means)$emmean[which(summary(adjusted_means)$env_col == control_env)]
  trait_treatment_mean = mean(summary(adjusted_means)$emmean[summary(adjusted_means)$env_col != control_env], na.rm = TRUE)
  
  # Calculate Plasticity (Pi) using the adjusted means
  Pi = abs((trait_treatment_mean - trait_control_mean) / trait_control_mean)
  
  return(Pi)
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
  Column1 = c(2, 2, 3, 3),      # Control (2) and Treatment (3)
  Column0 = c(10, 12, 20, 22),  # Response variable (trait)
  Column2 = c(1, 1, 2, 2)       # Covariate (e.g., biomass)
)

calculate_Pi_with_adjusted_means(df_test4,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)
calculate_Pi_with_adjusted_means(df_test5,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)
calculate_Pi_with_adjusted_means(df_test_simple,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)

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



# Example usage with synthetic data
external_light = rep(c(0.4, 0.6, 0.8), each = 100)
external_water = sample(rep(c("Low", "High"), each = 150))

combined_factors= combine_factors(synthetic_data1, factors_not_in_dataframe = list(external_light,external_water), factors=1)

##########################################

#' Calculate the Phenotypic Plasticity Index (PPF) Based on Least Square Means
#'
#' This function calculates the Phenotypic Plasticity Index (PPF), which quantifies the phenotypic plasticity 
#' of a trait between two environments or environmental groupings. The PPF is calculated using the least square means (LSMs) 
#' of the trait in each environment, with optional adjustment for covariates that may influence the trait.
#'
#' Covariates can be included to adjust the model for additional environmental influences or biological factors that might affect the trait.
#'
#' @param data A data frame containing the environmental indicators and trait data.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number, name, or a vector representing the environmental conditions.
#' @param covariates (Optional) A vector of column names or indices to include as covariates in the model.
#' @return The Phenotypic Plasticity Index (PPF) value.
#' @examples
#' df = data.frame(
#'   Environment = rep(c("Env1", "Env2"), each = 10),
#'   Height = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 20, 22, 21, 23, 24, 25, 23, 24, 22, 23),
#'   SoilQuality = c(3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3)
#' )
#' 
#' # Calculate PPF with a covariate
#' PPF = calculate_PPF(df, trait_col = "Height", env_col = "Environment", covariates = "SoilQuality")
#' print(PPF)
#' @export
calculate_PPF = function(data, trait_col, env_col, covariates = NULL) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  env_col=as.factor(env_col)
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Fit the linear model
  if (is.null(covariates)) {
    model = lm(trait_col ~ env_col, data = data)
  } else {
    covariates = if (is.numeric(covariates)) names(data)[covariates] else covariates
    formula = as.formula(paste("trait_col ~ env_col +", paste(covariates, collapse = " + ")))
    model = lm(formula, data = data)
  }
  
  # Calculate least square means (LSMs) for each environment
  lsm = as.data.frame(emmeans::emmeans(model, ~ env_col))
  
  # Ensure there are exactly two environments
  if (nrow(lsm) != 2) {
    stop("PPF calculation requires exactly two environments.")
  }
  
  # Calculate PPF using the formula: 100 x ((LSM in one environment - LSM in the other) / LSM in the first environment)
  PPF = 100 * abs((lsm[1, "emmean"] - lsm[2, "emmean"]) / lsm[1, "emmean"])
  
  return(PPF)
}

synthetic_data2=combine_factors(synthetic_data1,factors=NULL, factors_not_in_dataframe=list(external_water))

## test - passed on synthetic dataset

df_test6 = data.frame(
  Column0 = c(rep(3, 15), rep(2, 15)),   # Response variable
  Column1 = c(rep(4, 15), rep(2, 15)),   # Control (2) and Treatment (3)
  Column2 = c(rep(3, 10), rep(2, 10), rep(1, 10))    # Covariate (matches values of Column0)
)
df_test6$Column0 = as.factor(df_test6$Column0)
calculate_PPF(df_test6,env_col = 1, trait_col = 2,covariates = NULL)
model = lm(Column1 ~ Column0 , data = df_test6)
lsmeans_env = emmeans(model, ~ Column0)
summary_lsmeans = summary(lsmeans_env)
# Extract only the LSMeans (adjusted means)
lsmeans_values = summary_lsmeans$emmean
100*abs((lsmeans_values[[1]]-lsmeans_values[[2]])/lsmeans_values[[1]])



###########################################


#' Calculate the Phenotypic Plasticity Index (Pi)
#'
#' This function calculates the Phenotypic Plasticity Index (Pi), defined as the 
#' difference between the maximum and minimum values of a trait divided by the maximum value.
#' For the calculation of the Pi of all traits, the function needs to be called repeatedly.
#'
#' @param data A numeric vector, data frame, or matrix containing the trait data. 
#' If a data frame or matrix is provided, each column represents a trait.
#' @param trait_col (Optional) The column number or name of the trait to analyze if `data` is a data frame or matrix.
#' @return The Phenotypic Plasticity Index (Pi) value(s).

#' @export
calculate_Phenotypic_Plasticity_Index = function(data, trait_col = NULL) {
  if (is.null(trait_col)) {
    # If trait_col is not provided, assume data is a numeric vector
    max_value = max(data, na.rm = TRUE)
    min_value = min(data, na.rm = TRUE)
  } else {
    # If trait_col is provided, assume data is a data frame or matrix
    max_value = max(data[[trait_col]], na.rm = TRUE)
    min_value = min(data[[trait_col]], na.rm = TRUE)
  }
  
  # Calculate Pi as (Max - Min) / Max
  Pi = (max_value - min_value) / max_value
  
  return(Pi)
}


##test - passed on a synthetic dataset (look at the dataset for confirmation)

calculate_Phenotypic_Plasticity_Index(df_test6,trait_col = 2)




####################################

#' Calculate the Proportional Inter-Median Difference (PImd)
#'
#' This function calculates the Proportional Inter-Median Difference (PImd), 
#' defined as the difference between the maximum and minimum medians of a trait 
#' across different environmental conditions divided by the maximum median. 
#' It provides a measure of the relative variability in the trait across the specified conditions.
#'
#' @param data A data frame containing the trait and environmental condition data.
#' @param trait_col The column number or name of the trait to analyze. It can be a column index (integer) 
#' or a column name (string) within the data frame.
#' @param env_col The column number, name, or a vector representing the environmental conditions. 
#' It can be a column index (integer), a column name (string), or a vector of values.
#' @return The Proportional Inter-Median Difference (PImd) value.
#' @examples
#' # Example usage with a data frame
#' data = data.frame(
#'   trait_col = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'   env_col = factor(c("A", "A", "B", "B", "C", "C", "A", "B", "C", "A"))
#' )
#' result = calculate_PImd(data, "trait_col", "env_col")
#' print(result)
#' @export
calculate_PImd = function(data, trait_col, env_col) {
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    # If env_col is a vector, keep it as is
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Handle trait_col
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  
  # Ensure env_col is a factor to use levels
  if (!is.factor(env_col)) {
    env_col = factor(env_col)
  }
  # Get levels of env_col
  levels = levels(env_col)
  medians=c()
  for (level in levels){
    medians=c(medians,median(trait_col[env_col==level]))
  }
  PImd=(max(medians)-min(medians))/max(medians)
  return(PImd)
}

##test - passed on synthetic dataset (check  dataset for confirmation)

calculate_PImd(df_test6,trait_col = 2,env_col = 1)

###############################################



#' Calculate the Proportional Inter-Least Square Mean Difference (PILSM)
#' and Plot LSMs Against the Trait Data
#'
#' This function calculates the Proportional Inter-Least Square Mean Difference (PILSM), 
#' defined as the difference between the maximum and minimum least square means (LSMs) 
#' of a trait across different environments, divided by the maximum LSM.
#' Additionally, it plots the fitted LSMs against the original trait data.
#'
#' The function assumes normality in the data and can optionally adjust for covariates.
#'
#' @param data A data frame containing the trait data and environmental indicators.
#' @param trait_col The column number or name of the trait to analyze.
#' @param env_col The column number, name, or a vector representing the environmental conditions.
#' @param covariates (Optional) A vector of column names or indices to include as covariates in the model.
#' @return A list with the PILSM value and the LSM plot.
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
#' result = calculate_PILSM(df, trait_col = "Height", env_col = "Environment", covariates = "SoilQuality")
#' print(result$PILSM)
#' print(result$plot)
#' @export
calculate_PILSM = function(data, trait_col, env_col, covariates = NULL) {
  # List of required packages
  required_packages = c("emmeans", "ggplot2","ggplot")
  
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
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  trait_col = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]
  env_col=as.factor(env_col)
  # Fit the linear model
  if (is.null(covariates)) {
    model = lm(trait_col ~ env_col, data = data)
  } else {
    covariates = if (is.numeric(covariates)) names(data)[covariates] else covariates
    formula = as.formula(paste("trait_col ~ env_col +", paste(covariates, collapse = " + ")))
    model = lm(formula, data = data)
  }
  
  # Calculate least square means (LSMs) for each environment using emmeans
  lsm = as.data.frame(emmeans::emmeans(model, ~ env_col))
  print(lsm)
  # Calculate PILSM
  max_lsm = max(lsm$emmean, na.rm = TRUE)
  min_lsm = min(lsm$emmean, na.rm = TRUE)
  print(max_lsm)
  print(min_lsm)
  PILSM = (max_lsm - min_lsm) / max_lsm
  
  # Plot the LSMs against the original data
  plot_data = data.frame(Environment = env_col, Trait = trait_col)
  lsm_plot = ggplot2::ggplot(plot_data, ggplot2::aes(x = Environment, y = Trait)) +
    ggplot2::geom_point(color = "blue", alpha = 0.5) +
    ggplot2::geom_point(data = lsm, ggplot2::aes(x = env_col, y = emmean), color = "red", size = 3) +
    ggplot2::geom_line(data = lsm, ggplot2::aes(x = env_col, y = emmean), color = "red", linetype = "dashed") +
    ggplot2::labs(title = "Least Square Means (LSMs) by Environment", x = "Environment", y = "Trait") +
    ggplot2::theme_minimal()
  
  
  return(list(PILSM = PILSM, plot = lsm_plot))
}


## test - passed on synthetic dataset 


calculate_PILSM(df_test6,trait_col=2,env_col=1)

model = lm(Column1 ~ Column0 , data = df_test6)
lsmeans_env = emmeans(model, ~ Column0)
summary_lsmeans = summary(lsmeans_env)
# Extract only the LSMeans (adjusted means)
lsmeans_values = summary_lsmeans$emmean
(max(lsmeans_values)-min(lsmeans_values))/max(lsmeans_values)

################################################


#' Calculate the Relative Trait Response (RTR)
#'
#' This function calculates the Relative Trait Response (RTR) score, defined as the 
#' difference between the mean trait value at one end of an environmental gradient 
#' and the mean trait value at the opposite end, divided by the absolute maximum value of the trait.
#'
#' @param data A data frame containing the trait data and environmental conditions.
#' @param trait_col The column name or number for the trait to analyze.
#' @param env_col The column name or number for the environmental conditions.
#' @param env_low The value of the environmental condition representing one end of the gradient.
#' @param env_high The value of the environmental condition representing the opposite end of the gradient.
#' @return The RTR value.
#' @examples
#' df = data.frame(
#'   Environment = rep(c("Low", "High"), each = 10),
#'   Height = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'             20, 22, 21, 23, 24, 25, 23, 24, 22, 23)
#' )
#' 
#' RTR = calculate_RTR(df, trait_col = "Height", env_col = "Environment", env_low = "Low", env_high = "High")
#' print(RTR)
#' @export
calculate_RTR = function(data, trait_col, env_col, env_low, env_high) {
  
  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }
  
  # Subset data for the specified environmental conditions
  data_low = data[env_col == env_low, ]
  data_high = data[env_col == env_high, ]
  # Calculate mean trait values for each end of the gradient
  mean_low = mean(data_low[[trait_col]], na.rm = TRUE)
  mean_high = mean(data_high[[trait_col]], na.rm = TRUE)
  print(mean_high)
  print(mean_low)
  # Calculate the RTR value
  RTR = (mean_high - mean_low) / max(abs(trait_col), na.rm = TRUE)
  
  return(RTR)
}



## test - passed on synthetic dataset
calculate_RTR(df_test2,trait_col=2,env_col=1,env_low=1,env_high=2)


######################################

#' @title Calculate Phenotypic Instability Ratio (PIR)
#'
#' @description This function calculates the Phenotypic Instability Ratio (PIR) for a given trait across different environments,
#' following the method described by Robinson (1989). PIR is calculated as the ratio of the difference between the maximum and minimum
#' mean trait values across environments to the mean trait value in the environment with the maximum relative growth rate (RGR).
#' The PIR metric provides an intermediate measure of plasticity and assumes normality in the data. It also requires prior knowledge of the
#' relative growth rates (RGRs) for each environment. Note that the method has statistical limitations.
#'
#' @param data A data frame containing the data for trait values and corresponding environment and RGR information.
#' @param trait_col A column name or numeric vector indicating the trait values. Can be a numeric index of the column in the data frame or a vector.
#' @param env_col A column name or numeric vector indicating the environment identifiers. Can be a numeric index of the column in the data frame or a vector.
#' @param rgr_col A column name or numeric vector indicating the relative growth rate (RGR) values. Can be a numeric index of the column in the data frame or a vector.
#'
#' @details The function first converts the environment column into a factor and computes the mean trait values for each environment.
#' It then identifies the environment with the maximum relative growth rate and uses the trait's mean value in that environment
#' to compute the PIR score. The calculation follows the formula:
#' \deqn{PIR = \frac{MaxMean - MinMean}{Mean_{MaxRGR}}}
#' where \code{MaxMean} is the maximum mean trait value, \code{MinMean} is the minimum mean trait value,
#' and \code{Mean_{MaxRGR}} is the mean trait value in the environment with the maximum RGR.
#'
#' @return A numeric value representing the Phenotypic Instability Ratio (PIR).
#'
#' @references
#' Robinson, D. (1989). \emph{Plasticity in plant growth and resource use as a trait for crop breeding}. Field Crops Research, 11(2-3), 153-159.
#'
#' @examples
#' \dontrun{
#' # Sample usage with a data frame 'df' with columns 'Trait', 'Env', and 'RGR':
#' calculate_PIR(data = df, trait_col = "Trait", env_col = "Env", rgr_col = "RGR")
#' }
#' 
#' @export
calculate_PIR = function(data, trait_col, env_col, rgr_col) {
  
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
  
  # Handle trait_col
  trait_col = if (is.numeric(trait_col) && length(trait_col) == 1) data[[trait_col]] else trait_col
  
  # Calculate mean trait values for each environment
  means = tapply(trait_col, env_col, mean, na.rm = TRUE)

  # Identify maximum and minimum means
  max_mean = max(means, na.rm = TRUE)
  min_mean = min(means, na.rm = TRUE)

  
  # Identify the environment with the maximum growth rate
  max_rgr_env = levels(env_col)[which.max(tapply(rgr_col, env_col, mean, na.rm = TRUE))]
  
  # Find the mean of the trait at the environment where the growth rate is maximum
  mean_at_max_rgr = max(means[max_rgr_env])
  
  
  # Calculate PIR
  PIR = (max_mean - min_mean) / mean_at_max_rgr
  
  return(PIR)
}

## test - passed on synthetic dataset
specific_growthrate=c(rep(10,10),rep(20,10))

calculate_PIR(df_test2 , trait_col = 2 , env_col = 1, rgr_col = specific_growthrate)




