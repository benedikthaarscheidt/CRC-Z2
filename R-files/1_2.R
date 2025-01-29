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
library(purrr)
library(dplyr)

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
#' This function calculates the Coefficient of Variation (CVt) for a single vector of trait values.
#' The user is responsible for pre-processing the data to generate the input vector.
#'
#' @param trait_values A numeric vector containing the trait values for a single genotype.
#' @return A numeric value representing the CVt, or NA if the input vector length is less than 2.
#'
#' @examples
#' # Example usage
#' trait_values = c(100.2, 92, 83.5)
#' CVt = calculate_CVt_single(trait_values)
#' print(CVt)
#'
#' @export
calculate_CVt = function(trait_values) {
  if (length(trait_values) < 2) {
    return(NA)  # Avoid division by zero for single values
  }
  
  if (any(is.na(trait_values))) {
    stop("The input vector contains missing values. Please handle missing data before calling this function.")
  }
  
  # Calculate and return the CVt
  return(sd(trait_values) / mean(trait_values))
}




### test - passed in by hand calculation with this dataset

df_test1 = data.frame(genotype=c(rep(1,10),rep(1,10)),Column1 = c(rep(4, 10), rep(2, 10)), Column2 = c(rep(10, 10), rep(1, 10)))
#
##traitwise
#cv=calculate_CVt(df_test1,genotype_col = 1, trait_col = 2)
#print(sd(df_test1[,2])/mean(df_test1[,2]))
##print(sd(df_test1[,2])/mean(df_test1[,2]))
#print(cv)
##total
#df_test1_flattened=as.numeric(unlist(df_test1))
#calculate_CVt(df_test1,traitwise = F)
#print(sd(df_test1_flattened)/mean(df_test1_flattened))



################################


calculate_reaction_norm_slope = function(trait_values, environments=NULL) {
  if (length(trait_values) < 2) {
    return(NA)  
  }
  
  if(is.null(environments)){
    environments = seq_along(trait_values)
  }

  lm_fit = lm(trait_values ~ environments)
  
  return(coef(lm_fit)["environments"])
}


##############################

#' Calculate Nonlinear Reaction Norm for a Single Genotype
#'
#' This function fits a polynomial model to calculate the nonlinear reaction norm
#' for a single genotype, given a vector of trait values. The environmental values
#' are assumed to be equidistant.
#'
#' @param trait_values A numeric vector of trait values for a single genotype.
#' @param degree An integer specifying the degree of the polynomial model to fit (default is 2).
#' @return A list containing:
#'   - `Coefficients`: A numeric vector of coefficients for the fitted polynomial model.
#'   - `P_Value`: The p-value for the overall model fit (F-test from ANOVA).
#'   - `Linear_Slope`: The linear term (coefficient of degree 1) if present.
#'   - `Quadratic_Term`: The quadratic term (coefficient of degree 2) if present.
#'   - `Valid`: Logical indicating whether the model was successfully fitted.
#'
#' @examples
#' trait_values = c(10, 12, 15, 20, 25)
#' result = calculate_reaction_norm_nonlinear_single(trait_values, degree = 2)
#' print(result)
#'
#' @export
calculate_reaction_norm_non_linear = function(trait_values, degree = 2) {
  environments = seq_along(trait_values)
  
  # Ensure enough data points for the model
  if (length(trait_values) < degree + 1) {
    return(data.frame(
      Intercept = NA,
      Coefficients = rep(NA, degree),
      P_Value = NA,
      Valid = FALSE
    ))
  }
  
  # Try fitting the polynomial model
  tryCatch({
    model = lm(trait_values ~ poly(environments, degree))
    
    # Extract coefficients and p-value
    coefficients = coef(model)
    p_value = anova(model)$`Pr(>F)`[1]
    
    # Create a structured output
    output = data.frame(
      Intercept = coefficients[1],
      P_Value = p_value,
      Valid = TRUE
    )
    
    # Add polynomial coefficients as separate columns
    for (i in 2:length(coefficients)) {
      output[[paste0("Coefficient_", i - 1)]] = coefficients[i]
    }
    
    return(output)
  }, error = function(e) {
    # Return a structured NA output if an error occurs
    output = data.frame(
      Intercept = NA,
      P_Value = NA,
      Valid = FALSE
    )
    
    for (i in 1:degree) {
      output[[paste0("Coefficient_", i)]] = NA
    }
    
    return(output)
  })
}



#test - passed with synthetic dataset

#df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 10), rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
#calculate_reaction_norm_slope(df_test2,env_col = 1,trait_cols = c(2,3),plot = T)

################################


#NOTE: if the resource availability is an actual measurement of a metabolite which is being used by the plant then the grouping of the plants into high vs low resource availability should be done by clustering.
# Sadly I am missing the insight into the common practices in the field.


#' Calculate the D Slope from Trait Values
#'
#' This function calculates the D slope, which quantifies the scope of plasticity 
#' by computing the difference between the mean of the highest and lowest fractions 
#' of a sorted trait value vector. The lower and upper fractions are defined by the 
#' user.
#'
#' @param trait_values A numeric vector of trait values.
#' @param lower_fraction A numeric value (0 < lower_fraction < 1) specifying the fraction of the data 
#' to consider as the "lower" group. Defaults to 0.2 (the first 20% of sorted values).
#' @param upper_fraction A numeric value (0 < upper_fraction <= 1) specifying the fraction of the data 
#' to consider as the "upper" group. Defaults to 0.9 (the last 10% of sorted values).
#' @return A numeric value representing the D slope (difference between the means of 
#' the upper and lower groups).
#' @examples
#' # Example trait values
#' trait_values = c(10, 12, 20, 22, 25, 28, 30, 32, 35, 40)
#'
#' # Calculate D slope with 20% lower and 90% upper fractions
#' calculate_D_slope(trait_values, lower_fraction = 0.2, upper_fraction = 0.9)
#' 
#' # Calculate D slope with 10% lower and 80% upper fractions
#' calculate_D_slope(trait_values, lower_fraction = 0.1, upper_fraction = 0.8)
#' @export
calculate_D_slope = function(trait_values, lower_fraction = 0.2, upper_fraction = 0.8) {
  # Ensure trait values are sorted and valid
  sorted_values = sort(trait_values, na.last = TRUE)
  
  # Calculate boundaries
  lower_boundary_index = floor(length(sorted_values) * lower_fraction)
  upper_boundary_index = ceiling(length(sorted_values) * upper_fraction)
  
  # Extract lower and upper segments
  lower_values = sorted_values[1:lower_boundary_index]
  upper_values = sorted_values[upper_boundary_index:length(sorted_values)]
  
  # Calculate D slope
  D_slope = mean(upper_values, na.rm = TRUE) - mean(lower_values, na.rm = TRUE)
  
  return(D_slope)
}

#test - passed with synthetic 

#calculate_D_slope(df_test2,env_col = 1,trait_cols = 2)
#df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
#calculate_D_slope(df_test3,env_col = 1,trait_cols = 3)
#mean(df_test3[1:15,3])-mean(df_test3[16:nrow(df_test3),3])




################################

#' Calculate the Response Coefficient (RC) for Trait Values
#'
#' This function calculates the Response Coefficient (RC), which is the ratio of the mean trait values
#' between the highest and lowest fractions of a sorted trait value vector.
#'
#' @param trait_values A numeric vector of trait values.
#' @param lower_fraction A numeric value (0 < lower_fraction < 1) specifying the fraction of the data 
#' to consider as the "low" group. Defaults to 0.5 (the first 50% of sorted values).
#' @param upper_fraction A numeric value (0 < upper_fraction <= 1) specifying the fraction of the data 
#' to consider as the "high" group. Defaults to 0.5 (the last 50% of sorted values).
#' @return A numeric value representing the Response Coefficient (RC), calculated as the ratio of the 
#' mean of the "high" group to the mean of the "low" group.
#' @examples
#' # Example trait values
#' trait_values = c(10, 12, 20, 22, 25, 28, 30, 32, 35, 40)
#'
#' # Calculate RC with 50% lower and 50% upper fractions
#' calculate_RC(trait_values, lower_fraction = 0.5, upper_fraction = 0.5)
#'
#' # Calculate RC with 20% lower and 20% upper fractions
#' calculate_RC(trait_values, lower_fraction = 0.2, upper_fraction = 0.8)
#' @export
calculate_RC = function(trait_values, lower_fraction = 0.5, upper_fraction = 0.5) {
  # Ensure trait values are sorted
  sorted_values = sort(trait_values, na.last = TRUE)
  
  # Calculate indices for splitting
  lower_boundary_index = floor(length(sorted_values) * lower_fraction)
  upper_boundary_index = ceiling(length(sorted_values) * (1 - upper_fraction))
  
  # Extract lower and upper segments
  lower_values = sorted_values[1:lower_boundary_index]
  upper_values = sorted_values[(upper_boundary_index + 1):length(sorted_values)]
  
  # Calculate means
  mean_low = mean(lower_values, na.rm = TRUE)
  mean_high = mean(upper_values, na.rm = TRUE)
  
  # Calculate Response Coefficient
  RC = mean_high / mean_low
  
  return(RC)
}
# test - passed with synthetic dataset

#calculate_RC(df_test3,env_col = 1,trait_cols = 3)
#
#mean(df_test3[1:15,3])/mean(df_test3[16:nrow(df_test3),3])


################################

#' Calculate the Coefficient of Variation of Means (CVm) for Trait Values
#'
#' This function calculates the Coefficient of Variation of Means (CVm),
#' which is the ratio of the standard deviation of group means to the mean of the group means.
#' It can accept either a list of vectors (one vector per group) or a single vector of trait values 
#' along with a vector of group labels.
#'
#' @param trait_values Either a numeric vector of trait values or a list of numeric vectors (one per group).
#' @param group_labels An optional vector of group labels. Required if `trait_values` is a single vector.
#' @return A numeric value representing the CVm for the given trait values.
#' @examples
#' # Example 1: Using a list of vectors
#' trait_groups = list(
#'   Env1 = c(10, 12),
#'   Env2 = c(20, 22),
#'   Env3 = c(15, 18)
#' )
#' calculate_CVm(trait_values = trait_groups)
#'
#' # Example 2: Using a single vector with group labels
#' trait_values = c(10, 12, 20, 22, 15, 18)
#' group_labels = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3")
#' calculate_CVm(trait_values = trait_values, group_labels = group_labels)
#' @export
calculate_CVm = function(trait_values, group_labels = NULL) {
  # Handle case where trait_values is a list of vectors
  if (is.list(trait_values)) {
    # Calculate group means directly
    group_means = sapply(trait_values, mean, na.rm = TRUE)
  } else if (!is.null(group_labels)) {
    # Calculate group means based on group labels
    group_means = tapply(trait_values, group_labels, mean, na.rm = TRUE)
  } else {
    stop("If trait_values is not a list, group_labels must be provided.")
  }
  
  # Calculate standard deviation and mean of group means
  sd_of_means = sd(group_means, na.rm = TRUE)
  mean_of_means = mean(group_means, na.rm = TRUE)
  
  # Calculate and return CVm
  CVm = sd_of_means / mean_of_means
  return(CVm)
}

## test - passed with synthetic dataset

#calculate_CVm(df_test3,env_col = 1,trait_col = 3)
#
#sd(c(mean(df_test3[1:10,3]),mean(df_test3[11:20,3]),mean(df_test3[21:30,3])))/mean(mean(df_test3[,3]))


###############################


#' Calculate the Coefficient of Variation of Medians (CVmd) for Trait Values
#'
#' This function calculates the Coefficient of Variation of Medians (CVmd),
#' which is the ratio of the standard deviation of group medians to the mean of the group medians.
#' It can accept either a list of vectors (one vector per group) or a single vector of trait values
#' along with a vector of group labels.
#'
#' @param trait_values Either a numeric vector of trait values or a list of numeric vectors (one per group).
#' @param group_labels An optional vector of group labels. Required if `trait_values` is a single vector.
#' @return A numeric value representing the CVmd for the given trait values.
#' @examples
#' # Example 1: Using a list of vectors
#' trait_groups = list(
#'   Env1 = c(10, 12),
#'   Env2 = c(20, 22),
#'   Env3 = c(15, 18)
#' )
#' calculate_CVmd(trait_values = trait_groups)
#'
#' # Example 2: Using a single vector with group labels
#' trait_values = c(10, 12, 20, 22, 15, 18)
#' group_labels = c("Env1", "Env1", "Env2", "Env2", "Env3", "Env3")
#' calculate_CVmd(trait_values = trait_values, group_labels = group_labels)
#' @export
calculate_CVmd = function(trait_values, group_labels = NULL) {
  # Handle case where trait_values is a list of vectors
  if (is.list(trait_values)) {
    # Calculate group medians directly
    group_medians = sapply(trait_values, median, na.rm = TRUE)
  } else if (!is.null(group_labels)) {
    # Calculate group medians based on group labels
    group_medians = tapply(trait_values, group_labels, median, na.rm = TRUE)
  } else {
    stop("If trait_values is not a list, group_labels must be provided.")
  }
  
  # Calculate standard deviation and mean of group medians
  sd_of_medians = sd(group_medians, na.rm = TRUE)
  mean_of_medians = mean(group_medians, na.rm = TRUE)
  
  # Calculate and return CVmd
  CVmd = sd_of_medians / mean_of_medians
  return(CVmd)
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
#' @param trait_data A named list of numeric vectors representing the trait values, one vector per trait.
#' @param env_data A factor or vector representing the environmental conditions (must include control).
#' @param covariate_data A numeric vector representing the covariate (e.g., biomass) values.
#' @param control_env The label or value in `env_data` that represents the control environment.
#' @return A named numeric vector where each element represents the grand plasticity (Pi) for a single trait.
#' @examples
#' # Example data
#' trait_data = list(
#'   Height = c(10, 12, 20, 22),
#'   Weight = c(5, 6, 7, 8)
#' )
#' env_data = factor(c("Control", "Control", "Treatment1", "Treatment2"))
#' covariate_data = c(5, 6, 7, 8)  # Covariate (e.g., biomass)
#' control_env = "Control"
#'
#' grand_plasticity = calculate_grand_plasticity(
#'   trait_data = trait_data,
#'   env_data = env_data,
#'   covariate_data = covariate_data,
#'   control_env = control_env
#' )
#' print(grand_plasticity)
#' @export
calculate_grand_plasticity = function(trait_values, env_data, covariate_data = NULL, control_env = NULL) {

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("The 'emmeans' package is required but not installed.")
  }

  env_data = factor(env_data)

  if (!is.null(covariate_data)) {
    model = lm(trait_values ~ covariate_data + env_data)
  } else {
    model = lm(trait_values ~ env_data)
  }
  adjusted_means = emmeans::emmeans(model, ~ env_data)
  adjusted_means_summary = summary(adjusted_means)
  
  if (!is.null(control_env)) {
    treatment_means = adjusted_means_summary$emmean[adjusted_means_summary$env_data != control_env]
  } else {
    treatment_means = adjusted_means_summary$emmean
  }
  
  if (length(treatment_means) == 0) {
    stop("Could not extract treatment means. Check your input data and control environment.")
  }
  
  sd_means = sd(treatment_means, na.rm = TRUE)
  grand_mean = mean(treatment_means, na.rm = TRUE)
  
  grand_plasticity = sd_means / grand_mean
  
  return(grand_plasticity)
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

#calculate_grand_plasticity(trait_values=df_test4[,2], env_data=df_test4[,1], covariate_data = NULL, control_env = NULL) this should produce 0.5
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

#' @title Calculate the Phenotypic Plasticity Factor (PPF) for One or Multiple Traits
#'
#' @description This function calculates the Phenotypic Plasticity Factor (PPF) for one or multiple traits across environments.
#' PPF is defined as the proportional difference between least square means (LSMs) of a trait across two environments.
#' Covariates can be included to adjust for additional factors. The function dynamically handles single or multiple covariates.
#'
#' @param trait_values A numeric vector containing the trait values.
#' @param env_groups A vector (numeric, character, or factor) representing the environmental groups corresponding to `trait_values`.
#' @param covariate_values (Optional) A numeric vector (for one covariate) or a data frame (for multiple covariates) corresponding to `trait_values`.
#' @param env_pairs (Optional) A list of environment pairs to calculate the PPF for. Each pair should be a vector of two elements.
#' If `NULL`, the function calculates the PPF for all possible combinations of environments.
#'
#' @return A named numeric vector where each element represents the PPF value for a specific environment pair.
#'
#' @examples
#' trait_values = c(10, 12, 14, 16, 11, 13, 15, 17, 20, 22)
#' env_groups = rep(c("Env1", "Env2"), each = 5)
#' covariate_values = data.frame(Covariate1 = c(5, 5.1, 5.2, 5.3, 5.4, 6, 6.1, 6.2, 6.3, 6.4))
#' 
#' # Calculate PPF for all pairs
#' PPF_all = calculate_PPF(trait_values, env_groups, covariate_values)
#' print(PPF_all)
#'
#' # Calculate PPF for specific pairs
#' PPF_specific = calculate_PPF(trait_values, env_groups, covariate_values, env_pairs = list(c("Env1", "Env2")))
#' print(PPF_specific)
#' @export
calculate_PPF = function(trait_values, env_groups, covariate_values = NULL, env_pairs = NULL) {
  # Ensure the environment groups are treated as factors
  env_groups = as.factor(env_groups)
  
  # Default: Calculate PPF for all possible combinations of environment pairs
  if (is.null(env_pairs)) {
    env_pairs = combn(levels(env_groups), 2, simplify = FALSE)
  } else if (is.vector(env_pairs) && length(env_pairs) == 2) {
    env_pairs = list(env_pairs)  # Convert a single pair into a list of one pair
  }
  
  # Initialize a vector to store PPF values for each environment pair
  PPF_values = numeric(length(env_pairs))
  names(PPF_values) = sapply(env_pairs, function(pair) paste(pair, collapse = "-"))
  
  # Loop through each environment pair
  for (i in seq_along(env_pairs)) {
    env_pair = env_pairs[[i]]
    
    # Subset the data to include only the current environment pair
    subset_idx = which(env_groups %in% env_pair)
    subset_trait_values = trait_values[subset_idx]
    subset_env_groups = env_groups[subset_idx]
    
    if (!is.null(covariate_values)) {
      # Handle single or multiple covariates
      if (is.vector(covariate_values)) {
        subset_covariates = covariate_values[subset_idx]
        data = data.frame(trait_values = subset_trait_values, env_groups = subset_env_groups, Covariate = subset_covariates)
      } else if (is.data.frame(covariate_values)) {
        subset_covariates = covariate_values[subset_idx, , drop = FALSE]
        data = data.frame(trait_values = subset_trait_values, env_groups = subset_env_groups, subset_covariates)
      }
      formula = as.formula("trait_values ~ env_groups + .")
      model = lm(formula, data = data)
    } else {
      model = lm(subset_trait_values ~ subset_env_groups)
    }
    
    # Calculate least square means (LSMs) for the current environment pair
    lsm = as.data.frame(emmeans::emmeans(model, ~ env_groups))
    
    # Calculate PPF for the current pair: 100 * (|LSM1 - LSM2| / LSM1)
    PPF_values[i] = 100 * abs((lsm[1, "emmean"] - lsm[2, "emmean"]) / lsm[1, "emmean"])
  }
  
  return(PPF_values)
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
#' The function accepts vectors for traits and environmental conditions as input.
#'
#' @param trait_values A numeric vector containing the trait values.
#' @param env A vector (numeric, character, or factor) representing the environmental conditions.
#' @return A numeric value representing the PImd for the trait.
#' @examples
#' # Example usage
#' trait_values = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' env = factor(c("A", "A", "B", "B", "C", "C", "A", "B", "C", "A"))
#' result = calculate_PImd(trait_values, env)
#' print(result)
#' @export
calculate_PImd = function(trait_values, env) {
  # Ensure env is treated as a factor
  if (!is.factor(env)) {
    env = factor(env)
  }
  
  # Get the levels of the environment
  levels = levels(env)
  medians = numeric(length(levels))
  
  # Calculate medians for each environment level
  for (i in seq_along(levels)) {
    level = levels[i]
    medians[i] = median(trait_values[env == level], na.rm = TRUE)
  }
  
  # Check if maximum median is zero to avoid division by zero
  if (max(medians) == 0) {
    warning("Maximum median is zero; PImd cannot be calculated.")
    return(NA)
  }
  
  # Calculate PImd as (Max - Min) / Max
  PImd = (max(medians) - min(medians)) / max(medians)
  
  return(PImd)
}


##test - passed on synthetic dataset (check  dataset for confirmation)

#calculate_PImd(df_test6,trait_cols = c(2,3),env_col = 1)

###############################################



#' Calculate the Proportional Inter-Least Square Mean Difference (PILSM)
#' for One or Multiple Traits
#'
#' This function calculates the Proportional Inter-Least Square Mean Difference (PILSM), 
#' defined as the difference between the maximum and minimum least square means (LSMs)  
#' of a trait across different environments, divided by the maximum LSM.
#'
#' @param trait_values A numeric vector containing the trait values.
#' @param env A vector (numeric, character, or factor) representing the environmental conditions.
#' @param covariates (Optional) A numeric vector or data frame for covariate values.
#' @return A numeric value representing the PILSM for the trait.
#' @examples
#' trait_values = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'                  20, 22, 21, 23, 24, 25, 23, 24, 22, 23, 
#'                  30, 32, 31, 33, 34, 35, 33, 34, 32, 33)
#' env = factor(rep(c("Env1", "Env2", "Env3"), each = 10))
#' covariates = data.frame(SoilQuality = c(rep(3, 10), rep(4, 10), rep(5, 10)))
#' result = calculate_PILSM(trait_values, env, covariates)
#' print(result)
#' @export
calculate_PILSM = function(trait_values, env, covariates = NULL) {
  # Ensure env is treated as a factor
  if (!is.factor(env)) {
    env = factor(env)
  }
  
  # Fit the linear model
  if (is.null(covariates)) {
    model = lm(trait_values ~ env)
  } else {
    # If covariates are a vector, convert to a data frame
    if (is.vector(covariates)) {
      covariates = data.frame(Covariate = covariates)
    }
    # Combine env and covariates into a single data frame
    data = data.frame(trait_values = trait_values, env = env, covariates)
    formula = as.formula("trait_values ~ env + .")
    model = lm(formula, data = data)
  }
  
  # Calculate least square means (LSMs) for each environment
  lsm = as.data.frame(emmeans::emmeans(model, ~ env))
  
  # Calculate PILSM
  max_lsm = max(lsm$emmean, na.rm = TRUE)
  min_lsm = min(lsm$emmean, na.rm = TRUE)
  PILSM = (max_lsm - min_lsm) / max_lsm
  
  return(PILSM)
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
#' @param trait_values A numeric vector containing the trait measurements.
#' @param env_values A numeric vector representing the environmental gradient corresponding to `trait_values`.
#' @param env_low Either a specific value or a fraction (e.g., 0.2 for the bottom 20%) representing the low end of the environmental gradient.
#' @param env_high Either a specific value or a fraction (e.g., 0.2 for the top 20%) representing the high end of the environmental gradient.
#' @return A numeric value representing the RTR score for the trait.
#' @examples
#' trait_values = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'                  20, 22, 21, 23, 24, 25, 23, 24, 22, 23)
#' env_values = seq(1, 10, length.out = 20)
#' RTR_value = calculate_RTR(trait_values, env_values, env_low = 0.3, env_high = 0.3)
#' print(RTR_value)
#' @export
calculate_RTR = function(trait_values, env_values, env_low = 0.2, env_high = 0.2) {
  # Validate input lengths
  if (length(trait_values) != length(env_values)) {
    stop("`trait_values` and `env_values` must have the same length.")
  }
  
  
  if (env_low < 1) {
    low_threshold = quantile(env_values, probs = env_low, na.rm = TRUE)
  } else {
    low_threshold = env_low
  }
  
  if (env_high < 1) {
    high_threshold = quantile(env_values, probs = 1 - env_high, na.rm = TRUE)
  } else {
    high_threshold = env_high
  }
  
  # Subset data based on thresholds
  low_idx = which(env_values <= low_threshold)
  high_idx = which(env_values >= high_threshold)
  
  # Calculate mean trait values for each end of the gradient
  mean_low = mean(trait_values[low_idx], na.rm = TRUE)
  mean_high = mean(trait_values[high_idx], na.rm = TRUE)
  
  # Calculate the RTR value
  RTR_value = (mean_high - mean_low) / max(abs(trait_values), na.rm = TRUE)
  
  return(RTR_value)
}



## test - passed on synthetic dataset
#calculate_RTR(df_test2,trait_col=2,env_col=1,env_low=1,env_high=2)


######################################

#' @title Calculate Phenotypic Instability Ratio (PIR) for One or Multiple Traits
#'
#' @description This function calculates the Phenotypic Instability Ratio (PIR) for a given trait across different environments,
#' following the method described by Robinson (1989). PIR is calculated as the ratio of the difference between the maximum and minimum
#' mean trait values across environments to the mean trait value in the environment with the maximum relative growth rate (RGR).
#' The PIR metric provides an intermediate measure of plasticity and assumes normality in the data. If relative growth rates (RGRs) are
#' not provided, they are dynamically estimated from the mean trait values across environments. If environment values are not provided,
#' they are assumed to be equidistant.
#'
#' @param trait_values A numeric vector containing the trait values.
#' @param env_values (Optional) A vector (numeric, character, or factor) representing the environmental identifiers corresponding to `trait_values`.
#' If `NULL`, environments are assumed to be equidistant.
#' @param rgr_values (Optional) A numeric vector of relative growth rate (RGR) values corresponding to `trait_values`. If `NULL`, RGR is estimated dynamically.
#'
#' @details The function first converts the environment column into a factor (or generates equidistant environments if `env_values` is `NULL`) 
#' and computes the mean trait values for each environment. It then identifies the environment with the maximum relative growth rate and uses the
#' trait's mean value in that environment to compute the PIR score. If RGR is not provided, it is estimated as the relative change in mean trait values
#' between consecutive environments.
#' The calculation follows the formula:
#' \deqn{PIR = \frac{MaxMean - MinMean}{Mean_{MaxRGR}}}
#' where \code{MaxMean} is the maximum mean trait value, \code{MinMean} is the minimum mean trait value,
#' and \code{Mean_{MaxRGR}} is the mean trait value in the environment with the maximum RGR.
#'
#' @return A numeric value representing the PIR for the trait.
#'
#' @references
#' Robinson, D. (1989). \emph{Plasticity in plant growth and resource use as a trait for crop breeding}. Field Crops Research, 11(2-3), 153-159.
#'
#' @examples
#' trait_values = c(10, 12, 11, 13, 14, 15, 13, 14, 12, 13, 
#'                  20, 22, 21, 23, 24, 25, 23, 24, 22, 23, 
#'                  30, 32, 31, 33, 34, 35, 33, 34, 32, 33)
#' env_values = rep(c("Env1", "Env2", "Env3"), each = 10)
#' PIR_value = calculate_PIR(trait_values, env_values)
#' print(PIR_value)
#'
#' # Example with equidistant environments
#' PIR_value_equidistant = calculate_PIR(trait_values)
#' print(PIR_value_equidistant)
#' @export
calculate_PIR = function(trait_values, env_values = NULL, rgr_values = NULL) {
  
  if (is.null(env_values)) {
    env_values = factor(seq_along(trait_values))
  } else {
    env_values = as.factor(env_values)
  }
  

  means = tapply(trait_values, env_values, mean, na.rm = TRUE)
  

  if (is.null(means) || length(means) != length(levels(env_values))) {
    stop("Mismatch between environments and trait values. Ensure proper input alignment.")
  }
  

  max_mean = max(means, na.rm = TRUE)
  min_mean = min(means, na.rm = TRUE)
  

  if (is.null(rgr_values)) {
    rgr_values = numeric(length(means))  
    rgr_values[2:length(means)] = diff(means) / means[-length(means)]
    rgr_values[1] = NA  
  } else {
    
    rgr_values = tapply(rgr_values, env_values, mean, na.rm = TRUE)
  }
  

  if (length(rgr_values) != length(means)) {
    stop("Mismatch between calculated RGR values and environment means.")
  }
  

  max_rgr_index = which.max(rgr_values)
  
  
  mean_at_max_rgr = unname(means[max_rgr_index])  
  
  
  PIR_value = (max_mean - min_mean) / mean_at_max_rgr
  
  return(PIR_value)  
}
