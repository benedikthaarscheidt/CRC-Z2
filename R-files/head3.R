source("~/CRC 1622 - Z2/R-files/1_2.R")
#coefficient-of variation total (calculate_CVt) - tested,
#slope of norm reaction (calculate_reaction_norm_slope) - tested,
#slope of plastic response (D) (calculate_D_slope)- tested,
#response coefficient (RC) (calculate_RC)- tested,
#Standard deviation of means (CVm) (calculate_CVm)- tested,
#Standard deviation of medians (CVmd)(calculate_CVmd)- tested,
#Grand plasticity (calculate_GPi)- tested,
#Phenotypic Plasticity Index (calculate_PPF)- tested,
#Phenotypic Plasticity Index (calculate_Phenotypic_Plasticity_Index)- tested,
#PImd (calculate_PImd)- tested,
#PILSM (calculate_PILSM)- tested,
#RTR (calculate_RTR)- tested,
#PIR (calculate_PIR) - tested
source("~/CRC 1622 - Z2/R-files/2.R")
#RDPI	(rdpi_calculation) - tested,
#RDPIs (rdpi_mean_calculation) - tested,
#ESPI (calculate_ESPI) - tested,
#ESPIid (espiid_calculation) - tested,
#evwpi_calculation (idea from Benedikt)
source("~/CRC 1622 - Z2/R-files/3.R")
#Phenotypic Stability Index (calculate_PSI),
#Relative Plasticity Index (calculate_RPI) - tested,
#Plasticity Quotient (calculate_PQ) - tested,
#Phenotypic Range (PR) (calculate_PR) - tested,
#Norm of reaction width (calculate_NRW) - tested,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
#Calculate Plasticity Differential (PD) (calculate_PD) - tested,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,
source("~/CRC 1622 - Z2/R-files/4.R")
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
source("~/CRC 1622 - Z2/R-files/norms_generator2_nonlinear.R")
source("~/CRC 1622 - Z2/R-files/norms_generator2_linear.R")

n_datasets=3
gaussian_data=gaussian_data
sinusoidal_data=sinusoidal_data
wave_data=wave_data
linear_data=individual_norms

combined=rbind(gaussian_data, wave_data, sinusoidal_data, linear_data)
############################################################

#' Preprocess Dataset for Modular Index Functions
#'
#' This function preprocesses a dataset by averaging trait values over replicates
#' for each genotype and environmental condition. It outputs a list of trait vectors,
#' one for each genotype, with the length corresponding to the number of environmental conditions.
#'
#' @param data A data frame containing the dataset.
#' @param trait_col The column name or index containing the trait values.
#' @param genotype_col The column name or index containing the genotype information.
#' @param replicate_col (Optional) The column name or index containing replicate information.
#'   If provided, the function averages trait values across replicates for each genotype-environment combination.
#' @return A named list where each element corresponds to a genotype and contains:
#'   - A numeric vector of averaged trait values for that genotype.
#'
#' @examples
#' # Example dataset
#' data = data.frame(
#'   Genotype = rep(c("G1", "G2"), each = 10),
#'   Environment = rep(1:5, times = 4),
#'   Replicate = rep(1:2, each = 5, times = 2),
#'   Trait = c(10, 12, 15, 20, 22, 15, 17, 20, 18, 20,
#'             23, 25, 28, 30, 32, 25, 27, 30, 28, 30)
#' )
#'
#' # Preprocess the dataset
#' preprocessed = preprocess_data(
#'   data, trait_col = "Trait", genotype_col = "Genotype",
#'   replicate_col = "Replicate"
#' )
#'
#' print(preprocessed)
#'
#' @export
preprocess_data = function(data, trait_col, genotype_col, replicate_col = NULL) {
  if (is.numeric(trait_col)) trait_col = colnames(data)[trait_col]
  if (is.numeric(genotype_col)) genotype_col = colnames(data)[genotype_col]
  if (!is.null(replicate_col) && is.numeric(replicate_col)) replicate_col = colnames(data)[replicate_col]
  
  genotype_split = split(data, data[[genotype_col]])
  result = list()
  
  for (genotype in names(genotype_split)) {
    genotype_data = genotype_split[[genotype]]
    
    if (is.null(replicate_col)) {
      result[[genotype]] = genotype_data[[trait_col]][seq(1, length(genotype_data[[trait_col]]), 5)]
    } else {
      unique_replicates = unique(genotype_data[[replicate_col]])#unique(genotype_data[[replicate_col]][genotype_data[[replicate_col]] != 0])#
      
      temp = do.call(
        cbind,
        lapply(unique_replicates, function(rep) {
          genotype_data[genotype_data[[replicate_col]] == rep, trait_col]
        })
      )
  
      averaged_vector = rowMeans(temp, na.rm = TRUE)
      result[[genotype]] = averaged_vector[seq(1, length(averaged_vector), 5)] # for more or less samples of the same norm for a different environmental value just change this value here
    }
  }
  
  return(result)
}


call_function = function(data_list, score, additional_1 = NULL, additional_2 = NULL) {
  results = list()  # Store results for each item
  
  for (i in names(data_list)) {
    # Construct arguments dynamically
    args = list(data_list[[i]])
    if (!is.null(additional_1)) args = c(args, list(additional_1))
    if (!is.null(additional_2)) args = c(args, list(additional_2))
    
    # Call the score function
    result = do.call(score, args)
    
    # Store the result in the list
    results[[i]] = result
  }
  
  # Combine results into a data frame or keep as list
  if (all(sapply(results, is.list))) {
    # Extract all unique names from the output lists
    output_names = unique(unlist(lapply(results, names)))
    
    # Create a data frame with one column per output name
    results_df = do.call(rbind, lapply(results, function(res) {
      # Fill missing outputs with NA
      sapply(output_names, function(name) if (name %in% names(res)) res[[name]] else NA)
    }))
    
    rownames(results_df) = names(results)
    return(as.data.frame(results_df))
  } else {
    # Return results as-is if not all outputs are lists
    return(results)
  }
}

post_process = function(data_list, score_function, type_labels, ...) {
  result_list = list()
  
  for (i in seq_along(data_list)) {
    # Calculate scores using the provided score function and additional arguments
    score_data = call_function(data_list[[i]], score_function, ...)
    
    # Add the Type label for the current dataset
    score_data = cbind(Score = score_data, Type = type_labels[i])
    
    # Store in the result list
    result_list[[i]] = score_data
  }
  
  # Combine all results into a single data frame
  final_result = do.call(rbind, result_list)
  
  return(final_result)
}

#############################################################################

combine_traits_with_groups = function(processed_data_list) {
  combined_results = list()
  
  # Loop through each dataset type
  for (dataset_name in names(processed_data_list)) {
    dataset = processed_data_list[[dataset_name]]
    
    # Combine all sublists within the dataset into a single vector
    traits = unlist(dataset, use.names = FALSE)
    
    env=seq(1:10)
    
    
    # Create a data frame for this dataset
    combined_results[[dataset_name]] = data.frame(
      Trait = traits,
      Env =env
    )
  }
  
  # Combine all datasets into one large data frame
  final_result = do.call(rbind, combined_results)
  
  # Add a dataset identifier to distinguish between datasets
  final_result$Dataset = rep(names(combined_results), 
                              times = sapply(combined_results, nrow))
  
  return(final_result)
}


##############################################################################

processed_data_list = list(
  gaussian = preprocess_data(gaussian_data, trait_col = 4, genotype_col = 1, replicate_col = 2),
  sinusoidal = preprocess_data(sinusoidal_data, trait_col = 4, genotype_col = 1, replicate_col = 2),
  wave = preprocess_data(wave_data, trait_col = 4, genotype_col = 1, replicate_col = 2),
  linear = preprocess_data(linear_data, trait_col = 4, genotype_col = 1, replicate_col = 2)
)
linear = preprocess_data(linear_data, trait_col = 4, genotype_col = 1, replicate_col = 2)
gaussian = preprocess_data(gaussian_data, trait_col = 4, genotype_col = 1, replicate_col = 2)
sinusoidal = preprocess_data(sinusoidal_data, trait_col = 4, genotype_col = 1, replicate_col = 2)
wave = preprocess_data(wave_data, trait_col = 4, genotype_col = 1, replicate_col = 2)
type_labels = c("gaussian", "sinusoidal", "wave", "linear")

gaussian_long=combine_traits_with_groups(gaussian)
wave_long=combine_traits_with_groups(wave)
sinus_long=combine_traits_with_groups(sinusoidal)
linear_long=combine_traits_with_groups(linear)



test=data.frame(genotype=rep(1,40),replicate=c(rep(0,10),rep(1,10),rep(2,10),rep(3,10)),Trait=c(rep(1,20),rep(3,20)))
test_processed=preprocess_data(test,trait_col = 3,genotype_col = 1,replicate_col = 2)
print(test_processed)
CV_t_test=call_function(test_processed,calculate_CVt)
RN_test=call_function(test_processed,calculate_reaction_norm_slope)
#################################################################################################################################### 1.R
####################################################################################################################################


CV_t=post_process(processed_data_list,calculate_CVt,type_labels)



RN=post_process(processed_data_list,calculate_reaction_norm_slope,type_labels)



RNN=post_process(processed_data_list,calculate_reaction_norm_non_linear,type_labels,3)



D_slope=post_process(processed_data_list,calculate_D_slope,type_labels)



RC=post_process(processed_data_list,calculate_RC,type_labels)



#comparative score
CVm_gaussian=calculate_CVm(trait_values = gaussian_long[,1], group_labels = gaussian_long[,2])
CVm_wave=calculate_CVm(trait_values = wave_long[,1], group_labels = wave_long[,2])
CVm_sinus=calculate_CVm(trait_values = sinus_long[,1], group_labels = sinus_long[,2])
CVm_linear=calculate_CVm(trait_values = linear_long[,1], group_labels = linear_long[,2])



CVmd_gaussian=calculate_CVmd(trait_values = gaussian_long[,1], group_labels = gaussian_long[,2])
CVmd_wave=calculate_CVmd(trait_values = wave_long[,1], group_labels = wave_long[,2])
CVmd_sinus=calculate_CVmd(trait_values = sinus_long[,1], group_labels = sinus_long[,2])
CVmd_linear=calculate_CVmd(trait_values = linear_long[,1], group_labels = linear_long[,2])



Covariate=rep(1,10) #static covariate to minimize influence on the realisation of the score 
env=seq(1:10)# sequential environment
gPi=post_process(processed_data_list,calculate_grand_plasticity,type_labels,env,Covariate)


PPF=post_process(processed_data_list,calculate_PPF,type_labels,env,Covariate)



PPi=post_process(processed_data_list,calculate_Phenotypic_Plasticity_Index,type_labels)



PImd=post_process(processed_data_list,calculate_PImd,type_labels,env)



PILSM=post_process(processed_data_list,calculate_PILSM,type_labels,env,Covariate)



RTR=post_process(processed_data_list,calculate_RTR,type_labels,env)



PIR=post_process(processed_data_list,calculate_PIR,type_labels,env)



###############


RDPI=post_process(processed_data_list,calculate_rdpi,type_labels)



ESPI=post_process(processed_data_list,calculate_ESPI,type_labels)



ESPIID=post_process(processed_data_list,calculate_espiid,type_labels)






















