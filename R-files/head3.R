# Script:    head3.R
# Purpose:   Load all plasticity‐score functions, set sampling parameters (range,
#            interval, indices, covariate, type labels, etc.), preprocess reaction‐
#            norm data into per‐genotype trait vectors, and compute every plasticity
#            score across the four genotype forms (linear, gaussian, sinusoidal, wave).
#
# Steps:
#   1. Source the five score‐definition files (1_2.R, 2.R, 3.R, 4.R, 5.R) plus
#      the linear and nonlinear norms generators.
#   2. Define default globals (indices, interval, env, Covariate, type_labels) and
#      only override them if already defined in the calling environment (which happens when the script is launched from the resolution.R file).
#   3. Ensure the raw reaction‐norm data frames (gaussian_data, sinusoidal_data,
#      wave_data, individual_norms) are present, assigning defaults if needed.
#   4. Combine the four raw data sets into one data frame for any combined‐data ops.
#   5. Preprocess each form’s data via `preprocess_data()` to:
#        • average replicates per genotype × environment,
#        • subsample at `indices`,
#        • return a named list of genotype‐vectors.
#   6. For each of the plasticity functions :
#        • Call them over each genotype list via `call_function()`,
#        • Post‐process into a unified data frame with `Type` labels.
#   7. Build four reaction‐norm matrices (`linear`, `gaussian`, `sinusoidal`, `wave`)
#      by row‐binding the per‐genotype vectors.
#   8. Make “long” trait + environment data frames for plotting or summary stats.
#
# Outputs (in memory):
#   • `processed_data_list`: list of four preprocessed genotype‐lists
#   • Score tables
#   • Reaction‐norm matrices: `linear_matrix`, `gaussian_matrix`, `sinusoidal_matrix`,
#     `wave_matrix`
#
# Globals used/produced:
#   indices, interval, env, Covariate, type_labels,
#   gaussian_data, sinusoidal_data, wave_data, linear_data

source("~/CRC_1644_Z2/R-files/Plasticity_scores/1_2.R")
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
source("~/CRC_1644_Z2/R-files/Plasticity_scores/2.R")
#RDPI	(rdpi_calculation) - tested,
#RDPIs (rdpi_mean_calculation) - tested,
#ESPI (calculate_ESPI) - tested,
#ESPIid (espiid_calculation) - tested,
#evwpi_calculation (idea from Benedikt)
source("~/CRC_1644_Z2/R-files/Plasticity_scores/3.R")
#Phenotypic Stability Index (calculate_PSI),
#Relative Plasticity Index (calculate_RPI) - tested,
#Plasticity Quotient (calculate_PQ) - tested,
#Phenotypic Range (PR) (calculate_PR) - tested,
#Norm of reaction width (calculate_NRW) - tested,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
#Calculate Plasticity Differential (PD) (calculate_PD) - tested,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,
source("~/CRC_1644_Z2/R-files/Plasticity_scores/4.R")
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

source("~/CRC_1644_Z2/R-files/Plasticity_scores/5.R")
source("~/CRC_1644_Z2/R-files/norms_generator2_nonlinear.R")
source("~/CRC_1644_Z2/R-files/norms_generator2_linear.R")


## ──────────────────────────────────────────────────────────────────────────────
##  1) Define defaults for any globals you need
## ──────────────────────────────────────────────────────────────────────────────
default_indices   <- seq(1, 50)   #sampling over the whole interval every index
default_env       <- seq_along(default_indices)
default_Covariate <- rep(1, length(default_indices))
default_type_labels <- c("gaussian", "sinusoidal", "wave", "linear")

## ──────────────────────────────────────────────────────────────────────────────
##  2) Only assign them if they don't already exist
## ──────────────────────────────────────────────────────────────────────────────
if (!exists("indices",   envir = .GlobalEnv)) indices      <- default_indices
if (!exists("env",       envir = .GlobalEnv)) env          <- default_env
if (!exists("Covariate", envir = .GlobalEnv)) Covariate    <- default_Covariate
if (!exists("type_labels", envir = .GlobalEnv)) type_labels <- default_type_labels

## ──────────────────────────────────────────────────────────────────────────────
##  3) And likewise for your data‐frames, before you try to rbind them
## ──────────────────────────────────────────────────────────────────────────────
if (!exists("gaussian_data",  envir = .GlobalEnv))   gaussian_rn  = gaussian_data  
if (!exists("sinusoidal_data",envir = .GlobalEnv))   sinusoidal_rn = sinusoidal_data
if (!exists("wave_data",      envir = .GlobalEnv))   wave_rn       = wave_data       
if (!exists("linear_data",envir = .GlobalEnv))      linear_data     = individual_norms




combined=rbind(linear_data,gaussian_data, sinusoidal_data,wave_data )

combined_data=list(linear=linear_data,
                   gaussian=gaussian_data,
                   sinusoidal=sinusoidal_data,
                   wave=wave_data)


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
preprocess_data = function(data, trait_col, genotype_col, replicate_col = NULL, indices = NULL) {
  if (is.numeric(trait_col)) trait_col = colnames(data)[trait_col]
  if (is.numeric(genotype_col)) genotype_col = colnames(data)[genotype_col]
  if (!is.null(replicate_col) && is.numeric(replicate_col)) replicate_col = colnames(data)[replicate_col]
  
  genotype_split = split(data, data[[genotype_col]])
  result = list()
  
  for (genotype in names(genotype_split)) {
    genotype_data = genotype_split[[genotype]]
    
    if (is.null(replicate_col)) {
      result[[genotype]] = genotype_data[[trait_col]][seq(1, length(genotype_data[[trait_col]]), 1)]
    } else {
      unique_replicates = unique(genotype_data[[replicate_col]])
      
      temp = do.call(
        cbind,
        lapply(unique_replicates, function(rep) {
          genotype_data[genotype_data[[replicate_col]] == rep, trait_col]
        })
      )
      
      averaged_vector = rowMeans(temp, na.rm = TRUE)
      
      # If indices are provided externally, use them.
      # Otherwise, compute default indices over the entire vector.
      if (is.null(indices)) {
        indices = unique(c(seq(1, length(averaged_vector), by = interval),
                           length(averaged_vector)))
      }
    
      result[[genotype]] = averaged_vector[indices]
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
    
    env = seq(1, length(traits))
  
    
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

# Now call preprocess_data with the indices parameter.
processed_data_list <- list(
  gaussian = preprocess_data(gaussian_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices),
  sinusoidal = preprocess_data(sinusoidal_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices),
  wave = preprocess_data(wave_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices),
  linear = preprocess_data(linear_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
)

linear <- preprocess_data(linear_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
gaussian <- preprocess_data(gaussian_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
sinusoidal <- preprocess_data(sinusoidal_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
wave <- preprocess_data(wave_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
type_labels <- c("gaussian", "sinusoidal", "wave", "linear")


gaussian_long <- combine_traits_with_groups(gaussian)
wave_long <- combine_traits_with_groups(wave)
sinus_long <- combine_traits_with_groups(sinusoidal)
linear_long <- combine_traits_with_groups(linear)

#test <- data.frame(genotype = rep(1, 40),
#                   replicate = c(rep(0, 10), rep(1, 10), rep(2, 10), rep(3, 10)),
#                   Trait = c(rep(1, 20), rep(3, 20)))
#test_processed <- preprocess_data(test, trait_col = 3, genotype_col = 1, replicate_col = 2, indices = indices)
#
#CV_t_test <- call_function(test_processed, calculate_CVt)
#RN_test <- call_function(test_processed, calculate_reaction_norm_slope)
#################################################################################################################################### 1.R
####################################################################################################################################

CV_t=post_process(processed_data_list,calculate_CVt,type_labels)

RN=post_process(processed_data_list,calculate_reaction_norm_slope,type_labels)

RNN=post_process(processed_data_list,calculate_reaction_norm_non_linear,type_labels,3)

D_slope=post_process(processed_data_list,calculate_D_slope,type_labels)

RC=post_process(processed_data_list,calculate_RC,type_labels)

#comparative score -- not one score for a single genotype but for a population of different genotypes. Therefore not suitable for direct comparison 
CVm_gaussian=calculate_CVm(trait_values = gaussian)
CVm_wave=calculate_CVm(trait_values = wave)
CVm_sinus=calculate_CVm(trait_values = sinusoidal)
CVm_linear=calculate_CVm(trait_values = linear)

#comparative score -- not one score for a single genotype but for a population of different genotypes. Therefore not suitable for direct comparison 
CVmd_gaussian=calculate_CVmd(trait_values = gaussian)
CVmd_wave=calculate_CVmd(trait_values = wave)
CVmd_sinus=calculate_CVmd(trait_values = sinusoidal)
CVmd_linear=calculate_CVmd(trait_values = linear)

Covariate=rep(1,length(indices)) #static covariate to minimize influence on the realisation of the score 
# sequential environment


gPi=post_process(processed_data_list,calculate_grand_plasticity,type_labels,env,Covariate)

PPF=post_process(processed_data_list,calculate_PPF,type_labels,env,Covariate)

PPi=post_process(processed_data_list,calculate_Phenotypic_Plasticity_Index,type_labels)

PImd=post_process(processed_data_list,calculate_PImd,type_labels,env)

PILSM=post_process(processed_data_list,calculate_PILSM,type_labels,env,Covariate)

RTR=post_process(processed_data_list,calculate_RTR,type_labels,env)

PIR=post_process(processed_data_list,calculate_PIR,type_labels,env)



############### 2.R
###############

RDPI=post_process(processed_data_list,calculate_rdpi,type_labels)

ESPI=post_process(processed_data_list,calculate_ESPI,type_labels)

ESPIID=post_process(processed_data_list,calculate_espiid,type_labels)



################# 3.R
#################


PSI=post_process(processed_data_list,calculate_PSI,type_labels)

RPI=post_process(processed_data_list,calculate_RPI,type_labels)

PQ=post_process(processed_data_list,calculate_PQ,type_labels)

PR=post_process(processed_data_list,calculate_PR,type_labels) #not really comparable as it required individual level data

NRW=post_process(processed_data_list,calculate_NRW,type_labels)

ESP=post_process(processed_data_list,calculate_ESP,type_labels)

PD=post_process(processed_data_list,calculate_general_PD,type_labels) # not comparable as it required stress env and control env

FPI=post_process(processed_data_list,calculate_FPI,type_labels) # not comparable as it required stress env and control env

################## 4.R
##################


add_arg=rep(1,length(linear[1][[1]]))

#DPI=post_process(processed_data_list,calculate_DPI,type_labels,add_arg) # not comparable as this requires time resolved data

CEV=post_process(processed_data_list,calculate_CEV,type_labels) 

#PRI=post_process(processed_data_list,calculate_PRI,type_labels) not comparable as extreme env needs to be specified 

PFI=post_process(processed_data_list,calculate_PFI,type_labels) 

APC=post_process(processed_data_list,calculate_APC,type_labels) 

SI=post_process(processed_data_list,calculate_SI,type_labels) 

RSI=post_process(processed_data_list,calculate_RSI,type_labels) 

EVS=post_process(processed_data_list,calculate_EVS,type_labels) 

MVPi=post_process(processed_data_list,calculate_MVPi,type_labels) 

################### 5.R
###################

#Plasticity=post_process(processed_data_list,calculate_plasticity,type_labels)
#
#env_cov=post_process(processed_data_list,cross_env_cov,type_labels)
#

keys <- names(processed_data_list)
fw_matrices <- setNames(vector("list", length(keys)), keys)

for (k in keys) {
  gl <- processed_data_list[[k]]
  G <- names(gl)
  L <- max(lengths(gl))
  env_ids <- seq_len(L)
  M <- matrix(NA_real_, nrow=length(G), ncol=L,
              dimnames=list(G, paste0("E", env_ids)))
  for (g in G) {
    v <- as.numeric(gl[[g]])
    M[g, seq_along(v)] <- v
  }
  fw_matrices[[k]] <- as.data.frame(M)
}

linear_matrix     <- fw_matrices$linear
gaussian_matrix   <- fw_matrices$gaussian
sinusoidal_matrix <- fw_matrices$sinusoidal
wave_matrix       <- fw_matrices$wave

fw_long_all <- data.frame(Dataset=character(), genotype=character(),
                          environment=integer(), y=double(), 
                          stringsAsFactors=FALSE)

for (k in keys) {
  gl <- processed_data_list[[k]]
  G <- names(gl)
  for (g in G) {
    v <- as.numeric(gl[[g]])
    n <- length(v)
    df <- data.frame(Dataset=k, genotype=g, environment=seq_len(n), y=v)
    fw_long_all <- rbind(fw_long_all, df)
  }
}









keys <- names(processed_data_list)
row_list <- list()
fw_long_all <- data.frame(Dataset=character(), genotype=character(), environment=integer(), y=double(), stringsAsFactors=FALSE)
Lmax <- max(unlist(lapply(processed_data_list, function(gl) max(lengths(gl)))))

for (k in keys) {
  gl <- processed_data_list[[k]]
  G <- names(gl)
  for (g in G) {
    uid <- paste(k, g, sep="_")
    v <- as.numeric(gl[[g]])
    n <- length(v)
    row_list[[uid]] <- c(v, rep(NA_real_, Lmax - n))
    fw_long_all <- rbind(fw_long_all, data.frame(Dataset=k, genotype=uid, environment=seq_len(n), y=v))
  }
}

fw_matrix <- do.call(rbind, row_list)
rownames(fw_matrix) <- names(row_list)
colnames(fw_matrix) <- paste0("E", seq_len(ncol(fw_matrix)))
fw_df <- as.data.frame(fw_matrix)



FW=calculate_finlay_wilkinson(fw_df,plot=T)


library(ggplot2)
Y <- as.matrix(fw_df[28:35,])
gnames <- rownames(Y); if (is.null(gnames)) gnames <- paste0("G", seq_len(nrow(Y)))
env_values <- if (exists("env_values") && !is.null(env_values)) as.numeric(env_values) else seq_len(ncol(Y))
df <- data.frame(genotype=rep(gnames, each=ncol(Y)),
                 environment=rep(env_values, times=nrow(Y)),
                 y=as.vector(t(Y)))
df <- df[is.finite(df$y) & is.finite(df$environment), ]
ggplot(df, aes(environment, y, color=genotype, group=genotype)) +
  geom_line() + geom_point() +
  labs(x="Environment", y="Trait", title="Trait by Environment") +
  theme_bw()




