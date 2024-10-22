source("~/CRC 1622 - Z2/R-files/1.R")
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

############## Dataset generation - n=20

n_plants = 100
environmental_factor=c(1,1,1)

generate_random_baseline_and_variance = function(n_environments, n_traits) {
  baseline_values = matrix(
    rnorm(n_environments * n_traits, mean = sample(50:150, 1), sd = sample(5:15, 1)), 
    nrow = n_environments, ncol = n_traits
  )
  
  within_variance = matrix(
    rnorm(n_environments * n_traits, mean = sample(5:10, 1), sd = sample(1:3, 1)), 
    nrow = n_environments, ncol = n_traits
  )
  
  return(list(baseline_values = baseline_values, within_variance = within_variance))
}
########################################

generate_synthetic_data_with_skewness_type = function(n_plants, baseline_values, within_variance, skewed_traits = NULL, skewness_type = NULL, covariate_values = NULL) {
  n_traits = ncol(baseline_values)
  n_environments = nrow(baseline_values)
  
  # Ensure variances are strictly positive to avoid issues with log transformations
  within_variance[within_variance <= 0] = 1e-5
  
  # Create an additional column for covariate (or whatever column should be placed after the env_col)
  if (is.null(covariate_values)) {
    covariate_values = rep(1, n_plants * n_environments)  # Default covariate if not provided
  } else if (length(covariate_values) != n_plants * n_environments) {
    stop("Length of covariate_values must match the number of rows in the dataset")
  }
  
  # Initialize data frame with an additional column for the covariate
  synthetic_data = data.frame(matrix(nrow = n_plants * n_environments, ncol = n_traits + 2))
  rownames(synthetic_data) = rep(1:n_plants, times = n_environments) + rep(seq(0.1, (n_environments - 1) / 10 + 0.1, by = 0.1), each = n_plants)
  
  # Update column names to include the covariate after Env_indicator
  colnames(synthetic_data) = c("Identificator", "Environmental Factor", paste("Trait", 1:n_traits, sep = "_"))
  
  # If only one skewness type is provided, apply it to all skewed traits
  if (!is.null(skewness_type) && length(skewness_type) == 1) {
    skewness_type = rep(skewness_type, length(skewed_traits))
  }
  
  for (env in 1:n_environments) {
    start_idx = (env - 1) * n_plants + 1
    end_idx = env * n_plants
    
    # Assign the environment indicator and covariate values
    synthetic_data[start_idx:end_idx, 1] = env  # Env_indicator
    synthetic_data[start_idx:end_idx, 2] = covariate_values[start_idx:end_idx]  # Covariate
    
    for (i in 1:n_traits) {
      # Check if the trait should be skewed
      if (!is.null(skewed_traits) && i %in% skewed_traits) {
        trait_idx = which(skewed_traits == i)
        skew_type = skewness_type[trait_idx]
        
        # Apply the chosen skewness type
        if (!is.na(skew_type) && skew_type == "lognormal") {
          baseline_values[env, i] = max(baseline_values[env, i], 1e-5)
          within_variance[env, i] = max(within_variance[env, i], 1e-5)
          synthetic_data[start_idx:end_idx, i + 2] = rlnorm(n_plants, meanlog = log(baseline_values[env, i]), sdlog = within_variance[env, i])
        } else if (!is.na(skew_type) && skew_type == "gamma") {
          shape = (baseline_values[env, i] / within_variance[env, i])
          rate = baseline_values[env, i] / (within_variance[env, i]^2)
          synthetic_data[start_idx:end_idx, i + 2] = rgamma(n_plants, shape = shape, rate = rate)
        } else if (!is.na(skew_type) && skew_type == "weibull") {
          shape = 2  # You can adjust the shape parameter for Weibull as needed
          scale = baseline_values[env, i]
          synthetic_data[start_idx:end_idx, i + 2] = rweibull(n_plants, shape = shape, scale = scale)
        } else {
          stop("Unsupported skewness type specified.")
        }
      } else {
        # Default to normal distribution
        synthetic_data[start_idx:end_idx, i + 2] = rnorm(n_plants, mean = baseline_values[env, i], sd = within_variance[env, i])
      }
    }
  }
  
  return(synthetic_data)
}
########################################

# Generate 20 different datasets, each with unique baseline values and variances
n_datasets = 20
datasets = list()
datasets_log_skewed=list()
datasets_log_log=list()
datasets_gamma_skewed=list()
datasets_gamma_gamma=list()
datasets_weibull_skewed=list()
datasets_weibull_weibull=list()
datasets_log_gamma=list()
datasets_log_weibull=list()
datasets_gamma_weibull=list()
skewness_type = c("lognormal","gamma","weibull")



for (i in 1:n_datasets) {
  set.seed(12345 + i)  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets[[i]] = generate_synthetic_data(n_plants, params$baseline_values, params$within_variance,environmental_factor)
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_log_skewed[[i]] = generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=1, skewness_type[1])
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_log_log[[i]] = generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=c(1,2), skewness_type[1])
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_log_gamma[[i]] = generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=c(1,2), c(skewness_type[1],skewness_type[2]))
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_log_weibull[[i]] = generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=c(1,2), c(skewness_type[1],skewness_type[3]))
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_gamma_skewed[[i]]=generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=1, skewness_type[2])
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_gamma_gamma[[i]]=generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=c(1,2), skewness_type[2])
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_gamma_weibull[[i]]=generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=c(1,2), c(skewness_type[2],skewness_type[3]))
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_weibull_skewed[[i]]=generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=1, skewness_type[3])
  
  params = generate_random_baseline_and_variance(n_environments = 3, n_traits = 3)
  datasets_weibull_weibull[[i]] = generate_synthetic_data_with_skewness_type(n_plants, params$baseline_values, params$within_variance, skewed_traits=c(1,2), c(skewness_type[3],skewness_type[3]))
  
}


#####################################



CV_t=list(
  non_traitwise=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  
  CV_t[[1]][[i]]=calculate_CVt(datasets[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[2]][[i]]=calculate_CVt(datasets_log_skewed[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[3]][[i]]=calculate_CVt(datasets_log_log[[i]],exclude_col=c(1,2),traitwise=T) #note that for datasets with differing skewness combinations it doesnt make sense to calculate the score traitwise since otherise the score loses its info on the skewness
  CV_t[[4]][[i]]=calculate_CVt(datasets_log_gamma[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[5]][[i]]=calculate_CVt(datasets_log_weibull[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[6]][[i]]=calculate_CVt(datasets_gamma_skewed[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[7]][[i]]=calculate_CVt(datasets_gamma_gamma[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[8]][[i]]=calculate_CVt(datasets_gamma_weibull[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[9]][[i]]=calculate_CVt(datasets_weibull_skewed[[i]],exclude_col=c(1,2),traitwise=T)
  CV_t[[10]][[i]]=calculate_CVt(datasets_weibull_weibull[[i]],exclude_col=c(1,2),traitwise=T)
}

####################################


#reaction norm slope
RNS=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  RNS[[1]][[i]]=calculate_reaction_norm_slope(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[2]][[i]]=calculate_reaction_norm_slope(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[3]][[i]]=calculate_reaction_norm_slope(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  RNS[[4]][[i]]=calculate_reaction_norm_slope(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[5]][[i]]=calculate_reaction_norm_slope(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[6]][[i]]=calculate_reaction_norm_slope(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[7]][[i]]=calculate_reaction_norm_slope(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[8]][[i]]=calculate_reaction_norm_slope(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[9]][[i]]=calculate_reaction_norm_slope(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RNS[[10]][[i]]=calculate_reaction_norm_slope(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}

####################################

D_slope=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)



for(i in 1:n_datasets){
  
  D_slope[[1]][[i]]=calculate_D_slope(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[2]][[i]]=calculate_D_slope(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[3]][[i]]=calculate_D_slope(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  D_slope[[4]][[i]]=calculate_D_slope(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[5]][[i]]=calculate_D_slope(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[6]][[i]]=calculate_D_slope(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[7]][[i]]=calculate_D_slope(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[8]][[i]]=calculate_D_slope(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[9]][[i]]=calculate_D_slope(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  D_slope[[10]][[i]]=calculate_D_slope(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}


#########################################


RC=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)



for(i in 1:n_datasets){
  RC[[1]][[i]]=calculate_RC(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[2]][[i]]=calculate_RC(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[3]][[i]]=calculate_RC(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  RC[[4]][[i]]=calculate_RC(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[5]][[i]]=calculate_RC(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[6]][[i]]=calculate_RC(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[7]][[i]]=calculate_RC(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[8]][[i]]=calculate_RC(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[9]][[i]]=calculate_RC(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RC[[10]][[i]]=calculate_RC(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}
########################################


CVm=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)



for(i in 1:n_datasets){
  CVm[[1]][[i]]=calculate_CVm(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[2]][[i]]=calculate_CVm(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[3]][[i]]=calculate_CVm(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  CVm[[4]][[i]]=calculate_CVm(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[5]][[i]]=calculate_CVm(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[6]][[i]]=calculate_CVm(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[7]][[i]]=calculate_CVm(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[8]][[i]]=calculate_CVm(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[9]][[i]]=calculate_CVm(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVm[[10]][[i]]=calculate_CVm(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}


#############################################



CVmd=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)



for(i in 1:n_datasets){
  CVmd[[1]][[i]]=calculate_CVmd(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[2]][[i]]=calculate_CVmd(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[3]][[i]]=calculate_CVmd(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  CVmd[[4]][[i]]=calculate_CVmd(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[5]][[i]]=calculate_CVmd(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[6]][[i]]=calculate_CVmd(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[7]][[i]]=calculate_CVmd(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[8]][[i]]=calculate_CVmd(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[9]][[i]]=calculate_CVmd(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  CVmd[[10]][[i]]=calculate_CVmd(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}


##############################################





Pi_adj_mean=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)



for(i in 1:n_datasets){
  Pi_adj_mean[[1]][[i]]=calculate_grand_plasticity(datasets[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[2]][[i]]=calculate_grand_plasticity(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[3]][[i]]=calculate_grand_plasticity(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2) 
  Pi_adj_mean[[4]][[i]]=calculate_grand_plasticity(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[5]][[i]]=calculate_grand_plasticity(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[6]][[i]]=calculate_grand_plasticity(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[7]][[i]]=calculate_grand_plasticity(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[8]][[i]]=calculate_grand_plasticity(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col = 2)
  Pi_adj_mean[[9]][[i]]=calculate_grand_plasticity(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col =2 )
  Pi_adj_mean[[10]][[i]]=calculate_grand_plasticity(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),control_env = 1,covariate_col =2 )
}


#############################################

#calculate_PPF could potentially cause issues in when comparing to other scores as the score is comparative bewteen 2 envs. In the case of multiple environmenrts the score needs to be calculated pairwise


PPF=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  PPF[[1]][[i]]=calculate_PPF(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[2]][[i]]=calculate_PPF(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[3]][[i]]=calculate_PPF(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  PPF[[4]][[i]]=calculate_PPF(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[5]][[i]]=calculate_PPF(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[6]][[i]]=calculate_PPF(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[7]][[i]]=calculate_PPF(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[8]][[i]]=calculate_PPF(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[9]][[i]]=calculate_PPF(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PPF[[10]][[i]]=calculate_PPF(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}


#############################################



Pi=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  Pi[[1]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets[[i]],trait_cols = c(3,4,5))
  Pi[[2]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_log_skewed[[i]],trait_cols = c(3,4,5))
  Pi[[3]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_log_log[[i]],trait_cols = c(3,4,5)) 
  Pi[[4]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_log_gamma[[i]],trait_cols = c(3,4,5))
  Pi[[5]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_log_weibull[[i]],trait_cols = c(3,4,5))
  Pi[[6]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_gamma_skewed[[i]],trait_cols = c(3,4,5))
  Pi[[7]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_gamma_gamma[[i]],trait_cols = c(3,4,5))
  Pi[[8]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_gamma_weibull[[i]],trait_cols = c(3,4,5))
  Pi[[9]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_weibull_skewed[[i]],trait_cols = c(3,4,5))
  Pi[[10]][[i]]=calculate_Phenotypic_Plasticity_Index(datasets_weibull_weibull[[i]],trait_cols = c(3,4,5))
  
}
  ########################################
  


PImd=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)

  
  for(i in 1:n_datasets){
    PImd[[1]][[i]]=calculate_PImd(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[2]][[i]]=calculate_PImd(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[3]][[i]]=calculate_PImd(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
    PImd[[4]][[i]]=calculate_PImd(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[5]][[i]]=calculate_PImd(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[6]][[i]]=calculate_PImd(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[7]][[i]]=calculate_PImd(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[8]][[i]]=calculate_PImd(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[9]][[i]]=calculate_PImd(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
    PImd[[10]][[i]]=calculate_PImd(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  }
  
############################################





PILSM=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  PILSM[[1]][[i]]=calculate_PILSM(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[2]][[i]]=calculate_PILSM(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[3]][[i]]=calculate_PILSM(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5)) 
  PILSM[[4]][[i]]=calculate_PILSM(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[5]][[i]]=calculate_PILSM(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[6]][[i]]=calculate_PILSM(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[7]][[i]]=calculate_PILSM(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[8]][[i]]=calculate_PILSM(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[9]][[i]]=calculate_PILSM(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PILSM[[10]][[i]]=calculate_PILSM(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
}


#################################################


RTR=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  RTR[[1]][[i]]=calculate_RTR(datasets[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[2]][[i]]=calculate_RTR(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[3]][[i]]=calculate_RTR(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2) 
  RTR[[4]][[i]]=calculate_RTR(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[5]][[i]]=calculate_RTR(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[6]][[i]]=calculate_RTR(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[7]][[i]]=calculate_RTR(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[8]][[i]]=calculate_RTR(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[9]][[i]]=calculate_RTR(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
  RTR[[10]][[i]]=calculate_RTR(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),env_low=1,env_high=2)
}

###############################################


PIR=list(
  non=list(),
  log=list(),
  log_log=list(),
  log_gamma=list(),
  log_weibull=list(),
  gamma=list(),
  gamma_gamma=list(),
  gamma_weibull=list(),
  weibull=list(),
  weibull_weibull=list()
)


for(i in 1:n_datasets){
  PIR[[1]][[i]]=calculate_PIR(datasets[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[2]][[i]]=calculate_PIR(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[3]][[i]]=calculate_PIR(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2) 
  PIR[[4]][[i]]=calculate_PIR(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[5]][[i]]=calculate_PIR(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[6]][[i]]=calculate_PIR(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[7]][[i]]=calculate_PIR(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[8]][[i]]=calculate_PIR(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[9]][[i]]=calculate_PIR(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
  PIR[[10]][[i]]=calculate_PIR(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),rgr_col=2)
}


############################ correlation
# Initialize list to store correlation matrices for each trait separately across all dataset types
correlation_matrices_all_datasets = list()

dataset_types=list(datasets,
                   datasets_log_skewed,
                   datasets_log_log,
                   datasets_gamma_skewed,
                   datasets_gamma_gamma,
                   datasets_weibull_skewed,
                   datasets_weibull_weibull,
                   datasets_log_gamma,
                   datasets_log_weibull,
                   datasets_gamma_weibull)


# Loop through the 10 different dataset types
for (dataset_type_index in 1:length(dataset_types)) {
  
  # Initialize list to store correlation matrices per trait for the current dataset type
  correlation_matrices_per_trait = list()
  
  # Loop through each trait (3 traits)
  for (trait_index in 1:3) {
    
    
    # Initialize a data frame to store scores for all datasets for the current trait
    scores_df = data.frame(
      CV_t = numeric(n_datasets),
      RNS = numeric(n_datasets),
      D_slope = numeric(n_datasets),
      RC = numeric(n_datasets),
      CVm = numeric(n_datasets),
      CVmd = numeric(n_datasets),
      Pi_adj_mean = numeric(n_datasets),
      #PPF = numeric(n_datasets),
      Pi = numeric(n_datasets),
      PImd = numeric(n_datasets),
      PILSM = numeric(n_datasets),
      RTR = numeric(n_datasets),
      PIR = numeric(n_datasets)
    )
    
    # Loop over the 20 datasets for the current dataset type
    for (i in 1:n_datasets) {
      
      # For each dataset, extract the scores for the current trait and store them in the scores_df
      scores_df$CV_t[i] = as.numeric(CV_t[[dataset_type_index]][[i]][trait_index])
      scores_df$RNS[i] = as.numeric(RNS[[dataset_type_index]][[i]][trait_index])
      scores_df$D_slope[i] = as.numeric(D_slope[[dataset_type_index]][[i]][trait_index])
      scores_df$RC[i] = as.numeric(RC[[dataset_type_index]][[i]][trait_index])
      scores_df$CVm[i] = as.numeric(CVm[[dataset_type_index]][[i]][trait_index])
      scores_df$CVmd[i] = as.numeric(CVmd[[dataset_type_index]][[i]][trait_index])
      scores_df$Pi_adj_mean[i] = as.numeric(Pi_adj_mean[[dataset_type_index]][[i]][trait_index])
      #scores_df$PPF[i] = as.numeric(PPF[[dataset_type_index]][[i]][trait_index])
      scores_df$Pi[i] = as.numeric(Pi[[dataset_type_index]][[i]][trait_index])
      scores_df$PImd[i] = as.numeric(PImd[[dataset_type_index]][[i]][trait_index])
      scores_df$PILSM[i] = as.numeric(PILSM[[dataset_type_index]][[i]][trait_index][[1]])
      scores_df$RTR[i] = as.numeric(RTR[[dataset_type_index]][[i]][trait_index])
      scores_df$PIR[i] = as.numeric(PIR[[dataset_type_index]][[i]][trait_index])
    }
    
    # Calculate the correlation matrix for the current trait across the 20 datasets
    correlation_matrix = cor(scores_df, use = "pairwise.complete.obs")
    
    # Store the correlation matrix for the current trait in the list
    correlation_matrices_per_trait[[trait_index]] = correlation_matrix
  }
  
  # Store the correlation matrices for all 3 traits for the current dataset type
  correlation_matrices_all_datasets[[dataset_type_index]] = correlation_matrices_per_trait
}

# Output: correlation_matrices_all_traits will be a list of correlation matrices
# - Outer list: dataset types (10 types)
# - Inner list: traits (3 traits)
# - Each element: a correlation matrix for the scores for that trait and dataset type
