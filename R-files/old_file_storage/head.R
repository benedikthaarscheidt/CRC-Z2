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

skewness_type_vector = c(
  # Dataset 1 (Non-skewed for all traits)
  "non", "non", "non",
  # Dataset 2 (Log skewed first trait)
  "lognormal", "non", "non",
  # Dataset 3 (Log skewed first two traits)
  "lognormal", "lognormal", "non",
  # Dataset 4 (Log first, Gamma second)
  "lognormal", "gamma", "non",
  # Dataset 5 (Log first, Weibull second)
  "lognormal", "weibull", "non",
  # Dataset 6 (Gamma skewed first trait)
  "gamma", "non", "non",
  # Dataset 7 (Gamma first two traits)
  "gamma", "gamma", "non",
  # Dataset 8 (Gamma first, Weibull second)
  "gamma", "weibull", "non",
  # Dataset 9 (Weibull first trait)
  "weibull", "non", "non",
  # Dataset 10 (Weibull first two traits)
  "weibull", "weibull", "non"
)









################################################################################################################################### 1.R
###################################################################################################################################




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



################################################################################################################################### 2.R
###################################################################################################################################
#dataframe, trait_cols, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL, stat_analysis = NULL
#this is a comparative score comparing species or environments
RDPI=list(
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
  RDPI[[1]][[i]]=rdpi_calculation(datasets[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[2]][[i]]=rdpi_calculation(datasets_log_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[3]][[i]]=rdpi_calculation(datasets_log_log[[i]],factors = 1,trait_cols = c(3,4,5)) 
  RDPI[[4]][[i]]=rdpi_calculation(datasets_log_gamma[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[5]][[i]]=rdpi_calculation(datasets_log_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[6]][[i]]=rdpi_calculation(datasets_gamma_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[7]][[i]]=rdpi_calculation(datasets_gamma_gamma[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[8]][[i]]=rdpi_calculation(datasets_gamma_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[9]][[i]]=rdpi_calculation(datasets_weibull_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI[[10]][[i]]=rdpi_calculation(datasets_weibull_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  print(i)
}



################################


#this is a comparative score comparing species or environments
RDPI_mean=list(
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
  RDPI_mean[[1]][[i]]=rdpi_mean_calculation(datasets[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[2]][[i]]=rdpi_mean_calculation(datasets_log_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[3]][[i]]=rdpi_mean_calculation(datasets_log_log[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[4]][[i]]=rdpi_mean_calculation(datasets_log_gamma[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[5]][[i]]=rdpi_mean_calculation(datasets_log_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[6]][[i]]=rdpi_mean_calculation(datasets_gamma_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[7]][[i]]=rdpi_mean_calculation(datasets_gamma_gamma[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[8]][[i]]=rdpi_mean_calculation(datasets_gamma_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[9]][[i]]=rdpi_mean_calculation(datasets_weibull_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  RDPI_mean[[10]][[i]]=rdpi_mean_calculation(datasets_weibull_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  print(i)
}


##############################



#this is a comparative score comparing species or environments
ESPI=list(
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
  ESPI[[1]][[i]]=calculate_ESPI(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[2]][[i]]=calculate_ESPI(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[3]][[i]]=calculate_ESPI(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[4]][[i]]=calculate_ESPI(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[5]][[i]]=calculate_ESPI(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[6]][[i]]=calculate_ESPI(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[7]][[i]]=calculate_ESPI(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[8]][[i]]=calculate_ESPI(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[9]][[i]]=calculate_ESPI(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  ESPI[[10]][[i]]=calculate_ESPI(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  print(i)
}

##################################

# this score might not be suitable for comparison since it compares envs by trait combinations and therefor for the same dataset yields more scores 
ESPIID=list(
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
  ESPIID[[1]][[i]]=espiid_calculation(datasets[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[2]][[i]]=espiid_calculation(datasets_log_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[3]][[i]]=espiid_calculation(datasets_log_log[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[4]][[i]]=espiid_calculation(datasets_log_gamma[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[5]][[i]]=espiid_calculation(datasets_log_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[6]][[i]]=espiid_calculation(datasets_gamma_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[7]][[i]]=espiid_calculation(datasets_gamma_gamma[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[8]][[i]]=espiid_calculation(datasets_gamma_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[9]][[i]]=espiid_calculation(datasets_weibull_skewed[[i]],factors = 1,trait_cols = c(3,4,5))
  ESPIID[[10]][[i]]=espiid_calculation(datasets_weibull_weibull[[i]],factors = 1,trait_cols = c(3,4,5))
  print(i)
}

################################################################################################################################### 3.R
###################################################################################################################################


PSI=list(
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
  PSI[[1]][[i]]=calculate_PSI(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[2]][[i]]=calculate_PSI(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[3]][[i]]=calculate_PSI(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[4]][[i]]=calculate_PSI(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[5]][[i]]=calculate_PSI(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[6]][[i]]=calculate_PSI(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[7]][[i]]=calculate_PSI(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[8]][[i]]=calculate_PSI(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[9]][[i]]=calculate_PSI(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  PSI[[10]][[i]]=calculate_PSI(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  print(i)
}



###############################







# this is to be compared to ESPIID and some other score which I still have to determine
RPI=list(
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
  RPI[[1]][[i]]=calculate_RPI(datasets[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[2]][[i]]=calculate_RPI(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[3]][[i]]=calculate_RPI(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[4]][[i]]=calculate_RPI(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[5]][[i]]=calculate_RPI(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[6]][[i]]=calculate_RPI(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[7]][[i]]=calculate_RPI(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[8]][[i]]=calculate_RPI(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[9]][[i]]=calculate_RPI(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5))
  RPI[[10]][[i]]=calculate_RPI(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5))
  print(i)
}

########################################




PQ=list(
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
  PQ[[1]][[i]]=calculate_PQ(datasets[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[2]][[i]]=calculate_PQ(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[3]][[i]]=calculate_PQ(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[4]][[i]]=calculate_PQ(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[5]][[i]]=calculate_PQ(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[6]][[i]]=calculate_PQ(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[7]][[i]]=calculate_PQ(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[8]][[i]]=calculate_PQ(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[9]][[i]]=calculate_PQ(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  PQ[[10]][[i]]=calculate_PQ(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),factor_col=1)
  print(i)
}









###############################################

#this is supposed to be compared with the ESPIID 
  
PR=list(
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
  PR[[1]][[i]]=calculate_PR(datasets[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[2]][[i]]=calculate_PR(datasets_log_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[3]][[i]]=calculate_PR(datasets_log_log[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[4]][[i]]=calculate_PR(datasets_log_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[5]][[i]]=calculate_PR(datasets_log_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[6]][[i]]=calculate_PR(datasets_gamma_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[7]][[i]]=calculate_PR(datasets_gamma_gamma[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[8]][[i]]=calculate_PR(datasets_gamma_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[9]][[i]]=calculate_PR(datasets_weibull_skewed[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  PR[[10]][[i]]=calculate_PR(datasets_weibull_weibull[[i]],env_col = 1,trait_cols = c(3,4,5),across=T)
  print(i)
}

###############################################



NRW=list(
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
  NRW[[1]][[i]]=calculate_NRW(datasets[[i]],trait_cols = c(3,4,5))
  NRW[[2]][[i]]=calculate_NRW(datasets_log_skewed[[i]],trait_cols = c(3,4,5))
  NRW[[3]][[i]]=calculate_NRW(datasets_log_log[[i]],trait_cols = c(3,4,5))
  NRW[[4]][[i]]=calculate_NRW(datasets_log_gamma[[i]],trait_cols = c(3,4,5))
  NRW[[5]][[i]]=calculate_NRW(datasets_log_weibull[[i]],trait_cols = c(3,4,5))
  NRW[[6]][[i]]=calculate_NRW(datasets_gamma_skewed[[i]],trait_cols = c(3,4,5))
  NRW[[7]][[i]]=calculate_NRW(datasets_gamma_gamma[[i]],trait_cols = c(3,4,5))
  NRW[[8]][[i]]=calculate_NRW(datasets_gamma_weibull[[i]],trait_cols = c(3,4,5))
  NRW[[9]][[i]]=calculate_NRW(datasets_weibull_skewed[[i]],trait_cols = c(3,4,5))
  NRW[[10]][[i]]=calculate_NRW(datasets_weibull_weibull[[i]],trait_cols = c(3,4,5))
  print(i)
}


###############################################



##this will go with the ESPIID and PR and some other scores 

ESP=list(
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
  ESP[[1]][[i]]=calculate_ESP(datasets[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[2]][[i]]=calculate_ESP(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[3]][[i]]=calculate_ESP(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[4]][[i]]=calculate_ESP(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[5]][[i]]=calculate_ESP(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[6]][[i]]=calculate_ESP(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[7]][[i]]=calculate_ESP(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[8]][[i]]=calculate_ESP(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[9]][[i]]=calculate_ESP(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  ESP[[10]][[i]]=calculate_ESP(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  print(i)
}


###############################################



##this is not to be calculated for comparison since it needs definition of control and stress environment

PD=list(
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
  PD[[1]][[i]]=calculate_PD(datasets[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[2]][[i]]=calculate_PD(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[3]][[i]]=calculate_PD(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[4]][[i]]=calculate_PD(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[5]][[i]]=calculate_PD(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[6]][[i]]=calculate_PD(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[7]][[i]]=calculate_PD(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[8]][[i]]=calculate_PD(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[9]][[i]]=calculate_PD(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  PD[[10]][[i]]=calculate_PD(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  print(i)
}

###########################################


##this is not to be calculated for comparison since it needs definition of control and stress environment

FPI=list(
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
  FPI[[1]][[i]]=calculate_FPI(datasets[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[2]][[i]]=calculate_FPI(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[3]][[i]]=calculate_FPI(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[4]][[i]]=calculate_FPI(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[5]][[i]]=calculate_FPI(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[6]][[i]]=calculate_FPI(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[7]][[i]]=calculate_FPI(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[8]][[i]]=calculate_FPI(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[9]][[i]]=calculate_FPI(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  FPI[[10]][[i]]=calculate_FPI(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5),control_env = 1, stress_env = 2)
  print(i)
}



###################################

##this is not to be calculated for comparison since it needs definition of control and stress environment

TSP=list(
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
  TSP[[1]][[i]]=calculate_TPS(datasets[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[2]][[i]]=calculate_TPS(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[3]][[i]]=calculate_TPS(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[4]][[i]]=calculate_TPS(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[5]][[i]]=calculate_TPS(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[6]][[i]]=calculate_TPS(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[7]][[i]]=calculate_TPS(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[8]][[i]]=calculate_TPS(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[9]][[i]]=calculate_TPS(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  TSP[[10]][[i]]=calculate_TPS(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5),native_env = 1, transplanted_env = 2)
  print(i)
}



################################################################################################################################### 4.R
###################################################################################################################################

#calculate_DPI cannot be used since it is used for time resolved data 



CEV=list(
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
  CEV[[1]][[i]]=calculate_CEV(datasets[[i]],trait_cols = c(3,4,5))
  CEV[[2]][[i]]=calculate_CEV(datasets_log_skewed[[i]],trait_cols = c(3,4,5))
  CEV[[3]][[i]]=calculate_CEV(datasets_log_log[[i]],trait_cols = c(3,4,5))
  CEV[[4]][[i]]=calculate_CEV(datasets_log_gamma[[i]],trait_cols = c(3,4,5))
  CEV[[5]][[i]]=calculate_CEV(datasets_log_weibull[[i]],trait_cols = c(3,4,5))
  CEV[[6]][[i]]=calculate_CEV(datasets_gamma_skewed[[i]],trait_cols = c(3,4,5))
  CEV[[7]][[i]]=calculate_CEV(datasets_gamma_gamma[[i]],trait_cols = c(3,4,5))
  CEV[[8]][[i]]=calculate_CEV(datasets_gamma_weibull[[i]],trait_cols = c(3,4,5))
  CEV[[9]][[i]]=calculate_CEV(datasets_weibull_skewed[[i]],trait_cols = c(3,4,5))
  CEV[[10]][[i]]=calculate_CEV(datasets_weibull_weibull[[i]],trait_cols = c(3,4,5))
  print(i)
}


#################################

#calculate_PRI cannot be used since it is used for the comparison of environments 

#################################




PFI=list(
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
  PFI[[1]][[i]]=calculate_PFI(datasets[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[2]][[i]]=calculate_PFI(datasets_log_skewed[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[3]][[i]]=calculate_PFI(datasets_log_log[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[4]][[i]]=calculate_PFI(datasets_log_gamma[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[5]][[i]]=calculate_PFI(datasets_log_weibull[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[6]][[i]]=calculate_PFI(datasets_gamma_skewed[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[7]][[i]]=calculate_PFI(datasets_gamma_gamma[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[8]][[i]]=calculate_PFI(datasets_gamma_weibull[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[9]][[i]]=calculate_PFI(datasets_weibull_skewed[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  PFI[[10]][[i]]=calculate_PFI(datasets_weibull_weibull[[i]],baseline_cols=c(1,1,1),trait_cols = c(3,4,5))
  print(i)
}


#################################

#calculate_SPI cannot be used for cluster analysis since it is used for the comparison of environments

#################################



APC=list(
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
  APC[[1]][[i]]=calculate_APC(datasets[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[2]][[i]]=calculate_APC(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[3]][[i]]=calculate_APC(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[4]][[i]]=calculate_APC(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[5]][[i]]=calculate_APC(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[6]][[i]]=calculate_APC(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[7]][[i]]=calculate_APC(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[8]][[i]]=calculate_APC(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[9]][[i]]=calculate_APC(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  APC[[10]][[i]]=calculate_APC(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  print(i)
}


##################################




SI=list(
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
  SI[[1]][[i]]=calculate_SI(datasets[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[2]][[i]]=calculate_SI(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[3]][[i]]=calculate_SI(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[4]][[i]]=calculate_SI(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[5]][[i]]=calculate_SI(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[6]][[i]]=calculate_SI(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[7]][[i]]=calculate_SI(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[8]][[i]]=calculate_SI(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[9]][[i]]=calculate_SI(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  SI[[10]][[i]]=calculate_SI(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  print(i)
}

##################################




RSI=list(
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
  RSI[[1]][[i]]=calculate_RSI(datasets[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[2]][[i]]=calculate_RSI(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[3]][[i]]=calculate_RSI(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[4]][[i]]=calculate_RSI(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[5]][[i]]=calculate_RSI(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[6]][[i]]=calculate_RSI(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[7]][[i]]=calculate_RSI(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[8]][[i]]=calculate_RSI(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[9]][[i]]=calculate_RSI(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  RSI[[10]][[i]]=calculate_RSI(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  print(i)
}

##################################





EVS=list(
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
  EVS[[1]][[i]]=calculate_EVS(datasets[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[2]][[i]]=calculate_EVS(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[3]][[i]]=calculate_EVS(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[4]][[i]]=calculate_EVS(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[5]][[i]]=calculate_EVS(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[6]][[i]]=calculate_EVS(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[7]][[i]]=calculate_EVS(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[8]][[i]]=calculate_EVS(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[9]][[i]]=calculate_EVS(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5))
  EVS[[10]][[i]]=calculate_EVS(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5))
  print(i)
}


###############################




MVPi=list(
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
  MVPi[[1]][[i]]=calculate_MVPi(datasets[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[2]][[i]]=calculate_MVPi(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[3]][[i]]=calculate_MVPi(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[4]][[i]]=calculate_MVPi(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[5]][[i]]=calculate_MVPi(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[6]][[i]]=calculate_MVPi(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[7]][[i]]=calculate_MVPi(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[8]][[i]]=calculate_MVPi(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[9]][[i]]=calculate_MVPi(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  MVPi[[10]][[i]]=calculate_MVPi(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5),n_axes = 2)
  print(i)
}


###############################

#calculate_SPM cannot be used for clustering analysis because it makes assumptioons about native and non native environments

###############################


Plasticity_Ratio=list(
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
  Plasticity_Ratio[[1]][[i]]=calculate_Plasticity_Ratio(datasets[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[2]][[i]]=calculate_Plasticity_Ratio(datasets_log_skewed[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[3]][[i]]=calculate_Plasticity_Ratio(datasets_log_log[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[4]][[i]]=calculate_Plasticity_Ratio(datasets_log_gamma[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[5]][[i]]=calculate_Plasticity_Ratio(datasets_log_weibull[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[6]][[i]]=calculate_Plasticity_Ratio(datasets_gamma_skewed[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[7]][[i]]=calculate_Plasticity_Ratio(datasets_gamma_gamma[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[8]][[i]]=calculate_Plasticity_Ratio(datasets_gamma_weibull[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[9]][[i]]=calculate_Plasticity_Ratio(datasets_weibull_skewed[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  Plasticity_Ratio[[10]][[i]]=calculate_Plasticity_Ratio(datasets_weibull_weibull[[i]],env_col=1,trait_cols = c(3,4,5),pop_col=1)
  print(i)
}


