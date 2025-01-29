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
source("~/CRC 1622 - Z2/R-files/norms_generator2_nonlinear.R")
source("~/CRC 1622 - Z2/R-files/norms_generator2_linear.R")

n_datasets=3
gaussian_data=gaussian_data
sinusoidal_data=sinusoidal_data
wave_data=wave_data
linear_data=individual_norms
################################################################################################################################### 1.R
###################################################################################################################################
test=data.frame(genotype=rep(1,40),replicate=c(rep(0,10),rep(1,10),rep(2,10),rep(3,10)),Trait=c(rep(1,20),rep(3,20)))
CV_t_test=calculate_CVt(test,genotype_col = 1, trait_col = 3)
print(CV_t_test)
CV_t_gaussian = calculate_CVt(gaussian_data, genotype_col = 1, trait_col = 4)
CV_t_gaussian = cbind(CV_t_gaussian, Type = "gaussian")  
print(CV_t_gaussian)
CV_t_sinus = calculate_CVt(sinusoidal_data, genotype_col = 1, trait_col = 4)
CV_t_sinus = cbind(CV_t_sinus, Type = "sinus")

CV_t_wave = calculate_CVt(wave_data, genotype_col = 1, trait_col = 4)
CV_t_wave = cbind(CV_t_wave, Type = "wave")

CV_t_linear = calculate_CVt(linear_data, genotype_col = 1, trait_col = 4)
CV_t_linear = cbind(CV_t_linear, Type = "linear")

CV_t = rbind(CV_t_gaussian, CV_t_sinus, CV_t_wave, CV_t_linear)
####################################

RN_linear = calculate_reaction_norm_slope(linear_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
RN_linear=cbind(RN_linear,Type="linear")

RN_gaussian = calculate_reaction_norm_slope(gaussian_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
RN_gaussian=cbind(RN_gaussian,Type="gaussian")

RN_sinus = calculate_reaction_norm_slope(sinusoidal_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
RN_sinus=cbind(RN_sinus,Type="sinus")

RN_wave = calculate_reaction_norm_slope(wave_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
RN_wave=cbind(RN_wave,Type="wave")

RN=rbind(RN_gaussian,RN_sinus,RN_wave,RN_linear)
#################

RNN_linear=calculate_reaction_norm_nonlinear(linear_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2, degree=2)
RNN_linear=cbind(RNN_linear,Type="linear")

RNN_gaussian=calculate_reaction_norm_nonlinear(gaussian_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2, degree=2)
RNN_gaussian=cbind(RNN_gaussian,Type="gaussian")

RNN_sinus=calculate_reaction_norm_nonlinear(sinusoidal_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2, degree=2)
RNN_sinus=cbind(RNN_sinus,Type="sinus")

RNN_wave=calculate_reaction_norm_nonlinear(wave_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2, degree=2)
RNN_wave=cbind(RNN_wave,Type="wave")
RNN=rbind(RNN_gaussian,RNN_sinus,RNN_wave,RNN_linear)



#####################


Dslope_linear = calculate_D_slope(linear_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
Dslope_linear=cbind(Dslope_linear,Type="linear")

Dslope_gaussian = calculate_D_slope(gaussian_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
Dslope_gaussian=cbind(Dslope_gaussian,Type="gaussian")

Dslope_sinus = calculate_D_slope(sinusoidal_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
Dslope_sinus=cbind(Dslope_sinus,Type="sinus")

Dslope_wave = calculate_D_slope(wave_data,env_col = 3,trait_cols = 4,genotype_col = 1,replicate_col = 2)
Dslope_wave=cbind(Dslope_wave,Type="wave")

Dslope=rbind(Dslope_gaussian,Dslope_sinus,Dslope_wave,Dslope_linear)

# Example test dataset
test_data <- data.frame(
  Genotype = rep(1:2, each = 6),                # Two genotypes
  Replicate = rep(1:3, times = 4),             # Three replicates per genotype
  Environment = rep(c(1, 2, 5, 6, 9, 10), 2),  # Six environment levels
  Trait = c(10, 12, 15, 17, 20, 22,            # Linear increase for Genotype 1
            8, 10, 12, 14, 16, 18)             # Linear increase for Genotype 2
)

# Print test data
print(test_data)


D_slope_test <- calculate_D_slope(
  data = test_data,
  env_col = 3,               # Column for Environment
  trait_cols = 4,            # Column for Trait
  genotype_col = 1,          # Column for Genotype
  replicate_col = 2,         # Column for Replicate
  low_boundary = 2,          # Low environment threshold
  high_boundary = 9          # High environment threshold
)

# Print the output
print(D_slope_test)
