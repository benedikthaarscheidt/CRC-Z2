#source("~/CRC 1622 - Z2/R-files/1.R")
##coefficient-of variation total (calculate_CVt) - tested,
##slope of norm reaction (calculate_reaction_norm_slope) - tested,
##slope of plastic response (D) (calculate_D_slope)- tested,
##response coefficient (RC) (calculate_RC)- tested,
##Standard deviation of means (CVm) (calculate_CVm)- tested,
##Standard deviation of medians (CVmd)(calculate_CVmd)- tested,
##Grand plasticity (calculate_GPi)- tested,
##Phenotypic Plasticity Index (calculate_PPF)- tested,
##Phenotypic Plasticity Index (calculate_Phenotypic_Plasticity_Index)- tested,
##PImd (calculate_PImd)- tested,
##PILSM (calculate_PILSM)- tested,
##RTR (calculate_RTR)- tested,
##PIR (calculate_PIR) - tested
#source("~/CRC 1622 - Z2/R-files/2.R")
##RDPI	(rdpi_calculation) - tested,
##RDPIs (rdpi_mean_calculation) - tested,
##ESPI (calculate_ESPI) - tested,
##ESPIid (espiid_calculation) - tested,
##evwpi_calculation (idea from Benedikt)
#source("~/CRC 1622 - Z2/R-files/3.R")
##Phenotypic Stability Index (calculate_PSI),
##Relative Plasticity Index (calculate_RPI) - tested,
##Plasticity Quotient (calculate_PQ) - tested,
##Phenotypic Range (PR) (calculate_PR) - tested,
##Norm of reaction width (calculate_NRW) - tested,
##Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
##Calculate Plasticity Differential (PD) (calculate_PD) - tested,
##Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
##Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,
#source("~/CRC 1622 - Z2/R-files/4.R")
##Calculate Developmental Plasticity Index (DPI)(calculate_DPI) - tested,
##Calculate Coefficient of Environmental Variation (CEV)(calculate_CEV) - tested,
##Calculate Plasticity Response Index (PRI)(calculate_PRI) - tested,
##Calculate Phenotypic Flexibility Index (PFI)(calculate_PFI) - tested,
##Calculate Standardized Plasticity Index (SPI)(calculate_SPI) - tested,
##Calculate Absolute Plasticity Coefficient (APC)(calculate_APC) - tested,
##Calculate Stability Index (SI)(calculate_SI) - tested,
##Calculate Relative Stability Index (RSI)(calculate_RSI),
##Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
##Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
##Calculate Multivariate Plasticity Index (MVPi)(calculate_MVPi) NEEDS TO BE TESTED WITH THE REQUESTED DATASET FROM PROF. BARBOSA,
##Calculate Standardized Plasticity Metric (SPM)(calculate_SPM) - tested,
##Calculate SSpop/SStotal Plasticity Ratio(calculate_Plasticity_Ratio) - tested

source("~/CRC 1622 - Z2/R-files/norms_generator.R")

library(reshape2)

# Function to generate synthetic data with three replicates per genotype and environmental factor for multiple traits
generate_synthetic_data = function(baseline_list, noise_interval, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Number of genotypes and environments
  n_genotypes = nrow(baseline_list[[1]])
  n_environments = ncol(baseline_list[[1]])
  n_replicates = 3
  n_traits = length(baseline_list)
  
  # Initialize an empty data frame to store the results
  synthetic_data = data.frame(
    Genotype = rep(1:n_genotypes, each = n_environments * n_replicates),
    Environment = rep(1:n_environments, each = n_replicates, times = n_genotypes)
  )
  
  # Add columns for each trait
  for (t in 1:n_traits) {
    trait_name = paste("Trait", t, sep = "_")
    synthetic_data[[trait_name]] = NA
    
    baseline_values = baseline_list[[t]]
    
    # Generate synthetic data with noise for each trait
    for (genotype in 1:n_genotypes) {
      for (env in 1:n_environments) {
        mean_value = baseline_values[genotype, env]
        for (rep in 1:n_replicates) {
          row_index = (genotype - 1) * n_environments * n_replicates + (env - 1) * n_replicates + rep
          synthetic_data[row_index, trait_name] = rnorm(1, mean = mean_value, sd = runif(1, noise_interval[1], noise_interval[2]))
        }
      }
    }
  }
  
  return(synthetic_data)
}




# Prepare baseline values by removing the 'Row' column from each dataset
baseline_RN_linear_NC = as.matrix(RN_linear_NC[, 1:(ncol(RN_linear_NC) - 1)])
baseline_RN_linear_C = as.matrix(RN_linear_C[, 1:(ncol(RN_linear_C) - 1)])
RN_nonlinear_NC_smoothed = as.matrix(RN_nonlinear_NC_smoothed[, 1:(ncol(RN_nonlinear_NC_smoothed) - 1)])
RN_nonlinear_C_smoothed = as.matrix(RN_nonlinear_C_smoothed[, 1:(ncol(RN_nonlinear_C_smoothed) - 1)])

# Combine all baselines into a list for multiple traits
baseline_list = list(baseline_RN_linear_NC, baseline_RN_linear_C, RN_nonlinear_NC_smoothed, RN_nonlinear_C_smoothed)
# Specify noise interval 
noise_interval = c(0.3,0.4)
# Generate synthetic data for all traits
synthetic_data = generate_synthetic_data(baseline_list, noise_interval, seed = 123)
# Convert synthetic data to long format for plotting
synthetic_data_long = melt(synthetic_data, id.vars = c("Genotype", "Environment"), variable.name = "Trait", value.name = "Value")
synthetic_data_long$Genotype = as.factor(synthetic_data_long$Genotype)
synthetic_data_long$Environment = as.numeric(as.character(synthetic_data_long$Environment))
# Prepare individual synthetic data for each trait
synthetic_data_linear_NC = subset(synthetic_data_long, Trait == "Trait_1")
synthetic_data_linear_C = subset(synthetic_data_long, Trait == "Trait_2")
synthetic_data_nonlinear_NC = subset(synthetic_data_long, Trait == "Trait_3")
synthetic_data_nonlinear_C = subset(synthetic_data_long, Trait == "Trait_4")
# Ensure that only the sampled genotypes are used in each plot
synthetic_data_linear_NC = synthetic_data_linear_NC[synthetic_data_linear_NC$Genotype %in% genotypes, ]
synthetic_data_linear_C = synthetic_data_linear_C[synthetic_data_linear_C$Genotype %in% genotypes, ]
synthetic_data_nonlinear_NC = synthetic_data_nonlinear_NC[synthetic_data_nonlinear_NC$Genotype %in% genotypes, ]
synthetic_data_nonlinear_C = synthetic_data_nonlinear_C[synthetic_data_nonlinear_C$Genotype %in% genotypes, ]


library(dplyr)

# Function to create a data frame with min and max for each Env and Genotype
create_min_max_df <- function(data) {
  data %>%
    group_by(Environment, Genotype) %>%
    summarise(min_value = min(Value), max_value = max(Value)) %>%
    ungroup()
}

# Create data frames for each dataset with min and max regions
min_max_linear_NC <- create_min_max_df(synthetic_data_linear_NC)
min_max_linear_C <- create_min_max_df(synthetic_data_linear_C)
min_max_nonlinear_NC <- create_min_max_df(synthetic_data_nonlinear_NC)
min_max_nonlinear_C <- create_min_max_df(synthetic_data_nonlinear_C)



pdf("~/CRC 1622 - Z2/R-files/reaction_norms.pdf", width = 10, height = 6)

# Function to add a version without numeric scales
remove_numeric_scales <- function(plot) {
  plot + theme(axis.text = element_blank(), axis.ticks = element_blank())
}

# Linear Non-Crossing with Min-Max Regions
plot1_gradient = ggplot() +
  geom_line(data = RN_linear_NC_long_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_ribbon(data = min_max_linear_NC, aes(x = Environment, ymin = min_value, ymax = max_value, fill = Genotype), alpha = 0.3) +
  labs(title = "Linear Non-Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  guides(fill = "none")
print(plot1_gradient)
print(remove_numeric_scales(plot1_gradient))

# Linear Non-Crossing with Min-Max Error Bars
plot1_error = ggplot() +
  geom_line(data = RN_linear_NC_long_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_errorbar(data = min_max_linear_NC, aes(x = Environment, ymin = min_value, ymax = max_value, color = Genotype, group = Genotype), width = 0.2) +
  labs(title = "Linear Non-Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank())
print(plot1_error)
print(remove_numeric_scales(plot1_error))

# Linear Crossing with Min-Max Regions
plot2_gradient = ggplot() +
  geom_line(data = RN_linear_C_long_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_ribbon(data = min_max_linear_C, aes(x = Environment, ymin = min_value, ymax = max_value, fill = Genotype), alpha = 0.3) +
  labs(title = "Linear Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  guides(fill = "none")
print(plot2_gradient)
print(remove_numeric_scales(plot2_gradient))

# Linear Crossing with Min-Max Error Bars
plot2_error = ggplot() +
  geom_line(data = RN_linear_C_long_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_errorbar(data = min_max_linear_C, aes(x = Environment, ymin = min_value, ymax = max_value, color = Genotype, group = Genotype), width = 0.2) +
  labs(title = "Linear Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank())
print(plot2_error)
print(remove_numeric_scales(plot2_error))

# Non-Linear Non-Crossing with Min-Max Regions
plot3_gradient = ggplot() +
  geom_line(data = RN_nonlinear_NC_smoothed_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_ribbon(data = min_max_nonlinear_NC, aes(x = Environment, ymin = min_value, ymax = max_value, fill = Genotype), alpha = 0.3) +
  labs(title = "Non-Linear Non-Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  guides(fill = "none")
print(plot3_gradient)
print(remove_numeric_scales(plot3_gradient))

# Non-Linear Non-Crossing with Min-Max Error Bars
plot3_error = ggplot() +
  geom_line(data = RN_nonlinear_NC_smoothed_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_errorbar(data = min_max_nonlinear_NC, aes(x = Environment, ymin = min_value, ymax = max_value, color = Genotype, group = Genotype), width = 0.2) +
  labs(title = "Non-Linear Non-Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank())
print(plot3_error)
print(remove_numeric_scales(plot3_error))

# Non-Linear Crossing with Min-Max Regions
plot4_gradient = ggplot() +
  geom_line(data = RN_nonlinear_C_smoothed_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_ribbon(data = min_max_nonlinear_C, aes(x = Environment, ymin = min_value, ymax = max_value, fill = Genotype), alpha = 0.3) +
  labs(title = "Non-Linear Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  guides(fill = "none")
print(plot4_gradient)
print(remove_numeric_scales(plot4_gradient))

# Non-Linear Crossing with Min-Max Error Bars
plot4_error = ggplot() +
  geom_line(data = RN_nonlinear_C_smoothed_subset, aes(x = Env, y = Value, color = Row, group = Row), size = 1) +
  geom_errorbar(data = min_max_nonlinear_C, aes(x = Environment, ymin = min_value, ymax = max_value, color = Genotype, group = Genotype), width = 0.2) +
  labs(title = "Non-Linear Crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Phenotype Value") +
  theme_minimal() +
  theme(legend.title = element_blank())
print(plot4_error)
print(remove_numeric_scales(plot4_error))

# Close the PDF device
dev.off()









