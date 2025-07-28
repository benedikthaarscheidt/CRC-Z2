library(MASS)    # For multivariate normal distributions
library(ggplot2) # For plotting
library(dplyr)   # For data manipulation

set.seed(12)
environmental_factors = seq(1, 10, length.out = 50)  # 50 environments
num_genotypes = 20   # Total number of genotypes
num_replicates = 3   # Number of replicates for each genotype

# NON-LINEAR REACTION NORM FUNCTIONS
generate_gaussian_reaction_norm = function(env, base_shift, slope, amplitude = 4, width = 0.2, center = 5) {
  base_shift + slope * amplitude * exp(-width * (env - center)^2)
}

generate_sinusoidal_reaction_norm = function(env, base_shift, slope, amplitude = 4, frequency = 0.5, phase = 0) {
  base_shift + slope * amplitude * sin(frequency * env + phase)
}

generate_wave_reaction_norm = function(env, base_shift, slope, amplitude = 7, frequency = 1.2) {
  base_shift + slope * amplitude * sin(frequency * (env - 5)) * exp(-0.1 * (env - 5)^2)
}

# Function to build covariance matrix
build_cov_matrix = function(var_offset, var_slope, correlation) {
  cov = correlation * sqrt(var_offset * var_slope)
  matrix(c(var_offset, cov, cov, var_slope), nrow = 2)
}

# New function to generate a fixed non-linear norm with outlier modifications
generate_fixed_nonlin_norm = function(env, base_shift, slope, norm_function, outliers = "None") {
  base_curve = norm_function(env, base_shift, slope)
  
  if (outliers == "None") {
    return(base_curve)
  } else if (outliers == "positive") {
    i = sample(1:length(env), 1)
    cat("Original value at index", i, ":", base_curve[i], "\n")
    base_curve[i] = base_curve[i] * runif(1, 1, 2)
    cat("Modified value at index", i, ":", base_curve[i], "\n")
  } else if (outliers == "negative") {
    i = sample(1:length(env), 1)
    cat("Original value at index", i, ":", base_curve[i], "\n")
    base_curve[i] = base_curve[i] * runif(1, 0, 1)
    cat("Modified value at index", i, ":", base_curve[i], "\n")
  } else if (outliers == "both") {
    i = sample(1:length(env), 1)
    cat("Original value at index", i, ":", base_curve[i], "\n")
    base_curve[i] = base_curve[i] * runif(1, 0, 4)
    cat("Modified value at index", i, ":", base_curve[i], "\n")
  }
  
  return(base_curve)
}

# Modified function to generate data for a single genotype that inherits the mother's outlier
generate_genotype_data = function(base_shift, slope, env, norm_function, num_replicates, 
                                  var_offset, var_slope, correlation, genotype_id, norm_type,
                                  outliers = "None") {
  cov_matrix = build_cov_matrix(var_offset, var_slope, correlation)
  
  # Compute the mother (fixed) curve with outlier modifications
  mother_curve = generate_fixed_nonlin_norm(env, base_shift, slope, norm_function, outliers)
  # Also compute the base curve (without outliers) for calculating deviations
  base_curve = norm_function(env, base_shift, slope)
  
  mother_data = data.frame(
    Genotype = genotype_id,
    Replicate = 0,
    Environment = env,
    Trait = mother_curve,
    ReactionNorm = norm_type,
    BaseShift = base_shift,
    Slope = slope,
    VarianceOffset = NA,
    VarianceSlope = NA,
    Covariance = NA
  )
  
  # Generate replicates that inherit the mother's outlier pattern
  replicate_data = do.call(rbind, lapply(1:num_replicates, function(rep) {
    perturbations = MASS::mvrnorm(1, mu = c(0, 0), Sigma = cov_matrix)
    offset_perturbation = perturbations[1]
    slope_perturbation = perturbations[2]
    
    # Compute the perturbed curve without outlier modifications
    perturbed_curve = norm_function(env, base_shift + offset_perturbation, slope + slope_perturbation)
    # Calculate the deviation relative to the base (non-outlier) curve
    deviation = perturbed_curve - base_curve
    # Add the deviation to the mother curve so that the outlier is inherited
    replicate_curve = mother_curve + deviation
    
    data.frame(
      Genotype = genotype_id,
      Replicate = rep,
      Environment = env,
      Trait = replicate_curve,
      ReactionNorm = norm_type,
      BaseShift = base_shift + offset_perturbation,
      Slope = slope + slope_perturbation,
      VarianceOffset = var_offset,
      VarianceSlope = var_slope,
      Covariance = correlation * sqrt(var_offset * var_slope)
    )
  }))
  
  return(rbind(mother_data, replicate_data))
}

# Modified function to generate data for all genotypes for a specific norm type
generate_reaction_norm_data = function(norm_function, norm_type, base_shifts, slopes, 
                                       env, num_replicates, var_offset, var_slope, correlation,
                                       outliers = "None") {
  do.call(rbind, lapply(1:num_genotypes, function(genotype_id) {
    base_shift = base_shifts[genotype_id]
    slope = slopes[genotype_id]
    generate_genotype_data(base_shift, slope, env, norm_function, num_replicates, 
                           var_offset, var_slope, correlation, genotype_id, norm_type, outliers)
  }))
}

base_shifts  <- runif(num_genotypes, 10, 30)

# Randomize both magnitude and sign of each genotype’s slope
slope_mags  <- runif(num_genotypes, 0.5, 1.5)
slope_signs <- sample(c(-1, 1), num_genotypes, replace = TRUE)
slopes      <- slope_mags * slope_signs

# Now generate the data as before…
gaussian_data = generate_reaction_norm_data(
  generate_gaussian_reaction_norm, "Gaussian", base_shifts, slopes, 
  environmental_factors, num_replicates, var_offset = 0.5, var_slope = 0.1, correlation = -0.5,
  outliers = "none"
)
sinusoidal_data = generate_reaction_norm_data(
  generate_sinusoidal_reaction_norm, "Sinusoidal", base_shifts, slopes, 
  environmental_factors, num_replicates, var_offset = 0.3, var_slope = 0.05, correlation = -0.3,
  outliers = "none"
)
wave_data = generate_reaction_norm_data(
  generate_wave_reaction_norm, "Wave", base_shifts, slopes, 
  environmental_factors, num_replicates, var_offset = 0.7, var_slope = 0.15, correlation = -0.7,
  outliers = "none"
)
all_data = rbind(gaussian_data, sinusoidal_data, wave_data)

# Plotting functions (unchanged)
create_plot = function(data, genotypes_per_plot = 3, show_replicates = TRUE) {
  set.seed(44)
  selected_genotypes = sample(unique(data$Genotype), genotypes_per_plot, replace = TRUE)
  plot_data = data[data$Genotype %in% selected_genotypes, ]
  
  ggplot(plot_data, aes(x = Environment, y = Trait, color = factor(Genotype))) +
    geom_line(aes(group = interaction(Genotype, Replicate), linetype = as.factor(Replicate))) +
    facet_wrap(~ ReactionNorm) +
    labs(
      title = "Reaction Norms",
      x = "Environment",
      y = "Trait Value",
      color = "Genotype",
      linetype = "Replicate"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
}

create_plot2 = function(gaussian_data, sinusoidal_data, wave_data, show_replicates = TRUE) {
  set.seed(44)
  
  gaussian_sample = gaussian_data %>% 
    filter(Genotype == sample(unique(Genotype), 1))
  
  sinusoidal_sample = sinusoidal_data %>% 
    filter(Genotype == sample(unique(Genotype), 1))
  
  wave_sample = wave_data %>% 
    filter(Genotype == sample(unique(Genotype), 1))
  
  sampled_data = rbind(gaussian_sample, sinusoidal_sample, wave_sample)
  
  fixed_data = sampled_data[sampled_data$Replicate == 0, ]
  replicate_data = sampled_data[sampled_data$Replicate != 0, ]
  
  p = ggplot(fixed_data, aes(x = Environment, y = Trait, color = ReactionNorm, group = interaction(ReactionNorm, Genotype))) +
    geom_line(size = 1, linetype = "solid") +
    labs(
      title = "Randomly Sampled Reaction Norms (One Per Type)",
      x = "Environment",
      y = "Trait Value",
      color = "Reaction Norm Type"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  if (show_replicates) {
    p = p + geom_line(
      data = replicate_data,
      aes(group = interaction(ReactionNorm, Genotype, Replicate)),
      linetype = "dashed", size = 0.5
    )
  }
  
  return(p)
}

pdf("~/CRC_1644_Z2/plots/reaction_norms_plots_combined.pdf", width = 12, height = 8)
print(create_plot(gaussian_data))   # Plot for Gaussian norms
print(create_plot(sinusoidal_data))   # Plot for Sinusoidal norms
print(create_plot(wave_data))         # Plot for Wave norms
print(create_plot2(gaussian_data, sinusoidal_data, wave_data)) # Combined plot with 3 randomly sampled norms
dev.off()

ensure_dir_exists <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

ensure_dir_exists("~/CRC_1644_Z2/synthetic_data/fixed_full")
write.csv(gaussian_data, "~/CRC_1644_Z2/synthetic_data/fixed_full/gaussian_data.csv", row.names = FALSE)
write.csv(sinusoidal_data, "~/CRC_1644_Z2/synthetic_data/fixed_full/sinusoidal_data.csv", row.names = FALSE)
write.csv(wave_data, "~/CRC_1644_Z2/synthetic_data/fixed_full/wave_data.csv", row.names = FALSE)
write.csv(all_data, "~/CRC_1644_Z2/synthetic_data/fixed_full/all_combined_data.csv", row.names = FALSE)
