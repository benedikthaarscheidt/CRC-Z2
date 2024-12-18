# Load required libraries
library(MASS)  # For multivariate normal distributions
library(ggplot2)  # For plotting
library(grid)  # For adding axis arrows

# Set environmental factors and number of genotypes
set.seed(13)
environmental_factors = seq(1, 10, length.out = 50)  # 50 environments
num_genotypes = 50  # Total number of genotypes
num_replicates = 3  # Number of replicates for each genotype

# Function to generate Gaussian-like reaction norm
generate_gaussian_reaction_norm = function(env, base_shift, slope, amplitude = 4, width = 0.2, center = 5) {
  curve = base_shift + slope * amplitude * exp(-width * (env - center)^2)
  return(curve)
}

# Function to generate sinusoidal reaction norm
generate_sinusoidal_reaction_norm = function(env, base_shift, slope, amplitude = 4, frequency = 0.5, phase = 0) {
  curve = base_shift + slope * amplitude * sin(frequency * env + phase)
  return(curve)
}

# Function to generate wave-like reaction norm
generate_wave_reaction_norm = function(env, base_shift, slope, amplitude = 7, frequency = 1.2) {
  curve = base_shift + slope * amplitude * sin(frequency * (env - 5)) * exp(-0.1 * (env - 5)^2)
  return(curve) 
}

# Function to build covariance matrix for random effects
build_cov_matrix = function(var_offset, var_slope, correlation) {
  cov = correlation * sqrt(var_offset * var_slope)
  return(matrix(c(var_offset, cov, cov, var_slope), nrow = 2))
}

# Function to generate replicates using correlated random effects for offset and slope
generate_replicates = function(base_shift, slope, env, norm_function, num_replicates = 3, var_offset = 0.5, var_slope = 0.1, correlation = -1) {
  cov_matrix = build_cov_matrix(var_offset, var_slope, correlation)
  
  replicate_data = do.call(rbind, lapply(1:num_replicates, function(rep) {
    perturbations = MASS::mvrnorm(1, mu = c(0, 0), Sigma = cov_matrix)
    offset_perturbation = perturbations[1]
    slope_perturbation = perturbations[2]
    perturbed_curve = norm_function(env, base_shift + offset_perturbation, slope + slope_perturbation)
    data.frame(
      Environment = env,
      Trait = perturbed_curve,
      Replicate = rep
    )
  }))
  
  return(replicate_data)
}

# Generate fixed effects for all genotypes
base_shifts = runif(num_genotypes, 10, 30)  # Random fixed offset for each genotype
slopes = runif(num_genotypes, 0.5, 1.5)  # Random fixed slope for each genotype
set.seed(123) 
reaction_norm_types = sample(c("Gaussian", "Sinusoidal", "Wave"), num_genotypes, replace = TRUE)

# Function to generate reaction norms for a single genotype
generate_single_reaction_norm = function(genotype_id, norm_type) {
  env = environmental_factors
  base_shift = base_shifts[genotype_id]
  slope = slopes[genotype_id]
  
  # Generate reaction norms based on the assigned type
  if (norm_type == "Gaussian") {
    fixed_curve = generate_gaussian_reaction_norm(env, base_shift, slope)
    replicate_data = generate_replicates(base_shift, slope, env, generate_gaussian_reaction_norm)
  } else if (norm_type == "Sinusoidal") {
    fixed_curve = generate_sinusoidal_reaction_norm(env, base_shift, slope)
    replicate_data = generate_replicates(base_shift, slope, env, generate_sinusoidal_reaction_norm)
  } else {
    fixed_curve = generate_wave_reaction_norm(env, base_shift, slope)
    replicate_data = generate_replicates(base_shift, slope, env, generate_wave_reaction_norm)
  }
  
  # Fixed effect data
  fixed_data = data.frame(
    Environment = env,
    Trait = fixed_curve,
    Genotype = genotype_id,
    ReactionNorm = norm_type
  )
  
  # Replicate data
  replicate_data$Genotype = genotype_id
  replicate_data$ReactionNorm = norm_type
  
  return(list(fixed = fixed_data, replicates = replicate_data))
}


all_reaction_norms = lapply(1:num_genotypes, function(i) {
  generate_single_reaction_norm(i, reaction_norm_types[i])
})

# Combine fixed and replicate data
fixed_data = do.call(rbind, lapply(all_reaction_norms, function(x) x$fixed))
replicate_data = do.call(rbind, lapply(all_reaction_norms, function(x) x$replicates))



create_plot = function(fixed_data, replicate_data, genotypes_per_plot = 3, var_offset = 0.5, var_slope = 0.1, correlation = -1,  show_replicates = TRUE) {
  
  selected_genotypes = lapply(unique(fixed_data$ReactionNorm), function(norm) {
    norm_data = fixed_data[fixed_data$ReactionNorm == norm, ]
    sample(unique(norm_data$Genotype), genotypes_per_plot, replace = FALSE)
  })
  selected_genotypes = unlist(selected_genotypes)
  
  # Filter fixed and replicate data for selected genotypes
  plot_fixed_data = fixed_data[fixed_data$Genotype %in% selected_genotypes, ]
  plot_replicate_data = replicate_data[replicate_data$Genotype %in% selected_genotypes, ]
  
  # Create the plot title
  plot_title = paste0(
    "Reaction Norms of ", genotypes_per_plot, " Genotypes per Reaction Norm Type\n",
    "Offset Variance: ", var_offset, 
    ", Slope Variance: ", var_slope, 
    ", Correlation: ", correlation
  )
  
  # Base plot with fixed effects 
  p = ggplot() +
    geom_line(data = plot_fixed_data, 
              aes(x = Environment, y = Trait, color = factor(Genotype), 
                  group = interaction(Genotype, ReactionNorm)), 
              size = 1) +
    facet_wrap(~ ReactionNorm) +  # Separate plots for each ReactionNorm type
    labs(title = plot_title, x = "Environment", y = "Trait") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Add replicates as dashed lines if show_replicates = TRUE
  if (show_replicates) {
    p = p + geom_line(data = plot_replicate_data, 
                       aes(x = Environment, y = Trait, color = factor(Genotype), 
                           group = interaction(Genotype, ReactionNorm, Replicate)), 
                       linetype = "dashed", size = 0.5)
  }
  
  return(p)
}

create_plot2 = function(fixed_data, replicate_data, genotypes_per_plot = 3, 
                         var_offset = 0.5, var_slope = 0.1, correlation = -1, 
                         show_replicates = TRUE) {
  # Randomly sample the specified number of genotypes
  selected_genotypes = sample(unique(fixed_data$Genotype), genotypes_per_plot, replace = FALSE)
  
  # Filter the fixed and replicate data
  plot_fixed = fixed_data[fixed_data$Genotype %in% selected_genotypes, ]
  plot_replicates = replicate_data[replicate_data$Genotype %in% selected_genotypes, ]
  
  # Create the plot title
  plot_title = paste0(
    "Reaction Norms of ", genotypes_per_plot, " Genotypes per ReactionNorm Type\n",
    "Offset Variance: ", var_offset, 
    ", Slope Variance: ", var_slope, 
    ", Correlation: ", correlation
  )
  
  # Plot fixed effects as solid lines
  p = ggplot(plot_fixed, aes(x = Environment, y = Trait, color = factor(Genotype))) +
    geom_line(aes(group = Genotype), size = 1, linetype = "solid") +
    labs(title = plot_title, x = "Environment", y = "Trait", color = "Genotype") +  # Use plot_title here
    theme_minimal()
  
  # Add replicates if enabled
  if (show_replicates) {
    p = p + geom_line(data = plot_replicates,
                       aes(group = interaction(Genotype, Replicate), 
                           color = factor(Genotype)), 
                       size = 0.5, linetype = "dashed")
  }
  
  return(p)
}



# Generate and save plots
pdf("reaction_norms_separated_and_mixed.pdf", width = 12, height = 8)
print(create_plot(fixed_data, replicate_data))
print(create_plot2(fixed_data, replicate_data))
dev.off()
