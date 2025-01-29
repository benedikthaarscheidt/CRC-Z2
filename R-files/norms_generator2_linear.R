# Function to build covariance matrix
build_cov_matrix = function(VE, VS, rES) {
  matrix(c(VE, rES * sqrt(VE * VS),
           rES * sqrt(VE * VS), VS), nrow = 2)
}


generate_fixed_norm = function(intercept, slope, env) {
  curve = intercept + slope * env
  curve[curve < 0] = 0  
  return(curve)
}


generate_reaction_norm_with_covariance = function(intercept, slope, env, cov_matrix) {
  random_effects = MASS::mvrnorm(3, mu = c(0, 0), Sigma = cov_matrix)
  deviations = t(sapply(1:3, function(i) {
    u0 = random_effects[i, 1]
    u1 = random_effects[i, 2]
    intercept_indiv = intercept + u0
    slope_indiv = slope + u1
    curve = intercept_indiv + slope_indiv * env
    return(c(intercept_indiv, slope_indiv, curve))  
  }))
  return(deviations)
}


set.seed(42)
environmental_factors = seq(1, 10, length.out = 50)  # 50 environments
num_genotypes = 20  # Total number of genotypes
num_individuals_per_genotype = 3  # Number of individuals per genotype

fixed_intercepts = runif(num_genotypes, 5, 30)  # Random intercepts for 20 genotypes
fixed_slopes = runif(num_genotypes, -1, 1)  # Random slopes for 20 genotypes

# Covariance matrix setup
VE = 1
VS = 0.2
rES = 1
cov_matrix = build_cov_matrix(VE, VS, rES)

# Generate individual (replicate) norms for each genotype
individual_norms = do.call(rbind, lapply(1:num_genotypes, function(i) {
  random_effects = MASS::mvrnorm(3, mu = c(0, 0), Sigma = cov_matrix)
  
  
  fixed_effect = generate_fixed_norm(fixed_intercepts[i], fixed_slopes[i], environmental_factors)
  
  # Generate individual replicates
  replicate_data = do.call(rbind, lapply(1:3, function(j) {
    u0 = random_effects[j, 1]
    u1 = random_effects[j, 2]
    base_shift_indiv = fixed_intercepts[i] + u0
    slope_indiv = fixed_slopes[i] + u1
    curve = generate_fixed_norm(base_shift_indiv, slope_indiv, environmental_factors)
    
    data.frame(
      Genotype = i,
      Replicate = j,
      Environment = environmental_factors,
      Trait = curve,
      ReactionNorm = "Linear",
      BaseShift = base_shift_indiv,
      Slope = slope_indiv,
      VarianceOffset = VE,
      VarianceSlope = VS,
      Covariance = rES * sqrt(VE * VS)
    )
  }))
  
  # Combine fixed effect with replicates
  fixed_data = data.frame(
    Genotype = i,
    Replicate = 0,
    Environment = environmental_factors,
    Trait = fixed_effect,
    ReactionNorm = "Linear",
    BaseShift = fixed_intercepts[i],
    Slope = fixed_slopes[i],
    VarianceOffset = NA,
    VarianceSlope = NA,
    Covariance = NA
  )
  
  return(rbind(fixed_data, replicate_data))
}))


create_plot = function(individual_data, genotypes_per_plot = 5, title = "Reaction Norms with Individual Deviations") {
  selected_genotypes = sample(unique(individual_data$Genotype), genotypes_per_plot)
  
  plot_data = individual_data[individual_data$Genotype %in% selected_genotypes, ]
  
  ggplot() +
    geom_line(data = plot_data[plot_data$Replicate == 0, ],
              aes(x = Environment, y = Trait, color = factor(Genotype), group = Genotype),
              size = 1) +
    geom_line(data = plot_data[plot_data$Replicate != 0, ],
              aes(x = Environment, y = Trait, color = factor(Genotype), group = interaction(Genotype, Replicate)),
              linetype = "dashed", size = 0.5) +
    labs(title = title, x = "Environment", y = "Trait Value", color = "Genotype") + # Add a legend title for Genotype
    theme_minimal() +
    theme(legend.position = "right") # Ensure the legend is displayed on the right
}


p = create_plot(individual_norms, genotypes_per_plot = 3, title = "Random Sample of 3 Genotypes")
pdf("linear_reaction_norms2.pdf", width = 8, height = 6)
print(p)
dev.off()
