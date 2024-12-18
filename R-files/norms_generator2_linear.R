library(MASS)  
library(ggplot2) 
library(reshape2) 


build_cov_matrix = function(VE, VS, rES) {
  matrix(c(VE, rES * sqrt(VE * VS), 
           rES * sqrt(VE * VS), VS), nrow = 2)
}


generate_fixed_norm = function(intercept, slope, env) {
  curve = intercept + slope * env
  curve[curve < 0] = 0  # Ensure no negative values
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
    return(curve)
  }))
  return(deviations)
}


set.seed(42)
environmental_factors = seq(1, 10, length.out = 50)  # 50 environments
num_genotypes = 20  # Total number of genotypes
num_individuals_per_genotype = 3  # Number of individuals per genotype 
sampled_genotypes = sample(seq_len(num_genotypes), num_genotypes)  


fixed_intercepts = runif(num_genotypes, 5, 30)  # Random intercepts for 20 genotypes
fixed_slopes = runif(num_genotypes, -1, 1)  # Random slopes for 20 genotypes

fixed_norms = data.frame(
  t(sapply(seq_len(num_genotypes), function(i) 
    generate_fixed_norm(fixed_intercepts[i], fixed_slopes[i], environmental_factors)
  ))
)
colnames(fixed_norms) = environmental_factors
fixed_norms$Genotype = 1:num_genotypes


VE = 1  
VS = 0.2 
rES = 1   
cov_matrix = build_cov_matrix(VE, VS, rES)

individual_norms = do.call(rbind, lapply(1:num_genotypes, function(i) {
  
  indiv_devs = generate_reaction_norm_with_covariance(fixed_intercepts[i], fixed_slopes[i], environmental_factors, cov_matrix)
  
  
  fixed_effect = generate_fixed_norm(fixed_intercepts[i], fixed_slopes[i], environmental_factors)
  
  
  indiv_devs_combined = rbind(fixed_effect, indiv_devs)  # 4 x 50 (1 fixed + 3 deviates)
  
  indiv_devs_long = data.frame(
    Trait = as.vector(t(indiv_devs_combined)),  # Flatten the 4x50 matrix into 200-length vector
    Genotype = i,  # Repeat the genotype ID for 200 rows (4 curves x 50 points)
    Individual = rep(0:num_individuals_per_genotype, each = length(environmental_factors)),  # 0 (fixed) and 1, 2, 3 for deviated
    Environment = rep(environmental_factors, times = num_individuals_per_genotype + 1)  # Repeat for each individual (4 total)
  )
  
  return(indiv_devs_long)
}))


create_plot = function(individual_data, genotypes_per_plot = 5, title = "Reaction Norms with Individual Deviations") {
  # Select a subset of genotypes to plot
  selected_genotypes = sample(unique(individual_data$Genotype), genotypes_per_plot)
  
  # Filter the data for the selected genotypes
  plot_data = individual_data[individual_data$Genotype %in% selected_genotypes,]
  
  ggplot() +
    geom_line(data = plot_data[plot_data$Individual == 0,], 
              aes(x = Environment, y = Trait, color = factor(Genotype), group = Genotype), 
              size = 1) +  # Solid line for fixed effect (Individual 0)
    geom_line(data = plot_data[plot_data$Individual != 0,], 
              aes(x = Environment, y = Trait, color = factor(Genotype), group = interaction(Genotype, Individual)), 
              linetype = "dashed", size = 0.5) +  # Dashed lines for deviated norms (Individuals 1, 2, 3)
    labs(title = title, x = "Environment", y = "Trait Value") +
    theme_minimal() +
    theme(legend.position = "none")
}


p = create_plot(individual_norms, genotypes_per_plot = 3, title = "Random Sample of 5 Genotypes")
pdf("linear_reaction_norms2.pdf", width = 8, height = 6)
print(p)
dev.off()
