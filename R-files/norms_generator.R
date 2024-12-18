# Load required libraries
library(ggplot2)
library(grid)

# Set environmental factors and number of genotypes
set.seed(13)
environmental_factors = seq(1, 10, length.out = 50)
num_genotypes = 50
num_genotypes_per_plot = 3

# Select genotypes for sampling
sampled_genotypes = sample(seq_len(num_genotypes)[1:10], num_genotypes_per_plot)

# Function to generate Gaussian-like reaction norm
generate_gaussian_reaction_norm = function(env, base_shift, amplitude = 4, width = 0.2, center = (max(environmental_factors)/2), slope_variation = 0.6) {
  slope = runif(1, 1 - slope_variation, 1 + slope_variation)
  curve = base_shift + slope * amplitude * exp(-width * (env - center)^2)
  return(curve)
}

# Function to generate sinusoidal reaction norm
generate_sinusoidal_reaction_norm = function(env, base_shift, amplitude = 4, frequency = 0.5, phase = 0, slope_variation = 0.4) {
  slope = runif(1, 1 - slope_variation, 1 + slope_variation)
  curve = base_shift + slope * amplitude * sin(frequency * env + phase)
  return(curve)
}

# Function to generate wave-like reaction norm
generate_wave_reaction_norm = function(env, base_shift, amplitude = 7, frequency = 1.2, slope_variation = 0.5) {
  slope = runif(1, 1 - slope_variation, 1 + slope_variation)
  curve = base_shift + slope * amplitude * sin(frequency * (env - 5)) * exp(-0.1 * (env - 5)^2)
  return(curve)
}

# Function to generate replicates with noise
generate_replicates = function(reaction_norms, num_replicates = 3, noise_interval = c(0.3, 0.4)) {
  n_genotypes = nrow(reaction_norms)
  n_environments = ncol(reaction_norms)
  replicates = data.frame()
  for (genotype in 1:n_genotypes) {
    for (env in seq_len(n_environments)) {  # Generate replicates for every environment position
      mean_value = reaction_norms[genotype, env]
      replicate_values = data.frame()
      for (rep in 1:num_replicates) {
        replicate_value = rnorm(1, mean = mean_value, sd = runif(1, noise_interval[1], noise_interval[2]))
        replicate_values = rbind(replicate_values, data.frame(Genotype = genotype, Environment = environmental_factors[env], Replicate = rep, Value = replicate_value))
      }
      replicates = rbind(replicates, replicate_values)
    }
  }
  return(replicates)
}

# Define base shifts for non-linear reaction norms
base_shifts = seq(10, 300, length.out = num_genotypes)

# Generate Gaussian-like reaction norms
reaction_norms_gaussian = data.frame(
  t(sapply(seq_along(base_shifts), function(i) generate_gaussian_reaction_norm(environmental_factors, base_shift = base_shifts[i])))
)
colnames(reaction_norms_gaussian) = environmental_factors
rownames(reaction_norms_gaussian) = paste0("Gaussian_", seq_along(base_shifts))

# Generate Sinusoidal reaction norms
reaction_norms_sinusoidal = data.frame(
  t(sapply(seq_along(base_shifts), function(i) generate_sinusoidal_reaction_norm(environmental_factors, base_shift = base_shifts[i])))
)
colnames(reaction_norms_sinusoidal) = environmental_factors
rownames(reaction_norms_sinusoidal) = paste0("Sinusoidal_", seq_along(base_shifts))

# Generate Wave-like reaction norms
reaction_norms_wave = data.frame(
  t(sapply(seq_along(base_shifts), function(i) generate_wave_reaction_norm(environmental_factors, base_shift = base_shifts[i])))
)
colnames(reaction_norms_wave) = environmental_factors
rownames(reaction_norms_wave) = paste0("Wave_", seq_along(base_shifts))

# Generate replicates for non-linear reaction norms
replicates_gaussian = generate_replicates(reaction_norms_gaussian)
replicates_sinusoidal = generate_replicates(reaction_norms_sinusoidal)
replicates_wave = generate_replicates(reaction_norms_wave)

# Filter replicates for sampled genotypes only and select every 5th environment for plotting
replicates_gaussian_selected = replicates_gaussian[replicates_gaussian$Genotype %in% sampled_genotypes & replicates_gaussian$Environment %in% environmental_factors[seq(1, length(environmental_factors), by = 5)], ]
replicates_sinusoidal_selected = replicates_sinusoidal[replicates_sinusoidal$Genotype %in% sampled_genotypes & replicates_sinusoidal$Environment %in% environmental_factors[seq(1, length(environmental_factors), by = 5)], ]
replicates_wave_selected = replicates_wave[replicates_wave$Genotype %in% sampled_genotypes & replicates_wave$Environment %in% environmental_factors[seq(1, length(environmental_factors), by = 5)], ]

# Reshape data for plotting
reshape_for_plot = function(data, type) {
  df = t(data)
  df = as.data.frame(df)
  colnames(df) = paste0(type, "_", seq_len(ncol(df)))
  df$Environment = environmental_factors
  return(df)
}

# Reshape sampled curves for plotting
gaussian_sampled = reaction_norms_gaussian[sampled_genotypes, ]
sinusoidal_sampled = reaction_norms_sinusoidal[sampled_genotypes, ]
wave_sampled = reaction_norms_wave[sampled_genotypes, ]

gaussian_reshaped = reshape_for_plot(gaussian_sampled, "Gaussian")
sinusoidal_reshaped = reshape_for_plot(sinusoidal_sampled, "Sinusoidal")
wave_reshaped = reshape_for_plot(wave_sampled, "Wave")

# Function to add axis lines with arrows
add_axis_arrows = function() {
  list(
    annotation_custom(
      grid::linesGrob(
        x = c(0, 1), y = c(0, 0), 
        arrow = grid::arrow(length = unit(0.2, "cm")), 
        gp = grid::gpar(col = "black")
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ),
    annotation_custom(
      grid::linesGrob(
        x = c(0, 0), y = c(0, 1), 
        arrow = grid::arrow(length = unit(0.2, "cm")), 
        gp = grid::gpar(col = "black")
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  )
}

# Function to create individual plots with and without axis labels
create_plot = function(data, title, line_colors, replicates_data) {
  plot_with_labels = ggplot(data, aes(x = Environment))
  
  for (i in 1:(ncol(data) - 1)) {
    plot_with_labels = plot_with_labels + geom_line(aes_string(y = colnames(data)[i]), color = line_colors[i])
  }
  
  plot_with_labels = plot_with_labels +
    geom_point(data = replicates_data, aes(x = Environment, y = Value, color = as.factor(Genotype)), size = 2, shape = 21) +
    scale_color_manual(values = setNames(line_colors, as.character(sampled_genotypes))) +
    labs(title = title) +
    ylab("Trait") +
    theme_minimal() +
    theme(legend.position = "none") +
    add_axis_arrows()
  
  plot_without_labels = plot_with_labels + 
    ylab("Trait") +
    xlab("Environment") +
    theme(axis.text = element_blank())
  
  return(list(with_labels = plot_with_labels, without_labels = plot_without_labels))
}

# Creating the individual plots with replicates for sampled genotypes
line_colors = scales::hue_pal()(num_genotypes_per_plot)  # Generate distinct colors for each genotype
p1_list = create_plot(gaussian_reshaped, "Gaussian Reaction Norm", line_colors, replicates_gaussian_selected)
p2_list = create_plot(sinusoidal_reshaped, "Sinusoidal Reaction Norm", line_colors, replicates_sinusoidal_selected)
p3_list = create_plot(wave_reshaped, "Wave Reaction Norm", line_colors, replicates_wave_selected)

# Create combined plot with one curve per type and replicates
combined_sampled = data.frame(
  Gaussian = gaussian_reshaped[, 2],
  Sinusoidal = sinusoidal_reshaped[, 2],
  Wave = wave_reshaped[, 2],
  Environment = environmental_factors
)

combined_replicates_selected = rbind(
  replicates_gaussian_selected[replicates_gaussian_selected$Genotype == sampled_genotypes[2], ],
  replicates_sinusoidal_selected[replicates_sinusoidal_selected$Genotype == sampled_genotypes[2], ],
  replicates_wave_selected[replicates_wave_selected$Genotype == sampled_genotypes[2], ]
)
combined_replicates_selected$Type = rep(c("Gaussian", "Sinusoidal", "Wave"), each = nrow(combined_replicates_selected) / 3)
combined_replicates_selected$Type = factor(combined_replicates_selected$Type, levels = c("Gaussian", "Sinusoidal", "Wave"))

p_combined_with_arrows = ggplot(combined_sampled, aes(x = Environment)) +
  geom_line(aes(y = Gaussian, color = "Gaussian")) +
  geom_line(aes(y = Sinusoidal, color = "Sinusoidal")) +
  geom_line(aes(y = Wave, color = "Wave")) +
  geom_point(data = combined_replicates_selected, aes(x = Environment, y = Value, color = Type), size = 2, shape=21) +
  scale_color_manual(values = scales::hue_pal()(3)) +
  labs(title = "Non-linear Reaction Norms Crossing with Replicates") +
  ylab("Trait") +
  theme_minimal() +
  theme(legend.position = "right") +
  add_axis_arrows()

# Version without axis labels for the combined plot
p_combined_without_labels = p_combined_with_arrows +
  theme(axis.text = element_blank(), axis.title = element_blank())
# Save all plots into a single PDF file
pdf("non_linear_reaction_norms.pdf", width = 8, height = 6)
print(p1_list$with_labels)
print(p1_list$without_labels)
print(p2_list$with_labels)
print(p2_list$without_labels)
print(p3_list$with_labels)
print(p3_list$without_labels)
print(p_combined_with_arrows)
print(p_combined_without_labels)
dev.off()
