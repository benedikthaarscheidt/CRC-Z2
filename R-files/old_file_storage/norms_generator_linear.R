# Set environmental factors and number of genotypes
set.seed(13)
environmental_factors = seq(1, 10, length.out = 50)
num_genotypes = 50
num_genotypes_per_plot = 3

# Select genotypes for sampling
sampled_genotypes = sample(seq_len(num_genotypes), num_genotypes_per_plot)

# Function to generate linear non-crossing reaction norm
generate_linear_non_crossing_reaction_norm = function(env, base_shift, slope) {
  curve = slope * env + base_shift
  return(curve)
}

# Function to generate linear crossing reaction norm
generate_linear_crossing_reaction_norm = function(env) {
  slope = runif(1, min = -1.5, max = 2)
  if (slope < 0) {
    base_shift = runif(1, min = 10, max = 15)
  } else {
    base_shift = runif(1, min = 0, max = 10)
  }
  curve = slope * env + base_shift
  curve[curve < 0] = 0  # Ensure no negative values
  return(curve)
}

# Define base shifts and slopes for non-crossing reaction norms
base_shifts = seq(1, 10, length.out = num_genotypes)
slope_values = seq(0, 5, length.out = num_genotypes)

# Generate Linear Non-Crossing reaction norms
reaction_norms_linear_nc = data.frame(
  t(sapply(seq_along(base_shifts), function(i) generate_linear_non_crossing_reaction_norm(environmental_factors, base_shifts[i], slope_values[i])))
)
colnames(reaction_norms_linear_nc) = environmental_factors
rownames(reaction_norms_linear_nc) = paste0("Linear_NC_", seq_along(base_shifts))

# Generate Linear Crossing reaction norms
reaction_norms_linear_c = data.frame(
  t(sapply(seq_len(num_genotypes), function(i) generate_linear_crossing_reaction_norm(environmental_factors)))
)
colnames(reaction_norms_linear_c) = environmental_factors
rownames(reaction_norms_linear_c) = paste0("Linear_C_", seq_len(num_genotypes))

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

# Generate replicates for Linear Non-Crossing and Linear Crossing reaction norms
replicates_linear_nc = generate_replicates(reaction_norms_linear_nc)
replicates_linear_c = generate_replicates(reaction_norms_linear_c)

# Filter replicates for sampled genotypes only and select every 5th environment for plotting
replicates_linear_nc_sampled = replicates_linear_nc[replicates_linear_nc$Genotype %in% sampled_genotypes & replicates_linear_nc$Environment %in% environmental_factors[seq(1, length(environmental_factors), by = 5)], ]
replicates_linear_c_sampled = replicates_linear_c[replicates_linear_c$Genotype %in% sampled_genotypes & replicates_linear_c$Environment %in% environmental_factors[seq(1, length(environmental_factors), by = 5)], ]

# Reshape data for plotting
reshape_for_plot = function(data, type) {
  df = t(data)
  df = as.data.frame(df)
  colnames(df) = paste0(type, "_", seq_len(ncol(df)))
  df$Environment = environmental_factors
  return(df)
}

# Reshape sampled curves for plotting
linear_nc_sampled = reaction_norms_linear_nc[sampled_genotypes, ]
linear_c_sampled = reaction_norms_linear_c[sampled_genotypes, ]

linear_nc_reshaped = reshape_for_plot(linear_nc_sampled, "Linear_NC")
linear_c_reshaped = reshape_for_plot(linear_c_sampled, "Linear_C")

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
p1_list = create_plot(linear_nc_reshaped, "Linear non-crossing reaction norms", line_colors, replicates_linear_nc_sampled)
p2_list = create_plot(linear_c_reshaped, "Linear crossing reaction norms", line_colors, replicates_linear_c_sampled)

# Save all plots into a single PDF file
pdf("linear_reaction_norms.pdf", width = 8, height = 6)
print(p1_list$with_labels)
print(p1_list$without_labels)
print(p2_list$with_labels)
print(p2_list$without_labels)
dev.off()
