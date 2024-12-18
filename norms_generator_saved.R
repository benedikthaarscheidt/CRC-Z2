library(ggplot2)

# Function to generate a Gaussian-like reaction norm
generate_gaussian_reaction_norm <- function(env, base_shift, amplitude = 7, width = 0.2, center = 5, slope_variation = 0.6) {
  slope <- runif(1, 1 - slope_variation, 1 + slope_variation)
  curve <- base_shift + slope * amplitude * exp(-width * (env - center)^2)
  return(curve)
}

# Function to generate a sinusoidal reaction norm
generate_sinusoidal_reaction_norm <- function(env, base_shift, amplitude = 8, frequency = 0.5, phase = 0, slope_variation = 0.4) {
  slope <- runif(1, 1 - slope_variation, 1 + slope_variation)
  curve <- base_shift + slope * amplitude * sin(frequency * env + phase)
  return(curve)
}

# Function to generate a wave-like reaction norm
generate_wave_reaction_norm <- function(env, base_shift, amplitude = 10, frequency = 1, slope_variation = 0.4) {
  slope <- runif(1, 1 - slope_variation, 1 + slope_variation)
  curve <- base_shift + slope * amplitude * sin(frequency * (env - 5)) * exp(-0.1 * (env - 5)^2)
  return(curve)
}

# Set environmental factors
environmental_factors <- seq(1, 10, length.out = 100)

# Define base shifts for reaction norms
base_shifts <- seq(10, 300, length.out = 50)

# Generate Gaussian-like reaction norms
reaction_norms_gaussian <- data.frame(Environment = environmental_factors)
for (i in seq_along(base_shifts)) {
  reaction_norms_gaussian[[paste0("Gaussian_", i)]] <- generate_gaussian_reaction_norm(environmental_factors, base_shifts[i])
}

# Generate Sinusoidal reaction norms
reaction_norms_sinusoidal <- data.frame(Environment = environmental_factors)
for (i in seq_along(base_shifts)) {
  reaction_norms_sinusoidal[[paste0("Sinusoidal_", i)]] <- generate_sinusoidal_reaction_norm(environmental_factors, base_shifts[i])
}

# Generate Wave-like reaction norms
reaction_norms_wave <- data.frame(Environment = environmental_factors)
for (i in seq_along(base_shifts)) {
  reaction_norms_wave[[paste0("Wave_", i)]] <- generate_wave_reaction_norm(environmental_factors, base_shifts[i])
}

# Reshape the data for plotting
library(reshape2)
reaction_norms_gaussian_melted <- melt(reaction_norms_gaussian, id.vars = "Environment", 
                                       variable.name = "Curve", value.name = "TraitValue")
reaction_norms_sinusoidal_melted <- melt(reaction_norms_sinusoidal, id.vars = "Environment", 
                                         variable.name = "Curve", value.name = "TraitValue")
reaction_norms_wave_melted <- melt(reaction_norms_wave, id.vars = "Environment", 
                                   variable.name = "Curve", value.name = "TraitValue")

# Create separate plots for each reaction norm type
p1 <- ggplot(reaction_norms_gaussian_melted, aes(x = Environment, y = TraitValue, color = Curve)) +
  geom_line() +
  labs(title = "Gaussian-like Reaction Norms", x = "Environmental Factor", y = "Trait Values") +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- ggplot(reaction_norms_sinusoidal_melted, aes(x = Environment, y = TraitValue, color = Curve)) +
  geom_line() +
  labs(title = "Sinusoidal Reaction Norms", x = "Environmental Factor", y = "Trait Values") +
  theme_minimal() +
  theme(legend.position = "none")

p3 <- ggplot(reaction_norms_wave_melted, aes(x = Environment, y = TraitValue, color = Curve)) +
  geom_line() +
  labs(title = "Wave-like Reaction Norms", x = "Environmental Factor", y = "Trait Values") +
  theme_minimal() +
  theme(legend.position = "none")

# Combine all reaction norms into a single data frame
combined_curves <- rbind(
  cbind(reaction_norms_gaussian_melted, Type = "Gaussian"),
  cbind(reaction_norms_sinusoidal_melted, Type = "Sinusoidal"),
  cbind(reaction_norms_wave_melted, Type = "Wave")
)

# Combined plot for all reaction norms
p_combined <- ggplot(combined_curves, aes(x = Environment, y = TraitValue, color = Type, group = Curve)) +
  geom_line() +
  labs(title = "Combined Reaction Norms", x = "Environmental Factor", y = "Trait Values") +
  theme_minimal() +
  theme(legend.position = "right")

# Display plots
print(p1)
print(p2)
print(p3)
print(p_combined)
