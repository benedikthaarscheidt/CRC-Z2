library(ggplot2)

# Function to generate diverse non-linear crossing functions with visible traits
generate_diverse_functions <- function(num_functions = 5, x_range = c(0, 10), length_out = 1000) {
  # Generate x values
  x <- seq(x_range[1], x_range[2], length.out = length_out)
  
  # Create a data frame to hold all function data
  data <- data.frame()
  
  for (i in 1:num_functions) {
    # Randomly choose a function type
    chosen_type <- sample(c("quadratic", "polynomial", "sigmoid", "logistic"),1)#, "sinusoidal"), 1)
    
    # Generate y-values based on the chosen type
    if (chosen_type == "quadratic") {
      a <- runif(1, -0.00000000005, 0.000000000005)
      b <- runif(1, -0.001, 0.001)
      c <- runif(1, 10, 50)
      shift <- runif(1, 0, 50)
      y <- a * (x - shift)^2 + b * (x - shift) + c + i *10 # Vertical offset added
      print("quadratic")
    } else if (chosen_type == "polynomial") {
      coeff <- runif(4, -3, 3)
      c <- runif(1, 10, 50)
      y <- coeff[1] * x^3 + coeff[2] * x^2 + coeff[3] * x + coeff[4] + i * 10 + c # Vertical offset added
      print("polynomial")
    } else if (chosen_type == "sigmoid") {
      scale <- runif(1, 2, 100)
      offset <- runif(1, 0, 10)
      center <- runif(1, 10, 70)
      y <- scale / (1 + exp(-0.1 * (x - center))) + offset + i * 10  # Vertical offset added
      print("sigmoid")
    } else if (chosen_type == "logistic") {
      scale <- runif(1, 2, 5)
      offset <- runif(1, -5, 5)
      y <- scale / (1 + exp(-0.1 * x)) + offset + i * 10  # Vertical offset added
      print("logistic")
    } else if (chosen_type == "sinusoidal") {
      amplitude <- runif(1, 1, 4)
      phase_shift <- runif(1, 0, 2 * pi)
      y <- amplitude * sin(0.1 * x + phase_shift) + runif(1, -5, 5) + i * 10  # Vertical offset added
      print("sinusoidal")
    }
    
    
    # Store data in a data frame
    df <- data.frame(x = x, y = y, Genotype = paste("Genotype", i), Type = chosen_type)
    data <- rbind(data, df)
  }
  
  return(data)
}

# Generate and plot non-linear crossing functions
plot_data <- generate_diverse_functions(num_functions = 5)

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = Genotype)) +
  geom_line(size = 1.2) +
  labs(title = "Diverse Non-Linear Functions with Distinct Patterns",
       x = "X-axis",
       y = "Trait Value") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom")
