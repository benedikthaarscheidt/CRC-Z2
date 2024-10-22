# Generate log-normally distributed data
set.seed(123)
n <- 1000
meanlog <- 1  # Mean in log-space
sdlog <- 0.5  # Standard deviation in log-space

log_normal_data <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)

# Plot the histogram to visualize the distribution
hist(log_normal_data, breaks = 30, col = "lightblue", main = "Proper Log-Normal Distribution", xlab = "Values", ylab = "Frequency")

hist(log(log_normal_data+1),breaks = 30, col = "lightblue", main = "pups", xlab = "Values", ylab = "Frequency")