library(ggplot2)
library(reshape2)

#13 not bad but scale is bad for seeing sd 
set.seed(2)
rownames = seq(1, 50)
env = seq(1, 10)
num_genotypes_per_plot=3
genotypes=sample(seq(1:50),num_genotypes_per_plot)


set.seed(NULL)
########################## non crossing linear reaction norms


RN_linear_NC = data.frame(matrix(NA, nrow = length(rownames), ncol = length(env)))

colnames(RN_linear_NC) = env

RN_linear_NC$Row = rownames

# Apply linear values with increasing offsets
offsets = seq(1, 10, length.out = nrow(RN_linear_NC))  # Customize the range as needed
slope = seq(0,5,length.out=nrow(RN_linear_NC))
for (i in 1:nrow(RN_linear_NC)) {
  RN_linear_NC[i, 1:length(env)] = slope[i]*env + offsets[i]  # Add increasing offset to each row
}

# Transform to long format for ggplot
RN_linear_NC_long = melt(RN_linear_NC, id.vars = "Row", variable.name = "Env", value.name = "Value")
RN_linear_NC_long$Row = as.factor(RN_linear_NC_long$Row)
RN_linear_NC_long$Env = as.numeric(as.character(RN_linear_NC_long$Env))

# Subset for a clearer plot
RN_linear_NC_long_subset = RN_linear_NC_long[RN_linear_NC_long$Row %in% genotypes, ]

# Plot using ggplot2
ggplot(data = RN_linear_NC_long_subset, aes(x = Env, y = Value, color = Row, group = Row)) +
  geom_line(size = 1) +
  labs(title = "Linear non-crossing reaction norms with increasing offsets",
       x = "Environmental Factor",
       y = "Trait Values") +
  theme_minimal() +
  theme(legend.title = element_blank())


###################################### Linear crossing reaction norms 


RN_linear_C = data.frame(matrix(NA, nrow = length(rownames), ncol = length(env)))

colnames(RN_linear_C) = env

RN_linear_C$Row = rownames
set.seed(19)

for (i in 1:nrow(RN_linear_C)) {
  slope = runif(1, min = -1.5, max = 2)   
  if (slope < 0) {
    offset = runif(1, min = 10, max = 15) 
  } else {
    offset = runif(1, min = 0, max = 10)  
  }
  
  RN_linear_C[i, 1:length(env)] = slope * env + offset
}

# Ensure no negative values (just a safety net)
RN_linear_C[RN_linear_C < 0] = 0

RN_linear_C_long = melt(RN_linear_C, id.vars = "Row", variable.name = "Env", value.name = "Value")

RN_linear_C_long$Row = as.factor(RN_linear_C_long$Row)
RN_linear_C_long$Env = as.numeric(as.character(RN_linear_C_long$Env))

RN_linear_C_long_subset = RN_linear_C_long[RN_linear_C_long$Row %in% genotypes, ]

ggplot(data = RN_linear_C_long_subset, aes(x = Env, y = Value, color = Row, group = Row)) +
  geom_line(size = 1) +
  labs(title = "Linear crossing reaction norms",
       x = "Environmental Factor",
       y = "Trair Values") +
  theme_minimal() +
  theme(legend.title = element_blank())


#################################### Non-linear non-crossing reaction norms

RN_nonlinear_NC = data.frame(matrix(NA, nrow = length(rownames), ncol = length(env)))

colnames(RN_nonlinear_NC) = env

set.seed(31)

offset=runif(1,5,10)
RN_nonlinear_NC[1,1]=offset
for(i in 2:ncol(RN_nonlinear_NC)){
  random_increment = rnorm(1, mean=0, sd=1)
  RN_nonlinear_NC[1,i]=RN_nonlinear_NC[1,i-1]+random_increment
}

for (i in 2:nrow(RN_nonlinear_NC)) {
  prev_row = (i - 1) 
  
  random_increments1 = runif(length(env), min = 1, max = 3)
  random_increments2 = runif(length(env), min = 0.001, max = 1)
  
  new_row = numeric(length(env))
  for (j in 1:length(env)) {
    
    selected_increment = ifelse(runif(1) < 0.5, random_increments1[j], random_increments2[j])
    
    new_row[j] = RN_nonlinear_NC[prev_row, j] + selected_increment
  }
  
  new_row = pmax(0, new_row)
  
  RN_nonlinear_NC[i, 1:length(env)] = new_row
}



# Create an empty data frame for the smoothed data
RN_nonlinear_NC_smoothed = RN_nonlinear_NC

# Loop through each row and apply smoothing
for (i in 1:nrow(RN_nonlinear_NC_smoothed)) {
  env_values = env  
  trait_values = as.numeric(RN_nonlinear_NC_smoothed[i, ])    
  # Check if there are enough data points for smoothing
  if (length(env_values) > 1) {
    spline_fit = smooth.spline(env_values, trait_values, spar = 0.4)
    
    # Replace the current row's data  with smoothed values
    RN_nonlinear_NC_smoothed[i, -1] = approx(spline_fit$x, spline_fit$y, xout = env_values)$y
  } else {
    warning(paste("Not enough data points for smoothing for row", i))
  }
}

RN_nonlinear_NC_smoothed$Row=rownames
RN_nonlinear_NC_smoothed_long = melt(RN_nonlinear_NC_smoothed, id.vars = "Row", variable.name = "Env", value.name = "Value")


RN_nonlinear_NC_smoothed_long$Env = as.numeric(as.character(RN_nonlinear_NC_smoothed_long$Env))
RN_nonlinear_NC_smoothed_long$Row = as.character(RN_nonlinear_NC_smoothed_long$Row)
RN_nonlinear_NC_smoothed_subset = RN_nonlinear_NC_smoothed_long[RN_nonlinear_NC_smoothed_long$Row %in% genotypes , ]



ggplot(data = RN_nonlinear_NC_smoothed_subset, aes(x = Env, y = Value, color = Row, group = Row)) +
  geom_line(size = 1) +
  labs(title = "Non-linear Non-crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Trait Values") +
  theme_minimal() +
  theme(legend.title = element_blank())
################################################# Non-linear crossing reaction norms

RN_nonlinear_C = data.frame(matrix(NA, nrow = length(rownames), ncol = length(env)))

colnames(RN_nonlinear_C) = env

set.seed(292)

for (i in 1:nrow(RN_nonlinear_C)) {
  
  offset = runif(1, min = 0, max = 30)
  
  RN_nonlinear_C[i, 1] = offset
  
  for (j in 2:ncol(RN_nonlinear_C)) {
    
    random_change = rnorm(1, mean = 0, sd = 1.5)
    
    RN_nonlinear_C[i, j] = RN_nonlinear_C[i, j - 1] + random_change
    
    RN_nonlinear_C[i, j] = pmax(0, RN_nonlinear_C[i, j])
  }
}

# Create an empty data frame for the smoothed data
RN_nonlinear_C_smoothed = RN_nonlinear_C

# Loop through each row and apply smoothing
for (i in 1:nrow(RN_nonlinear_C_smoothed)) {
  env_values = env  
  trait_values = as.numeric(RN_nonlinear_C_smoothed[i, ])    
  # Check if there are enough data points for smoothing
  if (length(env_values) > 1) {
    spline_fit = smooth.spline(env_values, trait_values, spar = 0.45)
    
    # Replace the current row's data  with smoothed values
    RN_nonlinear_C_smoothed[i, -1] = approx(spline_fit$x, spline_fit$y, xout = env_values)$y
  } else {
    warning(paste("Not enough data points for smoothing for row", i))
  }
}

RN_nonlinear_C_smoothed$Row=rownames
RN_nonlinear_C_smoothed_long = melt(RN_nonlinear_C_smoothed, id.vars = "Row", variable.name = "Env", value.name = "Value")

RN_nonlinear_C_smoothed_long$Env = as.numeric(as.character(RN_nonlinear_C_smoothed_long$Env))
RN_nonlinear_C_smoothed_long$Row = as.character(RN_nonlinear_C_smoothed_long$Row)
RN_nonlinear_C_smoothed_subset = RN_nonlinear_C_smoothed_long[RN_nonlinear_C_smoothed_long$Row %in% genotypes , ]

ggplot(data = RN_nonlinear_C_smoothed_subset, aes(x = Env, y = Value, color = Row, group = Row)) +
  geom_line(size = 1) +
  labs(title = "Non-linear crossing Reaction Norms",
       x = "Environmental Factor",
       y = "Trait Values") +
  theme_minimal() +
  theme(legend.title = element_blank())


##############################


