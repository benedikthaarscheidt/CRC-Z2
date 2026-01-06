# Script: generate_reaction_norms.R
# Purpose: Generate raw reaction norm data (phenotypes) from genetic parameters.
#          This separates the "Simulation" step from the "Score Calculation" step.

library(dplyr)
library(MASS)

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------
set.seed(42)
output_dir <- "~/CRC_1644_Z2/synthetic_data/fixed_full"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

params_file <- "~/CRC_1644_Z2/synthetic_data/genetics/genotypic_parameters.csv"
if (!file.exists(params_file)) stop("Genetic parameters not found: ", params_file)

message("Loading genetic parameters...")
params_df <- read.csv(params_file)

# Source generators (they now have if(sys.nframe()==0) blocks so they won't auto-run)
source("~/CRC_1644_Z2/R-files/norms_generator2_linear.R")
source("~/CRC_1644_Z2/R-files/norms_generator2_nonlinear.R")

# -----------------------------------------------------------------------------
# 2. Define Generation Logic
# -----------------------------------------------------------------------------
environmental_factors <- seq(1, 10, length.out = 50)
num_replicates <- 3

# Helper to extract parameters for a specific genotype row
get_params_for_row <- function(row_idx) {
  p <- params_df[row_idx, ]
  list(
    BaseShift = p$BaseShift,
    Slope = p$Slope,
    VE = if("VE" %in% colnames(p)) p$VE else 1,
    VS = if("VS" %in% colnames(p)) p$VS else 0.2,
    Amplitude = if("Amplitude" %in% colnames(p)) p$Amplitude else 4,
    Width = if("Width" %in% colnames(p)) p$Width else 0.2,
    Center = if("Center" %in% colnames(p)) p$Center else 5,
    Frequency = if("Frequency" %in% colnames(p)) p$Frequency else 0.5,
    Phase = if("Phase" %in% colnames(p)) p$Phase else 0
  )
}

# -----------------------------------------------------------------------------
# 3. Generate Data
# -----------------------------------------------------------------------------
all_data_list <- list()

message("Generating reaction norms for 80 genotypes...")

for (i in 1:nrow(params_df)) {
  p <- get_params_for_row(i)
  form <- params_df$Form[i]
  geno_id <- params_df$Genotype_Name[i]
  
  # Determine which generator to use
  if (form == "Linear") {
    # Linear Generator Logic
    # We need to adapt the logic from norms_generator2_linear.R to work for a single genotype
    
    cov_matrix <- build_cov_matrix(p$VE, p$VS, rES = 1) # Assuming rES=1
    random_effects <- MASS::mvrnorm(num_replicates, mu = c(0, 0), Sigma = cov_matrix)
    
    fixed_effect <- generate_fixed_norm(p$BaseShift, p$Slope, environmental_factors, outliers = "none")
    
    reps <- do.call(rbind, lapply(1:num_replicates, function(r) {
      u0 <- random_effects[r, 1]
      u1 <- random_effects[r, 2]
      curve <- fixed_effect + (u0 + u1 * environmental_factors)
      data.frame(Genotype = geno_id, Replicate = r, Environment = environmental_factors, Trait = curve, ReactionNorm = "Linear")
    }))
    
    # Add Replicate 0 (Mean)
    rep0 <- data.frame(Genotype = geno_id, Replicate = 0, Environment = environmental_factors, Trait = fixed_effect, ReactionNorm = "Linear")
    all_data_list[[i]] <- rbind(rep0, reps)
    
  } else if (form == "Gaussian") {
    # Gaussian Logic
    norm_func <- function(e, b, s) generate_gaussian_reaction_norm(e, b, s, amplitude=p$Amplitude, width=p$Width, center=p$Center)
    
    dat <- generate_genotype_data(
      base_shift = p$BaseShift, slope = p$Slope, env = environmental_factors, norm_function = norm_func,
      num_replicates = num_replicates, var_offset = p$VE, var_slope = p$VS, correlation = -0.5,
      genotype_id = geno_id, norm_type = "Gaussian", outliers = "none"
    )
    # Fix columns to match Linear (remove extra cols for now if needed, or keep them)
    # We keep minimal columns for consistency
    all_data_list[[i]] <- dat[, c("Genotype", "Replicate", "Environment", "Trait", "ReactionNorm")]
    
  } else if (form == "Sinusoidal") {
    norm_func <- function(e, b, s) generate_sinusoidal_reaction_norm(e, b, s, amplitude=p$Amplitude, frequency=p$Frequency, phase=p$Phase)
    
    dat <- generate_genotype_data(
      base_shift = p$BaseShift, slope = p$Slope, env = environmental_factors, norm_function = norm_func,
      num_replicates = num_replicates, var_offset = p$VE, var_slope = p$VS, correlation = -0.3,
      genotype_id = geno_id, norm_type = "Sinusoidal", outliers = "none"
    )
    all_data_list[[i]] <- dat[, c("Genotype", "Replicate", "Environment", "Trait", "ReactionNorm")]
    
  } else if (form == "Wave") {
    norm_func <- function(e, b, s) generate_wave_reaction_norm(e, b, s, amplitude=p$Amplitude, frequency=p$Frequency)
    
    dat <- generate_genotype_data(
      base_shift = p$BaseShift, slope = p$Slope, env = environmental_factors, norm_function = norm_func,
      num_replicates = num_replicates, var_offset = p$VE, var_slope = p$VS, correlation = -0.7,
      genotype_id = geno_id, norm_type = "Wave", outliers = "none"
    )
    all_data_list[[i]] <- dat[, c("Genotype", "Replicate", "Environment", "Trait", "ReactionNorm")]
  }
}

combined_data <- do.call(rbind, all_data_list)

# -----------------------------------------------------------------------------
# 4. Save Output
# -----------------------------------------------------------------------------
output_file <- file.path(output_dir, "reaction_norms_full.rds")
message("Saving reaction norms to: ", output_file)
saveRDS(combined_data, output_file)

# Also save as CSV for inspection
write.csv(combined_data, file.path(output_dir, "reaction_norms_full.csv"), row.names = FALSE)

message("Done.")
