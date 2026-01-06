# Script: evaluate_gwas_accuracy.R
# Purpose: Compare GWAS results against simulated ground truth to validate plasticity scores.
source("~/CRC_1644_Z2/R-files/run_gwas_simulation.R")
library(dplyr)
library(ggplot2)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
results_dir <- latest_run
plots_dir <- paste(latest_run, "plots", sep = "/")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
# Try to find the results file (MLM preferred)
gwas_file <- file.path(results_dir, "gwas_results_mlm_80.csv")
if (!file.exists(gwas_file)) {
  gwas_file <- file.path(results_dir, "gwas_results_lm_80.csv")
}

if (!file.exists(gwas_file)) {
  stop("GWAS results file not found. Please wait for the simulation to finish.")
}

message("Loading GWAS results from: ", gwas_file)
gwas_res <- read.csv(gwas_file)

truth_file <- file.path(results_dir, "causal_snps_truth.csv")
message("Loading Truth from: ", truth_file)
truth <- read.csv(truth_file)

# -----------------------------------------------------------------------------
# 2. Define Significance Threshold
# -----------------------------------------------------------------------------
alpha <- 0.05
n_snps <- length(unique(gwas_res$SNP))
threshold <- alpha / n_snps
message(paste("Significance Threshold (Bonferroni): P <", format(threshold, scientific = TRUE)))

# === ADDED: Capture Genotype Count from Global SNP Matrix ===
if (exists("snp_matrix") && is.matrix(snp_matrix)) {
  n_genotypes <- nrow(snp_matrix)
} else {
  # Fallback: Infer from Kinship matrix or SNP column names if matrix is not global
  if (exists("kinship_matrix") && is.matrix(kinship_matrix)) {
    n_genotypes <- nrow(kinship_matrix)
  } else {
    # Estimate from data (less reliable, but better than nothing)
    n_genotypes <- length(unique(gsub("_Rep_.*", "", unique(gwas_res$Genotype_ID))))
  }
}

# -----------------------------------------------------------------------------
# 3. Evaluate Performance
# -----------------------------------------------------------------------------
# Filter significant hits
sig_hits <- gwas_res %>% filter(P_Value < threshold)

# Initialize Summary Table
# We want to know: For each Score, how many of EACH Parameter's SNPs did it find?
scores <- unique(gwas_res$Score)
params <- unique(truth$Parameter)

evaluation_list <- list()

for (s in scores) {
  # Get hits for this score (across all types/forms)
  # Note: You might want to analyze per Form (Type) or aggregated.
  # Let's aggregate for now, or check if specific forms detected specific things.
  # Since the truth is "Global" (SNPs affect the generator parameters),
  # and we ran GWAS on "Type" subsets, we should look at the union of hits?
  # Or maybe the score is calculated per form.

  # Let's look at unique SNPs found by this score across ANY form
  score_hits <- unique(sig_hits$SNP[sig_hits$Score == s])

  for (p in params) {
    # True SNPs for this parameter
    true_snps <- truth$SNP[truth$Parameter == p]
    n_true <- length(true_snps)

    # Matches
    matches <- intersect(score_hits, true_snps)
    n_found <- length(matches)

    # False Positives (SNPs significant for this score but NOT in this parameter's truth)
    # Note: This is tricky because a SNP might be causal for *another* parameter.
    # Strict FP: Significant but not causal for ANY parameter?
    # Or Specific FP: Significant but not causal for THIS parameter?
    # Usually, if it detects a Slope SNP, that's a "True Positive" for Slope.
    # If it detects a VE SNP, is that a False Positive for Slope?
    # Yes, if we interpret the score as "Slope".
    # But here we are asking "What does this score detect?".

    evaluation_list[[paste(s, p)]] <- data.frame(
      Score = s,
      Parameter = p,
      True_SNPs = n_true,
      Found = n_found,
      Power = n_found / n_true
    )
  }

  # Global False Positives for this Score
  # SNPs found that are NOT in the truth table at all
  all_causal <- unique(truth$SNP)
  n_fp_global <- length(setdiff(score_hits, all_causal))

  evaluation_list[[paste(s, "Global_FP")]] <- data.frame(
    Score = s,
    Parameter = "False_Positives_Global",
    True_SNPs = NA,
    Found = n_fp_global,
    Power = NA
  )
}

eval_df <- do.call(rbind, evaluation_list)

# -----------------------------------------------------------------------------
# 4. Save and Plot
# -----------------------------------------------------------------------------
write.csv(eval_df, file.path(results_dir, "gwas_accuracy_summary.csv"), row.names = FALSE)

# Heatmap of Power
power_df <- eval_df %>%
  filter(Parameter != "False_Positives_Global") %>%
  mutate(Power = as.numeric(Power))

subtitle_text <- paste0(
  "Simulation: N_Genotypes = ", n_genotypes, 
  ", N_SNPs = ", n_snps, 
  " (Bonferroni P < ", format(threshold, digits=2, scientific=TRUE), ")"
)

p <- ggplot(power_df, aes(x = Parameter, y = Score, fill = Power)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  geom_text(aes(label = round(Power, 2)), color = "black", size = 3) +
  labs(
    title = "GWAS Power: Which Scores Detect Which Parameters?",
    subtitle = subtitle_text, # <<< USING THE NEW SUBTITLE
    x = "Causal Parameter (Truth)", y = "Plasticity Score"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file.path(plots_dir, "gwas_power_heatmap.pdf"), width = 10, height = 8)
print(p)
dev.off()

message("Evaluation Complete!")
message("Summary: ", file.path(results_dir, "gwas_accuracy_summary.csv"))
message("Plot: ", file.path(plots_dir, "gwas_power_heatmap.pdf"))
