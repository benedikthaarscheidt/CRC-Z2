
if (!exists("data_loaded") || !data_loaded) {
  source("~/CRC 1622 - Z2/R-files/head3.R")
  data_loaded = TRUE
}

library(factoextra)
library(dbscan)
library(cluster)
library(mclust)
library(pheatmap)
library(corrplot)
library(vegan) 
library(energy)
library(knitr)
# Helper function: prepare a score vector for a given genotype form from a score dataframe.
# Assumes each score dataframe has two columns: first = score value, second = genotype form.
prepare_score_vector = function(score_df, genotype_form) {
  # Subset the rows corresponding to the given genotype form and return the score column
  subset(score_df, score_df[,2] == genotype_form)[,1]
}


permTestARI <- function(clusters1, clusters2, n_perm = 1000, alternative = "greater") {
  observed_ari <- adjustedRandIndex(clusters1, clusters2)
  permuted_ari <- numeric(n_perm)
  for (k in 1:n_perm) {
    permuted_labels <- sample(clusters2)
    permuted_ari[k] <- adjustedRandIndex(clusters1, permuted_labels)
  }
  if (alternative == "greater") {
    p_val <- sum(permuted_ari >= observed_ari) / n_perm
  } else if (alternative == "less") {
    p_val <- sum(permuted_ari <= observed_ari) / n_perm
  } else if (alternative == "two.sided") {
    null_mean <- mean(permuted_ari)
    p_val <- sum(abs(permuted_ari - null_mean) >= abs(observed_ari - null_mean)) / n_perm
  } else {
    stop("alternative must be one of 'greater', 'less', or 'two.sided'")
  }
  return(list(observed_ARI = observed_ari, p_value = p_val))
}

analyze_scores_by_genotype_correlation <- function(scores_list, genotype_forms, n_cluster, rn_matrix_list, eps = 0.1) {
  
  ## Initialize lists to store results
  correlation_distances              <- list()  # Distance matrices from correlations
  correlation_distances_per_geno     <- list()  # Correlation matrices computed per genotype form
  clustering_results                 <- list()  # Clustering assignments per genotype form
  silhouette_scores                  <- list()  # Average silhouette scores per genotype form
  matching_index_matrix_genotype_forms <- list() # ARI (with p-values) across genotype forms
  rn_comparison                      <- list()  # Reaction norm vs. plasticity score comparisons
  
  ## Define correlation and clustering methods
  corr_methods       <- c("kendall")
  clustering_methods <- c("hclust", "pam", "dbscan")
  
  pdf(file = "correlation_plots.pdf", width = 12, height = 12)
  tryCatch({
    ############################################
    ## 1. Per-Genotype Form Analysis
    ############################################
    for (form in genotype_forms) {
      # Extract score vectors for the current genotype form
      score_vectors <- lapply(scores_list, function(x) {
        prepare_score_vector(x, form)
      })
      common_length <- min(sapply(score_vectors, length))
      score_matrix  <- as.data.frame(sapply(score_vectors, function(v) as.numeric(v[1:common_length])))
      score_names   <- names(scores_list)
      
      ## Initialize storage for this genotype form
      correlation_distances[[form]] <- list()
      clustering_results[[form]]    <- list()
      silhouette_scores[[form]]     <- list()
      
      for (cm in clustering_methods) {
        clustering_results[[form]][[cm]] <- list()
        silhouette_scores[[form]][[cm]]  <- list()
      }
      
      # Loop over correlation methods for per-genotype analysis
      for (cor_method in corr_methods) {
        # Compute correlation matrix (transpose so that correlation is computed across columns)
        corr_mat <- cor(t(score_matrix), method = cor_method)
        heatmap(corr_mat,
                main = paste("Correlation Matrix (", cor_method, ") for", form),
                col = heat.colors(256), scale = "none",
                labRow = score_names, labCol = score_names)
        
        # Compute distance matrix (using get_dist if available)
        dist_matrix <- as.data.frame(as.matrix(get_dist(t(score_matrix), method = cor_method)))
        correlation_distances[[form]][[cor_method]] <- dist_matrix
        correlation_distances_per_geno[[form]][[cor_method]] <- as.data.frame(
          as.matrix(cor(t(score_matrix), method = cor_method))
        )
        
        ## Hierarchical clustering
        hc <- hclust(as.dist(dist_matrix), method = "complete")
        clusters_hc <- cutree(hc, k = n_cluster)
        clustering_results[[form]][["hclust"]][[cor_method]] <- clusters_hc
        sil_hc <- mean(silhouette(clusters_hc, as.dist(dist_matrix))[, "sil_width"])
        silhouette_scores[[form]][["hclust"]][[cor_method]] <- sil_hc
        plot(hc, labels = score_names,
             main = paste("Hierarchical Clustering (", cor_method, ") for", form),
             xlab = "", sub = "")
        
        ## PAM clustering
        pam_res <- pam(as.dist(dist_matrix), k = n_cluster)
        clusters_pam <- pam_res$clustering
        clustering_results[[form]][["pam"]][[cor_method]] <- clusters_pam
        sil_pam <- mean(silhouette(clusters_pam, as.dist(dist_matrix))[, "sil_width"])
        silhouette_scores[[form]][["pam"]][[cor_method]] <- sil_pam
        invisible(capture.output({
          fviz_silhouette(silhouette(pam_res),
                          main = paste("PAM Silhouette (", cor_method, ") for", form))
        }))
        
        ## DBSCAN clustering
        dbscan_res <- dbscan(as.dist(dist_matrix), eps = eps)
        clusters_db <- dbscan_res$cluster
        clustering_results[[form]][["dbscan"]][[cor_method]] <- clusters_db
        sil_db <- if (length(unique(clusters_db)) > 1)
          mean(silhouette(clusters_db, as.dist(dist_matrix))[, "sil_width"])
        else 0
        silhouette_scores[[form]][["dbscan"]][[cor_method]] <- sil_db
        mds <- cmdscale(as.dist(dist_matrix))
        plot(mds, col = clusters_db + 1, pch = 19,
             main = paste("DBSCAN Clustering (", cor_method, ") for", form),
             xlab = "MDS Dimension 1", ylab = "MDS Dimension 2")
        legend("topright", legend = unique(clusters_db),
               col = unique(clusters_db) + 1, pch = 19)
      }
    }
    
    ############################################
    ## Compute ARI Across Genotype Forms with Permutation Test
    ############################################
    # For each clustering method and each correlation method, we compare the clustering
    # results across the different genotype forms using a permutation test to obtain p-values.
    for (cm in clustering_methods) {
      for (cor_method in corr_methods) {
        forms <- genotype_forms
        n_forms <- length(forms)
        ari_obs <- matrix(NA, nrow = n_forms, ncol = n_forms)
        ari_p   <- matrix(NA, nrow = n_forms, ncol = n_forms)
        rownames(ari_obs) <- colnames(ari_obs) <- forms
        rownames(ari_p) <- colnames(ari_p) <- forms
        
        if (n_forms > 1) {
          for (i in 1:(n_forms - 1)) {
            for (j in (i + 1):n_forms) {
              clusters_i <- clustering_results[[forms[i]]][[cm]][[cor_method]]
              clusters_j <- clustering_results[[forms[j]]][[cm]][[cor_method]]
              test_res <- permTestARI(clusters_i, clusters_j, n_perm = 1000, alternative = "greater")
              ari_obs[i, j] <- test_res$observed_ARI
              ari_obs[j, i] <- test_res$observed_ARI
              ari_p[i, j]   <- test_res$p_value
              ari_p[j, i]   <- test_res$p_value
            }
          }
        }
        matching_index_matrix_genotype_forms[[paste(cm, cor_method, sep = "_")]] <- list(observed = ari_obs, p_value = ari_p)
      }
    }
    
    ######################################################
    ## 2. Combined Analysis (Separate from per-genotype)
    ######################################################
    combined_form <- "combined"
    correlation_distances[[combined_form]] <- list()
    clustering_results[[combined_form]]  <- list()
    silhouette_scores[[combined_form]]   <- list()
    for (cm in clustering_methods) {
      clustering_results[[combined_form]][[cm]] <- list()
      silhouette_scores[[combined_form]][[cm]]  <- list()
    }
    
    # Create a matrix where each column is one score from the raw data:
    score_matrix_combined <- do.call(cbind, lapply(scores_list, function(x) as.numeric(x[, 1])))
    # Trim to a common length across all columns:
    common_length_combined <- min(sapply(as.data.frame(score_matrix_combined), length))
    score_matrix_combined <- as.data.frame(sapply(as.data.frame(score_matrix_combined), function(v) v[1:common_length_combined]))
    score_names <- names(scores_list)
    
    # Now compute correlation across columns (i.e., between scores).
    for (cor_method in corr_methods) {
      corr_mat <- cor(score_matrix_combined, method = cor_method)
      heatmap(corr_mat,
              main = paste("Correlation Matrix (", cor_method, ") for Combined Analysis"),
              col = heat.colors(256), scale = "none",
              labRow = score_names, labCol = score_names)
      
      # Compute a distance matrix from the correlation matrix (e.g., 1 - correlation)
      dist_matrix <- as.dist(1 - corr_mat)
      
      # Store the combined analysis distance matrix in the return object:
      correlation_distances[[combined_form]][[cor_method]] <- as.matrix(dist_matrix)
      
      ## Hierarchical clustering on combined data.
      hc <- hclust(dist_matrix, method = "complete")
      clusters_hc <- cutree(hc, k = n_cluster)
      clustering_results[[combined_form]][["hclust"]][[cor_method]] <- clusters_hc
      sil_hc <- mean(silhouette(clusters_hc, dist_matrix)[, "sil_width"])
      silhouette_scores[[combined_form]][["hclust"]][[cor_method]] <- sil_hc
      plot(hc, labels = score_names,
           main = paste("Hierarchical Clustering (", cor_method, ") for Combined Analysis"),
           xlab = "", sub = "")
      
      ## PAM clustering on combined data.
      pam_res <- pam(dist_matrix, k = n_cluster)
      clusters_pam <- pam_res$clustering
      clustering_results[[combined_form]][["pam"]][[cor_method]] <- clusters_pam
      sil_pam <- mean(silhouette(clusters_pam, dist_matrix)[, "sil_width"])
      silhouette_scores[[combined_form]][["pam"]][[cor_method]] <- sil_pam
      fviz_silhouette(silhouette(pam_res),
                      main = paste("PAM Silhouette (", cor_method, ") for Combined Analysis"))
      
      ## DBSCAN clustering on combined data.
      dbscan_res <- dbscan(dist_matrix, eps = eps)
      clusters_db <- dbscan_res$cluster
      clustering_results[[combined_form]][["dbscan"]][[cor_method]] <- clusters_db
      sil_db <- if (length(unique(clusters_db)) > 1)
        mean(silhouette(clusters_db, dist_matrix)[, "sil_width"])
      else 0
      silhouette_scores[[combined_form]][["dbscan"]][[cor_method]] <- sil_db
      mds <- cmdscale(dist_matrix)
      plot(mds, col = clusters_db + 1, pch = 19,
           main = paste("DBSCAN Clustering (", cor_method, ") for Combined Analysis"),
           xlab = "MDS Dimension 1", ylab = "MDS Dimension 2")
      legend("topright", legend = unique(clusters_db),
             col = unique(clusters_db) + 1, pch = 19)
    }
    
  }, error = function(e) {
    message("An error occurred: ", e$message)
  }, finally = {
    dev.off()
  })
  
  return(list(
    correlation_distances             = correlation_distances,
    clustering_results                = clustering_results,
    silhouette_scores                 = silhouette_scores,
    matching_index_matrix_genotype_forms = matching_index_matrix_genotype_forms,
    reaction_norm_comparison          = rn_comparison
  ))
}

############################

analyze_direct_clusterings = function(scores_list, genotype_forms, rn_matrix_list, n_cluster = 3, eps = 0.1) {
  # Load required packages (install if needed)
  if (!require(cluster)) stop("Package 'cluster' required")
  if (!require(dbscan)) stop("Package 'dbscan' required")
  if (!require(factoextra)) stop("Package 'factoextra' required")
  
  # --- Helper functions for distance measures ---
  cosine_distance <- function(mat) {
    # Normalize rows (each row becomes a unit vector)
    norm_mat <- mat / sqrt(rowSums(mat^2))
    # Compute cosine similarity: dot product between rows
    cosine_sim <- norm_mat %*% t(norm_mat)
    # Cosine distance = 1 - cosine similarity
    as.dist(1 - cosine_sim)
  }
  
  correlation_distance <- function(mat) {
    # Compute Pearson correlation between rows (compare each instance)
    corr_mat <- cor(t(mat))
    # Define distance as 1 - correlation
    as.dist(1 - corr_mat)
  }
  
  # Define a vector with distance measure names to loop over
  dist_methods <- c("euclidean", "manhattan", "cosine", "correlation")
  
  # Prepare lists to store results for each analysis
  rn_cluster_results    = list()  # Reaction Norm clustering results
  score_cluster_results = list()  # Score clustering results
  
  # Open one PDF device for all plots.
  pdf(file = "direct_clustering_plots.pdf", width = 12, height = 12)
  tryCatch({
    
    ###############################
    ## 1. Reaction Norm Clustering
    ###############################
    # Assume each rn_matrix_list[[form]] is a numeric matrix with rows = genotypes for that form.
    # We combine them (stack rows) so that each row is one genotype.
    rn_matrix_combined = do.call(rbind, rn_matrix_list)
    rn_names = rownames(rn_matrix_combined)
    
    rn_cluster_results = list()
    for (method in dist_methods) {
      if (method == "euclidean") {
        dmat <- dist(rn_matrix_combined, method = "euclidean")
      } else if (method == "manhattan") {
        dmat <- dist(rn_matrix_combined, method = "manhattan")
      } else if (method == "cosine") {
        dmat <- cosine_distance(rn_matrix_combined)
      } else if (method == "correlation") {
        dmat <- correlation_distance(rn_matrix_combined)
      } else {
        stop("Unknown method: ", method)
      }
      
      rn_cluster_results[[method]] <- list(distance = dmat)
      
      # Plot heatmap of distance matrix
      heatmap(as.matrix(dmat),
              main = paste("Reaction Norm Heatmap (", method, ")", sep = ""),
              col = heat.colors(256), scale = "none",
              labRow = rn_names, labCol = rn_names)
      
      # Hierarchical clustering (complete linkage)
      hc <- hclust(dmat, method = "complete")
      rn_cluster_results[[method]]$hclust <- hc
      plot(hc, labels = rn_names,
           main = paste("Reaction Norm hclust (", method, ")", sep = ""),
           xlab = "", sub = "")
      
      # PAM clustering
      pam_res <- pam(dmat, k = n_cluster)
      rn_cluster_results[[method]]$pam <- pam_res$clustering
      sil <- silhouette(pam_res)
      fviz_silhouette(sil, main = paste("Reaction Norm PAM Silhouette (", method, ")", sep = ""))
      
      # DBSCAN clustering
      db_res <- dbscan(dmat, eps = eps)
      rn_cluster_results[[method]]$dbscan <- db_res$cluster
      mds <- cmdscale(dmat)
      plot(mds, col = db_res$cluster + 1, pch = 19,
           main = paste("Reaction Norm DBSCAN (MDS) (", method, ")", sep = ""),
           xlab = "MDS Dimension 1", ylab = "MDS Dimension 2")
      legend("topright", legend = unique(db_res$cluster),
             col = unique(db_res$cluster) + 1, pch = 19)
    }
    
    ###############################
    ## 2. Score Clustering
    ###############################
    # For each genotype form, use prepare_score_vector to extract the raw scores.
    # Each element in scores_list is a data frame with two columns (Score, Type).
    # We'll loop over genotype_forms and then combine by stacking rows.
    score_matrices_list = lapply(genotype_forms, function(form) {
      score_vectors = lapply(scores_list, function(x) {
        prepare_score_vector(x, form)
      })
      common_length = min(sapply(score_vectors, length))
      as.data.frame(sapply(score_vectors, function(v) as.numeric(v[1:common_length])))
    })
    # Combine all genotype forms by stacking rows.
    score_matrix_combined = do.call(rbind, score_matrices_list)
    # Set column names from scores_list names.
    score_names = names(scores_list)
    colnames(score_matrix_combined) <- score_names
    score_geno_names = rownames(score_matrix_combined)
    
    score_cluster_results = list()
    for (method in dist_methods) {
      if (method == "euclidean") {
        dmat <- dist(score_matrix_combined, method = "euclidean")
      } else if (method == "manhattan") {
        dmat <- dist(score_matrix_combined, method = "manhattan")
      } else if (method == "cosine") {
        dmat <- cosine_distance(score_matrix_combined)
      } else if (method == "correlation") {
        dmat <- correlation_distance(score_matrix_combined)
      } else {
        stop("Unknown method: ", method)
      }
      
      score_cluster_results[[method]] <- list(distance = dmat)
      
      # Plot heatmap of the distance matrix for scores.
      heatmap(as.matrix(dmat),
              main = paste("Score Heatmap (", method, ")", sep = ""),
              col = heat.colors(256), scale = "none",
              labRow = score_geno_names, labCol = score_names)
      
      # Hierarchical clustering (complete linkage) on scores.
      hc <- hclust(dmat, method = "complete")
      score_cluster_results[[method]]$hclust <- hc
      plot(hc, labels = score_geno_names,
           main = paste("Score hclust (", method, ")", sep = ""),
           xlab = "", sub = "")
      
      # PAM clustering on scores.
      pam_res <- pam(dmat, k = n_cluster)
      score_cluster_results[[method]]$pam <- pam_res$clustering
      sil <- silhouette(pam_res)
      fviz_silhouette(sil, main = paste("Score PAM Silhouette (", method, ")", sep = ""))
      
      # DBSCAN clustering on scores.
      db_res <- dbscan(dmat, eps = eps)
      score_cluster_results[[method]]$dbscan <- db_res$cluster
      mds <- cmdscale(dmat)
      plot(mds, col = db_res$cluster + 1, pch = 19,
           main = paste("Score DBSCAN (MDS) (", method, ")", sep = ""),
           xlab = "MDS Dimension 1", ylab = "MDS Dimension 2")
      legend("topright", legend = unique(db_res$cluster),
             col = unique(db_res$cluster) + 1, pch = 19)
    }
    
  }, error = function(e) {
    message("An error occurred: ", e$message)
  }, finally = {
    dev.off()
  })
  
  # Return a list with results for both parts.
  return(list(
    reaction_norm_clustering = rn_cluster_results,
    score_clustering         = score_cluster_results
  ))
}

########################### Main analysis



scores_list = list(
  CV_t = CV_t,
  RN= RN,
  RNN = RNN,
  D_slope = D_slope,
  RC = RC,
  #CVm = CVm,comparative
  #CVmd = CVmd,comparative
  gPi = gPi,
  PPF = PPF,
  PPi = PPi,
  PImd = PImd,
  PILSM = PILSM,
  RTR = RTR,
  PIR = PIR,
  RDPI = RDPI,
  #RDPI_mean = RDPI_mean, comparative 
  ESPI = ESPI,
  ESPIID = ESPIID,
  
  PSI = PSI,
  RPI = RPI, 
  PQ = PQ,
  PR = PR, 
  NRW = NRW,
  ESP = ESP, 
  #PD = PD, control stress
  #FPI = FPI,control stress
  #TPS = TSP,control stress
  CEV = CEV,
  PFI = PFI,
  APC = APC,
  SI = SI,
  RSI = RSI,
  EVS = EVS
  #MVPi = MVPi, multivariate
  #Plasticity_Ratio = Plasticity_Ratio
)

scores_list = list(
  CV_t = CV_t,
  RN= RN,
  RNN = RNN,
  D_slope = D_slope,
  RC = RC,
 
  gPi = gPi,
  PPF = PPF,
  PPi = PPi,
  PImd = PImd,
  PILSM = PILSM,
  RTR = RTR,
  PIR = PIR,
  RDPI = RDPI,
  
  ESPI = ESPI,
  ESPIID = ESPIID,
  PSI = PSI,
  RPI = RPI, 
  PQ = PQ,
  PR = PR, 
  NRW = NRW,
  ESP = ESP, 
  
  CEV = CEV,
  PFI = PFI,
  APC = APC,
  SI = SI,
  RSI = RSI,
  EVS = EVS
 
)

linear_matrix=do.call(rbind,linear)
gaussian_matrix=do.call(rbind,gaussian)
sinusoidal_matrix=do.call(rbind,sinusoidal)
wave_matrix=do.call(rbind,wave)

rn_matrix_list = list(
  linear = linear_matrix,    
  gaussian = gaussian_matrix, 
  sinusoidal = sinusoidal_matrix,
  wave = wave_matrix
)

# genotype forms 
genotype_forms = c("linear", "gaussian", "sinusoidal", "wave")
#correlation and clustering methods
corr_methods = c("kendall")
clustering_methods = c("hclust")

# ------------------------------
# 1. Determine Optimal n_clusters for SEPARATED Analysis
# ------------------------------
max_mean_sep <- -Inf
best_n_clusters_sep <- NA

for (i in 2:10) {
  res <- analyze_scores_by_genotype_correlation(
    scores_list, 
    genotype_forms = genotype_forms, 
    n_cluster = i, 
    rn_matrix_list = rn_matrix_list,
    eps = 0.1
  )
  # Extract silhouette scores only from the separated analyses
  sep_sils <- unlist(res$silhouette_scores[genotype_forms])
  mean_sep <- mean(sep_sils, na.rm = TRUE)
  cat("Separated analysis: Mean silhouette for n_clusters =", i, ":", mean_sep, "\n")
  if (mean_sep > max_mean_sep) {
    max_mean_sep <- mean_sep
    best_n_clusters_sep <- i
  }
}
cat("Best n_clusters for SEPARATED analysis:", best_n_clusters_sep, "\n")

res <- analyze_scores_by_genotype_correlation(
  scores_list, 
  genotype_forms = genotype_forms, 
  n_cluster = best_n_clusters_sep, 
  rn_matrix_list = rn_matrix_list,
  eps = 0.1)

cat("### ARI Matrices Across Genotype Forms\n\n")
for (method in names(res$matching_index_matrix_genotype_forms)) {
  cat("#### ", method, "\n\n")
  
  # Extract the observed ARI matrix and the p-value matrix for the current method
  ari_obs <- res$matching_index_matrix_genotype_forms[[method]]$observed
  ari_p   <- res$matching_index_matrix_genotype_forms[[method]]$p_value
  
  cat("**Observed ARI Matrix:**\n")
  print(kable(ari_obs, digits = 3, caption = paste("Observed ARI Matrix for", method)))
  
  cat("\n**Permutation Test p-value Matrix:**\n")
  print(kable(format(ari_p, scientific = TRUE, digits = 6), caption = paste("ARI p-value Matrix for", method)))
  
  cat("\n\n")
}

# ------------------------------
# 2. Determine Optimal n_clusters for COMBINED Analysis
# ------------------------------
max_mean_comb <- -Inf
best_n_clusters_comb <- NA

for (i in 2:10) {
  res <- analyze_scores_by_genotype_correlation(
    scores_list, 
    genotype_forms = genotype_forms, 
    n_cluster = i, 
    rn_matrix_list = rn_matrix_list,
    eps = 0.1
  )
  # Extract silhouette scores only from the combined analysis:
  if (!is.null(res$silhouette_scores[["combined"]])) {
    comb_sils <- unlist(res$silhouette_scores[["combined"]])
    mean_comb <- mean(comb_sils, na.rm = TRUE)
    cat("Combined analysis: Mean silhouette for n_clusters =", i, ":", mean_comb, "\n")
    if (mean_comb > max_mean_comb) {
      max_mean_comb <- mean_comb
      best_n_clusters_comb <- i
    }
  } else {
    cat("No combined analysis results for n_clusters =", i, "\n")
  }
}
cat("Best n_clusters for COMBINED analysis:", best_n_clusters_comb, "\n")



complete_res <- analyze_scores_by_genotype_correlation(
  scores_list, 
  genotype_forms = genotype_forms, 
  n_cluster = best_n_clusters_sep,  
  rn_matrix_list = rn_matrix_list,
  eps = 0.1
)

# --- Heatmaps ---
pdf("heatmaps_all_best_clustering.pdf", width = 10, height = 10)
# Separated genotype forms (using Kendall distance)
for (form in genotype_forms) {
  best_distance_matrix <- complete_res$correlation_distances[[form]][["kendall"]]
  similarity_matrix <- 1 / (1 + as.matrix(best_distance_matrix))
  pheatmap(similarity_matrix,
           clustering_distance_rows = as.dist(best_distance_matrix),
           clustering_distance_cols = as.dist(best_distance_matrix),
           clustering_method = "complete",
           cutree_rows = best_n_clusters_sep,
           cutree_cols = best_n_clusters_sep,
           main = paste("Heatmap - Genotype Form", form))
}
# Combined analysis heatmap (using Kendall distance and best_n_clusters_comb)
if (!is.null(complete_res$correlation_distances[["combined"]]) &&
    !is.null(complete_res$correlation_distances[["combined"]][["kendall"]])) {
  best_distance_matrix <- complete_res$correlation_distances[["combined"]][["kendall"]]
  similarity_matrix <- 1 / (1 + as.matrix(best_distance_matrix))
  pheatmap(similarity_matrix,
           clustering_distance_rows = as.dist(best_distance_matrix),
           clustering_distance_cols = as.dist(best_distance_matrix),
           clustering_method = "complete",
           cutree_rows = best_n_clusters_comb,
           cutree_cols = best_n_clusters_comb,
           main = "Heatmap - Combined Analysis")
} else {
  cat("Combined analysis Kendall distance matrix is missing.\n")
}
dev.off()

# --- Silhouette Plots ---
pdf("silhouettes_all_best_clustering.pdf", width = 10, height = 10)
# Separated genotype forms
for (form in genotype_forms) {
  best_distance_matrix <- complete_res$correlation_distances[[form]][["kendall"]]
  row_clusters <- cutree(hclust(as.dist(best_distance_matrix), method = "complete"), k = best_n_clusters_sep)
  sil_scores <- silhouette(row_clusters, as.dist(best_distance_matrix))
  rownames(sil_scores) <- names(scores_list)
  plot(sil_scores, main = paste("Silhouette - Genotype Form", form))
}
# Combined analysis silhouette using best_n_clusters_comb
if (!is.null(complete_res$correlation_distances[["combined"]]) &&
    !is.null(complete_res$correlation_distances[["combined"]][["kendall"]])) {
  best_distance_matrix <- complete_res$correlation_distances[["combined"]][["kendall"]]
  row_clusters <- cutree(hclust(as.dist(best_distance_matrix), method = "complete"), k = best_n_clusters_comb)
  sil_scores <- silhouette(row_clusters, as.dist(best_distance_matrix))
  rownames(sil_scores) <- names(scores_list)
  plot(sil_scores, main = "Silhouette - Combined Analysis")
} else {
  cat("Combined analysis Kendall distance matrix is missing.\n")
}
dev.off()

# --- MDS Plots ---
pdf("mds_all_best_clustering.pdf", width = 10, height = 10)
# Separated genotype forms
for (form in genotype_forms) {
  best_distance_matrix <- complete_res$correlation_distances[[form]][["kendall"]]
  mds_results <- cmdscale(as.dist(best_distance_matrix))
  row_clusters <- cutree(hclust(as.dist(best_distance_matrix), method = "complete"), k = best_n_clusters_sep)
  plot(mds_results, col = row_clusters, pch = 19, xlab = "MDS1", ylab = "MDS2",
       main = paste("MDS - Genotype Form", form))
}
# Combined analysis MDS using best_n_clusters_comb
if (!is.null(complete_res$correlation_distances[["combined"]]) &&
    !is.null(complete_res$correlation_distances[["combined"]][["kendall"]])) {
  best_distance_matrix <- complete_res$correlation_distances[["combined"]][["kendall"]]
  mds_results <- cmdscale(as.dist(best_distance_matrix))
  row_clusters <- cutree(hclust(as.dist(best_distance_matrix), method = "complete"), k = best_n_clusters_comb)
  plot(mds_results, col = row_clusters, pch = 19, xlab = "MDS1", ylab = "MDS2",
       main = "MDS - Combined Analysis")
} else {
  cat("Combined analysis Kendall distance matrix is missing.\n")
}
dev.off()

# --- Dendrograms ---
pdf("dendrograms_best_clustering.pdf", width = 10, height = 10)
# Separated genotype forms (using Pearson distance for dendrograms)
for (form in genotype_forms) {
  best_distance_matrix <- complete_res$correlation_distances[[form]][["kendall"]]
  hc <- hclust(as.dist(best_distance_matrix), method = "complete")
  plot(hc, labels = names(scores_list),
       main = paste("Dendrogram - Genotype Form", form))
}
# Combined analysis dendrogram (using Pearson distance)
if (!is.null(complete_res$correlation_distances[["combined"]]) &&
    !is.null(complete_res$correlation_distances[["combined"]][["kendall"]])) {
  best_distance_matrix <- complete_res$correlation_distances[["combined"]][["kendall"]]
  hc <- hclust(as.dist(best_distance_matrix), method = "complete")
  plot(hc, labels = names(scores_list),
       main = "Dendrogram - Combined Analysis")
} else {
  cat("Combined analysis Pearson distance matrix is missing for dendrograms.\n")
}
dev.off()

##################################### Co-occurrence analysis across clustering methods and correlation methods
n_scores = length(scores_list)
co_occurrence_matrix = matrix(0, nrow = n_scores, ncol = n_scores)
rownames(co_occurrence_matrix) = names(scores_list)
colnames(co_occurrence_matrix) = names(scores_list)

# Loop over correlation methods, clustering methods and genotype forms
for (cor_method in corr_methods) {
  for (clust_method in clustering_methods) {
    for (form in genotype_forms) {
      clustering = complete_res$clustering_results[[form]][[clust_method]][[cor_method]]
      for (i in 1:(n_scores - 1)) {
        for (j in (i + 1):n_scores) {
          if (!is.na(clustering[i]) && !is.na(clustering[j]) && clustering[i] == clustering[j]) {
            co_occurrence_matrix[i, j] = co_occurrence_matrix[i, j] + 1
            co_occurrence_matrix[j, i] = co_occurrence_matrix[j, i] + 1
          }
        }
      }
    }
  }
}
pdf("co-occurrence_scores.pdf")
co_occurrence_df = as.data.frame(as.table(co_occurrence_matrix))
colnames(co_occurrence_df) = c("Score1", "Score2", "Co_occurrences")
co_occurrence_df = co_occurrence_df[co_occurrence_df$Score1 != co_occurrence_df$Score2, ]
top_co_occurrences = co_occurrence_df[order(-co_occurrence_df$Co_occurrences), ]
print("Top co-occurring score pairs:")
print(head(top_co_occurrences, n = 10))

distance_matrix = as.dist(1 / (1 + co_occurrence_matrix))
pheatmap(co_occurrence_matrix,
         clustering_distance_rows = distance_matrix,
         clustering_distance_cols = distance_matrix,
         clustering_method = "complete",
         main = "Clustering of Scores Based on Co-occurrence (Similarity)")
dev.off()



###################### Covariance analysis across genotype forms and scores 


analysis2 = function(rn_matrix_list, scores_list, plot = TRUE, 
                     do_combined = TRUE, do_cov_combined = TRUE) {
  # Use names of rn_matrix_list as genotype forms
  genotype_forms = names(rn_matrix_list)
  
  # Storage containers for separate analyses
  correlation_distances_per_geno = list()
  rn_comparison = list()
  
  pdf(file = "correlation_plots2.pdf", width = 12, height = 12)
  
  # --- Separate Analyses (for each genotype form) ---
  for (form in genotype_forms) {
    # Prepare score vectors and build score matrix for the current form
    score_vectors = lapply(scores_list, function(df) {
      prepare_score_vector(df, form)
    })
    common_length = min(sapply(score_vectors, length))
    score_matrix = as.data.frame(sapply(score_vectors, function(v) as.numeric(v[1:common_length])))
    score_names = names(scores_list)
    
    # Compute plasticity score correlation matrix (using Kendall)
    corr_geno = as.matrix(cor(t(score_matrix), method = "kendall"))
    correlation_distances_per_geno[[form]] = corr_geno
    
    # Reaction norm data for the current form
    rn_matrix = rn_matrix_list[[form]]
    # Force a 2D matrix for covariance computation
    rn_cor = as.matrix(cov(t(rn_matrix), method = "kendall"))
    
    # Plot heatmaps for reaction norm covariance and plasticity score correlation
    par(mfrow = c(1, 2), mar = c(4,4,4,2))
    heatmap(rn_cor, symm = TRUE,
            main = paste("Reaction Norm Covariance", form),
            margins = c(10,10))
    heatmap(corr_geno, scale = "none", symm = TRUE,
            main = paste("Plasticity Score Correlation", form),
            margins = c(10,10))
    
    # Compute distance matrices for Mantel test (using Euclidean distances)
    rn_dist_l2 = dist(rn_matrix, method = "canberra")
    score_dist_l2 = dist(score_matrix, method = "canberra")
    
    # Run Mantel test on the Euclidean distances
    mantel_out = vegan::mantel(rn_dist_l2, score_dist_l2, permutations = 999)
    rn_comparison[[form]] = mantel_out
    
    # Hierarchical clustering based on correlation-based distances
    rn_dist = as.dist(1 - rn_cor)
    geno_dist = as.dist(1 - corr_geno)
    
    rn_hc = hclust(rn_dist)
    geno_hc = hclust(geno_dist)
    rn_dist_matrix = as.matrix(rn_dist)
    geno_dist_matrix = as.matrix(geno_dist)
    
    par(mfrow = c(1,2), mar = c(4,4,4,2))
    heatmap(rn_dist_matrix, symm = TRUE, scale = "none",
            main = paste("Reaction Norm Distance", form),
            margins = c(10,10))
    heatmap(geno_dist_matrix, symm = TRUE,
            main = paste("Plasticity Score Distance", form),
            margins = c(10,10))
    
    # L² (Euclidean) distance analysis with clustering
    rn_l2_hc = hclust(rn_dist_l2)
    score_l2_hc = hclust(score_dist_l2)
    rn_dist_l2_matrix = as.matrix(rn_dist_l2)
    score_dist_l2_matrix = as.matrix(score_dist_l2)
    
    par(mfrow = c(1,2), mar = c(4,4,4,2))
    heatmap(rn_dist_l2_matrix,
            Rowv = as.dendrogram(rn_l2_hc),
            Colv = as.dendrogram(rn_l2_hc),
            symm = TRUE,
            main = paste("Reaction Norm L2 Distance", form),
            margins = c(10,10))
    heatmap(score_dist_l2_matrix,
            Rowv = as.dendrogram(score_l2_hc),
            Colv = as.dendrogram(score_l2_hc),
            symm = TRUE,
            main = paste("Plasticity Score L2 Distance", form),
            margins = c(10,10))
    
    # Plot all reaction norms (overlay all curves)
    par(mfrow = c(1,1), mar = c(4,4,4,4))
    matplot(t(rn_matrix), type = "l", lty = 1,
            main = paste("All Reaction Norms", form),
            xlab = "Condition/Time", ylab = "Response",
            col = rainbow(nrow(rn_matrix)))
    legend("topright", legend = rownames(rn_matrix),
           col = rainbow(nrow(rn_matrix)), lty = 1, cex = 0.7)
    
    # Plot histograms for the separate analysis (trait values and scores)
    layout_matrix = rbind(
      c(1, 1, 2, 3, 4),
      c(1, 1, 5, 6, 7),
      c(8, 9, 10, 11, 12),
      c(13, 14, 15, 16, 17),
      c(18, 19, 20, 21, 22),
      c(23, 24, 25, 26, 27),
      c(28, 0, 0, 0, 0)
    )
    layout(layout_matrix)
    old_mar = par("mar")
    par(mar = c(4,4,4,2))
    
    hist(rn_matrix,
         breaks = 30,
         probability = TRUE,
         main = paste("Distribution of Direct Trait Values (", form, ")", sep = ""),
         xlab = "Trait Value",
         ylab = "Density",
         col = "lightblue",
         border = "white")
    lines(density(as.numeric(rn_matrix)), col = "red", lwd = 2)
    par(mar = old_mar)
    
    n_scores = ncol(score_matrix)
    for (i in 1:n_scores) {
      hist(score_matrix[[i]],
           breaks = 10,
           probability = TRUE,
           main = paste("Score:", score_names[i], "\n(", form, ")"),
           xlab = "Score",
           ylab = "Density",
           col = "lightblue",
           border = "white")
      if(length(score_matrix[[i]]) > 1 && length(unique(score_matrix[[i]])) > 1) {
        lines(density(score_matrix[[i]]), col = "red", lwd = 2)
      }
    }
    
    dev.flush()
  }  # end for separate analyses
  
  # --- Combined Analysis ---
  if (do_combined) {
    # Combine the separate reaction norm matrices into one matrix.
    combined_matrix = do.call(rbind, rn_matrix_list)
    rownames(combined_matrix) <- as.character(seq_len(nrow(combined_matrix)))
    score_matrix_combined = do.call(cbind, lapply(scores_list, function(x) as.numeric(x[,1])))
    rownames(score_matrix_combined) <- as.character(seq_len(nrow(score_matrix_combined)))
    score_matrix_combined_corr = as.matrix(cor(t(score_matrix_combined), method = "kendall"))
    
    if (do_cov_combined) {
      # Compute covariance on the combined matrix.
      rn_cor_combined = as.matrix(cov(t(combined_matrix), method = "kendall"))
      
      par(mfrow = c(1,2), mar = c(4,4,4,2))
      heatmap(rn_cor_combined, symm = TRUE,
              main = "Reaction Norm Covariance",
              margins = c(10,10))
      # Also compute a correlation matrix for the combined matrix:
      corr_combined = as.matrix(cor(t(combined_matrix), method = "kendall"))
      
      heatmap(corr_combined, scale = "none", symm = TRUE,
              main = "Reaction Norm Correlation",
              margins = c(10,10))
      
      heatmap(score_matrix_combined_corr, scale = "none", symm = TRUE,
              main = "Plasticity Score Correlation",
              margins = c(10,10))
      score_mat_dist = 1 - score_matrix_combined_corr
      rn_dist_combined = as.dist(1 - rn_cor_combined)
      geno_dist_combined = as.dist(1 - corr_combined)
      mantel_out_combined = vegan::mantel(geno_dist_combined, score_mat_dist, permutations = 999)
      rn_comparison[["combined"]] = mantel_out_combined
      
      rn_hc_combined = hclust(rn_dist_combined)
      geno_hc_combined = hclust(geno_dist_combined)
      rn_dist_matrix_combined = as.matrix(rn_dist_combined)
      geno_dist_matrix_combined = as.matrix(geno_dist_combined)
      
      rn_dist_l2 = dist(combined_matrix, method = "canberra")
      score_dist_l2 = dist(score_matrix_combined, method = "canberra")
      # Convert dist objects to matrices before heatmapping
      rn_dist_l2_matrix = as.matrix(rn_dist_l2)
      score_dist_l2_matrix = as.matrix(score_dist_l2)
      
      par(mfrow = c(1,2), mar = c(4,4,4,2))
      heatmap(rn_dist_l2_matrix, symm = TRUE, scale = "none",
              main = "Reaction Norm Distance",
              margins = c(10,10))
      heatmap(score_dist_l2_matrix, symm = TRUE,
              main = "Plasticity Score Distance",
              margins = c(10,10))
      
      rn_dist_l2_combined = dist(combined_matrix, method = "euclidean")
      rn_l2_hc_combined = hclust(rn_dist_l2_combined)
      rn_dist_l2_matrix_combined = as.matrix(rn_dist_l2_combined)
      
      par(mfrow = c(1,2), mar = c(4,4,4,2))
      heatmap(rn_dist_l2_matrix_combined,
              Rowv = as.dendrogram(rn_l2_hc_combined),
              Colv = as.dendrogram(rn_l2_hc_combined),
              symm = TRUE,
              main = "Combined Reaction Norm L2 Distance",
              margins = c(10,10))
    }
    
    layout_matrix = rbind(
      c(1, 1, 2, 3, 4),
      c(1, 1, 5, 6, 7),
      c(8, 9, 10, 11, 12),
      c(13, 14, 15, 16, 17),
      c(18, 19, 20, 21, 22),
      c(23, 24, 25, 26, 27),
      c(28, 0, 0, 0, 0)
    )
    layout(layout_matrix)
    old_mar = par("mar")
    par(mar = c(4,4,4,2))
    
    # Big histogram: Distribution of direct trait values (from the combined matrix)
    hist(as.numeric(combined_matrix),
         breaks = 50,
         probability = TRUE,
         main = "Distribution of Direct Trait Values (Combined)",
         xlab = "Trait Value",
         ylab = "Density",
         col = "lightblue",
         border = "white")
    lines(density(as.numeric(combined_matrix)), col = "red", lwd = 2)
    
    # Now plot the small histograms for each plasticity score.
    n_scores = ncol(score_matrix_combined)
    
    score_names_combined = colnames(score_matrix_combined)
    for (i in 1:n_scores) {
      hist(score_matrix_combined[,i],
           breaks = 30,
           probability = TRUE,
           main = paste("Score:", ifelse(is.null(score_names_combined),
                                         paste("Score", i), score_names_combined[i]), "\n(Combined)"),
           xlab = "Score",
           ylab = "Density",
           col = "lightblue",
           border = "white")
      if(length(score_matrix_combined[[i]]) > 1 &&
         length(unique(score_matrix_combined[[i]])) > 1) {
        lines(density(score_matrix_combined[[i]]), col = "red", lwd = 2)
      }
    }
    par(mar = old_mar)
  }
  
  dev.off()
  
  return(list(correlation_distances_per_geno = correlation_distances_per_geno,
              reaction_norm_comparison = rn_comparison))
}


# Example call:
result = analysis2(rn_matrix_list, scores_list, plot = TRUE, 
                   do_combined = TRUE, do_cov_combined = TRUE)

mantel_results <- data.frame(
  Genotype_Form = character(),
  Mantel_Statistic = numeric(),
  p_value = character(),
  stringsAsFactors = FALSE
)

# Loop over each genotype form in the reaction_norm_comparison list
for(form in names(result$reaction_norm_comparison)) {
  m_test <- result$reaction_norm_comparison[[form]]
  
  # Extract the Mantel statistic and p-value. The p-value is in m_test$signif
  # Use format() to display small values in scientific notation.
  mantel_results <- rbind(mantel_results, data.frame(
    Genotype_Form = form,
    Mantel_Statistic = round(m_test$statistic, 3),
    p_value = format(m_test$signif, scientific = TRUE, digits = 3),
    stringsAsFactors = FALSE
  ))
}

# Print the table
print(kable(mantel_results, caption = "Mantel Test Results by Genotype Form"))

####################### PCA analysis of Scores - I am not sure how sensical this is 


pca_results = list()


pdf("PCA_plasticity_scores.pdf", width = 8, height = 8)


for (form in genotype_forms) {
  message("Performing PCA for form: ", form)
  
  score_vectors = lapply(scores_list, function(df) {
    prepare_score_vector(df, form)
  })
  
  common_length = min(sapply(score_vectors, length))
  score_matrix = as.data.frame(sapply(score_vectors, function(v) as.numeric(v[1:common_length])))
  rownames(score_matrix) = 1:nrow(score_matrix)
  
  # --- Perform PCA ---
  # Each row corresponds to a genotype, and each column to a different plasticity score.
  pca_res = prcomp(score_matrix, scale. = TRUE)
  pca_results[[form]] = pca_res
  
  # --- Plot PCA Biplot ---
  # fviz_pca_biplot() plots both the individuals (genotypes) and variable loadings.
  p = fviz_pca_biplot(pca_res,
                       geom.ind = "point",    # show points for individuals
                       pointshape = 21,
                       pointsize = 2,
                       fill.ind = "lightblue",  # color for individuals
                       col.var = "red",         # color for variables (scores)
                       repel = TRUE,            # avoid label overlap
                       title = paste("PCA Biplot for", form, "Plasticity Scores"))
  
  print(p)
}

dev.off()



#################### Procrustes analysis 


analyze_procrustes_all = function(rn_matrix_list, scores_list , scale = TRUE) {

  genotype_forms = names(rn_matrix_list)
  pdf(file = "procrustes_plots.pdf", width = 8, height = 8)
  
  # List to store Procrustes results
  proc_results = list()
  for (form in genotype_forms) {
    message("Analyzing form: ", form)
    
    # Extract Trait Data 
    trait_data = rn_matrix_list[[form]]
    
    # Extract Plasticity Score Data 
    score_vectors = lapply(scores_list, function(df) {
      prepare_score_vector(df, form)
    })
    common_length = min(sapply(score_vectors, length))
    plasticity_scores = as.data.frame(sapply(score_vectors, function(v) as.numeric(v[1:common_length])))
    
    if (is.null(rownames(trait_data))) {
      rownames(trait_data) = 1:nrow(trait_data)
    }
    rownames(plasticity_scores) = rownames(trait_data)[1:common_length]
    
    pca_trait = prcomp(trait_data, scale. = scale)
    trait_scores = pca_trait$x[, 1:2]  # Use first two PCs
    
    pca_plasticity = prcomp(plasticity_scores, scale. = scale)
    plasticity_scores_pc = pca_plasticity$x[, 1:2]  # Use first two PCs
    
    proc = vegan::procrustes(trait_scores, plasticity_scores_pc)
    proc_results[[form]] = proc
    target = proc$Yrot  # The rotated configuration (aligned to the target)
    # Alternatively, the target configuration (the reference) is stored in proc$X
    reference = proc$X
    
    # Plot the reference points
    plot(reference, col = "blue", pch = 16, xlim = range(c(reference[,1], target[,1])), 
         ylim = range(c(reference[,2], target[,2])),
         main = "Procrustes Analysis: Separate Configurations",
         xlab = "PC1", ylab = "PC2")
    # Add the target configuration in a different symbol and color
    points(target, col = "red", pch = 17)
    legend("topright", legend = c("Reference (Trait Data)", "Transformed (Plasticity Scores)"),
           col = c("blue", "red"), pch = c(16, 17))
    mtext(paste("Procrustes sum-of-squares:", round(proc$ss, 3)), side = 3, line = 0.5)
    
    prot = vegan::protest(trait_scores, plasticity_scores_pc, permutations = 999)
    mtext(paste("Protest p-value:", round(prot$signif, 3)), side = 1, line = 0.5)
    
  }  
  
  dev.off()
  
  return(proc_results)
}

proc_results = analyze_procrustes_all(rn_matrix_list, scores_list, scale = TRUE)


