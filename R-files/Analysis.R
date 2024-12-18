############################ correlation
library(factoextra)
if (!exists("data_loaded") || !data_loaded) {
  source("~/CRC 1622 - Z2/R-files/head.R")
  data_loaded = T
}

# Initialize matrices to store results for Pearson, Spearman, and Kendall correlations
correlation_distance_matrices_pearson = list()
correlation_distance_matrices_spearman = list()
correlation_distance_matrices_kendall = list()

# Define dataset types
dataset_types = list(
  datasets, datasets_log_skewed, datasets_log_log, 
  datasets_gamma_skewed, datasets_gamma_gamma, 
  datasets_weibull_skewed, datasets_weibull_weibull, 
  datasets_log_gamma, datasets_log_weibull, datasets_gamma_weibull
)

# Loop through dataset types
for (dataset_type_index in 1:length(dataset_types)) {
  
  # Initialize lists for storing correlation distances for each trait
  pearson_distances_per_trait = list()
  spearman_distances_per_trait = list()
  kendall_distances_per_trait = list()
  
  # Loop through each trait (3 traits)
  for (trait_index in 1:3) {
    
    # Create a data frame to store scores for all datasets for the current trait
    scores_df = data.frame(
      CV_t = numeric(n_datasets),
      RNS = numeric(n_datasets),
      D_slope = numeric(n_datasets),
      RC = numeric(n_datasets),
      CVm = numeric(n_datasets),
      CVmd = numeric(n_datasets),
      Pi_adj_mean = numeric(n_datasets),
      Pi = numeric(n_datasets),
      PImd = numeric(n_datasets),
      PILSM = numeric(n_datasets),
      RTR = numeric(n_datasets),
      PIR = numeric(n_datasets)
    )
    
    # Loop over the 20 datasets for the current dataset type
    for (i in 1:n_datasets) {
      scores_df$CV_t[i] = as.numeric(CV_t[[dataset_type_index]][[i]][trait_index])
      scores_df$RNS[i] = as.numeric(RNS[[dataset_type_index]][[i]][trait_index])
      scores_df$D_slope[i] = as.numeric(D_slope[[dataset_type_index]][[i]][trait_index])
      scores_df$RC[i] = as.numeric(RC[[dataset_type_index]][[i]][trait_index])
      scores_df$CVm[i] = as.numeric(CVm[[dataset_type_index]][[i]][trait_index])
      scores_df$CVmd[i] = as.numeric(CVmd[[dataset_type_index]][[i]][trait_index])
      scores_df$Pi_adj_mean[i] = as.numeric(Pi_adj_mean[[dataset_type_index]][[i]][trait_index])
      scores_df$Pi[i] = as.numeric(Pi[[dataset_type_index]][[i]][trait_index])
      scores_df$PImd[i] = as.numeric(PImd[[dataset_type_index]][[i]][trait_index])
      scores_df$PILSM[i] = as.numeric(PILSM[[dataset_type_index]][[i]][trait_index][[1]])
      scores_df$RTR[i] = as.numeric(RTR[[dataset_type_index]][[i]][trait_index])
      scores_df$PIR[i] = as.numeric(PIR[[dataset_type_index]][[i]][trait_index])
    }
    
    # Calculate Pearson, Spearman, and Kendall correlation distances using get_dist()
    pearson_distances_per_trait[[trait_index]] = as.data.frame(as.matrix(get_dist(t(scores_df), method = "pearson")))
    spearman_distances_per_trait[[trait_index]] = as.data.frame(as.matrix(get_dist(t(scores_df), method = "spearman")))
    kendall_distances_per_trait[[trait_index]] = as.data.frame(as.matrix(get_dist(t(scores_df), method = "kendall")))
  }
  
  # Store the correlation distances for all 3 traits for the current dataset type
  correlation_distance_matrices_pearson[[dataset_type_index]] = pearson_distances_per_trait
  correlation_distance_matrices_spearman[[dataset_type_index]] = spearman_distances_per_trait
  correlation_distance_matrices_kendall[[dataset_type_index]] = kendall_distances_per_trait
}

# Combine the matrices into correlation_matrices_all_datasets
correlation_matrices_all_datasets = list()
for (dataset_type_index in 1:length(dataset_types)) {
  correlation_matrices_all_datasets[[dataset_type_index]] = list(
    pearson = correlation_distance_matrices_pearson[[dataset_type_index]],
    spearman = correlation_distance_matrices_spearman[[dataset_type_index]],
    kendall = correlation_distance_matrices_kendall[[dataset_type_index]]
  )
}

############################### clustering 
library(dbscan)
library(cluster)
library(mclust)

# Initialize data frames for all cluster assignments (hclust, DBSCAN, PAM) for each correlation method
hclust_df_pearson = data.frame(matrix(nrow = 30, ncol = 12))
hclust_df_spearman = data.frame(matrix(nrow = 30, ncol = 12))
hclust_df_kendall = data.frame(matrix(nrow = 30, ncol = 12))

dbscan_df_pearson = data.frame(matrix(nrow = 30, ncol = 12))
dbscan_df_spearman = data.frame(matrix(nrow = 30, ncol = 12))
dbscan_df_kendall = data.frame(matrix(nrow = 30, ncol = 12))

pam_df_pearson = data.frame(matrix(nrow = 30, ncol = 12))
pam_df_spearman = data.frame(matrix(nrow = 30, ncol = 12))
pam_df_kendall = data.frame(matrix(nrow = 30, ncol = 12))

# Define column names and row names
column_names = c("CV_t", "RNS", "D_slope", "RC", "CVm", "CVmd", "Pi_adj_mean", "Pi", "PImd", "PILSM", "RTR", "PIR")
row_names = paste0("Type", rep(1:10, each = 3), "_Trait", rep(1:3, 10))

# Assign column and row names to the data frames
for (df in list(hclust_df_pearson, hclust_df_spearman, hclust_df_kendall,
                dbscan_df_pearson, dbscan_df_spearman, dbscan_df_kendall,
                pam_df_pearson, pam_df_spearman, pam_df_kendall)) {
  colnames(df) = column_names
  rownames(df) = row_names
}

# Initialize vectors for silhouette scores
silhouette_scores_hclust_pearson = numeric(30)
silhouette_scores_hclust_spearman = numeric(30)
silhouette_scores_hclust_kendall = numeric(30)

silhouette_scores_pam_pearson = numeric(30)
silhouette_scores_pam_spearman = numeric(30)
silhouette_scores_pam_kendall = numeric(30)

silhouette_scores_dbscan_pearson = numeric(30)
silhouette_scores_dbscan_spearman = numeric(30)
silhouette_scores_dbscan_kendall = numeric(30)

# Initialize a matrix for ARI comparisons between clustering methods
ari_between_methods = matrix(0, nrow = 30, ncol = 3)  # Store ARI for hclust vs PAM, hclust vs DBSCAN, PAM vs DBSCAN

# Perform clustering and store results
for (dataset_type_index in 1:length(correlation_matrices_all_datasets)) {
  
  # Loop through each trait (3 traits)
  for (trait_index in 1:3) {
    row_idx = (dataset_type_index - 1) * 3 + trait_index
    
    # -------------------- Hierarchical Clustering --------------------
    # Pearson
    dist_pearson = as.dist(correlation_distance_matrices_pearson[[dataset_type_index]][[trait_index]])
    clusters_hc_pearson = cutree(hclust(dist_pearson, method = "complete"), k = 3)
    hclust_df_pearson[row_idx, ] = clusters_hc_pearson
    silhouette_scores_hclust_pearson[row_idx] = mean(silhouette(clusters_hc_pearson, dist_pearson)[, "sil_width"])
    
    # Spearman
    dist_spearman = as.dist(correlation_distance_matrices_spearman[[dataset_type_index]][[trait_index]])
    clusters_hc_spearman = cutree(hclust(dist_spearman, method = "complete"), k = 3)
    hclust_df_spearman[row_idx, ] = clusters_hc_spearman
    silhouette_scores_hclust_spearman[row_idx] = mean(silhouette(clusters_hc_spearman, dist_spearman)[, "sil_width"])
    
    # Kendall
    dist_kendall = as.dist(correlation_distance_matrices_kendall[[dataset_type_index]][[trait_index]])
    clusters_hc_kendall = cutree(hclust(dist_kendall, method = "complete"), k = 3)
    hclust_df_kendall[row_idx, ] = clusters_hc_kendall
    silhouette_scores_hclust_kendall[row_idx] = mean(silhouette(clusters_hc_kendall, dist_kendall)[, "sil_width"])
    
    # -------------------- DBSCAN Clustering --------------------
    # Pearson
    clusters_dbscan_pearson = dbscan(dist_pearson, eps = 0.1)$cluster
    dbscan_df_pearson[row_idx, ] = clusters_dbscan_pearson
    silhouette_scores_dbscan_pearson[row_idx] = if (length(unique(clusters_dbscan_pearson)) > 1) {
      mean(silhouette(clusters_dbscan_pearson, dist_pearson)[, "sil_width"])
    } else { NA }
    
    # Spearman
    clusters_dbscan_spearman = dbscan(dist_spearman, eps = 0.1)$cluster
    dbscan_df_spearman[row_idx, ] = clusters_dbscan_spearman
    silhouette_scores_dbscan_spearman[row_idx] = if (length(unique(clusters_dbscan_spearman)) > 1) {
      mean(silhouette(clusters_dbscan_spearman, dist_spearman)[, "sil_width"])
    } else { NA }
    
    # Kendall
    clusters_dbscan_kendall = dbscan(dist_kendall, eps = 0.1)$cluster
    dbscan_df_kendall[row_idx, ] = clusters_dbscan_kendall
    silhouette_scores_dbscan_kendall[row_idx] = if (length(unique(clusters_dbscan_kendall)) > 1) {
      mean(silhouette(clusters_dbscan_kendall, dist_kendall)[, "sil_width"])
    } else { NA }
    
    # -------------------- PAM Clustering --------------------
    # Pearson
    clusters_pam_pearson = pam(dist_pearson, k = 3)$clustering
    pam_df_pearson[row_idx, ] = clusters_pam_pearson
    silhouette_scores_pam_pearson[row_idx] = mean(silhouette(clusters_pam_pearson, dist_pearson)[, "sil_width"])
    
    # Spearman
    clusters_pam_spearman = pam(dist_spearman, k = 3)$clustering
    pam_df_spearman[row_idx, ] = clusters_pam_spearman
    silhouette_scores_pam_spearman[row_idx] = mean(silhouette(clusters_pam_spearman, dist_spearman)[, "sil_width"])
    
    # Kendall
    clusters_pam_kendall = pam(dist_kendall, k = 3)$clustering
    pam_df_kendall[row_idx, ] = clusters_pam_kendall
    silhouette_scores_pam_kendall[row_idx] = mean(silhouette(clusters_pam_kendall, dist_kendall)[, "sil_width"])
    
    # -------------------- Adjusted Rand Index (Comparing Clustering Methods) --------------------
    # ARI for Pearson
    ari_between_methods[row_idx, 1] = adjustedRandIndex(clusters_hc_pearson, clusters_pam_pearson)
    ari_between_methods[row_idx, 2] = adjustedRandIndex(clusters_hc_pearson, clusters_dbscan_pearson)
    ari_between_methods[row_idx, 3] = adjustedRandIndex(clusters_pam_pearson, clusters_dbscan_pearson)
  }
}

# Create summary data frames for silhouette and matching index scores
silhouette_df = data.frame(
  Hclust_Pearson = silhouette_scores_hclust_pearson,
  Hclust_Spearman = silhouette_scores_hclust_spearman,
  Hclust_Kendall = silhouette_scores_hclust_kendall,
  PAM_Pearson = silhouette_scores_pam_pearson,
  PAM_Spearman = silhouette_scores_pam_spearman,
  PAM_Kendall = silhouette_scores_pam_kendall,
  DBSCAN_Pearson = silhouette_scores_dbscan_pearson,
  DBSCAN_Spearman = silhouette_scores_dbscan_spearman,
  DBSCAN_Kendall = silhouette_scores_dbscan_kendall
)

matching_index_df = data.frame(
  ARI_Hclust_PAM = ari_between_methods[, 1],
  ARI_Hclust_DBSCAN = ari_between_methods[, 2],
  ARI_PAM_DBSCAN = ari_between_methods[, 3]
)

# Output silhouette and matching index scores
print("Silhouette Scores for Clustering Results:")
print(silhouette_df)

print("Adjusted Rand Index (Matching Index) for Clustering Method Comparisons:")
print(matching_index_df)
