
if (!exists("data_loaded") || !data_loaded) {
  source("~/CRC 1622 - Z2/R-files/head.R")
  data_loaded = T
}

# Main function to perform the entire analysis
analyze_scores = function(scores_list, dataset_types, n_datasets = 20, n_traits = 3, n_cluster , eps = 0.1) {
  # Required libraries
  library(factoextra)
  library(dbscan)
  library(cluster)
  library(mclust)
  
  # Initialize results containers
  correlation_distances = list(
    pearson = vector("list", length(dataset_types)),
    spearman = vector("list", length(dataset_types)),
    kendall = vector("list", length(dataset_types))
  )
  for (i in seq_along(correlation_distances$pearson)) {
    correlation_distances$pearson[[i]] = vector("list", n_traits)
    correlation_distances$spearman[[i]] = vector("list", n_traits)
    correlation_distances$kendall[[i]] = vector("list", n_traits)
  }
  
  clustering_results = list(
    hclust = list(
      pearson = vector("list", length(dataset_types)),
      spearman = vector("list", length(dataset_types)),
      kendall = vector("list", length(dataset_types))
    ),
    dbscan = list(
      pearson = vector("list", length(dataset_types)),
      spearman = vector("list", length(dataset_types)),
      kendall = vector("list", length(dataset_types))
    ),
    pam = list(
      pearson = vector("list", length(dataset_types)),
      spearman = vector("list", length(dataset_types)),
      kendall = vector("list", length(dataset_types))
    )
  )
  for (cor_method in names(clustering_results$hclust)) {
    for (i in seq_along(clustering_results$hclust[[cor_method]])) {
      clustering_results$hclust[[cor_method]][[i]] = vector("list", n_traits)
      clustering_results$dbscan[[cor_method]][[i]] = vector("list", n_traits)
      clustering_results$pam[[cor_method]][[i]] = vector("list", n_traits)
    }
  }
  
  # Initialize silhouette score and ARI matrices
  silhouette_scores = list(
    hclust = list(pearson = numeric(length(dataset_types) * n_traits), spearman = numeric(length(dataset_types) * n_traits), kendall = numeric(length(dataset_types) * n_traits)),
    pam = list(pearson = numeric(length(dataset_types) * n_traits), spearman = numeric(length(dataset_types) * n_traits), kendall = numeric(length(dataset_types) * n_traits)),
    dbscan = list(pearson = numeric(length(dataset_types) * n_traits), spearman = numeric(length(dataset_types) * n_traits), kendall = numeric(length(dataset_types) * n_traits))
  )
  
  matching_index_matrix_correlation = matrix(0, nrow = length(dataset_types) * n_traits, ncol = 9,
                                             dimnames = list(skewness_type_vector,
                                                             c("ARI_Hclust_PAM_Pearson", "ARI_Hclust_DBSCAN_Pearson", "ARI_PAM_DBSCAN_Pearson",
                                                               "ARI_Hclust_PAM_Spearman", "ARI_Hclust_DBSCAN_Spearman", "ARI_PAM_DBSCAN_Spearman",
                                                               "ARI_Hclust_PAM_Kendall", "ARI_Hclust_DBSCAN_Kendall", "ARI_PAM_DBSCAN_Kendall")))
  matching_index_matrix_across_corr = matrix(0, nrow = length(dataset_types) * n_traits, ncol = 3,
                                             dimnames = list(skewness_type_vector,
                                                             c("ARI_Hclust_across", "ARI_PAM_across", "ARI_DBSCAN_across")))
  row_counter = 1
  
  # Function to calculate correlation distance matrices
  calculate_correlation_distances = function(scores_list, dataset_type_index, trait_index) {
    scores_df = data.frame(matrix(ncol = length(scores_list), nrow = n_datasets))
    colnames(scores_df) = names(scores_list)
    
    for (i in 1:n_datasets) {
      for (score_name in 1:length(scores_list)) {
        if (!is.null(scores_list[[score_name]][[dataset_type_index]][[i]][[trait_index]])) {
          scores_df[i, score_name] = as.numeric(scores_list[[score_name]][[dataset_type_index]][[i]][[trait_index]])
        } else {
          scores_df[i, score_name] = NA
        }
      }
    }
    
    list(
      pearson = as.data.frame(as.matrix(get_dist(t(scores_df), method = "pearson"))),
      spearman = as.data.frame(as.matrix(get_dist(t(scores_df), method = "spearman"))),
      kendall = as.data.frame(as.matrix(get_dist(t(scores_df), method = "kendall")))
    )
  }
  
  # Loop over dataset types and traits
  for (dataset_type_index in 1:length(dataset_types)) {
    for (trait_index in 1:n_traits) {
      correlations = calculate_correlation_distances(scores_list, dataset_type_index, trait_index)
      
      correlation_distances$pearson[[dataset_type_index]][[trait_index]] = correlations$pearson
      correlation_distances$spearman[[dataset_type_index]][[trait_index]] = correlations$spearman
      correlation_distances$kendall[[dataset_type_index]][[trait_index]] = correlations$kendall
      
      for (cor_method in c("pearson", "spearman", "kendall")) {
        dist_matrix = correlation_distances[[cor_method]][[dataset_type_index]][[trait_index]]
        
        if (!is.null(dist_matrix)) {
          clustering_results$hclust[[cor_method]][[dataset_type_index]][[trait_index]] = 
            cutree(hclust(as.dist(dist_matrix), method = "complete"), k = n_cluster)
          
          clustering_results$pam[[cor_method]][[dataset_type_index]][[trait_index]] = 
            pam(as.dist(dist_matrix), k = n_cluster)$clustering
          
          clustering_results$dbscan[[cor_method]][[dataset_type_index]][[trait_index]] = 
            dbscan(as.dist(dist_matrix), eps = eps)$cluster
          
          # Calculate silhouette scores
          silhouette_scores$hclust[[cor_method]][row_counter] = mean(silhouette(clustering_results$hclust[[cor_method]][[dataset_type_index]][[trait_index]], dist_matrix)[, "sil_width"])
          silhouette_scores$pam[[cor_method]][row_counter] = mean(silhouette(clustering_results$pam[[cor_method]][[dataset_type_index]][[trait_index]], dist_matrix)[, "sil_width"])
          
          if (length(unique(clustering_results$dbscan[[cor_method]][[dataset_type_index]][[trait_index]])) > 1) {
            silhouette_scores$dbscan[[cor_method]][row_counter] = mean(silhouette(clustering_results$dbscan[[cor_method]][[dataset_type_index]][[trait_index]], dist_matrix)[, "sil_width"])
          } else {
            silhouette_scores$dbscan[[cor_method]][row_counter] = 0
          }
        }
      }
      
      # Matching index calculations within each correlation and across correlations
      hclust_labels_pearson = as.numeric(clustering_results$hclust$pearson[[dataset_type_index]][[trait_index]])
      pam_labels_pearson = as.numeric(clustering_results$pam$pearson[[dataset_type_index]][[trait_index]])
      dbscan_labels_pearson = as.numeric(clustering_results$dbscan$pearson[[dataset_type_index]][[trait_index]])
      
      hclust_labels_spearman = as.numeric(clustering_results$hclust$spearman[[dataset_type_index]][[trait_index]])
      pam_labels_spearman = as.numeric(clustering_results$pam$spearman[[dataset_type_index]][[trait_index]])
      dbscan_labels_spearman = as.numeric(clustering_results$dbscan$spearman[[dataset_type_index]][[trait_index]])
      
      hclust_labels_kendall = as.numeric(clustering_results$hclust$kendall[[dataset_type_index]][[trait_index]])
      pam_labels_kendall = as.numeric(clustering_results$pam$kendall[[dataset_type_index]][[trait_index]])
      dbscan_labels_kendall = as.numeric(clustering_results$dbscan$kendall[[dataset_type_index]][[trait_index]])
      
      # Within-correlation ARI
      matching_index_matrix_correlation[row_counter, 1:3] = c(
        adjustedRandIndex(hclust_labels_pearson, pam_labels_pearson),
        adjustedRandIndex(hclust_labels_pearson, dbscan_labels_pearson),
        adjustedRandIndex(pam_labels_pearson, dbscan_labels_pearson)
      )
      matching_index_matrix_correlation[row_counter, 4:6] = c(
        adjustedRandIndex(hclust_labels_spearman, pam_labels_spearman),
        adjustedRandIndex(hclust_labels_spearman, dbscan_labels_spearman),
        adjustedRandIndex(pam_labels_spearman, dbscan_labels_spearman)
      )
      matching_index_matrix_correlation[row_counter, 7:9] = c(
        adjustedRandIndex(hclust_labels_kendall, pam_labels_kendall),
        adjustedRandIndex(hclust_labels_kendall, dbscan_labels_kendall),
        adjustedRandIndex(pam_labels_kendall, dbscan_labels_kendall)
      )
      
      # Across-correlation ARI
      matching_index_matrix_across_corr[row_counter, ] = c(
        adjustedRandIndex(hclust_labels_pearson, hclust_labels_spearman) +
          adjustedRandIndex(hclust_labels_pearson, hclust_labels_kendall) +
          adjustedRandIndex(hclust_labels_spearman, hclust_labels_kendall),
        adjustedRandIndex(pam_labels_pearson, pam_labels_spearman) +
          adjustedRandIndex(pam_labels_pearson, pam_labels_kendall) +
          adjustedRandIndex(pam_labels_spearman, pam_labels_kendall),
        adjustedRandIndex(dbscan_labels_pearson, dbscan_labels_spearman) +
          adjustedRandIndex(dbscan_labels_pearson, dbscan_labels_kendall) +
          adjustedRandIndex(dbscan_labels_spearman, dbscan_labels_kendall)
      )
      
      row_counter = row_counter + 1
    }
  }
  
  list(
    correlation_distances = correlation_distances,
    clustering_results = clustering_results,
    silhouette_scores = silhouette_scores,
    matching_index_matrix_correlation = matching_index_matrix_correlation,
    matching_index_matrix_across_corr = matching_index_matrix_across_corr
  )
}
  
############################################


scores_list = list(
  CV_t = CV_t,
  RNS = RNS,
  D_slope = D_slope,
  RC = RC,
  CVm = CVm,
  CVmd = CVmd,
  Pi_adj_mean = Pi_adj_mean,
  #PPF = PPF, comparative
  Pi = Pi,
  PImd = PImd,
  PILSM = PILSM,
  RTR = RTR,
  PIR = PIR,
  RDPI = RDPI,
  #RDPI_mean = RDPI_mean, comparative 
  ESPI = ESPI,
  #ESPIID = ESPIID, comparative 
  PSI = PSI,
  #RPI = RPI, comparative 
  PQ = PQ,
  #PR = PR, comparative 
  NRW = NRW,
  #ESP = ESP, comparative 
  #PD = PD, control stress
  #FPI = FPI,control stress
  #TPS = TSP,control stress
  #CEV = CEV, time resolved
  PFI = PFI,
  APC = APC,
  SI = SI,
  RSI = RSI,
  EVS = EVS,
  #MVPi = MVPi, multivariat
  Plasticity_Ratio = Plasticity_Ratio
)



dataset_types = list(
  datasets, 
  datasets_log_skewed,
  datasets_log_log,
  datasets_gamma_skewed,
  datasets_gamma_gamma,
  datasets_weibull_skewed,
  datasets_weibull_weibull,
  datasets_log_gamma,
  datasets_log_weibull,
  datasets_gamma_weibull
)


library(cluster)
library(pheatmap)

# Initialize variables to find best n_clusters
max_mean_of_means = -Inf
overall_max_silhouette = -Inf
best_n_clusters = NULL
best_means = NULL


for (i in 2:10) {  
  complete_res = analyze_scores(scores_list, dataset_types = dataset_types, n_cluster = i)
  
  means = c(
    mean(complete_res$silhouette_scores$hclust$pearson),
    mean(complete_res$silhouette_scores$hclust$kendall),
    mean(complete_res$silhouette_scores$hclust$spearman),
    mean(complete_res$silhouette_scores$pam$pearson),
    mean(complete_res$silhouette_scores$pam$kendall),
    mean(complete_res$silhouette_scores$pam$spearman),
    mean(complete_res$silhouette_scores$dbscan$pearson),
    mean(complete_res$silhouette_scores$dbscan$kendall),
    mean(complete_res$silhouette_scores$dbscan$spearman)
  )
  
  mean_of_means = mean(means)
  print(paste("Mean of means for n_clusters =", i, ":", mean_of_means))
  
  
  if (mean_of_means > max_mean_of_means) {
    max_mean_of_means = mean_of_means
    best_n_clusters = i
    best_means = means
  }
  
  
  current_max_silhouette = max(means)
  if (current_max_silhouette > overall_max_silhouette) {
    overall_max_silhouette = current_max_silhouette
  }
}


print(paste("Best n_clusters:", best_n_clusters))
print("Best means across methods:")
print(best_means)
print(paste("Highest mean of means:", max_mean_of_means))
print(paste("Overall maximum silhouette score across all clusters and methods:", overall_max_silhouette))


pdf("heatmaps_all_traits_best_clustering.pdf")
for (trait_index in 1:3) {
  for (i in 1:length(dataset_types)) {
    best_distance_matrix = complete_res$correlation_distances$pearson[[i]][[trait_index]]
    similarity_matrix = 1 / (1 + best_distance_matrix)
    
    pheatmap(similarity_matrix,
             clustering_distance_rows = as.dist(best_distance_matrix),
             clustering_distance_cols = as.dist(best_distance_matrix),
             clustering_method = "complete",
             cutree_rows = best_n_clusters,
             cutree_cols = best_n_clusters,
             main = paste("Best Clustering Heatmap - Dataset", skewness_type_vector[i], "- Trait", trait_index))
  }
}
dev.off()


score_names = names(scores_list)

pdf("silhouettes_all_datasets_best_clustering.pdf")
for (trait_index in 1:3) {
  for (i in 1:length(dataset_types)) {
    best_distance_matrix = complete_res$correlation_distances$pearson[[i]][[trait_index]]
    row_clusters = cutree(hclust(as.dist(best_distance_matrix), method = "complete"), k = best_n_clusters)
    
    
    silhouette_scores = silhouette(row_clusters, as.dist(best_distance_matrix))
    rownames(silhouette_scores) = score_names
    plot(silhouette_scores, main = paste("Silhouette Plot for Best Clustering - Dataset", skewness_type_vector[i], "- Trait", trait_index))
  }
}
dev.off()


#########################
pdf("mds_all_datasets_best_clustering.pdf")
for (trait_index in 1:3) {
  for (i in 1:length(dataset_types)) {
    
    best_distance_matrix = complete_res$correlation_distances$pearson[[i]][[trait_index]]
    
    
    mds_results = cmdscale(as.dist(best_distance_matrix))
    
    
    row_clusters = cutree(hclust(as.dist(best_distance_matrix), method = "complete"), k = best_n_clusters)
    
    
    plot(mds_results, col = row_clusters, 
         pch = 19, xlab = "MDS1", ylab = "MDS2", 
         main = paste("MDS with Best Clustering - Dataset", skewness_type_vector[i], "- Trait", trait_index))
  }
}
dev.off()

##################################################### co-occurences


n_scores = length(scores_list)
co_occurrence_matrix = matrix(0, nrow = n_scores, ncol = n_scores)
rownames(co_occurrence_matrix) = names(scores_list)
colnames(co_occurrence_matrix) = names(scores_list)


for (cor_method in c("pearson", "spearman", "kendall")) {
  for (method in c("hclust", "pam", "dbscan")) {
    for (dataset_type_index in 1:length(dataset_types)) {
      for (trait_index in 1:3) {  # Adjust trait index range as needed
        
        clustering = complete_res$clustering_results[[method]][[cor_method]][[dataset_type_index]][[trait_index]]
        
        
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
}



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
         main = "Clustering of Scores Based on Co-occurrence in clusterings (Similarity)")



 