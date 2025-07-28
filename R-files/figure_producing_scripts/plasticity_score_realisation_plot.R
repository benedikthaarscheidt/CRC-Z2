
#this script is for the assembly of the second figure containing the boxplots for the distribution of the scores as well as the clustering results (hclust with kendalll correlation distance as distance metrics) in form of a dendrogram 
#and the adjusted rand index table comparing the clusterings across the different genotype forms. THe reszlts suggest (sampling interval=5) that the scores obey to a similar clustering for all different genotype forms

library(dplyr)
library(readr)
library(cluster)      # silhouette()
library(mclust)       # adjustedRandIndex()
library(factoextra)   # fviz_dend()
library(patchwork)    # layou
library(cluster)     # for silhouette()
library(pheatmap)
library(dendextend)  # for nicer dendrogram handling
library(ggbeeswarm)
library(tidyr)
library(ggplot2)
library(patchwork)
library(forcats)    # for fct_reorder()
library(ggridges)   # for geom_density_ridges()
library(pheatmap)   # for heatmaps
library(ggpubr)    # for ggtexttable & ggarrange
library(dplyr)     # for mutate_if
library(tibble)    # for rownames_to_column
library(aricode)
library(knitr)
library(scales)
# Read 
df <- read_csv("~/CRC_1644_Z2/R-files/regression_summary_stats/regression_data_full_interval_15_indices_1_16_31_46_50.csv")
na_counts <- colSums(is.na(df))

# keep only columns with zero NAs (there are some in the case of 1 and 2 samples across the environmental gradient as the RN and RNN cannot be calculated in this scenario)
df <- df[ , na_counts == 0]


df <- df %>% rename(Genotype = names(df)[1])
df=df[,1:28]

scores_long <- df %>%
  pivot_longer(
    cols       = -Genotype,       # all the metric columns
    names_to   = "Metric",
    values_to  = "Score"
  )

# Make one histogram + density per metric. This is a huge figure alone with all the distributions of all the plasticity scores. I suppose this can either go into the supplements or be left out entirely
make_plot <- function(metric_name) {
  sub <- scores_long %>% filter(Metric == metric_name)
  ggplot(sub, aes(x = Score)) +
    geom_histogram(aes(y = ..density..),
                   bins  = 30,
                   fill  = "grey70",
                   color = "white") +
    geom_density(size = 1, color = "black") +
    labs(title = metric_name, x = NULL, y = "Density") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.text.x      = element_text(size = 12),
      axis.text.y      = element_text(size = 11) )
    
}


plots   <- lapply(unique(scores_long$Metric), make_plot)
combined <- wrap_plots(plots, ncol = 4) 
print(combined)
ggsave("~/CRC_1644_Z2/plots/plasticity_scores_histograms.pdf", combined, width = 10, height = 8,
       dpi = 900, units = "in", device = "pdf")


##############################
pp=ggplot(scores_long, aes(x = Metric, y = Score)) +
  geom_boxplot(
    width         = 0.9,       # almost fill the space → boxes nearly touch
    outlier.size  = 0.5,
    fill          = "grey80",
    color         = "grey30"
  ) +
  # remove any extra padding around the categories
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(expand = FALSE) +
  labs(
    title = "All Metrics Boxplot",
    x     = NULL,
    y     = "Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),             # no horizontal lines
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.x        = element_text(
      angle = 90,                 # rotate labels if you need them legible
      hjust = 1,
      vjust = 0.5,
      size  = 15
    ),
    axis.text.y        = element_text(size = 12),
    axis.ticks.length  = unit(1, "pt"),
    plot.title         = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin        = margin(2, 2, 2, 2)
  )
print(pp)

ggsave("~/CRC_1644_Z2/plots/figures/plasticity_scores_boxplots.pdf", pp, width = 10, height = 8,
       dpi = 900, units = "in", device = "pdf")

##########################################################


# Zoran suggested that this might look good with violin plots but i think it doesnt 
pp_violin = ggplot(scores_long, aes(x = Metric, y = Score)) +
  geom_violin(
    trim   = FALSE,      # show full tails
    width  = 1.2,        # make the violins fatter
    adjust = 2,        # smooth out the density a bit more
    fill   = "grey80",
    color  = "grey30"
  ) +
  
  geom_jitter(
    aes(color = Metric),          # optional: color points by category
    position = position_jitter(
      width  = 0.2,   # jitter horizontally
      height = 0      # no vertical jitter
    ),
    size  = 0.7,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  # remove any extra padding around the categories
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(expand = FALSE) +
  labs(
    title = "All Metrics Violin Plot",
    x     = NULL,
    y     = "Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),             # no vertical grid lines
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.x        = element_text(
      angle = 90,                 # rotate labels for legibility
      hjust = 1,
      vjust = 0.5,
      size  = 15
    ),
    axis.text.y        = element_text(size = 10),
    axis.ticks.length  = unit(1, "pt"),
    plot.title         = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.margin        = margin(2, 2, 2, 2)
  )

print(pp_violin)
############################################################

df=df[,-1]



permTestARI <- function(cl1, cl2, n_perm = 1000, alternative = "greater") {
  #helper function for calculating significance values fpr the clusterings 
  obs  <- adjustedRandIndex(cl1, cl2)
  perm <- replicate(n_perm, adjustedRandIndex(cl1, sample(cl2)))
  if (alternative == "greater") {
    p <- mean(perm >= obs)
  } else if (alternative == "less") {
    p <- mean(perm <= obs)
  } else {
    m0 <- mean(perm)
    p  <- mean(abs(perm - m0) >= abs(obs - m0))
  }
  list(observed_ARI = obs, p_value = p)
}


form_ranges <- list(
  linear     = 1:20,
  gaussian   = 21:40,
  sinusoidal = 41:60,
  wave       = 61:80
)
forms      <- names(form_ranges)
all_labels <- c(forms, "combined")
k <-6 # number of clusters to cut the dendrogram into --> 3 was the number of clusters to cluster the data into that yielded the higherst agreement between the different genotype forms as per the adjusted rand index table
# k=6 also yields pretty good results though. This whole thing is just to show that the scores are similar across the different genotype forms and that the clustering is not just a random artefact of the data.

clusterings <- list()
hc_list     <- list()

for (f in forms) {
  sub_df <- df[ form_ranges[[f]], , drop = FALSE ]
  C      <- cor(sub_df, method = "kendall", use = "pairwise.complete.obs")
  D      <- as.dist(1 - C)
  hc     <- hclust(D, method = "complete")
  hc_list[[f]]     <- hc
  clusterings[[f]] <- cutree(hc, k = k)
}

# combined
C_all                  <- cor(df, method = "kendall", use = "pairwise.complete.obs")
D_all                  <- as.dist(1 - C_all)
hc_list[["combined"]]  <- hclust(D_all, method = "complete")
clusterings[["combined"]] <- cutree(hc_list[["combined"]], k = k)


nA      <- length(all_labels)
ARI_obs <- matrix(NA, nA, nA, dimnames = list(all_labels, all_labels))
ARI_p   <- ARI_obs

for (i in seq_len(nA - 1)) {
  for (j in (i + 1):nA) {
    a   <- all_labels[i]
    b   <- all_labels[j]
    res <- permTestARI(clusterings[[a]], clusterings[[b]], n_perm = 1000)
    ARI_obs[a, b] <- ARI_obs[b, a] <- res$observed_ARI
    ARI_p[  a, b] <- ARI_p[  b, a] <- res$p_value
  }
}


make_dend_plot <- function(hc, title, shrink = 0.5) {
  #this is a little helper function for visualising the dendrograms
  shrink_trans <- trans_new(
    name      = paste0("shrink_", shrink),
    transform = function(x) x * shrink,
    inverse   = function(x) x / shrink
  )
  p <- fviz_dend(hc,
                 k           = k,
                 rect        = FALSE,
                 show_labels = TRUE,
                 repel       = FALSE,    # use `cex` not ggrepel
                 k_colors    = "jco",
                 cex         = 0.8     # shrink tip labels
  ) +
    ggtitle(title) +
    scale_y_continuous(trans = shrink_trans) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 15) +
    theme(
      axis.ticks.x   = element_blank(),
      axis.title.x   = element_blank(),
      plot.margin    = margin(t = 5, r = 5, b = 15, l = 5)
    )
  return(p)
}







plots_sep <- lapply(forms, function(f) {
  make_dend_plot(hc_list[[f]], title = f)
})
plot_comb <- make_dend_plot(hc_list[["combined"]], title = "combined")
print(plot_comb)

sep_panel <- (plot_comb | plots_sep[[2]]) /
  (plots_sep[[3]] | plots_sep[[4]])

print(sep_panel)
ggsave(
  "~/CRC_1644_Z2/plots/figures/plasticity_scores_dendrograms.pdf",
  sep_panel,
  width = 10, height = 8, dpi = 900, units = "in", device = "pdf")

ggsave(
  "~/CRC_1644_Z2/plots/figures/dendrogram_linear.pdf",
  plots_sep[[1]],
  width = 10, height = 8, dpi = 900, units = "in", device = "pdf")

form_panel <- (plots_sep[[1]] | plots_sep[[2]]) /
  (plots_sep[[3]] | plots_sep[[4]])


########################################


# 1) Asterisks function (same as before)
sig_asterisks <- function(p) {
  case_when(
    p < 0.001           ~ "***",
    p < 0.01            ~ "**",
    p < 0.05            ~ "*",
    TRUE                ~ ""
  )
}

# 2) Turn your matrices into data frames
obs_df <- ARI_obs %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form")

p_df <- ARI_p %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form")

# 3) Build an “annotated ARI” data frame by pasting stars on the ARI values
obs_annot_df <- obs_df
for (col in names(obs_annot_df)[-1]) {
  obs_annot_df[[col]] <- sprintf(
    "%.3f%s",
    obs_df[[col]],
    sig_asterisks(p_df[[col]])
  )
}



# 1) Pivot to long form for ggplot
obs_long <- obs_annot_df %>%
  pivot_longer(-Form, names_to = "Comparison", values_to = "ARI_star")

ord <- levels(obs_long$Comparison)

tbl_obs_plot <- ggplot(obs_long, aes(x = Comparison, y = Form)) +
  # ensure x follows your Comparison order, labels on top
  scale_x_discrete(limits = ord, position = "top") +
  # use exactly the *reverse* of that same order on y
  scale_y_discrete(limits = ord) +
  
  # draw tiles + text
  geom_tile(fill = "white", color = "grey70") +
  geom_text(aes(label = ARI_star), size = 4) +
  
  # square cells
  coord_fixed() +
  
  # tidy up the theme
  theme_minimal(base_size = 15) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 0, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    axis.title       = element_blank(),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin      = margin(2, 2, 2, 2, "mm")
  )

print(tbl_obs_plot)


ggsave(
  "~/CRC_1644_Z2/plots/figures/plasticity_scores_ari_table.pdf",
  tbl_obs_plot,
  width = 6, height = 6, dpi = 900, units = "in", device = "pdf"
)


####################




final_fig3 <- wrap_plots(
  A = wrap_elements(full = pp),
  B = wrap_elements(full = plots_sep[[1]]),
  #C = wrap_elements(full = plot_comb),
  C = wrap_elements(full = tbl_obs_plot),  # now a ggplot
  design = "
AAAAAA
AAAAAA
AAAAAA
BBBBCC
BBBBCC
"
) +
  plot_annotation(
    tag_levels = "A",
    title      = "Plasticity Score distributions and hierarchical clustering of correlation distances",
    theme      = theme(
      plot.title        = element_text(size = 12, face = "bold"),
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )


print(final_fig3)

ggsave(
  "~/CRC_1644_Z2/plots/figures/figure2_2.pdf",
  final_fig3,
  width = 15, height = 18, dpi = 900, units = "in", device = "pdf"
)


################


#this is the same figure as the previous one but with the violin plots instead of the boxplots
final_fig4 <- wrap_plots(
  A = wrap_elements(full = pp_violin),
  B = wrap_elements(full = plots_sep[[1]]),
  #C = wrap_elements(full = plot_comb),
  C= wrap_elements(full = tbl_obs_plot),  # now a ggplot
  design = "
AAAAAA
AAAAAA
AAAAAA
AAAAAA
BBBBCC
BBBBCC

"
) +
  plot_annotation(
    tag_levels = "A",
    title      = "Plasticity Score distributions and hierarchical clustering of correlation distances",
    theme      = theme(
      plot.title        = element_text(size = 12, face = "bold"),
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

ggsave(
  "~/CRC_1644_Z2/plots/figures/figure2_3.pdf",
  final_fig4,
  width = 15, height = 18, dpi = 900, units = "in", device = "pdf"
)
print(final_fig4)
#########################################


#this figure is supposed to go in the supplements
final_fig7 <- wrap_plots(
  A = wrap_elements(full = plots_sep[[2]]),
  B = wrap_elements(full = plots_sep[[3]]),
  C= wrap_elements(full = plots_sep[[4]]),
  D = wrap_elements(full = plot_comb),
    # now a ggplot
  design = "
AABB
AABB
CCDD
CCDD
"
) +
  plot_annotation(
    tag_levels = "A",
    title      = "Hierarchical clustering of correlation distances",
    theme      = theme(
      plot.title        = element_text(size = 12, face = "bold"),
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

ggsave(
  "~/CRC_1644_Z2/plots/figures/dendro_supp.pdf",
  final_fig7,
  width = 15, height = 13, dpi = 900, units = "in", device = "pdf"
)

print(final_fig7)