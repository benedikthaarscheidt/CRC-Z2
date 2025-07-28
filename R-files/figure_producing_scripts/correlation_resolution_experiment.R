# here we do a similar robustness analysis as in correlation_resolution.R, however in the regression we additionally consider confounders (the summary statistics) 
#in order to distinguish effects which occur due to the sampling effect on the summary statistics and might unknowingly effect the slope of the sampling resolution

library(tidyverse)
library(lme4)
library(lmerTest)


directory <- "~/CRC_1644_Z2/R-files/regression_summary_stats"
files     <- c(
  "regression_data_full_interval_49_indices_1_50.csv",
  "regression_data_full_interval_25_indices_1_26_50.csv",
  "regression_data_full_interval_20_indices_1_21_41_50.csv",
  "regression_data_full_interval_15_indices_1_16_31_46_50.csv",
  "regression_data_full_interval_10_indices_1_11_21_31_41_50.csv",
  "regression_data_full_interval_5_indices_1_6_11_16_21_26_31_36_41_46_50.csv",
  "regression_data_full_interval_2_indices_1_3_5_7_9_11_13_15_17_19_21_23_25_27_29_31_33_35_37_39_41_43_45_47_49_50.csv",
  "regression_data_full_interval_1_indices_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40_41_42_43_44_45_46_47_48_49_50.csv"
)
n_samples <- c(2, 3, 4, 5, 6, 11, 26, 50)


first_dat    <- read_csv(file.path(directory, files[1]), show_col_types = FALSE)
score_cols   <- names(first_dat)[2:28]
summary_cols <- names(first_dat)[29:(ncol(first_dat) - 2)]


long <- tibble(
  variable    = character(),
  value       = numeric(),
  sample_size = numeric(),
  replicate   = integer()
)

for (i in seq_along(files)) {
  dat <- read_csv(file.path(directory, files[i]), show_col_types = FALSE)
  scores_df <- dat %>% select(all_of(score_cols))
  n   <- nrow(scores_df)
  nv  <- ncol(scores_df)
  
  long <- bind_rows(long,
                    tibble(
                      variable    = rep(score_cols, each = n),
                      value       = as.vector(as.matrix(scores_df)),
                      sample_size = n_samples[i],
                      replicate   = rep(seq_len(n), times = nv)
                    )
  )
}

long <- long %>%
  mutate(
    replicate   = factor(replicate),
    sample_size = as.numeric(sample_size)
  )


wide_summ <- tibble()

for (i in seq_along(files)) {
  dat <- read_csv(file.path(directory, files[i]), show_col_types = FALSE)
  
  sums_df <- dat %>%
    select(all_of(summary_cols)) %>%
    mutate(
      replicate   = seq_len(n()),
      sample_size = n_samples[i]
    )
  
  wide_summ <- bind_rows(wide_summ, sums_df)
}

wide_summ <- wide_summ %>%
  mutate(
    replicate   = factor(replicate),
    sample_size = as.numeric(sample_size)
  )


long2 <- left_join(long, wide_summ, by = c("replicate", "sample_size"))


vars <- unique(long2$variable)
out  <- tibble(
  variable   = vars,
  estimate   = NA_real_,
  std.error  = NA_real_,
  t.value    = NA_real_,
  p.value    = NA_real_
)


conf_terms <- paste(summary_cols, collapse = " + ")

for (j in seq_along(vars)) {
  vv  <- vars[j]
  sub <- filter(long2, variable == vv, !is.na(value))
  if (length(unique(sub$sample_size)) < 2) next
  
  # build formula: value ~ sample_size + <all summary_cols> + (1|replicate)
  fm <- as.formula(
    paste0("value ~ sample_size + ", conf_terms, " + (1 | replicate)")
  )
  
  mod <- lmer(fm, data = sub)
  cf  <- summary(mod)$coefficients
  
  out$estimate[j]  <- cf["sample_size", "Estimate"]
  out$std.error[j] <- cf["sample_size", "Std. Error"]
  out$t.value[j]   <- cf["sample_size", "t value"]
  out$p.value[j]   <- cf["sample_size", "Pr(>|t|)"]
}


out <- out %>%
  mutate(
    p.adj    = p.adjust(p.value, method = "BH"),
    ci_lower = estimate - qnorm(0.975) * std.error,
    ci_upper = estimate + qnorm(0.975) * std.error,
    class    = case_when(
      abs(estimate) < 0.001        ~ "robust",
      p.adj         < 0.05         ~ "sample-size-sensitive",
      TRUE                          ~ "robust"
    )
  )

print(out)

plot_tbl <- arrange(out, estimate)

old_mar <- par("mar")
par(mar = c(5, 12, 4, 2))

dotchart(
  plot_tbl$estimate,
  labels = plot_tbl$variable,
  pch    = 19,
  col    = ifelse(plot_tbl$class == "sample-size-sensitive",
                  "tomato", "skyblue"),
  xlab   = "Adjusted slope for sample_size",
  main   = "Effect of Sampling Effort on Plasticity Scores\n(adjusted for distribution confounders)",
  
)
for (i in seq_len(nrow(plot_tbl))) {
  arrows(
    plot_tbl$ci_lower[i], i,
    plot_tbl$ci_upper[i], i,
    angle = 90, code = 3, length = 0.05,
    col   = ifelse(plot_tbl$class[i] == "sample-size-sensitive",
                   "tomato", "skyblue")
  )
}

legend(
  "topleft",
  legend = c("robust", "sample-size-sensitive"),
  pch    = 19,
  col    = c("skyblue", "tomato"),
  bty    = "n"       # no box around the legend
)
par(mar = old_mar)

path="~/CRC_1644_Z2/plots/figures"
if (!dir.exists(dirname(path))) {
  dir.create(dirname(path), recursive = TRUE)
}

pdf(
  file   = file.path("~/CRC_1644_Z2/plots/figures/confoundermodel_sumstats.pdf"),
  width  = 8,
  height = 10
)
par(mar = c(5, 12, 4, 2))
dotchart(
  plot_tbl$estimate,
  labels = plot_tbl$variable,
  pch    = 19,
  col    = ifelse(plot_tbl$class == "sample-size-sensitive",
                  "tomato", "skyblue"),
  xlab   = "Adjusted slope for sample_size",
  main   = "Effect of Sampling Effort on Plasticity Scores\n(adjusted for distribution confounders)"
)
for (i in seq_len(nrow(plot_tbl))) {
  arrows(
    plot_tbl$ci_lower[i], i,
    plot_tbl$ci_upper[i], i,
    # angle = 90, code = 3, length = 0.05,
    col   = ifelse(plot_tbl$class[i] == "sample-size-sensitive",
                   "tomato", "skyblue")
  )
  
}
legend(
  "topleft",
  legend = c("robust", "sample-size-sensitive"),
  pch    = 19,
  col    = c("skyblue", "tomato"),
  bty    = "n"       # no box around the legend
)
dev.off()
