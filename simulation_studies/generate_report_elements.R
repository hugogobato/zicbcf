################################################################################
##  Generate Report Tables and Figures from Simulation Results (Updated Layout)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

RESULTS_DIR <- "/home/hugo_souto/Stuff/Research/zicbcf/simulation_studies/results"
OUTPUT_DIR  <- "/home/hugo_souto/Stuff/Research/zicbcf/simulation_studies"

cat("=== Initializing Report Elements Generator (Updated Layout) ===\n")

## ---- Helper to calculate metrics for standard BCF and ZIC-BCF-Smear ----------
process_file <- function(filepath, dgp_name, n_val, zi_level, zero_prop) {
  if (!file.exists(filepath)) {
    warning(sprintf("File not found: %s", filepath))
    return(NULL)
  }
  
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  
  # standard BCF metrics
  bcf_cate_rmse <- mean(df$Linear_CATE_RMSE, na.rm = TRUE)
  bcf_cate_bias <- mean(df$Linear_CATE_Abs_Bias, na.rm = TRUE)
  bcf_cate_cov  <- mean(df$Linear_CATE_Coverage, na.rm = TRUE)
  bcf_cate_corr <- mean(df$Linear_CATE_Correlation, na.rm = TRUE)
  bcf_cate_ci   <- mean(df$Linear_CATE_CI_Length, na.rm = TRUE)
  bcf_ate_rmse  <- sqrt(mean((df$Linear_Est_ATE - df$True_ATE)^2, na.rm = TRUE))
  bcf_ate_bias  <- abs(mean(df$Linear_Est_ATE - df$True_ATE, na.rm = TRUE))
  
  # ZIC-BCF-Smear metrics
  smear_cate_rmse <- mean(df$Smear_CATE_RMSE, na.rm = TRUE)
  smear_cate_bias <- mean(df$Smear_CATE_Abs_Bias, na.rm = TRUE)
  smear_cate_cov  <- mean(df$Smear_CATE_Coverage, na.rm = TRUE)
  smear_cate_corr <- mean(df$Smear_CATE_Correlation, na.rm = TRUE)
  smear_cate_ci   <- mean(df$Smear_CATE_CI_Length, na.rm = TRUE)
  smear_ate_rmse  <- sqrt(mean((df$Smear_Est_ATE - df$True_ATE)^2, na.rm = TRUE))
  smear_ate_bias  <- abs(mean(df$Smear_Est_ATE - df$True_ATE, na.rm = TRUE))
  
  # ZIC-BCF-Smear Hurdle metrics
  smear_hurdle_rmse <- mean(df$Smear_Hurdle_RMSE, na.rm = TRUE)
  smear_hurdle_bias <- mean(df$Smear_Hurdle_Abs_Bias, na.rm = TRUE)
  smear_hurdle_cov  <- mean(df$Smear_Hurdle_Coverage, na.rm = TRUE)
  smear_hurdle_corr <- mean(df$Smear_Hurdle_Correlation, na.rm = TRUE)
  smear_hurdle_ci   <- mean(df$Smear_Hurdle_CI_Length, na.rm = TRUE)
  smear_hurdle_ate_rmse <- sqrt(mean((df$Smear_Est_Hurdle_ATE - df$True_Hurdle_ATE)^2, na.rm = TRUE))
  smear_hurdle_ate_bias <- abs(mean(df$Smear_Est_Hurdle_ATE - df$True_Hurdle_ATE, na.rm = TRUE))
  
  list(
    bcf = data.frame(
      DGP = dgp_name,
      N = n_val,
      ZI_Level = zi_level,
      Zero_Proportion = zero_prop,
      Model = "BCF-Linear",
      CATE_RMSE = bcf_cate_rmse,
      CATE_Abs_Bias = bcf_cate_bias,
      CATE_Coverage = bcf_cate_cov,
      CATE_Correlation = bcf_cate_corr,
      CATE_CI_Length = bcf_cate_ci,
      ATE_RMSE = bcf_ate_rmse,
      ATE_Abs_Bias = bcf_ate_bias,
      Hurdle_CATE_RMSE = NA,
      Hurdle_CATE_Correlation = NA,
      Hurdle_CATE_Coverage = NA,
      Hurdle_CATE_CI_Length = NA,
      stringsAsFactors = FALSE
    ),
    smear = data.frame(
      DGP = dgp_name,
      N = n_val,
      ZI_Level = zi_level,
      Zero_Proportion = zero_prop,
      Model = "ZIC-BCF-Smear",
      CATE_RMSE = smear_cate_rmse,
      CATE_Abs_Bias = smear_cate_bias,
      CATE_Coverage = smear_cate_cov,
      CATE_Correlation = smear_cate_corr,
      CATE_CI_Length = smear_cate_ci,
      ATE_RMSE = smear_ate_rmse,
      ATE_Abs_Bias = smear_ate_bias,
      Hurdle_CATE_RMSE = smear_hurdle_rmse,
      Hurdle_CATE_Correlation = smear_hurdle_corr,
      Hurdle_CATE_Coverage = smear_hurdle_cov,
      Hurdle_CATE_CI_Length = smear_hurdle_ci,
      stringsAsFactors = FALSE
    )
  )
}

## ---- Mapping of scenarios and nominal zero proportions ----------------------
get_zero_prop <- function(dgp, zi_level) {
  if (dgp == "DGP A: Log-Normal Hurdle") {
    props <- c("1" = 0.854, "2" = 0.619, "3" = 0.367, "4" = 0.185, "5" = 0.066)
    return(props[as.character(zi_level)])
  } else if (dgp == "DGP B: Gamma Hurdle") {
    props <- c("1" = 0.851, "2" = 0.607, "3" = 0.395, "4" = 0.181, "5" = 0.049)
    return(props[as.character(zi_level)])
  } else if (dgp == "DGP C: Tweedie Semicontinuous") {
    props <- c("1" = 0.600, "2" = 0.397, "3" = 0.180, "4" = 0.115, "5" = 0.028)
    return(props[as.character(zi_level)])
  }
  return(0)
}

dgps <- c("DGP A: Log-Normal Hurdle", "DGP B: Gamma Hurdle", "DGP C: Tweedie Semicontinuous")
dgp_abbrev <- c("DGP A: Log-Normal Hurdle" = "dgp_a", "DGP B: Gamma Hurdle" = "dgp_b", "DGP C: Tweedie Semicontinuous" = "dgp_c")

## =============================================================================
##  1. AGGREGATE STANDARD MAIN RUNS (N = 500, ZI Level 3 / Normal)
## =============================================================================
cat("\nAggregating Standard Main Runs (N = 500)...\n")
std_results <- list()
for (dgp in dgps) {
  abbrev <- dgp_abbrev[dgp]
  file_path <- file.path(RESULTS_DIR, "Normal", sprintf("results_bcf_%s_N500.csv", abbrev))
  zero_p <- get_zero_prop(dgp, 3)
  res <- process_file(file_path, dgp, 500, 3, zero_p)
  if (!is.null(res)) {
    std_results[[length(std_results) + 1]] <- res$bcf
    std_results[[length(std_results) + 1]] <- res$smear
  }
}
std_df <- do.call(rbind, std_results)

# Create a markdown table of standard results
cat("\n--- STANDARD COMPARATIVE RESULTS TABLE (N=500, ZI Level 3) ---\n")
std_table_md <- std_df %>%
  select(DGP, Model, CATE_RMSE, CATE_Abs_Bias, CATE_Coverage, CATE_CI_Length, ATE_RMSE, ATE_Abs_Bias) %>%
  mutate(
    CATE_RMSE = round(CATE_RMSE, 3),
    CATE_Abs_Bias = round(CATE_Abs_Bias, 3),
    CATE_Coverage = percent(CATE_Coverage, accuracy = 0.1),
    CATE_CI_Length = round(CATE_CI_Length, 3),
    ATE_RMSE = round(ATE_RMSE, 3),
    ATE_Abs_Bias = round(ATE_Abs_Bias, 3)
  )
print(knitr::kable(std_table_md, format = "markdown"))

## =============================================================================
##  2. AGGREGATE SAMPLE SIZE N SENSITIVITY RUNS
## =============================================================================
cat("\nAggregating N-Sensitivity Runs...\n")
n_results <- list()
n_list <- c(100, 250, 500, 1000)

for (dgp in dgps) {
  abbrev <- dgp_abbrev[dgp]
  zero_p <- get_zero_prop(dgp, 3)
  
  for (n_val in n_list) {
    if (n_val == 500) {
      file_path <- file.path(RESULTS_DIR, "Normal", sprintf("results_bcf_%s_N500.csv", abbrev))
    } else {
      file_path <- file.path(RESULTS_DIR, sprintf("results_bcf_%s_N%d.csv", abbrev, n_val))
    }
    
    res <- process_file(file_path, dgp, n_val, 3, zero_p)
    if (!is.null(res)) {
      n_results[[length(n_results) + 1]] <- res$bcf
      n_results[[length(n_results) + 1]] <- res$smear
    }
  }
}
n_df <- do.call(rbind, n_results)

## =============================================================================
##  3. AGGREGATE ZERO-INFLATION SENSITIVITY RUNS
## =============================================================================
cat("\nAggregating Zero-Inflation Sensitivity Runs...\n")
zi_results <- list()
zi_levels <- 1:5

for (dgp in dgps) {
  abbrev <- dgp_abbrev[dgp]
  
  for (lvl in zi_levels) {
    zero_p <- get_zero_prop(dgp, lvl)
    
    if (lvl == 3) {
      file_path <- file.path(RESULTS_DIR, "Normal", sprintf("results_bcf_%s_N500.csv", abbrev))
    } else {
      file_path <- file.path(RESULTS_DIR, sprintf("ZI_%d", lvl), sprintf("results_bcf_%s_N500.csv", abbrev))
    }
    
    res <- process_file(file_path, dgp, 500, lvl, zero_p)
    if (!is.null(res)) {
      zi_results[[length(zi_results) + 1]] <- res$bcf
      zi_results[[length(zi_results) + 1]] <- res$smear
    }
  }
}
zi_df <- do.call(rbind, zi_results)

# Save the aggregated sensitivity outputs to CSVs
write.csv(n_df, file.path(RESULTS_DIR, "n_sensitivity_summary.csv"), row.names = FALSE)
write.csv(zi_df, file.path(RESULTS_DIR, "zi_sensitivity_summary.csv"), row.names = FALSE)
cat("[SUCCESS] Aggregated sensitivity CSVs saved to results/.\n")

## =============================================================================
##  4. PLOTTING SENSITIVITY ANALYSES (2x2 GRID, 4 PANELS PER IMAGE)
## =============================================================================
cat("\nGenerating 2x2 grid visualizations...\n")

premium_colors <- c(
  "BCF-Linear"     = "#64748b", # Cool slate grey
  "ZIC-BCF-Smear"  = "#8b5cf6"  # Premium purple
)

# --- 4.1 N-Sensitivity Visualizations ---
# Metrics to include: CATE RMSE, CATE Coverage, ATE RMSE, CATE CI Length (exactly 4 metrics)
n_long <- n_df %>%
  select(DGP, N, Model, CATE_RMSE, CATE_Coverage, ATE_RMSE, CATE_CI_Length) %>%
  pivot_longer(
    cols = c(CATE_RMSE, CATE_Coverage, ATE_RMSE, CATE_CI_Length),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "CATE_RMSE" ~ "CATE RMSE (Lower is Better)",
      Metric == "CATE_Coverage" ~ "CATE 95% Coverage (Nominal = 0.95)",
      Metric == "ATE_RMSE" ~ "ATE RMSE (Lower is Better)",
      Metric == "CATE_CI_Length" ~ "CATE Credible Interval Length"
    )
  )

for (dgp in dgps) {
  dgp_n_data <- n_long %>% filter(DGP == dgp)
  dgp_abbr <- dgp_abbrev[dgp]
  
  p <- ggplot(dgp_n_data, aes(x = N, y = Value, color = Model, group = Model)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = n_list, trans = "log10") +
    scale_color_manual(values = premium_colors) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "#1e293b"),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#64748b", margin = margin(b = 12)),
      strip.background = element_rect(fill = "#f8fafc", color = NA),
      strip.text = element_text(face = "bold", size = 11, color = "#334155"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f5f9"),
      panel.background = element_rect(fill = "white", color = "#cbd5e1"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold"),
      legend.margin = margin(t = 8),
      axis.title = element_text(face = "bold", size = 11),
      axis.title.y = element_blank()
    ) +
    labs(
      title = paste("Sample Size (N) Sensitivity Analysis:", dgp),
      subtitle = "Evaluating causal performance across sample sizes N in {100, 250, 500, 1000}",
      x = "Sample Size (N) - Log Scale"
    )
  
  # Add nominal 95% line for coverage facet
  coverage_nominal <- data.frame(
    Metric = "CATE 95% Coverage (Nominal = 0.95)",
    yintercept = 0.95
  )
  p <- p + geom_hline(data = coverage_nominal, aes(yintercept = yintercept), linetype = "dashed", color = "#475569", alpha = 0.7)
  
  file_name <- sprintf("n_sensitivity_%s.png", dgp_abbr)
  ggsave(file.path(RESULTS_DIR, file_name), plot = p, width = 10, height = 8, dpi = 300)
  cat(sprintf("[SUCCESS] N-Sensitivity Graph saved to: %s\n", file.path(RESULTS_DIR, file_name)))
}


# --- 4.2 Zero-Inflation Sensitivity Visualizations (Normalized RMSE) ---
# Metrics: Normalized CATE RMSE, CATE Coverage, Normalized ATE RMSE, CATE CI Length (exactly 4 metrics)
zi_norm <- zi_df %>%
  group_by(DGP, ZI_Level) %>%
  mutate(
    CATE_RMSE_Norm = CATE_RMSE / CATE_RMSE[Model == "BCF-Linear"],
    ATE_RMSE_Norm  = ATE_RMSE / ATE_RMSE[Model == "BCF-Linear"]
  ) %>%
  ungroup() %>%
  select(DGP, Zero_Proportion, Model, CATE_RMSE_Norm, CATE_Coverage, ATE_RMSE_Norm, CATE_CI_Length) %>%
  pivot_longer(
    cols = c(CATE_RMSE_Norm, CATE_Coverage, ATE_RMSE_Norm, CATE_CI_Length),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "CATE_RMSE_Norm" ~ "Normalized CATE RMSE (Linear = 1.0)",
      Metric == "CATE_Coverage" ~ "CATE 95% Coverage",
      Metric == "ATE_RMSE_Norm" ~ "Normalized ATE RMSE (Linear = 1.0)",
      Metric == "CATE_CI_Length" ~ "CATE Credible Interval Length"
    )
  )

for (dgp in dgps) {
  dgp_zi_data <- zi_norm %>% filter(DGP == dgp)
  dgp_abbr <- dgp_abbrev[dgp]
  
  p <- ggplot(dgp_zi_data, aes(x = Zero_Proportion, y = Value, color = Model, group = Model)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_color_manual(values = premium_colors) +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(color = "#1e293b"),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#64748b", margin = margin(b = 12)),
      strip.background = element_rect(fill = "#f8fafc", color = NA),
      strip.text = element_text(face = "bold", size = 11, color = "#334155"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f5f9"),
      panel.background = element_rect(fill = "white", color = "#cbd5e1"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold"),
      legend.margin = margin(t = 8),
      axis.title = element_text(face = "bold", size = 11),
      axis.title.y = element_blank()
    ) +
    labs(
      title = paste("Zero-Inflation Sensitivity Analysis:", dgp),
      subtitle = "RMSE is normalized using BCF-Linear as the baseline (1.0) at each level",
      x = "Zero Proportion (Zero Inflation Degree)"
    )
  
  # Add nominal 95% line for coverage facet
  coverage_nominal <- data.frame(
    Metric = "CATE 95% Coverage",
    yintercept = 0.95
  )
  p <- p + geom_hline(data = coverage_nominal, aes(yintercept = yintercept), linetype = "dashed", color = "#475569", alpha = 0.7)
  
  # Add baseline 1.0 horizontal line for the normalized facets
  baseline_one <- data.frame(
    Metric = c("Normalized CATE RMSE (Linear = 1.0)", "Normalized ATE RMSE (Linear = 1.0)"),
    yintercept = 1.0
  )
  p <- p + geom_hline(data = baseline_one, aes(yintercept = yintercept), linetype = "dotted", color = "#64748b", linewidth = 1.0, alpha = 0.8)
  
  file_name <- sprintf("zi_sensitivity_%s.png", dgp_abbr)
  ggsave(file.path(RESULTS_DIR, file_name), plot = p, width = 10, height = 8, dpi = 300)
  cat(sprintf("[SUCCESS] ZI-Sensitivity Graph saved to: %s\n", file.path(RESULTS_DIR, file_name)))
}

cat("\n=== All Report Elements Generated Successfully! ===\n")
