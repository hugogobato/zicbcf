################################################################################
##  Semicontinuous BCF Treatment Sensitivity: Master Orchestrator
##
##  This script orchestrates the treatment sensitivity analysis by running each
##  scenario (1 to 15) in a separate, isolated R process.
##
##  It normalizes CATE RMSE and ATE Absolute Bias using BCF-Linear as the
##  baseline (value = 1.0) for each treatment effect scale, enabling direct
##  relative comparisons across models as treatment effect magnitude changes.
################################################################################

library(ggplot2)
library(tidyr)
library(dplyr)

RESULTS_DIR <- "simulation_studies/reduced_simulations/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

cat("=== Starting Semicontinuous BCF Treatment Sensitivity Master Orchestrator ===\n")
cat("RAM Protection: Process-isolated sequential execution with checkpoints.\n\n")

t_start <- Sys.time()

# Loop through all 15 scenarios
for (i in 1:15) {
  checkpoint_file <- file.path(RESULTS_DIR, sprintf("scenario_treatment_results_%d.csv", i))
  
  if (file.exists(checkpoint_file)) {
    cat(sprintf("[Checkpoint] Scenario %d/15 already completed. Skipping.\n", i))
  } else {
    cat(sprintf("\n========================================================================\n"))
    cat(sprintf("=== Orchestrating Scenario %d/15 (Isolated Process) ===\n", i))
    cat(sprintf("========================================================================\n"))
    
    # Run the scenario in a completely separate R process
    cmd <- sprintf("Rscript simulation_studies/reduced_simulations/run_treatment_single.R %d", i)
    exit_code <- system(cmd)
    
    if (exit_code != 0) {
      stop(sprintf("\nError: Scenario %d failed with exit code %d. Stopped.\n", i, exit_code))
    }
  }
}

t_end <- Sys.time()
cat(sprintf("\n=== All Scenarios Finished in %.1f minutes ===\n", as.numeric(difftime(t_end, t_start, units = "mins"))))

## ---- Aggregate All Results --------------------------------------------------
cat("\nAggregating results...\n")
results_list <- list()

for (i in 1:15) {
  checkpoint_file <- file.path(RESULTS_DIR, sprintf("scenario_treatment_results_%d.csv", i))
  if (file.exists(checkpoint_file)) {
    results_list[[i]] <- read.csv(checkpoint_file, stringsAsFactors = FALSE)
  } else {
    stop(sprintf("Fatal Error: Missing results file for Scenario %d.", i))
  }
}

all_results <- do.call(rbind, results_list)
write.csv(all_results, file.path(RESULTS_DIR, "treatment_sensitivity_results.csv"), row.names = FALSE)
cat(sprintf("[SUCCESS] Master CSV saved to: %s\n", file.path(RESULTS_DIR, "treatment_sensitivity_results.csv")))

## ---- Generate Beautiful Visualizations with Normalization -------------------
cat("Generating ggplot2 sensitivity charts with BCF-Linear normalization...\n")

long_results <- all_results %>%
  pivot_longer(
    cols = starts_with(c("Linear_", "PathA_", "Gemini_")),
    names_to = c("Model", "Metric"),
    names_sep = "_"
  ) %>%
  mutate(
    Model = case_when(
      Model == "Linear" ~ "BCF-Linear",
      Model == "PathA"  ~ "ZIC-BCF (Path A)",
      Model == "Gemini" ~ "ZIC-BCF-Smear (Best_Path_Gemini)"
    )
  )

# Normalize RMSE and Bias by BCF-Linear (basis = 1.0) at each DGP and k level
norm_metrics <- long_results %>%
  filter(Metric %in% c("RMSE", "Bias")) %>%
  group_by(DGP, k, Metric) %>%
  mutate(
    basis_val = value[Model == "BCF-Linear"],
    # Handle division by zero or tiny values if ATE is exactly zero when k=0
    value_norm = ifelse(abs(basis_val) < 1e-6, 1.0, value / basis_val)
  ) %>%
  ungroup() %>%
  select(DGP, k, Model, Metric, value = value_norm) %>%
  mutate(Metric = ifelse(Metric == "RMSE", "Normalized CATE RMSE (Linear = 1.0)", "Normalized ATE Abs Bias (Linear = 1.0)"))

raw_metrics <- long_results %>%
  filter(Metric %in% c("Coverage", "Correlation")) %>%
  mutate(Metric = ifelse(Metric == "Coverage", "CATE 95% Coverage", "CATE Correlation")) %>%
  select(DGP, k, Model, Metric, value)

plot_data <- rbind(raw_metrics, norm_metrics)

premium_colors <- c(
  "BCF-Linear"                      = "#94a3b8", # Cool grey (slate)
  "ZIC-BCF (Path A)"                = "#10b981", # Emerald Green
  "ZIC-BCF-Smear (Best_Path_Gemini)" = "#8b5cf6"  # Purple (Gemini)
)

plot_dgp_sensitivity <- function(dgp_name, file_name) {
  dgp_data <- plot_data %>% filter(DGP == dgp_name)
  
  p <- ggplot(dgp_data, aes(x = k, y = value, color = Model, group = Model)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~Metric, scales = "free_y") +
    scale_color_manual(values = premium_colors) +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "sans", color = "#1e293b"),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 12)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#64748b", margin = margin(b = 16)),
      strip.background = element_rect(fill = "#f8fafc", color = NA),
      strip.text = element_text(face = "bold", size = 12, color = "#334155"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f5f9"),
      panel.background = element_rect(fill = "white", color = "#cbd5e1"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold"),
      legend.margin = margin(t = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_blank()
    ) +
    labs(
      title = paste("Treatment Sensitivity (Normalized):", dgp_name),
      subtitle = "Varying treatment magnitude multiplier k. RMSE and Bias normalized by BCF-Linear (1.0)",
      x = "Treatment Effect Multiplier (k)"
    )
  
  # Add nominal 95% line for coverage facet
  coverage_nominal <- data.frame(
    Metric = "CATE 95% Coverage",
    yintercept = 0.95
  )
  p <- p + geom_hline(data = coverage_nominal, aes(yintercept = yintercept), linetype = "dashed", color = "#475569", alpha = 0.7)
  
  # Add baseline 1.0 horizontal line for the normalized facets
  baseline_one <- data.frame(
    Metric = c("Normalized CATE RMSE (Linear = 1.0)", "Normalized ATE Abs Bias (Linear = 1.0)"),
    yintercept = 1.0
  )
  p <- p + geom_hline(data = baseline_one, aes(yintercept = yintercept), linetype = "dotted", color = "#94a3b8", linewidth = 1.0, alpha = 0.8)
  
  ggsave(file.path(RESULTS_DIR, file_name), plot = p, width = 11, height = 8.5, dpi = 300)
  cat(sprintf("[SUCCESS] Graph saved to: %s\n", file.path(RESULTS_DIR, file_name)))
}

# Generate plots
plot_dgp_sensitivity("DGP A: Log-Normal Hurdle", "dgp_a_treatment_sensitivity.png")
plot_dgp_sensitivity("DGP B: Gaussian Hurdle", "dgp_b_treatment_sensitivity.png")
plot_dgp_sensitivity("DGP C: Tweedie Compound", "dgp_c_treatment_sensitivity.png")

cat("\n[SUCCESS] Treatment sensitivity analysis complete. Plots generated successfully!\n")
