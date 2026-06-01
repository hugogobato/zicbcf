################################################################################
##  Semicontinuous BCF Extended Study: 5-Seed Orchestrator (Process Isolated)
##
##  This script runs the 5-seed comparative simulation study across the
##  three DGPs for BCF-Linear, ZIC-BCF (Path A), and ZIC-BCF-Smear.
##
##  Sequential execution in isolated R processes prevents memory buildup and
##  protects system RAM.
################################################################################

library(dplyr)
library(tidyr)

RESULTS_DIR <- "simulation_studies/reduced_simulations/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

cat("=== Starting 5-Seed Semicontinuous BCF Simulation Study ===\n")
cat("RAM Protection: Process-isolated sequential execution.\n\n")

dgps <- c("A", "B", "C")
seeds <- 1:5

t_start <- Sys.time()

# Loop through all DGP and seed combinations
for (dgp in dgps) {
  for (s in seeds) {
    checkpoint_file <- file.path(RESULTS_DIR, sprintf("checkpoint_5seed_DGP_%s_seed_%d.csv", dgp, s))
    
    if (file.exists(checkpoint_file)) {
      cat(sprintf("[Checkpoint] DGP %s Seed %d already completed. Skipping.\n", dgp, s))
    } else {
      cat(sprintf("\n========================================================================\n"))
      cat(sprintf("=== Orchestrating DGP %s Seed %d (Isolated Process) ===\n", dgp, s))
      cat(sprintf("========================================================================\n"))
      
      # Run in an isolated process
      cmd <- sprintf("Rscript simulation_studies/reduced_simulations/run_5seed_single.R %s %d", dgp, s)
      exit_code <- system(cmd)
      
      if (exit_code != 0) {
        stop(sprintf("\nError: DGP %s Seed %d failed with exit code %d. Stopped.\n", dgp, s, exit_code))
      }
    }
  }
}

t_end <- Sys.time()
cat(sprintf("\n=== All 15 Simulation Runs Finished in %.1f minutes ===\n", 
            as.numeric(difftime(t_end, t_start, units = "mins"))))

## ---- Aggregate results ------------------------------------------------------
cat("\nAggregating results...\n")
results_list <- list()
idx <- 1

for (dgp in dgps) {
  for (s in seeds) {
    checkpoint_file <- file.path(RESULTS_DIR, sprintf("checkpoint_5seed_DGP_%s_seed_%d.csv", dgp, s))
    if (file.exists(checkpoint_file)) {
      results_list[[idx]] <- read.csv(checkpoint_file, stringsAsFactors = FALSE)
      idx <- idx + 1
    } else {
      stop(sprintf("Fatal Error: Missing results file for DGP %s Seed %d.", dgp, s))
    }
  }
}

all_runs_df <- do.call(rbind, results_list)
write.csv(all_runs_df, file.path(RESULTS_DIR, "raw_5seed_runs.csv"), row.names = FALSE)
cat(sprintf("[SUCCESS] Raw runs CSV saved to: %s\n", file.path(RESULTS_DIR, "raw_5seed_runs.csv")))

## ---- Compute Aggregated Summary ---------------------------------------------
cat("Computing means, standard deviations, and aggregate RMSE...\n")

# We want to pivot longer to make summary calculations easy
# Columns prefix: Linear_, PathA_, Smear_
long_runs <- all_runs_df %>%
  pivot_longer(
    cols = starts_with(c("Linear_", "PathA_", "Smear_")),
    names_to = c("Model", "Type", "Metric"),
    names_pattern = "^(Linear|PathA|Smear)_(CATE|Hurdle|Est)_(.*)$"
  )

# Wait, let's look at the columns:
# Linear_CATE_RMSE, Linear_CATE_Abs_Bias, Linear_CATE_Coverage, Linear_CATE_Correlation, Linear_CATE_CI_Length, Linear_Est_ATE
# Linear_Hurdle_RMSE, Linear_Hurdle_Abs_Bias, Linear_Hurdle_Coverage, Linear_Hurdle_Correlation, Linear_Hurdle_CI_Length, Linear_Est_Hurdle_ATE
#
# Let's write a robust parser in R to aggregate each metric properly:

models <- c("Linear", "PathA", "Smear")
metrics_summary <- list()
sum_idx <- 1

for (dgp_name in unique(all_runs_df$DGP)) {
  dgp_subset <- all_runs_df %>% filter(DGP == dgp_name)
  
  for (model in models) {
    # Extract columns for this model
    m_rmse        <- dgp_subset[[paste0(model, "_CATE_RMSE")]]
    m_bias        <- dgp_subset[[paste0(model, "_CATE_Abs_Bias")]]
    m_cov         <- dgp_subset[[paste0(model, "_CATE_Coverage")]]
    m_cor         <- dgp_subset[[paste0(model, "_CATE_Correlation")]]
    m_ci_len      <- dgp_subset[[paste0(model, "_CATE_CI_Length")]]
    m_est_ate     <- dgp_subset[[paste0(model, "_Est_ATE")]]
    m_true_ate    <- dgp_subset$True_ATE
    
    # Hurdle columns
    h_rmse        <- dgp_subset[[paste0(model, "_Hurdle_RMSE")]]
    h_bias        <- dgp_subset[[paste0(model, "_Hurdle_Abs_Bias")]]
    h_cov         <- dgp_subset[[paste0(model, "_Hurdle_Coverage")]]
    h_cor         <- dgp_subset[[paste0(model, "_Hurdle_Correlation")]]
    h_ci_len      <- dgp_subset[[paste0(model, "_Hurdle_CI_Length")]]
    h_est_ate     <- dgp_subset[[paste0(model, "_Est_Hurdle_ATE")]]
    h_true_ate    <- dgp_subset$True_Hurdle_ATE
    
    # Standard ATE RMSE (single number across 5 seeds)
    ate_rmse <- sqrt(mean((m_est_ate - m_true_ate)^2))
    ate_bias_vec <- abs(m_est_ate - m_true_ate)
    
    # Hurdle ATE RMSE (single number across 5 seeds, if not NA)
    if (any(is.na(h_est_ate))) {
      hurdle_ate_rmse <- NA
      hurdle_ate_bias_vec <- NA
    } else {
      hurdle_ate_rmse <- sqrt(mean((h_est_ate - h_true_ate)^2))
      hurdle_ate_bias_vec <- abs(h_est_ate - h_true_ate)
    }
    
    metrics_summary[[sum_idx]] <- data.frame(
      DGP = dgp_name,
      Model = case_when(
        model == "Linear" ~ "BCF-Linear",
        model == "PathA"  ~ "ZIC-BCF (Path A)",
        model == "Smear"  ~ "ZIC-BCF-Smear (Best_Path_Gemini)"
      ),
      
      # CATE Metrics: mean (SD)
      CATE_RMSE_Mean = mean(m_rmse),
      CATE_RMSE_SD   = sd(m_rmse),
      
      CATE_Bias_Mean = mean(m_bias),
      CATE_Bias_SD   = sd(m_bias),
      
      CATE_Coverage_Mean = mean(m_cov),
      CATE_Coverage_SD   = sd(m_cov),
      
      CATE_Correlation_Mean = mean(m_cor),
      CATE_Correlation_SD   = sd(m_cor),
      
      CATE_CI_Length_Mean = mean(m_ci_len),
      CATE_CI_Length_SD   = sd(m_ci_len),
      
      # ATE Metrics
      Est_ATE_Mean = mean(m_est_ate),
      Est_ATE_SD   = sd(m_est_ate),
      
      ATE_Abs_Bias_Mean = mean(ate_bias_vec),
      ATE_Abs_Bias_SD   = sd(ate_bias_vec),
      
      ATE_RMSE = ate_rmse,
      
      # Hurdle CATE Metrics
      Hurdle_CATE_RMSE_Mean = mean(h_rmse),
      Hurdle_CATE_RMSE_SD   = sd(h_rmse),
      
      Hurdle_CATE_Bias_Mean = mean(h_bias),
      Hurdle_CATE_Bias_SD   = sd(h_bias),
      
      Hurdle_CATE_Coverage_Mean = mean(h_cov),
      Hurdle_CATE_Coverage_SD   = sd(h_cov),
      
      Hurdle_CATE_Correlation_Mean = mean(h_cor),
      Hurdle_CATE_Correlation_SD   = sd(h_cor),
      
      Hurdle_CATE_CI_Length_Mean = mean(h_ci_len),
      Hurdle_CATE_CI_Length_SD   = sd(h_ci_len),
      
      # Hurdle ATE Metrics
      Est_Hurdle_ATE_Mean = mean(h_est_ate),
      Est_Hurdle_ATE_SD   = sd(h_est_ate),
      
      Hurdle_ATE_Abs_Bias_Mean = mean(hurdle_ate_bias_vec),
      Hurdle_ATE_Abs_Bias_SD   = sd(hurdle_ate_bias_vec),
      
      Hurdle_ATE_RMSE = hurdle_ate_rmse,
      
      stringsAsFactors = FALSE
    )
    
    sum_idx <- sum_idx + 1
  }
}

summary_df <- do.call(rbind, metrics_summary)
write.csv(summary_df, file.path(RESULTS_DIR, "aggregated_5seed_results.csv"), row.names = FALSE)
cat(sprintf("[SUCCESS] Aggregated 5-seed results saved to: %s\n\n", file.path(RESULTS_DIR, "aggregated_5seed_results.csv")))

# Display summary beautifully
print(summary_df %>% select(DGP, Model, CATE_RMSE_Mean, CATE_Coverage_Mean, CATE_CI_Length_Mean, ATE_RMSE, Hurdle_CATE_RMSE_Mean, Hurdle_CATE_Coverage_Mean), digits = 4)
