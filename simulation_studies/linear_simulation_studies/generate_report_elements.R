################################################################################
##  Generate Report Tables and Figures from Simulation Results (All Models)
##
##  Models compared:
##    - BCF-Linear      (Linear_*      cols in results_bcf_*.csv)
##    - DPglm           (DPglm_*       cols in results_dpglm_*.csv)
##    - Gamma Hurdle    (GammaHurdle_* cols in results_gamma/folder*/results_gamma_hurdle_*.csv)
##    - Gamma +.01      (GammaP01_*    cols in results_gamma/folder*/results_gamma_plus01_*.csv)
##    - ZIC-BCF-Smear   (Smear_*       cols in results_bcf_*.csv)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

RESULTS_DIR <- "/home/hugo_souto/Stuff/Research/zicbcf/simulation_studies/results"
GAMMA_DIR   <- "/home/hugo_souto/Stuff/Research/zicbcf/simulation_studies/results_gamma"
OUTPUT_DIR  <- "/home/hugo_souto/Stuff/Research/zicbcf/simulation_studies"

cat("=== Initializing Report Elements Generator (All Benchmark Models) ===\n")

## ---- Model registry ---------------------------------------------------------
## Plot/legend/factor order. ZIC-BCF-Smear listed last (the proposed estimator).
MODELS <- c("BCF-Linear", "DPglm", "Gamma Hurdle", "Gamma +.01", "ZIC-BCF-Smear")

PREFIX_OF <- c(
  "BCF-Linear"    = "Linear",
  "ZIC-BCF-Smear" = "Smear",
  "DPglm"         = "DPglm",
  "Gamma Hurdle"  = "GammaHurdle",
  "Gamma +.01"    = "GammaP01"
)

## Verified results_gamma/ folder map (folder -> (kind, N, ZI level)):
##   hurdle: 33=N500/ZI3 34=N1000/ZI3 35=N100/ZI3 36=N250/ZI3 | 37=ZI1 38=ZI2 39=ZI4 40=ZI5 (all N500)
##   plus01: 41=N500/ZI3 42=N1000/ZI3 43=N100/ZI3 44=N250/ZI3 | 45=ZI1 46=ZI2 47=ZI4 48=ZI5 (all N500)
GAMMA_FOLDER <- list(
  hurdle = c("N3_500" = 33, "N3_1000" = 34, "N3_100" = 35, "N3_250" = 36,
             "ZI1" = 37, "ZI2" = 38, "ZI4" = 39, "ZI5" = 40),
  plus01 = c("N3_500" = 41, "N3_1000" = 42, "N3_100" = 43, "N3_250" = 44,
             "ZI1" = 45, "ZI2" = 46, "ZI4" = 47, "ZI5" = 48)
)

## ---- Path resolution --------------------------------------------------------
## A scenario cell is identified by (model, dgp abbrev, N, ZI level).
## N-sensitivity holds ZI level = 3; ZI-sensitivity holds N = 500.
base_loc <- function(stub, n_val, zi_level) {
  if (zi_level == 3) {
    if (n_val == 500) file.path(RESULTS_DIR, "Normal", stub)
    else              file.path(RESULTS_DIR, stub)
  } else {
    file.path(RESULTS_DIR, sprintf("ZI_%d", zi_level), stub)
  }
}

gamma_loc <- function(kind, abbrev, n_val, zi_level) {
  fol <- GAMMA_FOLDER[[kind]]
  key <- if (zi_level == 3) sprintf("N3_%d", n_val) else sprintf("ZI%d", zi_level)
  folder <- fol[[key]]
  n_in_name <- if (zi_level == 3) n_val else 500   # ZI cells are always N=500
  stub <- sprintf("results_gamma_%s_%s_N%d.csv", kind, abbrev, n_in_name)
  file.path(GAMMA_DIR, sprintf("folder%d", folder), stub)
}

get_path <- function(model, abbrev, n_val, zi_level) {
  if (model %in% c("BCF-Linear", "ZIC-BCF-Smear")) {
    base_loc(sprintf("results_bcf_%s_N%d.csv", abbrev, n_val), n_val, zi_level)
  } else if (model == "DPglm") {
    base_loc(sprintf("results_dpglm_%s_N%d.csv", abbrev, n_val), n_val, zi_level)
  } else if (model == "Gamma Hurdle") {
    gamma_loc("hurdle", abbrev, n_val, zi_level)
  } else if (model == "Gamma +.01") {
    gamma_loc("plus01", abbrev, n_val, zi_level)
  } else stop("unknown model: ", model)
}

## ---- Metric extraction for a single (model, cell) ---------------------------
extract_metrics <- function(model, dgp_name, n_val, zi_level, zero_prop) {
  abbrev <- dgp_abbrev[[dgp_name]]
  path   <- get_path(model, abbrev, n_val, zi_level)
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s", path)); return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if ("Status" %in% names(df)) df <- df[df$Status == "OK", , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  pre <- PREFIX_OF[[model]]
  col <- function(suffix) df[[paste0(pre, suffix)]]
  est_ate <- col("_Est_ATE")

  data.frame(
    DGP = dgp_name, N = n_val, ZI_Level = zi_level, Zero_Proportion = zero_prop,
    Model = model,
    CATE_RMSE        = mean(col("_CATE_RMSE"),        na.rm = TRUE),
    CATE_Abs_Bias    = mean(col("_CATE_Abs_Bias"),    na.rm = TRUE),
    CATE_Coverage    = mean(col("_CATE_Coverage"),    na.rm = TRUE),
    CATE_Correlation = mean(col("_CATE_Correlation"), na.rm = TRUE),
    CATE_CI_Length   = mean(col("_CATE_CI_Length"),   na.rm = TRUE),
    ATE_RMSE         = sqrt(mean((est_ate - df$True_ATE)^2, na.rm = TRUE)),
    ATE_Abs_Bias     = abs(mean(est_ate - df$True_ATE, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

## ---- Mapping of scenarios and nominal zero proportions ----------------------
get_zero_prop <- function(dgp, zi_level) {
  if (dgp == "DGP A: Log-Normal Hurdle") {
    props <- c("1" = 0.854, "2" = 0.619, "3" = 0.367, "4" = 0.185, "5" = 0.066)
  } else if (dgp == "DGP B: Gamma Hurdle") {
    props <- c("1" = 0.851, "2" = 0.607, "3" = 0.395, "4" = 0.181, "5" = 0.049)
  } else if (dgp == "DGP C: Tweedie Semicontinuous") {
    props <- c("1" = 0.600, "2" = 0.397, "3" = 0.180, "4" = 0.115, "5" = 0.028)
  } else return(0)
  props[as.character(zi_level)]
}

dgps <- c("DGP A: Log-Normal Hurdle", "DGP B: Gamma Hurdle", "DGP C: Tweedie Semicontinuous")
dgp_abbrev <- c("DGP A: Log-Normal Hurdle" = "dgp_a",
                "DGP B: Gamma Hurdle" = "dgp_b",
                "DGP C: Tweedie Semicontinuous" = "dgp_c")

collect <- function(cells) {
  out <- lapply(seq_len(nrow(cells)), function(i) {
    r <- cells[i, ]
    extract_metrics(r$model, r$dgp, r$n_val, r$zi_level, get_zero_prop(r$dgp, r$zi_level))
  })
  do.call(rbind, out)
}

## =============================================================================
##  1. STANDARD MAIN RUNS (N = 500, ZI Level 3 / Normal)
## =============================================================================
cat("\nAggregating Standard Main Runs (N = 500, ZI Level 3)...\n")
std_cells <- expand.grid(dgp = dgps, model = MODELS, n_val = 500, zi_level = 3,
                         stringsAsFactors = FALSE)
std_df <- collect(std_cells)
std_df$Model <- factor(std_df$Model, levels = MODELS)
std_df <- std_df %>% arrange(DGP, Model)

cat("\n--- STANDARD COMPARATIVE RESULTS TABLE (N=500, ZI Level 3) ---\n")
std_table_md <- std_df %>%
  select(DGP, Model, CATE_RMSE, CATE_Abs_Bias, CATE_Coverage,
         CATE_Correlation, CATE_CI_Length, ATE_RMSE, ATE_Abs_Bias) %>%
  mutate(
    CATE_RMSE        = round(CATE_RMSE, 3),
    CATE_Abs_Bias    = round(CATE_Abs_Bias, 3),
    CATE_Coverage    = percent(CATE_Coverage, accuracy = 0.1),
    CATE_Correlation = round(CATE_Correlation, 3),
    CATE_CI_Length   = round(CATE_CI_Length, 3),
    ATE_RMSE         = round(ATE_RMSE, 3),
    ATE_Abs_Bias     = round(ATE_Abs_Bias, 3)
  )
print(knitr::kable(std_table_md, format = "markdown"))

## =============================================================================
##  2. SAMPLE SIZE N SENSITIVITY (ZI Level 3; N in {100,250,500,1000})
## =============================================================================
cat("\n\nAggregating N-Sensitivity Runs...\n")
n_list <- c(100, 250, 500, 1000)
n_cells <- expand.grid(dgp = dgps, model = MODELS, n_val = n_list, zi_level = 3,
                       stringsAsFactors = FALSE)
n_df <- collect(n_cells)
n_df$Model <- factor(n_df$Model, levels = MODELS)

## =============================================================================
##  3. ZERO-INFLATION SENSITIVITY (N = 500; ZI levels 1..5)
## =============================================================================
cat("Aggregating Zero-Inflation Sensitivity Runs...\n")
zi_cells <- expand.grid(dgp = dgps, model = MODELS, n_val = 500, zi_level = 1:5,
                        stringsAsFactors = FALSE)
zi_df <- collect(zi_cells)
zi_df$Model <- factor(zi_df$Model, levels = MODELS)

## Save aggregated sensitivity outputs to CSVs (all models)
write.csv(n_df %>% arrange(DGP, Model, N),
          file.path(RESULTS_DIR, "n_sensitivity_summary.csv"), row.names = FALSE)
write.csv(zi_df %>% arrange(DGP, Model, ZI_Level),
          file.path(RESULTS_DIR, "zi_sensitivity_summary.csv"), row.names = FALSE)
cat("[SUCCESS] Aggregated sensitivity CSVs saved to results/.\n")

## Dump full long tables to stdout for the markdown write-up
cat("\n--- FULL N-SENSITIVITY TABLE ---\n")
print(knitr::kable(
  n_df %>% arrange(DGP, Model, N) %>%
    transmute(DGP, Model, N,
              CATE_RMSE = round(CATE_RMSE, 3),
              CATE_Cov = round(CATE_Coverage, 3),
              CATE_Corr = round(CATE_Correlation, 3),
              CATE_CI = round(CATE_CI_Length, 3),
              ATE_RMSE = round(ATE_RMSE, 3)),
  format = "markdown"))

cat("\n--- FULL ZI-SENSITIVITY TABLE ---\n")
print(knitr::kable(
  zi_df %>% arrange(DGP, Model, ZI_Level) %>%
    transmute(DGP, Model, ZI_Level, Zero_Proportion,
              CATE_RMSE = round(CATE_RMSE, 3),
              CATE_Cov = round(CATE_Coverage, 3),
              CATE_Corr = round(CATE_Correlation, 3),
              CATE_CI = round(CATE_CI_Length, 3),
              ATE_RMSE = round(ATE_RMSE, 3)),
  format = "markdown"))

## =============================================================================
##  4. PLOTTING SENSITIVITY ANALYSES (2x2 GRID, 4 PANELS PER IMAGE)
## =============================================================================
## Built with patchwork (one ggplot per metric) so each panel sets its own
## y-scale. Credible-interval-length panels use a log10 y-axis because DPglm's
## intervals are an order of magnitude wider than the other models'.
cat("\nGenerating 2x2 grid visualizations (patchwork)...\n")

premium_colors <- c(
  "BCF-Linear"    = "#64748b",  # cool slate grey
  "DPglm"         = "#0891b2",  # teal/cyan
  "Gamma Hurdle"  = "#f59e0b",  # amber
  "Gamma +.01"    = "#ef4444",  # red (naive)
  "ZIC-BCF-Smear" = "#8b5cf6"   # premium purple (proposed)
)

panel_theme <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "#1e293b"),
    plot.title = element_text(face = "bold", size = 11.5, color = "#334155",
                             hjust = 0.5, margin = margin(b = 6)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f1f5f9"),
    panel.background = element_rect(fill = "white", color = "#cbd5e1"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(face = "bold", size = 10.5)
  )

## A single metric panel.
##   xscale: a ggplot scale_x_* object (log N axis, or percent axis)
##   logy  : TRUE -> log10 y-axis (for credible-interval lengths)
##   hline : optional numeric (e.g. 0.95 nominal, 1.0 baseline) + linetype
metric_panel <- function(d, yvar, title, xvar, xscale, xlab,
                         logy = FALSE, hline = NA, hlt = "dashed") {
  p <- ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]],
                     color = Model, group = Model)) +
    geom_line(linewidth = 1.0, alpha = 0.85) +
    geom_point(size = 2.3, alpha = 0.9) +
    scale_color_manual(values = premium_colors, drop = FALSE) +
    xscale +
    labs(title = title, x = xlab, y = NULL) +
    panel_theme
  if (logy) p <- p + scale_y_log10(labels = label_number(accuracy = 0.1))
  if (!is.na(hline)) p <- p +
    geom_hline(yintercept = hline, linetype = hlt, color = "#475569", alpha = 0.7)
  p
}

assemble <- function(p1, p2, p3, p4, title, subtitle, file_name) {
  combined <- (p1 | p2) / (p3 | p4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                  margin = margin(b = 6)),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#64748b",
                                     margin = margin(b = 10))
      )
    ) & theme(legend.position = "bottom")
  ggsave(file.path(RESULTS_DIR, file_name), plot = combined,
         width = 10, height = 8.2, dpi = 300)
  cat(sprintf("[SUCCESS] Graph saved to: %s\n", file.path(RESULTS_DIR, file_name)))
}

## --- 4.0 Renderer (called once for all models, once excluding Gamma +.01) -----
## Gamma +.01's penalized-IRLS fit diverges in a handful of small-N / extreme-ZI
## cells (CATE RMSE up to ~1e9), which swamps every axis. The "_no_gp01" pass
## drops it so the remaining four estimators stay legible.
nx    <- scale_x_continuous(breaks = n_list, trans = "log10")
nxlab <- "Sample Size (N) - Log Scale"
zx    <- scale_x_continuous(labels = percent_format(accuracy = 1))
zxlab <- "Zero Proportion (Zero Inflation Degree)"

render_sensitivity <- function(models_keep, suffix, note = "") {
  ## ---- N-sensitivity ----
  nsub <- n_df %>% filter(Model %in% models_keep) %>% mutate(Model = droplevels(Model))
  for (dgp in dgps) {
    d <- nsub %>% filter(DGP == dgp)
    abbr <- dgp_abbrev[dgp]
    p1 <- metric_panel(d, "CATE_RMSE",      "CATE RMSE (lower is better)",            "N", nx, nxlab)
    p2 <- metric_panel(d, "CATE_Coverage",  "CATE 95% Coverage (Nominal = 0.95)",     "N", nx, nxlab, hline = 0.95)
    p3 <- metric_panel(d, "ATE_RMSE",       "ATE RMSE (lower is better)",             "N", nx, nxlab)
    p4 <- metric_panel(d, "CATE_CI_Length", "CATE Cred. Interval Length (log scale)", "N", nx, nxlab, logy = TRUE)
    assemble(p1, p2, p3, p4,
             paste("Sample Size (N) Sensitivity Analysis:", dgp),
             paste0("Causal performance across N in {100, 250, 500, 1000}", note),
             sprintf("n_sensitivity_%s%s.png", abbr, suffix))
  }
  ## ---- ZI-sensitivity (Normalized RMSE; BCF-Linear = 1.0) ----
  zi_norm <- zi_df %>%
    filter(Model %in% models_keep) %>% mutate(Model = droplevels(Model)) %>%
    group_by(DGP, ZI_Level) %>%
    mutate(
      CATE_RMSE_Norm = CATE_RMSE / CATE_RMSE[Model == "BCF-Linear"],
      ATE_RMSE_Norm  = ATE_RMSE  / ATE_RMSE[Model == "BCF-Linear"]
    ) %>%
    ungroup()
  for (dgp in dgps) {
    d <- zi_norm %>% filter(DGP == dgp)
    abbr <- dgp_abbrev[dgp]
    p1 <- metric_panel(d, "CATE_RMSE_Norm", "Normalized CATE RMSE (Linear = 1.0)",    "Zero_Proportion", zx, zxlab, hline = 1.0, hlt = "dotted")
    p2 <- metric_panel(d, "CATE_Coverage",  "CATE 95% Coverage",                      "Zero_Proportion", zx, zxlab, hline = 0.95)
    p3 <- metric_panel(d, "ATE_RMSE_Norm",  "Normalized ATE RMSE (Linear = 1.0)",     "Zero_Proportion", zx, zxlab, hline = 1.0, hlt = "dotted")
    p4 <- metric_panel(d, "CATE_CI_Length", "CATE Cred. Interval Length (log scale)", "Zero_Proportion", zx, zxlab, logy = TRUE)
    assemble(p1, p2, p3, p4,
             paste("Zero-Inflation Sensitivity Analysis:", dgp),
             paste0("RMSE normalized to BCF-Linear (=1.0) at each level", note),
             sprintf("zi_sensitivity_%s%s.png", abbr, suffix))
  }
}

## Full set (all 5 models)
render_sensitivity(MODELS, "")
## Legible set (drops the diverging Gamma +.01 benchmark)
render_sensitivity(setdiff(MODELS, "Gamma +.01"), "_no_gp01", note = "  (Gamma +.01 excluded)")

cat("\n=== All Report Elements Generated Successfully! ===\n")
