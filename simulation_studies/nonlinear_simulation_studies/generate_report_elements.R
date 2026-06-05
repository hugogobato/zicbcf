################################################################################
##  Generate Report Tables and Figures from NONLINEAR Simulation Results
##
##  Nonlinear-DGP stress test (see nonlinear_dgps.R / README.md). Five models,
##  same column schema as the linear study, but spread across two flat dirs:
##    results/                                  results_notebooks/
##      - BCF-Linear     (Linear_* in results_bcf_*)        - Gamma Hurdle (GammaHurdle_* in results_gamma_hurdle_*)
##      - DPglm          (DPglm_*  in results_dpglm_*)       - Gamma +.01   (GammaP01_*    in results_gamma_plus01_*)
##      - ZIC-BCF-Smear  (Smear_*  in results_zicbcf_smear_*)
##
##  Config naming:  N-sensitivity (ZI lvl 3) = "N<n>" for n in {100,250,500,1000};
##                  ZI-sensitivity (N=500)   = "N500_lvl_<k>" for k in {1,2,4,5}
##                  (the N=500 standard file is reserved for ZI level 3).
##
##  Mirrors simulation_studies/generate_report_elements.R but with the nonlinear
##  study's flat directory layout (no Normal/ / ZI_*/ / results_gamma/folderNN/).
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

RESULTS_DIR <- "/home/hugo_souto/Stuff/Research/zicbcf/nonlinear_simulation_studies/results"
NB_DIR      <- "/home/hugo_souto/Stuff/Research/zicbcf/nonlinear_simulation_studies/results_notebooks"

cat("=== Initializing Report Elements Generator (Nonlinear DGP Stress Test, All Models) ===\n")

## ---- Model registry ---------------------------------------------------------
## Plot/legend/factor order. ZIC-BCF-Smear listed last (the proposed estimator).
MODELS <- c("BCF-Linear", "DPglm", "Gamma Hurdle", "Gamma +.01", "ZIC-BCF-Smear")

PREFIX_OF <- c(
  "BCF-Linear"    = "Linear",
  "DPglm"         = "DPglm",
  "Gamma Hurdle"  = "GammaHurdle",
  "Gamma +.01"    = "GammaP01",
  "ZIC-BCF-Smear" = "Smear"
)

## file-name stub + directory per model
STUB_OF <- c(
  "BCF-Linear"    = "bcf",
  "DPglm"         = "dpglm",
  "Gamma Hurdle"  = "gamma_hurdle",
  "Gamma +.01"    = "gamma_plus01",
  "ZIC-BCF-Smear" = "zicbcf_smear"
)
DIR_OF <- c(
  "BCF-Linear"    = RESULTS_DIR,
  "DPglm"         = RESULTS_DIR,
  "Gamma Hurdle"  = NB_DIR,
  "Gamma +.01"    = NB_DIR,
  "ZIC-BCF-Smear" = RESULTS_DIR
)

## ---- Path resolution --------------------------------------------------------
## A scenario cell is identified by (model, dgp abbrev, N, ZI level).
## N-sensitivity holds ZI level = 3; ZI-sensitivity holds N = 500.
## ZI level 3 at N=500 is the shared "standard" cell ("N500", no lvl suffix).
config_tag <- function(n_val, zi_level) {
  if (zi_level == 3) sprintf("N%d", n_val) else sprintf("N500_lvl_%d", zi_level)
}

get_path <- function(model, abbrev, n_val, zi_level) {
  file.path(DIR_OF[[model]],
            sprintf("results_%s_%s_%s.csv",
                    STUB_OF[[model]], abbrev, config_tag(n_val, zi_level)))
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
  est_ate    <- col("_Est_ATE")
  cate_rmse  <- col("_CATE_RMSE")
  ate_abserr <- abs(est_ate - df$True_ATE)

  ## A seed is "diverged" if its per-seed CATE RMSE overflows the ~O(1-6) true
  ## CATE scale (penalized-IRLS / exp-link overflow in the parametric GLMs).
  DIVERGE_THRESH <- 100
  n_div <- sum(cate_rmse > DIVERGE_THRESH, na.rm = TRUE)

  data.frame(
    DGP = dgp_name, N = n_val, ZI_Level = zi_level, Zero_Proportion = zero_prop,
    Model = model,
    CATE_RMSE        = mean(cate_rmse,                na.rm = TRUE),
    CATE_Abs_Bias    = mean(col("_CATE_Abs_Bias"),    na.rm = TRUE),
    CATE_Coverage    = mean(col("_CATE_Coverage"),    na.rm = TRUE),
    CATE_Correlation = mean(col("_CATE_Correlation"), na.rm = TRUE),
    CATE_CI_Length   = mean(col("_CATE_CI_Length"),   na.rm = TRUE),
    ATE_RMSE         = sqrt(mean((est_ate - df$True_ATE)^2, na.rm = TRUE)),
    ATE_Abs_Bias     = abs(mean(est_ate - df$True_ATE, na.rm = TRUE)),
    ## Divergence-robust companions (median over seeds) + instability diagnostics
    CATE_RMSE_Median = median(cate_rmse,  na.rm = TRUE),
    ATE_AbsErr_Median = median(ate_abserr, na.rm = TRUE),
    N_Seeds          = sum(!is.na(cate_rmse)),
    N_Diverged       = n_div,
    Diverge_Rate     = n_div / sum(!is.na(cate_rmse)),
    stringsAsFactors = FALSE
  )
}

## ---- Mapping of scenarios and nominal zero proportions ----------------------
## Realized zero proportions for the recalibrated NONLINEAR ZI grid (README §7).
get_zero_prop <- function(dgp, zi_level) {
  if (dgp == "DGP A: Log-Normal Hurdle") {
    props <- c("1" = 0.85, "2" = 0.62, "3" = 0.40, "4" = 0.18, "5" = 0.07)
  } else if (dgp == "DGP B: Gamma Hurdle") {
    props <- c("1" = 0.85, "2" = 0.61, "3" = 0.40, "4" = 0.18, "5" = 0.05)
  } else if (dgp == "DGP C: Tweedie Semicontinuous") {
    props <- c("1" = 0.85, "2" = 0.60, "3" = 0.40, "4" = 0.11, "5" = 0.03)
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

## Pretty-printer: Gamma +.01 diverges (penalized-IRLS overflow) on the nonlinear
## DGPs, so its mean RMSE can be astronomically large. Show those in scientific
## notation; round everything else to 3 dp.
fmt_rmse <- function(x) ifelse(abs(x) >= 1e4, formatC(x, format = "e", digits = 2), as.character(round(x, 3)))

## =============================================================================
##  1. STANDARD MAIN RUNS (N = 500, ZI Level 3 / ~40% zeros)
## =============================================================================
cat("\nAggregating Standard Main Runs (N = 500, ZI Level 3, ~40% zeros)...\n")
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
    CATE_RMSE        = fmt_rmse(CATE_RMSE),
    CATE_Abs_Bias    = fmt_rmse(CATE_Abs_Bias),
    CATE_Coverage    = percent(CATE_Coverage, accuracy = 0.1),
    CATE_Correlation = round(CATE_Correlation, 3),
    CATE_CI_Length   = fmt_rmse(CATE_CI_Length),
    ATE_RMSE         = fmt_rmse(ATE_RMSE),
    ATE_Abs_Bias     = fmt_rmse(ATE_Abs_Bias)
  )
print(knitr::kable(std_table_md, format = "markdown"))
write.csv(std_df %>% arrange(DGP, Model),
          file.path(RESULTS_DIR, "standard_summary.csv"), row.names = FALSE)

## Numerical-stability view: the parametric Gamma benchmarks overflow on a
## minority of the 100 seeds (which is what inflates their MEAN RMSE above).
## The MEDIAN RMSE is the divergence-robust "typical" performance.
cat("\n--- NUMERICAL STABILITY (standard cell): divergence rate + median RMSE ---\n")
stab_md <- std_df %>%
  transmute(DGP, Model,
            Diverged = sprintf("%d / %d", N_Diverged, N_Seeds),
            Diverge_Rate = percent(Diverge_Rate, accuracy = 1),
            CATE_RMSE_Median = round(CATE_RMSE_Median, 3),
            ATE_AbsErr_Median = round(ATE_AbsErr_Median, 3),
            CATE_Correlation = round(CATE_Correlation, 3))
print(knitr::kable(stab_md, format = "markdown"))

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
              CATE_RMSE = fmt_rmse(CATE_RMSE),
              CATE_Cov = round(CATE_Coverage, 3),
              CATE_Corr = round(CATE_Correlation, 3),
              CATE_CI = fmt_rmse(CATE_CI_Length),
              ATE_RMSE = fmt_rmse(ATE_RMSE)),
  format = "markdown"))

cat("\n--- FULL ZI-SENSITIVITY TABLE ---\n")
print(knitr::kable(
  zi_df %>% arrange(DGP, Model, ZI_Level) %>%
    transmute(DGP, Model, ZI_Level, Zero_Proportion,
              CATE_RMSE = fmt_rmse(CATE_RMSE),
              CATE_Cov = round(CATE_Coverage, 3),
              CATE_Corr = round(CATE_Correlation, 3),
              CATE_CI = fmt_rmse(CATE_CI_Length),
              ATE_RMSE = fmt_rmse(ATE_RMSE)),
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

nx    <- scale_x_continuous(breaks = n_list, trans = "log10")
nxlab <- "Sample Size (N) - Log Scale"
zx    <- scale_x_continuous(labels = percent_format(accuracy = 1))
zxlab <- "Zero Proportion (Zero Inflation Degree)"

## --- 4.0 Renderer (called once for all models, once excluding Gamma +.01) -----
## Gamma +.01's penalized-IRLS fit diverges on the nonlinear DGPs (CATE RMSE up
## to ~1e13), which swamps every axis. The "_no_gp01" pass drops it so the
## remaining four estimators stay legible. Normalization is robust to a missing
## BCF-Linear baseline cell (the line just gaps there).
norm_to_linear <- function(x, model) {
  base <- x[model == "BCF-Linear"]
  if (length(base) == 0 || is.na(base[1])) return(rep(NA_real_, length(x)))
  x / base[1]
}

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
             paste0("Nonlinear DGP | Causal performance across N in {100, 250, 500, 1000}", note),
             sprintf("n_sensitivity_%s%s.png", abbr, suffix))
  }
  ## ---- ZI-sensitivity (Normalized RMSE; BCF-Linear = 1.0) ----
  zi_norm <- zi_df %>%
    filter(Model %in% models_keep) %>% mutate(Model = droplevels(Model)) %>%
    group_by(DGP, ZI_Level) %>%
    mutate(
      CATE_RMSE_Norm = norm_to_linear(CATE_RMSE, Model),
      ATE_RMSE_Norm  = norm_to_linear(ATE_RMSE,  Model)
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
             paste0("Nonlinear DGP | RMSE normalized to BCF-Linear (=1.0) at each level", note),
             sprintf("zi_sensitivity_%s%s.png", abbr, suffix))
  }
}

## Full set (all 5 models) -- archival; swamped because BOTH Gamma benchmarks
## diverge on the nonlinear DGPs (mean RMSE up to ~1e13).
render_sensitivity(MODELS, "")
## Legible set: drop BOTH diverging parametric benchmarks, leaving the three
## numerically-stable estimators. (Unlike the linear study, Gamma Hurdle also
## diverges here, so dropping only Gamma +.01 would not be enough.)
render_sensitivity(c("BCF-Linear", "DPglm", "ZIC-BCF-Smear"), "_stable",
                   note = "  (numerically-stable estimators; both Gamma benchmarks diverge)")

cat("\n=== All Report Elements Generated Successfully! ===\n")
