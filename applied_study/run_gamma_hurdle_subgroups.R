library(zicbcf)

# Subgroup analysis for the Gamma Hurdle (Oganisian et al. 2019) fit on the
# MEPS dental-insurance data. Mirrors run_cate_subgroups.R, but uses the
# posterior draws saved by run_gamma_hurdle.R and reports both the overall
# expenditure CATE and the hurdle-part CATE (probability of any expenditure).

cat("Loading dataset and saved draws...\n")
df <- read.csv("applied_study/h251.csv")

covars <- c("AGE23X", "SEX", "RACEV2X", "FAMINC23", "POVCAT23", "REGION23", "MARRY23X")
df_clean <- df[df$DNTINS23_M23 > 0, ]
for (col in covars) {
  df_clean <- df_clean[df_clean[[col]] >= 0, ]
}
df_clean <- df_clean[df_clean$DVTEXP23 >= 0, ]

draws <- readRDS("applied_study/gamma_hurdle_draws.rds")
cate_draws <- draws$cate_draws            # nsim x n  overall CATE (expected outcome)
hurdle_cate_draws <- draws$hurdle_cate_draws  # nsim x n  hurdle-part CATE (Pr(y>0))

# subgroup posterior: average unit-level effects within each draw
subgroup_pair <- function(draws_mat, idx) {
  sub <- draws_mat[, idx, drop = FALSE]
  m <- rowMeans(sub)
  c(N       = ncol(sub),
    Est     = mean(m),
    CI_low  = unname(quantile(m, 0.025)),
    CI_high = unname(quantile(m, 0.975)),
    P_gt_0  = mean(m > 0))
}

contrast_pair <- function(draws_mat, idx_a, idx_b) {
  da <- rowMeans(draws_mat[, idx_a, drop = FALSE])
  db <- rowMeans(draws_mat[, idx_b, drop = FALSE])
  d  <- da - db
  c(Diff    = mean(d),
    CI_low  = unname(quantile(d, 0.025)),
    CI_high = unname(quantile(d, 0.975)),
    P_gt_0  = mean(d > 0))
}

# subgroup labels (same scheme as run_cate_subgroups.R)
sex_lab <- factor(df_clean$SEX, levels = c(1, 2), labels = c("Male", "Female"))
race_lab <- factor(ifelse(df_clean$RACEV2X == 1, "White",
                   ifelse(df_clean$RACEV2X == 2, "Black",
                   ifelse(df_clean$RACEV2X == 4, "Asian", "Other/Multiple"))),
                   levels = c("White", "Black", "Asian", "Other/Multiple"))
pov_lab <- factor(df_clean$POVCAT23, levels = 1:5,
                  labels = c("Poor", "Near-poor", "Low-income",
                             "Middle-income", "High-income"))
region_lab <- factor(df_clean$REGION23, levels = 1:4,
                     labels = c("Northeast", "Midwest", "South", "West"))
marr_lab <- factor(ifelse(df_clean$MARRY23X == 1, "Married", "Not married"),
                   levels = c("Married", "Not married"))
age_lab <- cut(df_clean$AGE23X,
               breaks = c(-Inf, 17, 34, 49, 64, Inf),
               labels = c("0-17", "18-34", "35-49", "50-64", "65+"))
inc_lab <- cut(df_clean$FAMINC23,
               breaks = quantile(df_clean$FAMINC23, c(0, 1/3, 2/3, 1)),
               labels = c("Low-income tertile", "Middle-income tertile",
                          "High-income tertile"),
               include.lowest = TRUE)

group_vars <- list(
  Sex            = sex_lab,
  Race           = race_lab,
  `Age group`    = age_lab,
  `Poverty cat.` = pov_lab,
  `Family income (tertile)` = inc_lab,
  Region         = region_lab,
  `Marital status` = marr_lab
)

# subgroup table for both the overall CATE and the hurdle CATE
build_table <- function(draws_mat) {
  rows <- list()
  rows[[1]] <- data.frame(Covariate = "Overall", Level = "All",
                          t(subgroup_pair(draws_mat, rep(TRUE, ncol(draws_mat)))),
                          check.names = FALSE)
  for (gname in names(group_vars)) {
    lab <- group_vars[[gname]]
    for (lv in levels(lab)) {
      idx <- which(lab == lv)
      if (length(idx) < 30) next
      rows[[length(rows) + 1]] <- data.frame(
        Covariate = gname, Level = lv,
        t(subgroup_pair(draws_mat, idx)), check.names = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("Est", "CI_low", "CI_high")] <- round(out[, c("Est", "CI_low", "CI_high")], 4)
  out$P_gt_0 <- round(out$P_gt_0, 3)
  out
}

build_contrasts <- function(draws_mat) {
  crows <- list()
  for (gname in names(group_vars)) {
    lab <- group_vars[[gname]]
    lvs <- levels(lab)[table(lab) >= 30]
    if (length(lvs) < 2) next
    cmb <- combn(lvs, 2)
    for (k in seq_len(ncol(cmb))) {
      a <- cmb[1, k]; b <- cmb[2, k]
      ct <- contrast_pair(draws_mat, which(lab == a), which(lab == b))
      crows[[length(crows) + 1]] <- data.frame(
        Covariate = gname, Contrast = paste(a, "-", b),
        t(ct), check.names = FALSE)
    }
  }
  out <- do.call(rbind, crows)
  rownames(out) <- NULL
  out[, c("Diff", "CI_low", "CI_high")] <- round(out[, c("Diff", "CI_low", "CI_high")], 4)
  out$P_gt_0 <- round(out$P_gt_0, 3)
  out
}

cate_table <- build_table(cate_draws)
hurdle_table <- build_table(hurdle_cate_draws)
contrast_table <- build_contrasts(cate_draws)
hurdle_contrast_table <- build_contrasts(hurdle_cate_draws)

write.csv(cate_table, "applied_study/gamma_hurdle_cate_subgroups.csv", row.names = FALSE)
write.csv(hurdle_table, "applied_study/gamma_hurdle_hurdle_subgroups.csv", row.names = FALSE)
write.csv(contrast_table, "applied_study/gamma_hurdle_cate_contrasts.csv", row.names = FALSE)
write.csv(hurdle_contrast_table, "applied_study/gamma_hurdle_hurdle_contrasts.csv", row.names = FALSE)

cat("\n=== Subgroup CATE estimates (overall expenditure, $) ===\n")
print(cate_table)
cat("\n=== Subgroup CATE estimates (hurdle, probability points) ===\n")
print(hurdle_table)
cat("\n=== Between-group CATE contrasts (overall expenditure, $) ===\n")
print(contrast_table)
cat("\n=== Between-group CATE contrasts (hurdle, probability points) ===\n")
print(hurdle_contrast_table)

# Forest plot of subgroup overall CATEs
plot_df <- cate_table[cate_table$Covariate != "Overall", ]
plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
plot_df <- plot_df[nrow(plot_df):1, ]
overall_cate <- cate_table$Est[cate_table$Covariate == "Overall"]

png("applied_study/gamma_hurdle_cate_subgroups.png", width = 1000, height = 1100)
par(mar = c(5, 16, 4, 2))
yy <- seq_len(nrow(plot_df))
plot(plot_df$Est, yy, xlim = range(c(plot_df$CI_low, plot_df$CI_high)),
     yaxt = "n", pch = 19, col = "steelblue", cex = 1.3,
     xlab = "CATE: effect of dental insurance on annual dental expenditure ($)",
     ylab = "", main = "Gamma Hurdle - Subgroup CATE estimates with 95% CrI")
axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.9)
segments(plot_df$CI_low, yy, plot_df$CI_high, yy, col = "steelblue", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
abline(v = overall_cate, lty = 3, col = "firebrick", lwd = 2)
legend("bottomright",
       legend = c("Subgroup CATE (95% CrI)", "No effect",
                  sprintf("Overall ATE = $%.0f", overall_cate)),
       col = c("steelblue", "grey50", "firebrick"),
       pch = c(19, NA, NA), lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
dev.off()

cat("\nCompleted successfully.\n")
