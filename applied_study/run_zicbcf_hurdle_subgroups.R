library(zicbcf)

# ==========================================================================
# Hurdle-margin (participation) treatment-effect analysis for ZIC-BCF-Smear
# on the MEPS dental-insurance data.
#
# In addition to the overall dollar-scale CATE (fit$cate), ZIC-BCF-Smear's
# probit hurdle stage yields a posterior for the participation effect on the
# probability of any dental expenditure. Following the simulation-study code,
# the hurdle-margin CATE draws are:
#     p0 = pnorm(mu_b),  p1 = pnorm(mu_b + tau_b),  hurdle_cate = p1 - p0,
# where mu_b (= fit$mu_b) and tau_b (= fit$tau_b) are the probit-scale
# prognostic and treatment-effect draws (nsim x n). Subgroup CATEs and
# between-group contrasts are built draw-by-draw, exactly as for the overall
# dollar-scale CATE in run_cate_subgroups.R and for the Gamma Hurdle margin
# in run_gamma_hurdle_subgroups.R, so the three analyses are comparable.
# ==========================================================================

set.seed(1)

cat("Loading dataset...\n")
df <- read.csv("applied_study/h251.csv")

covars <- c("AGE23X", "SEX", "RACEV2X", "FAMINC23", "POVCAT23", "REGION23", "MARRY23X")
df_clean <- df[df$DNTINS23_M23 > 0, ]
for (col in covars) {
  df_clean <- df_clean[df_clean[[col]] >= 0, ]
}
df_clean <- df_clean[df_clean$DVTEXP23 >= 0, ]
cat("Cleaned sample size:", nrow(df_clean), "\n")

y <- df_clean$DVTEXP23
z <- ifelse(df_clean$DNTINS23_M23 == 1, 1, 0)   # 1 = insured, 0 = uninsured
X <- as.matrix(df_clean[, covars])

cat("Estimating Propensity Scores...\n")
ps_model <- glm(z ~ ., data = as.data.frame(X), family = binomial())
pihat <- predict(ps_model, type = "response")

cat("Fitting ZIC-BCF-Smear...\n")
fit <- zicbcf_smear(
  y = y, z = z,
  x_control = X, x_moderate = X,
  pihat = pihat,
  nburn = 500, nsim = 1000
)

# ---- hurdle-margin (participation) CATE draws --------------------------
p0_draws <- pnorm(fit$mu_b)                 # nsim x n  Pr(y>0 | z=0)
p1_draws <- pnorm(fit$mu_b + fit$tau_b)     # nsim x n  Pr(y>0 | z=1)
hurdle_cate_draws <- p1_draws - p0_draws    # nsim x n  participation CATE
hurdle_ate_draws  <- rowMeans(hurdle_cate_draws)

saveRDS(list(hurdle_cate_draws = hurdle_cate_draws,
             hurdle_ate_draws  = hurdle_ate_draws,
             hurdle_cate_mean  = colMeans(hurdle_cate_draws)),
        "applied_study/zicbcf_hurdle_draws.rds")

hurdle_ate_mean <- mean(hurdle_ate_draws)
hurdle_ate_ci   <- quantile(hurdle_ate_draws, c(0.025, 0.975))
cat("\nZIC-BCF hurdle-margin ATE:", round(hurdle_ate_mean, 4),
    " 95% CI [", round(hurdle_ate_ci[1], 4), ",",
    round(hurdle_ate_ci[2], 4), "]\n")

# ---- helper: subgroup posterior summary --------------------------------
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

# ---- subgroup labels (same scheme as the other scripts) ----------------
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

hurdle_table          <- build_table(hurdle_cate_draws)
hurdle_contrast_table <- build_contrasts(hurdle_cate_draws)

write.csv(hurdle_table, "applied_study/zicbcf_hurdle_subgroups.csv", row.names = FALSE)
write.csv(hurdle_contrast_table, "applied_study/zicbcf_hurdle_contrasts.csv", row.names = FALSE)

cat("\n=== ZIC-BCF hurdle-margin subgroup CATEs (probability points) ===\n")
print(hurdle_table)
cat("\n=== ZIC-BCF hurdle-margin between-group contrasts (probability points) ===\n")
print(hurdle_contrast_table)

# ---- histogram of unit-level hurdle CATEs ------------------------------
png("applied_study/zicbcf_hurdle_cate_histogram.png", width = 800, height = 600)
hist(colMeans(hurdle_cate_draws), breaks = 50, col = "lightgreen",
     main = "ZIC-BCF - Distribution of Hurdle CATE (Insurance on Any Dental Expenditure)",
     xlab = "Hurdle CATE (Probability Difference)")
dev.off()

# ---- forest plot of subgroup hurdle CATEs ------------------------------
plot_df <- hurdle_table[hurdle_table$Covariate != "Overall", ]
plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
plot_df <- plot_df[nrow(plot_df):1, ]
overall_h <- hurdle_table$Est[hurdle_table$Covariate == "Overall"]

png("applied_study/zicbcf_hurdle_subgroups.png", width = 1000, height = 1100)
par(mar = c(5, 16, 4, 2))
yy <- seq_len(nrow(plot_df))
plot(plot_df$Est, yy, xlim = range(c(plot_df$CI_low, plot_df$CI_high)),
     yaxt = "n", pch = 19, col = "seagreen", cex = 1.3,
     xlab = "Hurdle CATE: effect of dental insurance on Pr(any dental expenditure)",
     ylab = "", main = "ZIC-BCF - Subgroup hurdle-margin CATEs with 95% CrI")
axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.9)
segments(plot_df$CI_low, yy, plot_df$CI_high, yy, col = "seagreen", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
abline(v = overall_h, lty = 3, col = "firebrick", lwd = 2)
legend("bottomright",
       legend = c("Subgroup hurdle CATE (95% CrI)", "No effect",
                  sprintf("Overall hurdle ATE = %.3f", overall_h)),
       col = c("seagreen", "grey50", "firebrick"),
       pch = c(19, NA, NA), lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
dev.off()

cat("\nCompleted successfully.\n")
