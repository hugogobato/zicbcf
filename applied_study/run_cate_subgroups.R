library(zicbcf)

# ==========================================================================
# Proper CATE (subgroup) analysis for the MEPS dental-insurance study.
#
# The ZIC-BCF-Smear fit returns a full posterior of unit-level CATEs
# (cate_draws: nsim x n). A subgroup CATE is obtained by averaging the
# unit-level effects over the units in the subgroup *within each posterior
# draw*, which yields a posterior distribution for the subgroup effect and
# therefore a point estimate AND a 95% credible interval. Between-group
# contrasts (e.g. Male - Female) are computed the same way, giving a direct
# posterior test for treatment-effect heterogeneity.
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

cate_draws <- fit$cate            # nsim x n posterior of unit-level CATE
ate_draws  <- fit$ate
saveRDS(cate_draws, "applied_study/cate_draws.rds")

# ---- helper: subgroup posterior summary --------------------------------
# idx: logical/integer vector selecting the units (columns) in the subgroup.
subgroup_cate <- function(idx) {
  sub <- cate_draws[, idx, drop = FALSE]
  draw_means <- rowMeans(sub)              # posterior draws of the subgroup CATE
  c(N        = ncol(sub),
    CATE     = mean(draw_means),
    CI_low   = unname(quantile(draw_means, 0.025)),
    CI_high  = unname(quantile(draw_means, 0.975)),
    P_gt_0   = mean(draw_means > 0))
}

# posterior of the difference between two subgroups (heterogeneity test)
contrast_cate <- function(idx_a, idx_b) {
  da <- rowMeans(cate_draws[, idx_a, drop = FALSE])
  db <- rowMeans(cate_draws[, idx_b, drop = FALSE])
  d  <- da - db
  c(Diff    = mean(d),
    CI_low  = unname(quantile(d, 0.025)),
    CI_high = unname(quantile(d, 0.975)),
    P_gt_0  = mean(d > 0))
}

# ---- build interpretable subgroup labels -------------------------------
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

# ---- overall + per-subgroup CATE table ---------------------------------
rows <- list()
rows[[1]] <- data.frame(Covariate = "Overall", Level = "All",
                        t(subgroup_cate(rep(TRUE, ncol(cate_draws)))),
                        check.names = FALSE)

for (gname in names(group_vars)) {
  lab <- group_vars[[gname]]
  for (lv in levels(lab)) {
    idx <- which(lab == lv)
    if (length(idx) < 30) next          # skip tiny cells
    rows[[length(rows) + 1]] <- data.frame(
      Covariate = gname, Level = lv,
      t(subgroup_cate(idx)), check.names = FALSE)
  }
}
cate_table <- do.call(rbind, rows)
rownames(cate_table) <- NULL
cate_table[, c("CATE", "CI_low", "CI_high")] <-
  round(cate_table[, c("CATE", "CI_low", "CI_high")], 2)
cate_table$P_gt_0 <- round(cate_table$P_gt_0, 3)

write.csv(cate_table, "applied_study/cate_subgroups.csv", row.names = FALSE)
cat("\n=== Subgroup CATE estimates ===\n")
print(cate_table)

# ---- pairwise heterogeneity contrasts within each covariate -------------
crows <- list()
for (gname in names(group_vars)) {
  lab <- group_vars[[gname]]
  lvs <- levels(lab)[table(lab) >= 30]
  if (length(lvs) < 2) next
  cmb <- combn(lvs, 2)
  for (k in seq_len(ncol(cmb))) {
    a <- cmb[1, k]; b <- cmb[2, k]
    ct <- contrast_cate(which(lab == a), which(lab == b))
    crows[[length(crows) + 1]] <- data.frame(
      Covariate = gname, Contrast = paste(a, "-", b),
      t(ct), check.names = FALSE)
  }
}
contrast_table <- do.call(rbind, crows)
rownames(contrast_table) <- NULL
contrast_table[, c("Diff", "CI_low", "CI_high")] <-
  round(contrast_table[, c("Diff", "CI_low", "CI_high")], 2)
contrast_table$P_gt_0 <- round(contrast_table$P_gt_0, 3)

write.csv(contrast_table, "applied_study/cate_contrasts.csv", row.names = FALSE)
cat("\n=== Between-group CATE contrasts (heterogeneity) ===\n")
print(contrast_table)

# ---- plot: subgroup CATEs with 95% credible intervals -------------------
plot_df <- cate_table[cate_table$Covariate != "Overall", ]
plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
plot_df <- plot_df[nrow(plot_df):1, ]        # top-to-bottom order
overall_cate <- cate_table$CATE[cate_table$Covariate == "Overall"]

png("applied_study/11_MEPS_cate_subgroups.png", width = 1000, height = 1100)
par(mar = c(5, 16, 4, 2))
yy <- seq_len(nrow(plot_df))
plot(plot_df$CATE, yy, xlim = range(c(plot_df$CI_low, plot_df$CI_high)),
     yaxt = "n", pch = 19, col = "steelblue", cex = 1.3,
     xlab = "CATE: effect of dental insurance on annual dental expenditure ($)",
     ylab = "", main = "Subgroup CATE estimates with 95% credible intervals")
axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.9)
segments(plot_df$CI_low, yy, plot_df$CI_high, yy, col = "steelblue", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
abline(v = overall_cate, lty = 3, col = "firebrick", lwd = 2)
legend("bottomright", legend = c("Subgroup CATE (95% CrI)", "No effect",
       sprintf("Overall ATE = $%.0f", overall_cate)),
       col = c("steelblue", "grey50", "firebrick"),
       pch = c(19, NA, NA), lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
dev.off()

cat("\nATE:", round(mean(ate_draws), 2),
    " 95% CI [", round(quantile(ate_draws, .025), 2), ",",
    round(quantile(ate_draws, .975), 2), "]\n")
cat("Completed successfully.\n")
