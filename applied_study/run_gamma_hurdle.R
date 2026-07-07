library(zicbcf)

cat("Loading dataset...\n")
df <- read.csv("applied_study/h251.csv")

# Select columns and filter out missing/inapplicable (-1, -7, -8, -9)
covars <- c("AGE23X", "SEX", "RACEV2X", "FAMINC23", "POVCAT23", "REGION23", "MARRY23X")
df_clean <- df[df$DNTINS23_M23 > 0, ]
for (col in covars) {
  df_clean <- df_clean[df_clean[[col]] >= 0, ]
}
df_clean <- df_clean[df_clean$DVTEXP23 >= 0, ]

cat("Cleaned sample size:", nrow(df_clean), "\n")

# Prepare variables
y <- df_clean$DVTEXP23
z <- ifelse(df_clean$DNTINS23_M23 == 1, 1, 0)
X <- as.matrix(df_clean[, covars])

# Propensity Score Estimation
cat("Estimating Propensity Scores...\n")
ps_model <- glm(z ~ ., data = as.data.frame(X), family = binomial())
pihat <- predict(ps_model, type = "response")

# Fit Gamma Hurdle (Oganisian et al. 2019)
cat("Fitting Gamma Hurdle...\n")
fit <- gamma_hurdle(
  y = y,
  z = z,
  x = X,
  nburn = 1000,
  nsim = 1000,
  nthin = 1
)

# Analyze results
cate_draws <- fit$cate
ate_draws <- fit$ate
hurdle_cate_draws <- fit$p1 - fit$p0
hurdle_ate_draws <- rowMeans(hurdle_cate_draws)

cate_mean <- colMeans(cate_draws)
ate_mean <- mean(ate_draws)
ate_ci <- quantile(ate_draws, c(0.025, 0.975))

hurdle_ate_mean <- mean(hurdle_ate_draws)
hurdle_ate_ci <- quantile(hurdle_ate_draws, c(0.025, 0.975))

cat("ATE:", ate_mean, "\n")
cat("ATE 95% CI:", ate_ci[1], "-", ate_ci[2], "\n")
cat("Hurdle ATE:", hurdle_ate_mean, "\n")
cat("Hurdle ATE 95% CI:", hurdle_ate_ci[1], "-", hurdle_ate_ci[2], "\n")

# Save CATE histogram
png("applied_study/gamma_hurdle_cate_histogram.png", width=800, height=600)
hist(cate_mean, breaks=50, col="skyblue",
     main="Gamma Hurdle - Distribution of CATE (Insurance on Dental Expenditure, $)",
     xlab="CATE ($)")
dev.off()

# Save hurdle CATE histogram
png("applied_study/gamma_hurdle_hurdle_cate_histogram.png", width=800, height=600)
hist(colMeans(hurdle_cate_draws), breaks=50, col="lightgreen",
     main="Gamma Hurdle - Distribution of Hurdle CATE (Insurance on Any Dental Expenditure)",
     xlab="Hurdle CATE (Probability Difference)")
dev.off()

# Save results
res <- data.frame(
  Metric = c("ATE", "ATE Lower 95% CI", "ATE Upper 95% CI",
             "Mean CATE", "Median CATE",
             "Hurdle ATE", "Hurdle ATE Lower 95% CI", "Hurdle ATE Upper 95% CI",
             "Mean Hurdle CATE", "Median Hurdle CATE"),
  Value = c(ate_mean, ate_ci[1], ate_ci[2],
            mean(cate_mean), median(cate_mean),
            hurdle_ate_mean, hurdle_ate_ci[1], hurdle_ate_ci[2],
            mean(colMeans(hurdle_cate_draws)), median(colMeans(hurdle_cate_draws)))
)
write.csv(res, "applied_study/gamma_hurdle_results.csv", row.names=FALSE)

# Save CATE draws for later analysis
saveRDS(list(cate_draws = cate_draws, hurdle_cate_draws = hurdle_cate_draws,
             ate_draws = ate_draws, hurdle_ate_draws = hurdle_ate_draws,
             cate_mean = cate_mean, hurdle_cate_mean = colMeans(hurdle_cate_draws)),
        "applied_study/gamma_hurdle_draws.rds")

cat("Completed successfully.\n")
