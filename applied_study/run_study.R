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

# 1. Propensity Score Estimation
cat("Estimating Propensity Scores...\n")
ps_model <- glm(z ~ ., data = as.data.frame(X), family = binomial())
pihat <- predict(ps_model, type = "response")

# 2. Fit ZIC-BCF-Smear
cat("Fitting ZIC-BCF-Smear...\n")
fit <- zicbcf_smear(
  y = y,
  z = z,
  x_control = X,
  x_moderate = X,
  pihat = pihat,
  nburn = 500,
  nsim = 1000
)

# 3. Analyze Results
cate_draws <- fit$cate
ate_draws <- fit$ate

cate_mean <- colMeans(cate_draws)
ate_mean <- mean(ate_draws)
ate_ci <- quantile(ate_draws, c(0.025, 0.975))

cat("ATE:", ate_mean, "\n")
cat("ATE 95% CI:", ate_ci[1], "-", ate_ci[2], "\n")

# Identify heterogeneity (e.g. by Age and Income)
df_clean$CATE <- cate_mean

# Summarize heterogeneity
cat("\nSummary of CATE:\n")
print(summary(df_clean$CATE))

# Save a simple plot of CATE distribution
png("applied_study/cate_histogram.png", width=800, height=600)
hist(df_clean$CATE, breaks=50, col="skyblue", main="Distribution of Conditional Average Treatment Effects (CATE)", xlab="CATE (Impact of Insurance on Dental Expenditure, $)")
dev.off()

# Save ATE and summary to a text file
res <- data.frame(
  Metric = c("ATE", "Lower 95% CI", "Upper 95% CI", "Mean CATE", "Median CATE"),
  Value = c(ate_mean, ate_ci[1], ate_ci[2], mean(cate_mean), median(cate_mean))
)
write.csv(res, "applied_study/ate_results.csv", row.names=FALSE)
cat("Completed successfully.\n")
