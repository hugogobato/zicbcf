################################################################################
##  Nonlinear-DGP comparison: ZIC-BCF-Smear vs. Gamma Hurdle
##
##  Tests the conjecture that the Gamma-Hurdle benchmark's strong performance in
##  the main study is an artefact of *link-linear* DGPs. Here the three DGPs have
##  nonlinear conditional-mean structure (see nonlinear_dgps.R) that the
##  Gamma-Hurdle (1, z, X, z*X) design cannot represent, while BCF can.
##
##  Standard configuration only: N = 500, default c_shift, MCMC = 1000 burn-in +
##  1000 saved draws (identical budget to the linear simulation studies).
##
##  Usage:
##    Rscript run_nonlinear_comparison.R            # 3 seeds per DGP (default)
##    Rscript run_nonlinear_comparison.R 10         # 10 seeds per DGP
################################################################################
suppressMessages(library(zicbcf))

HERE <- tryCatch(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE))), error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
source(file.path(HERE, "nonlinear_dgps.R"))

RESULTS_DIR <- file.path(HERE, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Settings (standard configuration, matching the linear studies) ----------
args   <- commandArgs(trailingOnly = TRUE)
N_SIM  <- if (length(args) >= 1) as.integer(args[1]) else 3L
N      <- 500L
P      <- 5L
NBURN  <- 1000L
NSIM   <- 1000L
NTHIN  <- 1L

cat(sprintf("=== Nonlinear DGP comparison | %d seeds/DGP | N=%d | nburn=%d nsim=%d ===\n",
            N_SIM, N, NBURN, NSIM))

## ---- Per-(model, seed) fitter, returns a one-row data.frame ------------------
fit_one <- function(model, d, dgp_name, seed) {
  fit <- tryCatch({
    if (model == "ZIC-BCF-Smear") {
      f <- zicbcf_smear(y = d$y, z = d$z, x_control = d$x, pihat = d$pihat,
                        nburn = NBURN, nsim = NSIM, update_interval = 99999)
      # response-scale CATE draws + probit-hurdle contrast draws
      list(cate = f$cate, ate = f$ate,
           hurdle = pnorm(f$mu_b + f$tau_b) - pnorm(f$mu_b))
    } else { # Gamma Hurdle
      f <- gamma_hurdle(y = d$y, z = d$z, x = d$x, nburn = NBURN, nsim = NSIM, nthin = NTHIN)
      list(cate = f$cate, ate = f$ate, hurdle = f$p1 - f$p0)
    }
  }, error = function(e) { message(sprintf("  %s seed %d FAILED: %s", model, seed, conditionMessage(e))); NULL })

  if (is.null(fit)) {
    return(data.frame(DGP = dgp_name, Seed = seed, Model = model, Status = "FAILED",
                      True_ATE = d$true_ate, True_Hurdle_ATE = d$true_hurdle_ate,
                      CATE_RMSE = NA, CATE_Abs_Bias = NA, CATE_Coverage = NA,
                      CATE_Correlation = NA, CATE_CI_Length = NA,
                      Est_ATE = NA, ATE_AbsErr = NA,
                      Hurdle_RMSE = NA, Hurdle_Coverage = NA, Hurdle_Correlation = NA,
                      Hurdle_CI_Length = NA, Est_Hurdle_ATE = NA, stringsAsFactors = FALSE))
  }

  m  <- calc_cate_metrics(fit$cate, d$true_cate, fit$ate)
  mh <- calc_cate_metrics(fit$hurdle, d$true_hurdle_cate, rowMeans(fit$hurdle))
  data.frame(DGP = dgp_name, Seed = seed, Model = model, Status = "OK",
             True_ATE = d$true_ate, True_Hurdle_ATE = d$true_hurdle_ate,
             CATE_RMSE = m$rmse, CATE_Abs_Bias = abs(m$bias), CATE_Coverage = m$coverage,
             CATE_Correlation = m$correlation, CATE_CI_Length = m$ci_length,
             Est_ATE = m$est_ate_mean, ATE_AbsErr = abs(m$est_ate_mean - d$true_ate),
             Hurdle_RMSE = mh$rmse, Hurdle_Coverage = mh$coverage,
             Hurdle_Correlation = mh$correlation, Hurdle_CI_Length = mh$ci_length,
             Est_Hurdle_ATE = mh$est_ate_mean, stringsAsFactors = FALSE)
}

## ---- Main loop ---------------------------------------------------------------
all_rows <- list()
for (k in names(nl_dgps)) {
  g <- nl_dgps[[k]]
  cat(sprintf("\n--- %s ---\n", g$name))
  for (s in 1:N_SIM) {
    cat(sprintf("  [seed %d/%d] generate + fit both models...\n", s, N_SIM))
    d <- g$func(N, P, seed = s, c_shift = g$c_shift)
    for (model in c("ZIC-BCF-Smear", "Gamma Hurdle")) {
      all_rows[[length(all_rows) + 1]] <- fit_one(model, d, g$name, s)
    }
  }
}
res <- do.call(rbind, all_rows)
write.csv(res, file.path(RESULTS_DIR, "nonlinear_results_long.csv"), row.names = FALSE)

## ---- Aggregate across seeds --------------------------------------------------
ok <- res[res$Status == "OK", ]
agg <- do.call(rbind, lapply(split(ok, list(ok$DGP, ok$Model), drop = TRUE), function(df) {
  data.frame(
    DGP = df$DGP[1], Model = df$Model[1], Seeds = nrow(df),
    CATE_RMSE        = mean(df$CATE_RMSE),
    CATE_Correlation = mean(df$CATE_Correlation),
    CATE_Coverage    = mean(df$CATE_Coverage),
    CATE_CI_Length   = mean(df$CATE_CI_Length),
    ATE_RMSE         = sqrt(mean(df$ATE_AbsErr^2)),
    ATE_AbsErr_mean  = mean(df$ATE_AbsErr),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$DGP, agg$Model), ]
write.csv(agg, file.path(RESULTS_DIR, "nonlinear_summary.csv"), row.names = FALSE)

cat("\n\n================  AGGREGATED SUMMARY (mean over seeds)  ================\n")
print(within(agg, {
  CATE_RMSE        <- round(CATE_RMSE, 3)
  CATE_Correlation <- round(CATE_Correlation, 3)
  CATE_Coverage    <- round(CATE_Coverage, 3)
  CATE_CI_Length   <- round(CATE_CI_Length, 3)
  ATE_RMSE         <- round(ATE_RMSE, 3)
  ATE_AbsErr_mean  <- round(ATE_AbsErr_mean, 3)
}), row.names = FALSE)

cat(sprintf("\n[SAVED] %s\n[SAVED] %s\n",
            file.path(RESULTS_DIR, "nonlinear_results_long.csv"),
            file.path(RESULTS_DIR, "nonlinear_summary.csv")))
cat("\n=== Done ===\n")
