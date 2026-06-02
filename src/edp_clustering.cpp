#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// ---------------------------------------------------------------------------
// EDP+ZI cluster-membership sampler (Algorithm 8, Neal 2000) adapted from
// Kim et al. (2024) "Bayesian nonparametric model for heterogeneous treatment
// effects with zero-inflated data" (EDPcausal, sim2/clustering.cpp).
//
// ADAPTATION FOR zicbcf BENCHMARKS:
//   The original code hard-codes the covariate layout used in the paper's
//   simulation: design matrix xa = (1, A, X1[binary], X2, X3, X4, X5) where X1
//   is modelled with a Dirichlet/Bernoulli kernel and X2-X5 with Gaussian
//   kernels. The zicbcf simulation DGPs (A/B/C) generate FIVE continuous
//   covariates plus a binary treatment, so this version models the treatment
//   A with a 2-category Dirichlet kernel (xPAR_p0) and ALL FIVE covariates
//   X1-X5 with Gaussian kernels (xPAR_mu / xPAR_sig, now 5 rows). The binary
//   X1 kernel (xPAR_p1) of the original is removed.
//
//   The MODEL (EDP mixture, kernels, priors, the zero-inflation link
//   pi = expit(r), and the Neal-8 cluster weights) is exactly as in the paper.
//   Three implementation defects in the authors' clustering.cpp are corrected
//   WITHOUT altering the model:
//     (i)  the cluster-assignment zero term was written exp(r)/exp(1+r), which
//          collapses to the constant 1/e; it is restored to the model's
//          pi = expit(r), evaluated stably via -log1p(exp(-/+r));
//     (ii) the membership weights are accumulated in log-space and exp-
//          normalised, so an outlier subject whose kernel product underflows
//          no longer yields an all-zero P (the authors' code aborts there);
//     (iii)the truncated-normal outcome kernel is evaluated in closed form,
//          log N(y|b,sig) - log Phi(b/sig), removing a per-iteration R callback.
// ---------------------------------------------------------------------------

// Dirichlet sampler (from EDPcausal MCMC_utils.h)
inline arma::mat edp_rdirichlet_cpp(int num_samples, arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  arma::mat distribution = arma::zeros(num_samples, distribution_size);
  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    for (int j = 0; j < distribution_size; ++j) {
      double cur = R::rgamma(alpha_m[j], 1.0);
      distribution(i, j) = cur;
      sum_term += cur;
    }
    for (int j = 0; j < distribution_size; ++j) {
      distribution(i, j) = distribution(i, j) / sum_term;
    }
  }
  return (distribution);
}

// [[Rcpp::export]]
List edp_clustering_c(
    IntegerVector Sy,
    IntegerVector Sx,
    NumericMatrix xPAR_p0,
    NumericMatrix xPAR_mu,
    NumericMatrix xPAR_sig,
    NumericMatrix yPAR_beta,
    NumericMatrix yPAR_r,
    NumericVector yPAR_sig,
    double alpha_omega,
    double alpha_theta,
    NumericMatrix xa,
    NumericVector y
) {

  // Data preparation
  const int n = Sy.length();
  arma::mat ypar_beta = as<arma::mat>((yPAR_beta));
  arma::mat ypar_r = as<arma::mat>((yPAR_r));
  arma::vec ypar_sig = as<arma::vec>((yPAR_sig));

  arma::mat xpar_p0 = as<arma::mat>((xPAR_p0));
  arma::mat xpar_mu = as<arma::mat>((xPAR_mu));
  arma::mat xpar_sig = as<arma::mat>((xPAR_sig));
  int numY = yPAR_beta.ncol();
  int numX = xPAR_p0.nrow();
  int M = 1;                   // num of auxilary clusters = 1
  double mu_r = 0;
  double mu_beta = 0;
  double Sigma_r = 10;
  double Sigma_beta = 10;
  NumericMatrix ycoef_r(7, 1);    // num of col. in design mat = 7
  NumericMatrix ycoef_beta(7, 1); // num of col. in design mat = 7
  double ycoef_sig;
  NumericVector alpha00 = {0.1, 0.1}; // Dirichlet prior for binary treatment A
  NumericVector x_cat0(alpha00.length());
  double x_var1, x_var2, x_var3, x_var4, x_var5;
  double x_mu1, x_mu2, x_mu3, x_mu4, x_mu5;
  double nu0 = 2;
  double tau0 = 1;
  double c0 = 1;

  // R's order() is used to canonically sort the (Sy, Sx) cluster labels.
  Function order("order");

  // loop starts from here
  for (int h = 0; h < n; h++) {

  arma::mat ypar_beta_hold = (ypar_beta);
  arma::mat ypar_r_hold = (ypar_r);
  arma::vec ypar_sig_hold = (ypar_sig);

  arma::mat xpar_p0_hold = (xpar_p0);
  arma::mat xpar_mu_hold = (xpar_mu);
  arma::mat xpar_sig_hold = (xpar_sig);

  NumericMatrix uniqueYX(n, 2);
  uniqueYX(_, 0) = clone(Sy);
  uniqueYX(_, 1) = clone(Sx);

  arma::uvec ulmt = arma::zeros<arma::uvec>(uniqueYX.nrow());
  arma::mat uniqueyx = arma::mat(uniqueYX.begin(), uniqueYX.nrow(), uniqueYX.ncol(), false);

  for (arma::uword i = 0; i < uniqueYX.nrow(); i++) {
    for (arma::uword j = i + 1; j < uniqueYX.nrow(); j++) {
      if (arma::approx_equal(uniqueyx.row(i), uniqueyx.row(j), "absdiff", 0.00000001)) { ulmt(j) = 1; break; }
    }
  }
  arma::mat uniqueyxi_temp = uniqueyx.rows(find(ulmt == 0));
  NumericVector uniqueyxi_ind = order(uniqueyxi_temp.col(0), uniqueyxi_temp.col(1));
  arma::mat uniqueyxi = arma::zeros<arma::mat>(uniqueyxi_temp.n_rows, 2);
  for (int j = 0; j < uniqueyxi_ind.length(); j++) {
    uniqueyxi.row(j) = uniqueyxi_temp.row(uniqueyxi_ind(j) - 1);
  }

  IntegerVector Syid = Sy[Sy == Sy(h)];
  double lenY = Syid.length() - 1;

  IntegerVector Sxid = Sx[Sx == Sx(h) & Sy == Sy(h)];
  double lenX = Sxid.length() - 1;

  arma::vec indY_int = as<arma::vec>(wrap(match(clone(unique(Sy)).sort(), unique(Sy))));
  arma::uvec indY = arma::find(indY_int == Sy(h)); // position of Sy(h) in the sorted unique Sy

  arma::uvec indYX = arma::find(uniqueyxi.col(0) == Sy(h));

  arma::uvec indX = arma::find((uniqueyxi.col(0) == Sy(h)) && (uniqueyxi.col(1) == Sx(h)));

  arma::vec sy = as<arma::vec>((Sy));
  arma::vec sx = as<arma::vec>((Sx));

  if (lenX == 0)
  {
    if (lenY == 0)
    {
      arma::vec idy = as<arma::vec>(wrap(seq_len(numY)));
      ypar_beta = ypar_beta.cols(arma::find(idy != Sy(h)));
      ypar_r = ypar_r.cols(arma::find(idy != Sy(h)));
      ypar_sig = ypar_sig(arma::find(idy != Sy(h)));

      arma::vec idx = as<arma::vec>(wrap(seq_len(numX) - 1));
      xpar_p0 = xpar_p0.rows(arma::find(idx != as<int>(wrap(indX))));
      xpar_mu = xpar_mu.cols(arma::find(idx != as<int>(wrap(indX))));
      xpar_sig = xpar_sig.cols(arma::find(idx != as<int>(wrap(indX))));
    }
    else
    {
      arma::vec idx = as<arma::vec>(wrap(seq_len(numX) - 1));
      xpar_p0 = xpar_p0.rows(arma::find(idx != as<int>(wrap(indX))));
      xpar_mu = xpar_mu.cols(arma::find(idx != as<int>(wrap(indX))));
      xpar_sig = xpar_sig.cols(arma::find(idx != as<int>(wrap(indX))));

      arma::uvec indyy = arma::find(sy == Sy(h));

      arma::uvec indyx = indyy.elem(arma::find(indyy != h));

      arma::vec sxsub = sx.elem(indyx);

      arma::vec sxsub_unique = sort(unique(sxsub));

      arma::vec sxsub_len = as<arma::vec>(wrap(seq_len(sxsub_unique.n_elem)));

      for (int z = 0; z < indyx.n_elem; z++) {
        arma::uword idx = indyx[z];
        arma::uvec sx_idx = arma::find(sxsub_unique == sx(idx));
        sx(idx) = sxsub_len(sx_idx(0));
      }
    }

    arma::uvec sylar = arma::find(sy > sy[h]);
    if ((lenY == 0) && (sylar.n_elem > 0))
    {
      sy.elem(sylar) = sy.elem(sylar) - 1;
    }
  }

  IntegerVector newind = seq(0, n - 1);
  newind.erase(h);

  arma::vec newsx = sx.elem(as<arma::uvec>((newind)));
  arma::vec newsy = sy.elem(as<arma::uvec>((newind)));

  arma::mat newuniqueyx(newsx.n_elem, 2);
  newuniqueyx.col(0) = newsy;
  newuniqueyx.col(1) = newsx;
  arma::uvec newulmt = arma::zeros<arma::uvec>(newuniqueyx.n_rows);

  for (arma::uword i = 0; i < newuniqueyx.n_rows; i++) {
    for (arma::uword j = i + 1; j < newuniqueyx.n_rows; j++) {
      if (arma::approx_equal(newuniqueyx.row(i), newuniqueyx.row(j), "absdiff", 0.00000001)) { newulmt(j) = 1; break; }
    }
  }
  arma::mat newuniqueyxi = newuniqueyx.rows(find(newulmt == 0));

  // recalculate number of unique clusters
  numY = ypar_beta.n_cols;
  numX = xpar_p0.n_rows;

  // count # in each X-Y cluster excluding the h-th person
  NumericVector nji(numY);
  IntegerVector numXj(numY); // # of X clusters within each Y cluster
  IntegerMatrix nlji(numY, max(newsx));

  for (int i = 0; i < numY; i++) {
    nji(i) = (arma::find(newsy == (i + 1))).eval().n_elem;
    for (int j = 0; j < max(newsx); j++) {
      nlji(i, j) = (arma::find((newsy == (i + 1)) && (newsx == (j + 1)))).eval().n_elem;
    }
    numXj(i) = (arma::find(as<arma::mat>(wrap(nlji)).row(i) != 0)).eval().n_elem;
  }

  if (lenY * lenX != 0)
  {
    arma::mat x_mu_temp = arma::zeros<arma::mat>(xpar_mu.n_rows, xpar_mu.n_cols + numY * M + M);
    arma::mat x_sig_temp = arma::zeros<arma::mat>(xpar_sig.n_rows, xpar_sig.n_cols + numY * M + M);
    arma::mat x_p0_temp = arma::zeros<arma::mat>(xpar_p0.n_rows + numY * M + M, xpar_p0.n_cols);

    for (int i = 0; i < M; i++) {
      ycoef_r(_, 0) = Rcpp::rnorm(7, mu_r, sqrt(Sigma_r));
      ycoef_beta(_, 0) = Rcpp::rnorm(7, mu_beta, sqrt(Sigma_beta));
      ycoef_sig = Rcpp::rgamma(1, 0.1, 0.1)(0);

      x_cat0 = edp_rdirichlet_cpp(1, alpha00);
      x_var1 = (nu0 * tau0) / R::rchisq(nu0);
      x_var2 = (nu0 * tau0) / R::rchisq(nu0);
      x_var3 = (nu0 * tau0) / R::rchisq(nu0);
      x_var4 = (nu0 * tau0) / R::rchisq(nu0);
      x_var5 = (nu0 * tau0) / R::rchisq(nu0);
      x_mu1 = R::rnorm(0, sqrt(x_var1 / c0));
      x_mu2 = R::rnorm(0, sqrt(x_var2 / c0));
      x_mu3 = R::rnorm(0, sqrt(x_var3 / c0));
      x_mu4 = R::rnorm(0, sqrt(x_var4 / c0));
      x_mu5 = R::rnorm(0, sqrt(x_var5 / c0));

      x_mu_temp.col(xpar_mu.n_cols + numY * M + i) = arma::vec({x_mu1, x_mu2, x_mu3, x_mu4, x_mu5});
      x_sig_temp.col(xpar_sig.n_cols + numY * M + i) = arma::vec({x_var1, x_var2, x_var3, x_var4, x_var5});
      x_p0_temp.row(xpar_p0.n_rows + numY * M + i) = as<arma::mat>(x_cat0);

      arma::mat ypar_beta_temp(7, 1);
      ypar_beta_temp.col(0) = as<arma::mat>(ycoef_beta);
      ypar_beta = arma::join_rows(ypar_beta, ypar_beta_temp);

      arma::mat ypar_r_temp(7, 1);
      ypar_r_temp.col(0) = as<arma::mat>(ycoef_r);
      ypar_r = arma::join_rows(ypar_r, ypar_r_temp);

      int sz = ypar_sig.size();
      ypar_sig.resize(sz + 1);
      ypar_sig(sz) = ycoef_sig;

      int cumsum = 0;
      for (int j = 0; j < numY; j++) {

        cumsum += numXj(j) + M;
        x_cat0 = edp_rdirichlet_cpp(1, alpha00);
        x_var1 = (nu0 * tau0) / R::rchisq(nu0);
        x_var2 = (nu0 * tau0) / R::rchisq(nu0);
        x_var3 = (nu0 * tau0) / R::rchisq(nu0);
        x_var4 = (nu0 * tau0) / R::rchisq(nu0);
        x_var5 = (nu0 * tau0) / R::rchisq(nu0);
        x_mu1 = R::rnorm(0, sqrt(x_var1 / c0));
        x_mu2 = R::rnorm(0, sqrt(x_var2 / c0));
        x_mu3 = R::rnorm(0, sqrt(x_var3 / c0));
        x_mu4 = R::rnorm(0, sqrt(x_var4 / c0));
        x_mu5 = R::rnorm(0, sqrt(x_var5 / c0));

        x_mu_temp.col(cumsum - 1 - i) = arma::vec({x_mu1, x_mu2, x_mu3, x_mu4, x_mu5});
        x_sig_temp.col(cumsum - 1 - i) = arma::vec({x_var1, x_var2, x_var3, x_var4, x_var5});
        x_p0_temp.row(cumsum - 1 - i) = as<arma::mat>(x_cat0);
      }
    }

    int Cumsum = 0;
    int Cumsum1 = 0;
    for (int j = 0; j < numY; j++) {
      Cumsum1 += numXj(j) + M;
      for (int t = 0; t < numXj(j); t++) {
        x_mu_temp.col(Cumsum1 - numXj(j) - M + t) = xpar_mu.col(Cumsum);
        x_sig_temp.col(Cumsum1 - numXj(j) - M + t) = xpar_sig.col(Cumsum);
        x_p0_temp.row(Cumsum1 - numXj(j) - M + t) = xpar_p0.row(Cumsum);
        Cumsum++;
      }
    }

    xpar_mu = (x_mu_temp);
    xpar_sig = (x_sig_temp);
    xpar_p0 = (x_p0_temp);
  }

  // The current value of S[h] is in an existing y-cluster, but not in an existing x-subcluster
  if (lenY != 0 & lenX == 0) {
    arma::mat x_mu_temp = arma::zeros<arma::mat>(xpar_mu.n_rows, xpar_mu.n_cols + numY * M + M);
    arma::mat x_sig_temp = arma::zeros<arma::mat>(xpar_sig.n_rows, xpar_sig.n_cols + numY * M + M);
    arma::mat x_p0_temp = arma::zeros<arma::mat>(xpar_p0.n_rows + numY * M + M, xpar_p0.n_cols);
    for (int i = 0; i < M; i++) {
      ycoef_r(_, 0) = Rcpp::rnorm(7, mu_r, sqrt(Sigma_r));
      ycoef_beta(_, 0) = Rcpp::rnorm(7, mu_beta, sqrt(Sigma_beta));
      ycoef_sig = Rcpp::rgamma(1, 0.1, 0.1)(0);

      x_cat0 = edp_rdirichlet_cpp(1, alpha00);
      x_var1 = (nu0 * tau0) / R::rchisq(nu0);
      x_var2 = (nu0 * tau0) / R::rchisq(nu0);
      x_var3 = (nu0 * tau0) / R::rchisq(nu0);
      x_var4 = (nu0 * tau0) / R::rchisq(nu0);
      x_var5 = (nu0 * tau0) / R::rchisq(nu0);
      x_mu1 = R::rnorm(0, sqrt(x_var1 / c0));
      x_mu2 = R::rnorm(0, sqrt(x_var2 / c0));
      x_mu3 = R::rnorm(0, sqrt(x_var3 / c0));
      x_mu4 = R::rnorm(0, sqrt(x_var4 / c0));
      x_mu5 = R::rnorm(0, sqrt(x_var5 / c0));

      x_mu_temp.col(xpar_mu.n_cols + numY * M + i) = arma::vec({x_mu1, x_mu2, x_mu3, x_mu4, x_mu5});
      x_sig_temp.col(xpar_sig.n_cols + numY * M + i) = arma::vec({x_var1, x_var2, x_var3, x_var4, x_var5});
      x_p0_temp.row(xpar_p0.n_rows + numY * M + i) = as<arma::mat>(x_cat0);

      arma::mat ypar_beta_temp(7, 1);
      ypar_beta_temp.col(0) = as<arma::mat>(ycoef_beta);
      ypar_beta = arma::join_rows(ypar_beta, ypar_beta_temp);

      arma::mat ypar_r_temp(7, 1);
      ypar_r_temp.col(0) = as<arma::mat>(ycoef_r);
      ypar_r = arma::join_rows(ypar_r, ypar_r_temp);

      int sz = ypar_sig.size();
      ypar_sig.resize(sz + 1);
      ypar_sig(sz) = ycoef_sig;

      int cumsum = 0;
      for (int j = 0; j < numY; j++) {
        cumsum += numXj(j) + M;

        x_cat0 = edp_rdirichlet_cpp(1, alpha00);
        x_var1 = (nu0 * tau0) / R::rchisq(nu0);
        x_var2 = (nu0 * tau0) / R::rchisq(nu0);
        x_var3 = (nu0 * tau0) / R::rchisq(nu0);
        x_var4 = (nu0 * tau0) / R::rchisq(nu0);
        x_var5 = (nu0 * tau0) / R::rchisq(nu0);
        x_mu1 = R::rnorm(0, sqrt(x_var1 / c0));
        x_mu2 = R::rnorm(0, sqrt(x_var2 / c0));
        x_mu3 = R::rnorm(0, sqrt(x_var3 / c0));
        x_mu4 = R::rnorm(0, sqrt(x_var4 / c0));
        x_mu5 = R::rnorm(0, sqrt(x_var5 / c0));

        x_mu_temp.col(cumsum - 1 - i) = arma::vec({x_mu1, x_mu2, x_mu3, x_mu4, x_mu5});
        x_sig_temp.col(cumsum - 1 - i) = arma::vec({x_var1, x_var2, x_var3, x_var4, x_var5});
        x_p0_temp.row(cumsum - 1 - i) = as<arma::mat>(x_cat0);
      }
    }

    int Cumsum = 0;
    int Cumsum1 = 0;
    for (int j = 0; j < numY; j++) {
      Cumsum1 += numXj(j) + M;
      if (j == (Sy(h) - 1))
      {
        x_mu_temp.col(Cumsum1 - M) = xpar_mu_hold.col(indX(0));
        x_sig_temp.col(Cumsum1 - M) = xpar_sig_hold.col(indX(0));
        x_p0_temp.row(Cumsum1 - M) = xpar_p0_hold.row(indX(0));
      }
      for (int t = 0; t < numXj(j); t++) {

        x_mu_temp.col(Cumsum1 - numXj(j) - M + t) = xpar_mu.col(Cumsum);
        x_sig_temp.col(Cumsum1 - numXj(j) - M + t) = xpar_sig.col(Cumsum);
        x_p0_temp.row(Cumsum1 - numXj(j) - M + t) = xpar_p0.row(Cumsum);
        Cumsum++;
      }
    }

    xpar_mu = (x_mu_temp);
    xpar_sig = (x_sig_temp);
    xpar_p0 = (x_p0_temp);
  }

  // The current value of S[h] is not in an existing y-cluster
  if (lenY == 0) {
    arma::mat x_mu_temp = arma::zeros<arma::mat>(xpar_mu.n_rows, xpar_mu.n_cols + numY * M + M);
    arma::mat x_sig_temp = arma::zeros<arma::mat>(xpar_sig.n_rows, xpar_sig.n_cols + numY * M + M);
    arma::mat x_p0_temp = arma::zeros<arma::mat>(xpar_p0.n_rows + numY * M + M, xpar_p0.n_cols);
    for (int i = 0; i < M; i++) {

      if (i == 0)
      {
        ycoef_r(_, 0) = yPAR_r.column(indY(0));
        ycoef_beta(_, 0) = yPAR_beta.column(indY(0));
        ycoef_sig = yPAR_sig(indY(0));
      }
      else
      {
        ycoef_r(_, 0) = Rcpp::rnorm(7, mu_r, sqrt(Sigma_r));
        ycoef_beta(_, 0) = Rcpp::rnorm(7, mu_beta, sqrt(Sigma_beta));
        ycoef_sig = Rcpp::rgamma(1, 0.1, 0.1)(0);
      }

      x_cat0 = edp_rdirichlet_cpp(1, alpha00);
      x_var1 = (nu0 * tau0) / R::rchisq(nu0);
      x_var2 = (nu0 * tau0) / R::rchisq(nu0);
      x_var3 = (nu0 * tau0) / R::rchisq(nu0);
      x_var4 = (nu0 * tau0) / R::rchisq(nu0);
      x_var5 = (nu0 * tau0) / R::rchisq(nu0);
      x_mu1 = R::rnorm(0, sqrt(x_var1 / c0));
      x_mu2 = R::rnorm(0, sqrt(x_var2 / c0));
      x_mu3 = R::rnorm(0, sqrt(x_var3 / c0));
      x_mu4 = R::rnorm(0, sqrt(x_var4 / c0));
      x_mu5 = R::rnorm(0, sqrt(x_var5 / c0));

      x_mu_temp.col(xpar_mu.n_cols + numY * M + i) = arma::vec({x_mu1, x_mu2, x_mu3, x_mu4, x_mu5});
      x_sig_temp.col(xpar_sig.n_cols + numY * M + i) = arma::vec({x_var1, x_var2, x_var3, x_var4, x_var5});
      x_p0_temp.row(xpar_p0.n_rows + numY * M + i) = as<arma::mat>(x_cat0);

      arma::mat ypar_beta_temp(7, 1);
      ypar_beta_temp.col(0) = as<arma::mat>(ycoef_beta);
      ypar_beta = arma::join_rows(ypar_beta, ypar_beta_temp);

      arma::mat ypar_r_temp(7, 1);
      ypar_r_temp.col(0) = as<arma::mat>(ycoef_r);
      ypar_r = arma::join_rows(ypar_r, ypar_r_temp);

      int sz = ypar_sig.size();
      ypar_sig.resize(sz + 1);
      ypar_sig(sz) = ycoef_sig;

      int cumsum = 0;
      for (int j = 0; j < numY; j++) {
        cumsum += numXj(j) + M;

        x_cat0 = edp_rdirichlet_cpp(1, alpha00);
        x_var1 = (nu0 * tau0) / R::rchisq(nu0);
        x_var2 = (nu0 * tau0) / R::rchisq(nu0);
        x_var3 = (nu0 * tau0) / R::rchisq(nu0);
        x_var4 = (nu0 * tau0) / R::rchisq(nu0);
        x_var5 = (nu0 * tau0) / R::rchisq(nu0);
        x_mu1 = R::rnorm(0, sqrt(x_var1 / c0));
        x_mu2 = R::rnorm(0, sqrt(x_var2 / c0));
        x_mu3 = R::rnorm(0, sqrt(x_var3 / c0));
        x_mu4 = R::rnorm(0, sqrt(x_var4 / c0));
        x_mu5 = R::rnorm(0, sqrt(x_var5 / c0));

        if (j == numY - 1 & i == 0)
        {
          x_mu_temp.col(cumsum - 1 - i) = xpar_mu_hold.col(indX(0));
          x_sig_temp.col(cumsum - 1 - i) = xpar_sig_hold.col(indX(0));
          x_p0_temp.row(cumsum - 1 - i) = xpar_p0_hold.row(indX(0));
        }
        else
        {
          x_mu_temp.col(cumsum - 1 - i) = arma::vec({x_mu1, x_mu2, x_mu3, x_mu4, x_mu5});
          x_sig_temp.col(cumsum - 1 - i) = arma::vec({x_var1, x_var2, x_var3, x_var4, x_var5});
          x_p0_temp.row(cumsum - 1 - i) = as<arma::mat>(x_cat0);
        }
      }
    }

    int Cumsum = 0;
    int Cumsum1 = 0;
    for (int j = 0; j < numY; j++) {
      Cumsum1 += numXj(j) + M;
      for (int t = 0; t < numXj(j); t++) {
        x_mu_temp.col(Cumsum1 - numXj(j) - M + t) = xpar_mu.col(Cumsum);
        x_sig_temp.col(Cumsum1 - numXj(j) - M + t) = xpar_sig.col(Cumsum);
        x_p0_temp.row(Cumsum1 - numXj(j) - M + t) = xpar_p0.row(Cumsum);
        Cumsum++;
      }
    }

    xpar_mu = (x_mu_temp);
    xpar_sig = (x_sig_temp);
    xpar_p0 = (x_p0_temp);
  }

  int indbase = 0;
  // Membership log-probabilities (Eq. (3)/(4) Neal-8 weights). Evaluated in
  // log-space so the product of Gaussian/zero-inflated kernels cannot underflow
  // to an all-zero P; the zero-inflation probability is the model's
  // pi = expit(r) and the outcome kernel is the truncated-normal density on
  // (0, Inf) written in closed form, log N(y|b,sig) - log Phi(b/sig).
  NumericVector P = rep(R_NegInf, xpar_mu.n_cols);
  IntegerVector Py = rep(0, xpar_mu.n_cols);
  IntegerVector Px = rep(0, xpar_mu.n_cols);
  for (int j = 0; j < numY; j++) {
    for (int l = 0; l < numXj(j) + M; l++) {
      int idx = indbase + l;
      arma::mat r = as<arma::mat>(xa).row(h) * ypar_r.col(j);
      arma::mat b = as<arma::mat>(xa).row(h) * ypar_beta.col(j);
      double rr = r(0, 0), bb = b(0, 0), sd_y = sqrt(ypar_sig(j));
      double w = (l < numXj(j)) ? nji(j) * nlji(j, l) / (nji(j) + alpha_omega)
                                : nji(j) * (alpha_omega / M) / (nji(j) + alpha_omega);
      double log_kernel = (y(h) == 0)
          ? -log1p(exp(-rr))                                                          // log expit(rr) = log P(Y=0)
          : -log1p(exp(rr)) + R::dnorm(y(h), bb, sd_y, true) - R::pnorm(0.0, bb, sd_y, false, true); // log(1-pi)+log truncN
      P(idx) = log(w)
          + R::dbinom(xa(h, 1), 1, xpar_p0(idx, 1), true)
          + R::dnorm(xa(h, 2), xpar_mu(0, idx), sqrt(xpar_sig(0, idx)), true)
          + R::dnorm(xa(h, 3), xpar_mu(1, idx), sqrt(xpar_sig(1, idx)), true)
          + R::dnorm(xa(h, 4), xpar_mu(2, idx), sqrt(xpar_sig(2, idx)), true)
          + R::dnorm(xa(h, 5), xpar_mu(3, idx), sqrt(xpar_sig(3, idx)), true)
          + R::dnorm(xa(h, 6), xpar_mu(4, idx), sqrt(xpar_sig(4, idx)), true)
          + log_kernel;
      Py(idx) = j + 1;
      Px(idx) = l + 1;
    }
    indbase = indbase + numXj(j) + M;
  }
  for (int j = (numY); j < (numY + M); j++) {
    arma::mat r = as<arma::mat>(xa).row(h) * ypar_r.col(j);
    arma::mat b = as<arma::mat>(xa).row(h) * ypar_beta.col(j);
    double rr = r(0, 0), bb = b(0, 0), sd_y = sqrt(ypar_sig(j));
    Py(indbase) = j + 1;
    Px(indbase) = 1;
    double w = (alpha_theta / M) * (alpha_omega / M) / (alpha_omega);
    double log_kernel = (y(h) == 0)
        ? -log1p(exp(-rr))
        : -log1p(exp(rr)) + R::dnorm(y(h), bb, sd_y, true) - R::pnorm(0.0, bb, sd_y, false, true);
    P(indbase) = log(w)
        + R::dbinom(xa(h, 1), 1, xpar_p0(indbase, 1), true)
        + R::dnorm(xa(h, 2), xpar_mu(0, indbase), sqrt(xpar_sig(0, indbase)), true)
        + R::dnorm(xa(h, 3), xpar_mu(1, indbase), sqrt(xpar_sig(1, indbase)), true)
        + R::dnorm(xa(h, 4), xpar_mu(2, indbase), sqrt(xpar_sig(2, indbase)), true)
        + R::dnorm(xa(h, 5), xpar_mu(3, indbase), sqrt(xpar_sig(3, indbase)), true)
        + R::dnorm(xa(h, 6), xpar_mu(4, indbase), sqrt(xpar_sig(4, indbase)), true)
        + log_kernel;
    indbase++;
  }

  // exp-normalise (subtract max) so probabilities are valid and never all-zero
  double maxlp = R_NegInf;
  for (int i = 0; i < P.length(); i++) if (P(i) > maxlp) maxlp = P(i);
  for (int i = 0; i < P.length(); i++) P(i) = (P(i) == R_NegInf) ? 0.0 : exp(P(i) - maxlp);

  IntegerVector indsample = Rcpp::Range(0, P.length() - 1);
  int sSyx = sample(indsample, 1, false, P)(0);

  sy(h) = Py(sSyx);
  sx(h) = Px(sSyx);

  IntegerVector extraInd(numY + M);
  int extratemp = 0;
  for (int j = 0; j < numY; j++) {
    extratemp += numXj[j] + M;
    extraInd(j) = extratemp;
  }
  extraInd = extraInd - 1;
  extraInd(numY) = Py.length() - 1;

  if (sy(h) <= numY)
  {
    ypar_beta.shed_col(numY + M - 1);
    ypar_r.shed_col(numY + M - 1);
    ypar_sig = ypar_sig.elem(as<arma::uvec>(wrap(seq_len(numY))) - 1);

    if (sx(h) <= numXj((sy(h) - 1)))
    {
      xpar_mu.shed_cols(as<arma::uvec>(extraInd));
      xpar_sig.shed_cols(as<arma::uvec>(extraInd));
      xpar_p0.shed_rows(as<arma::uvec>(extraInd));
    }
    else
    {
      extraInd.erase((sy(h) - 1));
      xpar_mu.shed_cols(as<arma::uvec>(extraInd));
      xpar_sig.shed_cols(as<arma::uvec>(extraInd));
      xpar_p0.shed_rows(as<arma::uvec>(extraInd));
    }
  }
  else
  {
    extraInd.erase(numY);
    xpar_mu.shed_cols(as<arma::uvec>(extraInd));
    xpar_sig.shed_cols(as<arma::uvec>(extraInd));
    xpar_p0.shed_rows(as<arma::uvec>(extraInd));
  }

  Sy = Rcpp::IntegerVector(Rcpp::wrap((sy)));
  Sx = Rcpp::IntegerVector(Rcpp::wrap((sx)));

  numY = ypar_beta.n_cols;
  numX = xpar_p0.n_rows;

  yPAR_beta = Rcpp::NumericMatrix(Rcpp::wrap(ypar_beta));
  yPAR_r = Rcpp::NumericMatrix(Rcpp::wrap(ypar_r));
  yPAR_sig = Rcpp::NumericVector(Rcpp::wrap(ypar_sig));

  xPAR_p0 = Rcpp::NumericMatrix(Rcpp::wrap(xpar_p0));
  xPAR_mu = Rcpp::NumericMatrix(Rcpp::wrap(xpar_mu));
  xPAR_sig = Rcpp::NumericMatrix(Rcpp::wrap(xpar_sig));
  }

  List L = List::create(
    Named("Sx") = Sx,
    Named("Sy") = Sy,
    Named("xPAR_p0") = xPAR_p0,
    Named("xPAR_mu") = xPAR_mu,
    Named("xPAR_sig") = xPAR_sig,
    Named("yPAR_beta") = yPAR_beta,
    Named("yPAR_r") = yPAR_r,
    Named("yPAR_sig") = yPAR_sig
  );

  return L;
}
