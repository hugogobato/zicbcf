#include "arma_config.h"
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using namespace Rcpp;

// [[Rcpp::export]]
List pathc_bcfCore(arma::vec y_,
                   arma::vec I_,
                   arma::mat Omega_sel_con, arma::mat Omega_sel_mod,
                   arma::mat Omega_out_con, arma::mat Omega_out_mod,
                   arma::mat Omega_sel_con_est, arma::mat Omega_sel_mod_est,
                   arma::mat Omega_out_con_est, arma::mat Omega_out_mod_est,
                   NumericVector x_sel_con_, NumericVector x_sel_mod_,
                   NumericVector x_out_con_, NumericVector x_out_mod_,
                   NumericVector x_sel_con_est_, NumericVector x_sel_mod_est_,
                   NumericVector x_out_con_est_, NumericVector x_out_mod_est_,
                   List x_sel_con_info_list, List x_sel_mod_info_list,
                   List x_out_con_info_list, List x_out_mod_info_list,
                   int burn, int nd, int thin,
                   int ntree_sel_con, int ntree_sel_mod,
                   int ntree_out_con, int ntree_out_mod,
                   double lambda, double nu, // prior parameters for sigma0^2
                   arma::mat Sigma0_sel_con, arma::mat Sigma0_sel_mod,
                   arma::mat Sigma0_out_con, arma::mat Sigma0_out_mod,
                   double sel_con_alpha, double sel_con_beta,
                   double sel_mod_alpha, double sel_mod_beta,
                   double out_con_alpha, double out_con_beta,
                   double out_mod_alpha, double out_mod_beta,
                   double con_scale_df = 1, double mod_scale_df = -1,
                   int status_interval = 100)
{
  RNGScope scope;
  RNG gen;

  size_t n = y_.n_elem;

  //-------------------------------------------------------------------------//
  // Read and format y_ and I_                                               //
  //-------------------------------------------------------------------------//
  std::vector<double> y;
  std::vector<double> I;
  for(NumericVector::iterator it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);
  }
  for(NumericVector::iterator it = I_.begin(); it != I_.end(); ++it) {
    I.push_back(*it);
  }

  //-------------------------------------------------------------------------//
  // Read and format Selection Design and Covariate Matrices                 //
  //-------------------------------------------------------------------------//
  std::vector<double> x_sel_con;
  for(NumericVector::iterator it = x_sel_con_.begin(); it != x_sel_con_.end(); ++it) {
    x_sel_con.push_back(*it);
  }
  size_t p_sel_con = x_sel_con.size()/n;

  xinfo xi_sel_con;
  xi_sel_con.resize(p_sel_con);
  for(int i=0; i<p_sel_con; ++i) {
    NumericVector tmp = x_sel_con_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_sel_con[i] = tmp2;
  }

  std::vector<double> x_sel_mod;
  for(NumericVector::iterator it = x_sel_mod_.begin(); it != x_sel_mod_.end(); ++it) {
    x_sel_mod.push_back(*it);
  }
  size_t p_sel_mod = x_sel_mod.size()/n;

  xinfo xi_sel_mod;
  xi_sel_mod.resize(p_sel_mod);
  for(int i=0; i<p_sel_mod; ++i) {
    NumericVector tmp = x_sel_mod_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_sel_mod[i] = tmp2;
  }

  //-------------------------------------------------------------------------//
  // Read and format Outcome Design and Covariate Matrices                   //
  //-------------------------------------------------------------------------//
  std::vector<double> x_out_con;
  for(NumericVector::iterator it = x_out_con_.begin(); it != x_out_con_.end(); ++it) {
    x_out_con.push_back(*it);
  }
  size_t p_out_con = x_out_con.size()/n;

  xinfo xi_out_con;
  xi_out_con.resize(p_out_con);
  for(int i=0; i<p_out_con; ++i) {
    NumericVector tmp = x_out_con_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_out_con[i] = tmp2;
  }

  std::vector<double> x_out_mod;
  for(NumericVector::iterator it = x_out_mod_.begin(); it != x_out_mod_.end(); ++it) {
    x_out_mod.push_back(*it);
  }
  size_t p_out_mod = x_out_mod.size()/n;

  xinfo xi_out_mod;
  xi_out_mod.resize(p_out_mod);
  for(int i=0; i<p_out_mod; ++i) {
    NumericVector tmp = x_out_mod_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_out_mod[i] = tmp2;
  }

  //-------------------------------------------------------------------------//
  // Read and format out-of-sample prediction data                           //
  //-------------------------------------------------------------------------//
  std::vector<double> x_sel_con_est;
  for(NumericVector::iterator it = x_sel_con_est_.begin(); it != x_sel_con_est_.end(); ++it) {
    x_sel_con_est.push_back(*it);
  }
  size_t n_sel_con_est = x_sel_con_est.size()/p_sel_con;

  dinfo di_sel_con_est;
  di_sel_con_est.n = n_sel_con_est;
  di_sel_con_est.p = p_sel_con;
  if(n_sel_con_est > 0) {
    di_sel_con_est.x = &x_sel_con_est[0];
    di_sel_con_est.basis_dim = Omega_sel_con_est.n_rows;
    di_sel_con_est.omega = &Omega_sel_con_est[0];
  }

  std::vector<double> x_sel_mod_est;
  for(NumericVector::iterator it = x_sel_mod_est_.begin(); it != x_sel_mod_est_.end(); ++it) {
    x_sel_mod_est.push_back(*it);
  }
  size_t n_sel_mod_est = x_sel_mod_est.size()/p_sel_mod;

  dinfo di_sel_mod_est;
  di_sel_mod_est.n = n_sel_mod_est;
  di_sel_mod_est.p = p_sel_mod;
  if(n_sel_mod_est > 0) {
    di_sel_mod_est.x = &x_sel_mod_est[0];
    di_sel_mod_est.basis_dim = Omega_sel_mod_est.n_rows;
    di_sel_mod_est.omega = &Omega_sel_mod_est[0];
  }

  std::vector<double> x_out_con_est;
  for(NumericVector::iterator it = x_out_con_est_.begin(); it != x_out_con_est_.end(); ++it) {
    x_out_con_est.push_back(*it);
  }
  size_t n_out_con_est = x_out_con_est.size()/p_out_con;

  dinfo di_out_con_est;
  di_out_con_est.n = n_out_con_est;
  di_out_con_est.p = p_out_con;
  if(n_out_con_est > 0) {
    di_out_con_est.x = &x_out_con_est[0];
    di_out_con_est.basis_dim = Omega_out_con_est.n_rows;
    di_out_con_est.omega = &Omega_out_con_est[0];
  }

  std::vector<double> x_out_mod_est;
  for(NumericVector::iterator it = x_out_mod_est_.begin(); it != x_out_mod_est_.end(); ++it) {
    x_out_mod_est.push_back(*it);
  }
  size_t n_out_mod_est = x_out_mod_est.size()/p_out_mod;

  dinfo di_out_mod_est;
  di_out_mod_est.n = n_out_mod_est;
  di_out_mod_est.p = p_out_mod;
  if(n_out_mod_est > 0) {
    di_out_mod_est.x = &x_out_mod_est[0];
    di_out_mod_est.basis_dim = Omega_out_mod_est.n_rows;
    di_out_mod_est.omega = &Omega_out_mod_est[0];
  }

  //-------------------------------------------------------------------------//
  // Initialize Forests                                                      //
  //-------------------------------------------------------------------------//
  std::vector<tree> t_sel_con(ntree_sel_con);
  for(size_t i=0; i<ntree_sel_con; i++) t_sel_con[i].setm(zeros(Omega_sel_con.n_rows));

  std::vector<tree> t_sel_mod(ntree_sel_mod);
  for(size_t i=0; i<ntree_sel_mod; i++) t_sel_mod[i].setm(zeros(Omega_sel_mod.n_rows));

  std::vector<tree> t_out_con(ntree_out_con);
  for(size_t i=0; i<ntree_out_con; i++) t_out_con[i].setm(zeros(Omega_out_con.n_rows));

  std::vector<tree> t_out_mod(ntree_out_mod);
  for(size_t i=0; i<ntree_out_mod; i++) t_out_mod[i].setm(zeros(Omega_out_mod.n_rows));

  //-------------------------------------------------------------------------//
  // Prior info parameters                                                   //
  //-------------------------------------------------------------------------//
  pinfo pi_sel_con;
  pi_sel_con.pbd = 1.0; pi_sel_con.pb = 0.5;
  pi_sel_con.alpha = sel_con_alpha; pi_sel_con.beta = sel_con_beta;
  pi_sel_con.mu0 = zeros(Omega_sel_con.n_rows);
  pi_sel_con.Sigma0 = Sigma0_sel_con; pi_sel_con.Prec0 = pi_sel_con.Sigma0.i();
  pi_sel_con.logdetSigma0 = log(det(pi_sel_con.Sigma0));
  pi_sel_con.eta = 1.0; pi_sel_con.gamma = 1.0;
  pi_sel_con.scale_df = con_scale_df;
  pi_sel_con.vanilla = true;

  pinfo pi_sel_mod;
  pi_sel_mod.pbd = 1.0; pi_sel_mod.pb = 0.5;
  pi_sel_mod.alpha = sel_mod_alpha; pi_sel_mod.beta = sel_mod_beta;
  pi_sel_mod.mu0 = zeros(Omega_sel_mod.n_rows);
  pi_sel_mod.Sigma0 = Sigma0_sel_mod; pi_sel_mod.Prec0 = pi_sel_mod.Sigma0.i();
  pi_sel_mod.logdetSigma0 = log(det(pi_sel_mod.Sigma0));
  pi_sel_mod.eta = 1.0; pi_sel_mod.gamma = 1.0;
  pi_sel_mod.scale_df = mod_scale_df;
  pi_sel_mod.vanilla = false;

  pinfo pi_out_con;
  pi_out_con.pbd = 1.0; pi_out_con.pb = 0.5;
  pi_out_con.alpha = out_con_alpha; pi_out_con.beta = out_con_beta;
  pi_out_con.mu0 = zeros(Omega_out_con.n_rows);
  pi_out_con.Sigma0 = Sigma0_out_con; pi_out_con.Prec0 = pi_out_con.Sigma0.i();
  pi_out_con.logdetSigma0 = log(det(pi_out_con.Sigma0));
  pi_out_con.eta = 1.0; pi_out_con.gamma = 1.0;
  pi_out_con.scale_df = con_scale_df;
  pi_out_con.vanilla = true;

  pinfo pi_out_mod;
  pi_out_mod.pbd = 1.0; pi_out_mod.pb = 0.5;
  pi_out_mod.alpha = out_mod_alpha; pi_out_mod.beta = out_mod_beta;
  pi_out_mod.mu0 = zeros(Omega_out_mod.n_rows);
  pi_out_mod.Sigma0 = Sigma0_out_mod; pi_out_mod.Prec0 = pi_out_mod.Sigma0.i();
  pi_out_mod.logdetSigma0 = log(det(pi_out_mod.Sigma0));
  pi_out_mod.eta = 1.0; pi_out_mod.gamma = 1.0;
  pi_out_mod.scale_df = mod_scale_df;
  pi_out_mod.vanilla = false;

  //-------------------------------------------------------------------------//
  // Setup data wrappers (dinfo)                                             //
  //-------------------------------------------------------------------------//
  double* allfit_sel_con = new double[n];
  double* allfit_sel_mod = new double[n];
  double* allfit_out_con = new double[n];
  double* allfit_out_mod = new double[n];

  for(size_t i=0; i<n; i++) {
    allfit_sel_con[i] = 0.0;
    allfit_sel_mod[i] = 0.0;
    allfit_out_con[i] = 0.0;
    allfit_out_mod[i] = 0.0;
  }

  double* r_sel_con = new double[n];
  double* r_sel_mod = new double[n];
  double* r_out_con = new double[n];
  double* r_out_mod = new double[n];

  dinfo di_sel_con;
  di_sel_con.n = n; di_sel_con.p = p_sel_con; di_sel_con.x = &x_sel_con[0]; di_sel_con.y = r_sel_con;
  di_sel_con.basis_dim = Omega_sel_con.n_rows; di_sel_con.omega = &Omega_sel_con[0];

  dinfo di_sel_mod;
  di_sel_mod.n = n; di_sel_mod.p = p_sel_mod; di_sel_mod.x = &x_sel_mod[0]; di_sel_mod.y = r_sel_mod;
  di_sel_mod.basis_dim = Omega_sel_mod.n_rows; di_sel_mod.omega = &Omega_sel_mod[0];

  dinfo di_out_con;
  di_out_con.n = n; di_out_con.p = p_out_con; di_out_con.x = &x_out_con[0]; di_out_con.y = r_out_con;
  di_out_con.basis_dim = Omega_out_con.n_rows; di_out_con.omega = &Omega_out_con[0];

  dinfo di_out_mod;
  di_out_mod.n = n; di_out_mod.p = p_out_mod; di_out_mod.x = &x_out_mod[0]; di_out_mod.y = r_out_mod;
  di_out_mod.basis_dim = Omega_out_mod.n_rows; di_out_mod.omega = &Omega_out_mod[0];

  double* ftemp = new double[n];

  // Initialize node pointers
  std::vector<std::vector<tree::tree_cp> > node_pointers_sel_con(ntree_sel_con);
  std::vector<std::vector<tree::tree_cp> > node_pointers_sel_mod(ntree_sel_mod);
  std::vector<std::vector<tree::tree_cp> > node_pointers_out_con(ntree_out_con);
  std::vector<std::vector<tree::tree_cp> > node_pointers_out_mod(ntree_out_mod);

  for(size_t j=0; j<ntree_sel_con; ++j) {
    node_pointers_sel_con[j].resize(n);
    fit_basis(t_sel_con[j], xi_sel_con, di_sel_con, ftemp, node_pointers_sel_con[j], true, true);
  }
  for(size_t j=0; j<ntree_sel_mod; ++j) {
    node_pointers_sel_mod[j].resize(n);
    fit_basis(t_sel_mod[j], xi_sel_mod, di_sel_mod, ftemp, node_pointers_sel_mod[j], true, false);
  }
  for(size_t j=0; j<ntree_out_con; ++j) {
    node_pointers_out_con[j].resize(n);
    fit_basis(t_out_con[j], xi_out_con, di_out_con, ftemp, node_pointers_out_con[j], true, true);
  }
  for(size_t j=0; j<ntree_out_mod; ++j) {
    node_pointers_out_mod[j].resize(n);
    fit_basis(t_out_mod[j], xi_out_mod, di_out_mod, ftemp, node_pointers_out_mod[j], true, false);
  }

  //-------------------------------------------------------------------------//
  // Initialize Latent variables for selection and augmented outcome        //
  //-------------------------------------------------------------------------//
  std::vector<double> W(n);
  std::vector<double> V(n);

  for(size_t i=0; i<n; i++) {
    if(I[i] == 1.0) {
      W[i] = 1.0;
      V[i] = y[i];
    } else {
      W[i] = -1.0;
      V[i] = 0.0;
    }
  }

  // Initial joint covariance state
  double beta = 0.0;
  double sigma02 = 1.0;
  double sigma0 = 1.0;

  //-------------------------------------------------------------------------//
  // Storage for MCMC Output Draws                                           //
  //-------------------------------------------------------------------------//
  NumericMatrix sel_con_post(nd, n);
  NumericMatrix sel_mod_post(nd, n);
  NumericMatrix out_con_post(nd, n);
  NumericMatrix out_mod_post(nd, n);
  NumericMatrix yhat_post(nd, n);
  NumericMatrix sel_tau_post(nd, n);
  NumericMatrix out_tau_post(nd, n);

  NumericVector beta_post(nd);
  NumericVector sigma0_post(nd);
  NumericVector sigma_post(nd);
  NumericVector rho_post(nd);

  NumericMatrix W_post(nd, n);
  NumericMatrix V_post(nd, n);

  // Out of sample storage if needed
  NumericMatrix sel_con_est_post(nd, n_sel_con_est);
  NumericMatrix sel_mod_est_post(nd, n_sel_mod_est);
  NumericMatrix out_con_est_post(nd, n_out_con_est);
  NumericMatrix out_mod_est_post(nd, n_out_mod_est);

  //-------------------------------------------------------------------------//
  // MCMC Loop                                                               //
  //-------------------------------------------------------------------------//
  Rcout << "\nBeginning Path C Joint Copula-BCF MCMC:\n";
  time_t tp;
  int time1 = time(&tp);
  size_t save_ctr = 0;

  size_t total_iterations = nd * thin + burn;

  for(size_t i = 0; i < total_iterations; i++) {
    Rcpp::checkUserInterrupt();

    if(i % status_interval == 0) {
      Rcout << "iteration: " << i << " of " << total_iterations << " | beta: " << beta << " | sigma0: " << sigma0 << endl;
    }

    // Compute joint covariance values
    double rho = beta / sqrt(sigma02 + beta*beta);
    double sigma2 = sigma02 + beta*beta;
    double sigma = sqrt(sigma2);
    double var_sel_active = 1.0 - rho*rho;
    double sd_sel_active = sqrt(var_sel_active);

    // --- Step 1: Draw latent selection utilities W_i* ---
    for(size_t k=0; k<n; k++) {
      double eta_b_k = allfit_sel_con[k] + allfit_sel_mod[k];
      double eta_c_k = allfit_out_con[k] + allfit_out_mod[k];

      if(I[k] == 1.0) {
        // active: truncated to (0, infinity)
        double mean_active = eta_b_k + (beta / sigma2) * (V[k] - eta_c_k);
        W[k] = rtnormlo(mean_active, sd_sel_active, 0.0);
      } else {
        // inactive: truncated to (-infinity, 0)
        W[k] = -rtnormlo(-eta_b_k, 1.0, 0.0);
      }
    }

    // --- Step 2: Draw augmented log-intensity V_k for inactive units ---
    for(size_t k=0; k<n; k++) {
      if(I[k] != 1.0) {
        double eta_b_k = allfit_sel_con[k] + allfit_sel_mod[k];
        double eta_c_k = allfit_out_con[k] + allfit_out_mod[k];
        double mean_V = eta_c_k + beta * (W[k] - eta_b_k);
        V[k] = gen.normal(mean_V, sigma0);
      }
    }

    // --- Step 3: Update Selection forests (t_sel_con, t_sel_mod) ---
    double sd_sel = sqrt(1.0 - rho*rho);
    pi_sel_con.sigma = sd_sel / fabs(pi_sel_con.eta);
    pi_sel_mod.sigma = sd_sel / fabs(pi_sel_mod.eta);

    // Prognostic Selection Forest (vanilla = true)
    for(size_t j=0; j<ntree_sel_con; j++) {
      fit_basis(t_sel_con[j], xi_sel_con, di_sel_con, ftemp, node_pointers_sel_con[j], false, true);

      for(size_t k=0; k<n; k++) {
        double target_sel = W[k] - (beta / sigma2) * (V[k] - (allfit_out_con[k] + allfit_out_mod[k]));
        allfit_sel_con[k] = allfit_sel_con[k] - pi_sel_con.eta * ftemp[k];
        r_sel_con[k] = (target_sel - (allfit_sel_con[k] + allfit_sel_mod[k])) / pi_sel_con.eta;
      }

      bd_basis(t_sel_con[j], xi_sel_con, di_sel_con, pi_sel_con, gen, node_pointers_sel_con[j]);
      drmu_basis(t_sel_con[j], xi_sel_con, di_sel_con, pi_sel_con, gen);
      fit_basis(t_sel_con[j], xi_sel_con, di_sel_con, ftemp, node_pointers_sel_con[j], false, true);

      for(size_t k=0; k<n; k++) {
        allfit_sel_con[k] += pi_sel_con.eta * ftemp[k];
      }
    }

    // Moderating Selection Forest (vanilla = false)
    for(size_t j=0; j<ntree_sel_mod; j++) {
      fit_basis(t_sel_mod[j], xi_sel_mod, di_sel_mod, ftemp, node_pointers_sel_mod[j], false, false);

      for(size_t k=0; k<n; k++) {
        double target_sel = W[k] - (beta / sigma2) * (V[k] - (allfit_out_con[k] + allfit_out_mod[k]));
        allfit_sel_mod[k] = allfit_sel_mod[k] - pi_sel_mod.eta * ftemp[k];
        r_sel_mod[k] = (target_sel - (allfit_sel_con[k] + allfit_sel_mod[k])) / pi_sel_mod.eta;
      }

      bd_basis(t_sel_mod[j], xi_sel_mod, di_sel_mod, pi_sel_mod, gen, node_pointers_sel_mod[j]);
      drmu_basis(t_sel_mod[j], xi_sel_mod, di_sel_mod, pi_sel_mod, gen);
      fit_basis(t_sel_mod[j], xi_sel_mod, di_sel_mod, ftemp, node_pointers_sel_mod[j], false, false);

      for(size_t k=0; k<n; k++) {
        allfit_sel_mod[k] += pi_sel_mod.eta * ftemp[k];
      }
    }

    // --- Step 4: Update Outcome forests (t_out_con, t_out_mod) ---
    pi_out_con.sigma = sigma0 / fabs(pi_out_con.eta);
    pi_out_mod.sigma = sigma0 / fabs(pi_out_mod.eta);

    // Prognostic Outcome Forest (vanilla = true)
    for(size_t j=0; j<ntree_out_con; j++) {
      fit_basis(t_out_con[j], xi_out_con, di_out_con, ftemp, node_pointers_out_con[j], false, true);

      for(size_t k=0; k<n; k++) {
        double target_out = V[k] - beta * (W[k] - (allfit_sel_con[k] + allfit_sel_mod[k]));
        allfit_out_con[k] = allfit_out_con[k] - pi_out_con.eta * ftemp[k];
        r_out_con[k] = (target_out - (allfit_out_con[k] + allfit_out_mod[k])) / pi_out_con.eta;
      }

      bd_basis(t_out_con[j], xi_out_con, di_out_con, pi_out_con, gen, node_pointers_out_con[j]);
      drmu_basis(t_out_con[j], xi_out_con, di_out_con, pi_out_con, gen);
      fit_basis(t_out_con[j], xi_out_con, di_out_con, ftemp, node_pointers_out_con[j], false, true);

      for(size_t k=0; k<n; k++) {
        allfit_out_con[k] += pi_out_con.eta * ftemp[k];
      }
    }

    // Moderating Outcome Forest (vanilla = false)
    for(size_t j=0; j<ntree_out_mod; j++) {
      fit_basis(t_out_mod[j], xi_out_mod, di_out_mod, ftemp, node_pointers_out_mod[j], false, false);

      for(size_t k=0; k<n; k++) {
        double target_out = V[k] - beta * (W[k] - (allfit_sel_con[k] + allfit_sel_mod[k]));
        allfit_out_mod[k] = allfit_out_mod[k] - pi_out_mod.eta * ftemp[k];
        r_out_mod[k] = (target_out - (allfit_out_con[k] + allfit_out_mod[k])) / pi_out_mod.eta;
      }

      bd_basis(t_out_mod[j], xi_out_mod, di_out_mod, pi_out_mod, gen, node_pointers_out_mod[j]);
      drmu_basis(t_out_mod[j], xi_out_mod, di_out_mod, pi_out_mod, gen);
      fit_basis(t_out_mod[j], xi_out_mod, di_out_mod, ftemp, node_pointers_out_mod[j], false, false);

      for(size_t k=0; k<n; k++) {
        allfit_out_mod[k] += pi_out_mod.eta * ftemp[k];
      }
    }

    // --- Step 5: Update PX Scaling Factors (eta) ---
    double eta_old;

    // Selection Prognostic scale
    for(size_t k=0; k<n; k++) {
      double target_sel = W[k] - (beta / sigma2) * (V[k] - (allfit_out_con[k] + allfit_out_mod[k]));
      ftemp[k] = target_sel - allfit_sel_mod[k];
    }
    eta_old = pi_sel_con.eta;
    update_scale(ftemp, allfit_sel_con, n, sd_sel, pi_sel_con, gen);
    for(size_t k=0; k<n; ++k) {
      allfit_sel_con[k] = allfit_sel_con[k] * pi_sel_con.eta / eta_old;
    }
    pi_sel_con.sigma = sd_sel / fabs(pi_sel_con.eta);

    // Selection Moderating scale
    for(size_t k=0; k<n; k++) {
      double target_sel = W[k] - (beta / sigma2) * (V[k] - (allfit_out_con[k] + allfit_out_mod[k]));
      ftemp[k] = target_sel - allfit_sel_con[k];
    }
    eta_old = pi_sel_mod.eta;
    update_scale(ftemp, allfit_sel_mod, n, sd_sel, pi_sel_mod, gen);
    for(size_t k=0; k<n; ++k) {
      allfit_sel_mod[k] = allfit_sel_mod[k] * pi_sel_mod.eta / eta_old;
    }
    pi_sel_mod.sigma = sd_sel / fabs(pi_sel_mod.eta);

    // Outcome Prognostic scale
    for(size_t k=0; k<n; k++) {
      double target_out = V[k] - beta * (W[k] - (allfit_sel_con[k] + allfit_sel_mod[k]));
      ftemp[k] = target_out - allfit_out_mod[k];
    }
    eta_old = pi_out_con.eta;
    update_scale(ftemp, allfit_out_con, n, sigma0, pi_out_con, gen);
    for(size_t k=0; k<n; ++k) {
      allfit_out_con[k] = allfit_out_con[k] * pi_out_con.eta / eta_old;
    }
    pi_out_con.sigma = sigma0 / fabs(pi_out_con.eta);

    // Outcome Moderating scale
    for(size_t k=0; k<n; k++) {
      double target_out = V[k] - beta * (W[k] - (allfit_sel_con[k] + allfit_sel_mod[k]));
      ftemp[k] = target_out - allfit_out_con[k];
    }
    eta_old = pi_out_mod.eta;
    update_scale(ftemp, allfit_out_mod, n, sigma0, pi_out_mod, gen);
    for(size_t k=0; k<n; ++k) {
      allfit_out_mod[k] = allfit_out_mod[k] * pi_out_mod.eta / eta_old;
    }
    pi_out_mod.sigma = sigma0 / fabs(pi_out_mod.eta);

    // --- Step 6: Joint Conjugate Gibbs Update for (beta, sigma0^2) ---
    double S_xx = 0.0;
    double S_xy = 0.0;
    double S_yy = 0.0;

    for(size_t k=0; k<n; k++) {
      double eta_b_k = allfit_sel_con[k] + allfit_sel_mod[k];
      double eta_c_k = allfit_out_con[k] + allfit_out_mod[k];

      double delta_k = W[k] - eta_b_k;
      double epsilon_k = V[k] - eta_c_k;

      S_xx += delta_k * delta_k;
      S_xy += delta_k * epsilon_k;
      S_yy += epsilon_k * epsilon_k;
    }

    double v_beta = 100.0; // diffuse prior variance scale for beta
    double V_beta = 1.0 / (1.0 / v_beta + S_xx);
    double M_beta = V_beta * S_xy;

    double a0 = nu / 2.0;
    double b0 = nu * lambda / 2.0;

    double a_N = a0 + n / 2.0;
    double b_N = b0 + 0.5 * (S_yy - M_beta * M_beta / V_beta);

    if (b_N <= 0) b_N = 1e-6;

    // Draw sigma0^2 from Inverse-Gamma
    double gamma_draw = R::rgamma(a_N, 1.0);
    sigma02 = b_N / gamma_draw;
    sigma0 = sqrt(sigma02);

    // Draw beta from Normal
    beta = M_beta + gen.normal(0.0, 1.0) * sqrt(sigma02 * V_beta);

    //-----------------------------------------------------------------------//
    // Save posterior samples                                                //
    //-----------------------------------------------------------------------//
    if((i >= burn) && ((i - burn) % thin == 0)) {
      beta_post(save_ctr) = beta;
      sigma0_post(save_ctr) = sigma0;
      sigma_post(save_ctr) = sqrt(sigma02 + beta*beta);
      rho_post(save_ctr) = beta / sqrt(sigma02 + beta*beta);

      for(size_t k=0; k<n; k++) {
        sel_con_post(save_ctr, k) = allfit_sel_con[k];
        sel_mod_post(save_ctr, k) = allfit_sel_mod[k];
        out_con_post(save_ctr, k) = allfit_out_con[k];
        out_mod_post(save_ctr, k) = allfit_out_mod[k];
        yhat_post(save_ctr, k) = allfit_out_con[k] + allfit_out_mod[k];

        sel_tau_post(save_ctr, k) = pi_sel_mod.eta * fit_i_basis(k, t_sel_mod, xi_sel_mod, di_sel_mod, true);
        out_tau_post(save_ctr, k) = pi_out_mod.eta * fit_i_basis(k, t_out_mod, xi_out_mod, di_out_mod, true);

        W_post(save_ctr, k) = W[k];
        V_post(save_ctr, k) = V[k];
      }

      // Out of sample estimations
      if(n_sel_con_est > 0) {
        for(size_t k=0; k<n_sel_con_est; k++) {
          sel_con_est_post(save_ctr, k) = pi_sel_con.eta * fit_i_basis(k, t_sel_con, xi_sel_con, di_sel_con_est, true);
        }
      }
      if(n_sel_mod_est > 0) {
        for(size_t k=0; k<n_sel_mod_est; k++) {
          sel_mod_est_post(save_ctr, k) = pi_sel_mod.eta * fit_i_basis(k, t_sel_mod, xi_sel_mod, di_sel_mod_est, false);
        }
      }
      if(n_out_con_est > 0) {
        for(size_t k=0; k<n_out_con_est; k++) {
          out_con_est_post(save_ctr, k) = pi_out_con.eta * fit_i_basis(k, t_out_con, xi_out_con, di_out_con_est, true);
        }
      }
      if(n_out_mod_est > 0) {
        for(size_t k=0; k<n_out_mod_est; k++) {
          out_mod_est_post(save_ctr, k) = pi_out_mod.eta * fit_i_basis(k, t_out_mod, xi_out_mod, di_out_mod_est, false);
        }
      }

      save_ctr += 1;
    }
  }

  int time2 = time(&tp);
  Rcout << "MCMC Loop took " << time2 - time1 << " seconds.\n";

  // Clean memory
  delete[] allfit_sel_con;
  delete[] allfit_sel_mod;
  delete[] allfit_out_con;
  delete[] allfit_out_mod;
  delete[] r_sel_con;
  delete[] r_sel_mod;
  delete[] r_out_con;
  delete[] r_out_mod;
  delete[] ftemp;

  return(List::create(_["sel_con_post"] = sel_con_post,
                      _["sel_mod_post"] = sel_mod_post,
                      _["out_con_post"] = out_con_post,
                      _["out_mod_post"] = out_mod_post,
                      _["yhat_post"] = yhat_post,
                      _["beta_post"] = beta_post,
                      _["sigma0_post"] = sigma0_post,
                      _["sigma_post"] = sigma_post,
                      _["rho_post"] = rho_post,
                      _["W_post"] = W_post,
                      _["V_post"] = V_post,
                      _["sel_con_est_post"] = sel_con_est_post,
                      _["sel_mod_est_post"] = sel_mod_est_post,
                      _["out_con_est_post"] = out_con_est_post,
                      _["out_mod_est_post"] = out_mod_est_post,
                      _["sel_tau_post"] = sel_tau_post,
                      _["out_tau_post"] = out_tau_post
  ));
}
