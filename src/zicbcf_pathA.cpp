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
List zicbcfCore_pathA(
    arma::vec y_hurdle, arma::mat Omega_con_hurdle, arma::mat Omega_mod_hurdle,
    NumericVector x_con_hurdle_, NumericVector x_mod_hurdle_,
    List x_con_info_hurdle_list, List x_mod_info_hurdle_list,
    int ntree_con_hurdle, int ntree_mod_hurdle,
    arma::mat Sigma0_con_hurdle, arma::mat Sigma0_mod_hurdle,
    double con_alpha_hurdle, double con_beta_hurdle,
    double mod_alpha_hurdle, double mod_beta_hurdle,
    bool vanilla_hurdle,
    bool use_con_scale_hurdle, bool use_mod_scale_hurdle,
    double con_scale_df_hurdle, double mod_scale_df_hurdle,
    arma::vec y_continuous, arma::mat Omega_con_continuous, arma::mat Omega_mod_continuous,
    NumericVector x_con_continuous_, NumericVector x_mod_continuous_,
    List x_con_info_continuous_list, List x_mod_info_continuous_list,
    int ntree_con_continuous, int ntree_mod_continuous,
    arma::mat Sigma0_con_continuous, arma::mat Sigma0_mod_continuous,
    double con_alpha_continuous, double con_beta_continuous,
    double mod_alpha_continuous, double mod_beta_continuous,
    bool vanilla_continuous,
    bool use_con_scale_continuous, bool use_mod_scale_continuous,
    double con_scale_df_continuous, double mod_scale_df_continuous,
    NumericVector x_con_continuous_est_, NumericVector x_mod_continuous_est_,
    arma::mat Omega_con_continuous_est, arma::mat Omega_mod_continuous_est,
    int burn, int nd, int thin,
    double lambda, double nu,
    int status_interval = 100
) {

  RNGScope scope;
  RNG gen;

  //-------------------------------------------------------------------------//
  // 1. HURDLE PART VARIABLE SETUP                                           //
  //-------------------------------------------------------------------------//
  size_t n_b = y_hurdle.size();
  std::vector<double> y_lat_hurdle(n_b);
  for(size_t i=0; i<n_b; i++) {
    y_lat_hurdle[i] = y_hurdle[i] == 1 ? 1.0 : -1.0;
  }

  std::vector<double> x_con_hurdle;
  for(NumericVector::iterator it=x_con_hurdle_.begin(); it!=x_con_hurdle_.end(); ++it) {
    x_con_hurdle.push_back(*it);
  }
  size_t p_con_hurdle = x_con_hurdle.size()/n_b;

  xinfo xi_con_hurdle;
  xi_con_hurdle.resize(p_con_hurdle);
  for(int i=0; i<p_con_hurdle; ++i) {
    NumericVector tmp = x_con_info_hurdle_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_con_hurdle[i] = tmp2;
  }

  std::vector<double> x_mod_hurdle;
  for(NumericVector::iterator it=x_mod_hurdle_.begin(); it!=x_mod_hurdle_.end(); ++it) {
    x_mod_hurdle.push_back(*it);
  }
  size_t p_mod_hurdle = x_mod_hurdle.size()/n_b;

  xinfo xi_mod_hurdle;
  xi_mod_hurdle.resize(p_mod_hurdle);
  for(int i=0; i<p_mod_hurdle; ++i) {
    NumericVector tmp = x_mod_info_hurdle_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_mod_hurdle[i] = tmp2;
  }

  std::vector<tree> t_con_hurdle(ntree_con_hurdle);
  for(size_t i=0; i<ntree_con_hurdle; i++) t_con_hurdle[i].setm(zeros(Omega_con_hurdle.n_rows));

  std::vector<tree> t_mod_hurdle(ntree_mod_hurdle);
  for(size_t i=0; i<ntree_mod_hurdle; i++) t_mod_hurdle[i].setm(zeros(Omega_mod_hurdle.n_rows));

  pinfo pi_con_hurdle;
  pi_con_hurdle.pbd = 1.0; pi_con_hurdle.pb = 0.5;
  pi_con_hurdle.alpha = con_alpha_hurdle; pi_con_hurdle.beta = con_beta_hurdle;
  pi_con_hurdle.sigma = 1.0;
  pi_con_hurdle.mu0 = zeros(Omega_con_hurdle.n_rows);
  pi_con_hurdle.Sigma0 = Sigma0_con_hurdle;
  pi_con_hurdle.Prec0 = pi_con_hurdle.Sigma0.i();
  pi_con_hurdle.logdetSigma0 = log(det(pi_con_hurdle.Sigma0));
  pi_con_hurdle.eta = 1.0; pi_con_hurdle.gamma = 1.0;
  pi_con_hurdle.scale_df = con_scale_df_hurdle;
  pi_con_hurdle.dart = false;

  pinfo pi_mod_hurdle;
  pi_mod_hurdle.pbd = 1.0; pi_mod_hurdle.pb = 0.5;
  pi_mod_hurdle.alpha = mod_alpha_hurdle; pi_mod_hurdle.beta = mod_beta_hurdle;
  pi_mod_hurdle.sigma = 1.0;
  pi_mod_hurdle.mu0 = zeros(Omega_mod_hurdle.n_rows);
  pi_mod_hurdle.Sigma0 = Sigma0_mod_hurdle;
  pi_mod_hurdle.Prec0 = pi_mod_hurdle.Sigma0.i();
  pi_mod_hurdle.logdetSigma0 = log(det(pi_mod_hurdle.Sigma0));
  pi_mod_hurdle.eta = 1.0; pi_mod_hurdle.gamma = 1.0;
  pi_mod_hurdle.scale_df = mod_scale_df_hurdle;
  pi_mod_hurdle.dart = false;

  std::vector<double> allfit_con_hurdle(n_b, 0.0);
  std::vector<double> r_con_hurdle(n_b);

  dinfo di_con_hurdle;
  di_con_hurdle.n = n_b; di_con_hurdle.p = p_con_hurdle;
  di_con_hurdle.x = x_con_hurdle.data(); di_con_hurdle.y = r_con_hurdle.data();
  di_con_hurdle.basis_dim = Omega_con_hurdle.n_rows; di_con_hurdle.omega = &Omega_con_hurdle[0];

  std::vector<double> allfit_mod_hurdle(n_b, 0.0);
  std::vector<double> r_mod_hurdle(n_b);

  dinfo di_mod_hurdle;
  di_mod_hurdle.n = n_b; di_mod_hurdle.p = p_mod_hurdle;
  di_mod_hurdle.x = x_mod_hurdle.data(); di_mod_hurdle.y = r_mod_hurdle.data();
  di_mod_hurdle.basis_dim = Omega_mod_hurdle.n_rows; di_mod_hurdle.omega = &Omega_mod_hurdle[0];

  std::vector<double> ftemp_hurdle(n_b);
  std::vector<std::vector<tree::tree_cp> > node_pointers_con_hurdle(ntree_con_hurdle);
  std::vector<std::vector<tree::tree_cp> > node_pointers_mod_hurdle(ntree_mod_hurdle);
  for(size_t j=0; j<ntree_con_hurdle; ++j) {
    node_pointers_con_hurdle[j].resize(n_b);
    fit_basis(t_con_hurdle[j], xi_con_hurdle, di_con_hurdle, ftemp_hurdle.data(), node_pointers_con_hurdle[j], true, vanilla_hurdle);
  }
  for(size_t j=0; j<ntree_mod_hurdle; ++j) {
    node_pointers_mod_hurdle[j].resize(n_b);
    fit_basis(t_mod_hurdle[j], xi_mod_hurdle, di_mod_hurdle, ftemp_hurdle.data(), node_pointers_mod_hurdle[j], true, false);
  }

  std::vector<double> allfit_hurdle(n_b);
  for(size_t i=0; i<n_b; i++) {
    allfit_hurdle[i] = allfit_con_hurdle[i] + allfit_mod_hurdle[i];
  }

  //-------------------------------------------------------------------------//
  // 2. CONTINUOUS PART VARIABLE SETUP                                       //
  //-------------------------------------------------------------------------//
  size_t n_c = y_continuous.size();

  std::vector<double> x_con_continuous;
  for(NumericVector::iterator it=x_con_continuous_.begin(); it!=x_con_continuous_.end(); ++it) {
    x_con_continuous.push_back(*it);
  }
  size_t p_con_continuous = x_con_continuous.size()/n_c;

  xinfo xi_con_continuous;
  xi_con_continuous.resize(p_con_continuous);
  for(int i=0; i<p_con_continuous; ++i) {
    NumericVector tmp = x_con_info_continuous_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_con_continuous[i] = tmp2;
  }

  std::vector<double> x_mod_continuous;
  for(NumericVector::iterator it=x_mod_continuous_.begin(); it!=x_mod_continuous_.end(); ++it) {
    x_mod_continuous.push_back(*it);
  }
  size_t p_mod_continuous = x_mod_continuous.size()/n_c;

  xinfo xi_mod_continuous;
  xi_mod_continuous.resize(p_mod_continuous);
  for(int i=0; i<p_mod_continuous; ++i) {
    NumericVector tmp = x_mod_info_continuous_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
    xi_mod_continuous[i] = tmp2;
  }

  sinfo allys_continuous;
  for(size_t i=0; i<n_c; i++) {
    allys_continuous.sy += y_continuous[i];
    allys_continuous.sy2 += y_continuous[i]*y_continuous[i];
  }
  allys_continuous.n = n_c;
  double ybar_continuous = allys_continuous.sy/n_c;
  double shat_continuous = sqrt((allys_continuous.sy2 - n_c * ybar_continuous * ybar_continuous)/(n_c - 1));
  double sigma_continuous = shat_continuous;

  std::vector<tree> t_con_continuous(ntree_con_continuous);
  for(size_t i=0; i<ntree_con_continuous; i++) t_con_continuous[i].setm(zeros(Omega_con_continuous.n_rows));

  std::vector<tree> t_mod_continuous(ntree_mod_continuous);
  for(size_t i=0; i<ntree_mod_continuous; i++) t_mod_continuous[i].setm(zeros(Omega_mod_continuous.n_rows));

  pinfo pi_con_continuous;
  pi_con_continuous.pbd = 1.0; pi_con_continuous.pb = 0.5;
  pi_con_continuous.alpha = con_alpha_continuous; pi_con_continuous.beta = con_beta_continuous;
  pi_con_continuous.sigma = shat_continuous;
  pi_con_continuous.mu0 = zeros(Omega_con_continuous.n_rows);
  pi_con_continuous.Sigma0 = Sigma0_con_continuous;
  pi_con_continuous.Prec0 = pi_con_continuous.Sigma0.i();
  pi_con_continuous.logdetSigma0 = log(det(pi_con_continuous.Sigma0));
  pi_con_continuous.eta = 1.0; pi_con_continuous.gamma = 1.0;
  pi_con_continuous.scale_df = con_scale_df_continuous;
  pi_con_continuous.dart = false;

  pinfo pi_mod_continuous;
  pi_mod_continuous.pbd = 1.0; pi_mod_continuous.pb = 0.5;
  pi_mod_continuous.alpha = mod_alpha_continuous; pi_mod_continuous.beta = mod_beta_continuous;
  pi_mod_continuous.sigma = shat_continuous;
  pi_mod_continuous.mu0 = zeros(Omega_mod_continuous.n_rows);
  pi_mod_continuous.Sigma0 = Sigma0_mod_continuous;
  pi_mod_continuous.Prec0 = pi_mod_continuous.Sigma0.i();
  pi_mod_continuous.logdetSigma0 = log(det(pi_mod_continuous.Sigma0));
  pi_mod_continuous.eta = 1.0; pi_mod_continuous.gamma = 1.0;
  pi_mod_continuous.scale_df = mod_scale_df_continuous;
  pi_mod_continuous.dart = false;

  std::vector<double> allfit_con_continuous(n_c, 0.0);
  std::vector<double> r_con_continuous(n_c);

  dinfo di_con_continuous;
  di_con_continuous.n = n_c; di_con_continuous.p = p_con_continuous;
  di_con_continuous.x = x_con_continuous.data(); di_con_continuous.y = r_con_continuous.data();
  di_con_continuous.basis_dim = Omega_con_continuous.n_rows; di_con_continuous.omega = &Omega_con_continuous[0];

  std::vector<double> allfit_mod_continuous(n_c, 0.0);
  std::vector<double> r_mod_continuous(n_c);

  dinfo di_mod_continuous;
  di_mod_continuous.n = n_c; di_mod_continuous.p = p_mod_continuous;
  di_mod_continuous.x = x_mod_continuous.data(); di_mod_continuous.y = r_mod_continuous.data();
  di_mod_continuous.basis_dim = Omega_mod_continuous.n_rows; di_mod_continuous.omega = &Omega_mod_continuous[0];

  std::vector<double> ftemp_continuous(n_c);
  std::vector<std::vector<tree::tree_cp> > node_pointers_con_continuous(ntree_con_continuous);
  std::vector<std::vector<tree::tree_cp> > node_pointers_mod_continuous(ntree_mod_continuous);
  for(size_t j=0; j<ntree_con_continuous; ++j) {
    node_pointers_con_continuous[j].resize(n_c);
    fit_basis(t_con_continuous[j], xi_con_continuous, di_con_continuous, ftemp_continuous.data(), node_pointers_con_continuous[j], true, vanilla_continuous);
  }
  for(size_t j=0; j<ntree_mod_continuous; ++j) {
    node_pointers_mod_continuous[j].resize(n_c);
    fit_basis(t_mod_continuous[j], xi_mod_continuous, di_mod_continuous, ftemp_continuous.data(), node_pointers_mod_continuous[j], true, false);
  }

  std::vector<double> allfit_continuous(n_c);
  for(size_t i=0; i<n_c; i++) {
    allfit_continuous[i] = allfit_con_continuous[i] + allfit_mod_continuous[i];
  }

  //-------------------------------------------------------------------------//
  // 2B. CONTINUOUS OUT-OF-SAMPLE ESTIMATION SETUP                           //
  //-------------------------------------------------------------------------//
  std::vector<double> x_con_continuous_est;
  for(NumericVector::iterator it=x_con_continuous_est_.begin(); it!=x_con_continuous_est_.end(); ++it) {
    x_con_continuous_est.push_back(*it);
  }

  std::vector<double> x_mod_continuous_est;
  for(NumericVector::iterator it=x_mod_continuous_est_.begin(); it!=x_mod_continuous_est_.end(); ++it) {
    x_mod_continuous_est.push_back(*it);
  }

  dinfo di_con_continuous_est;
  di_con_continuous_est.n = n_b; di_con_continuous_est.p = p_con_continuous;
  di_con_continuous_est.x = &x_con_continuous_est[0]; di_con_continuous_est.y = 0;
  di_con_continuous_est.basis_dim = Omega_con_continuous_est.n_rows; di_con_continuous_est.omega = &Omega_con_continuous_est[0];

  dinfo di_mod_continuous_est;
  di_mod_continuous_est.n = n_b; di_mod_continuous_est.p = p_mod_continuous;
  di_mod_continuous_est.x = &x_mod_continuous_est[0]; di_mod_continuous_est.y = 0;
  di_mod_continuous_est.basis_dim = Omega_mod_continuous_est.n_rows; di_mod_continuous_est.omega = &Omega_mod_continuous_est[0];

  //-------------------------------------------------------------------------//
  // 3. STORAGE FOR OUTPUT DRAWS                                             //
  //-------------------------------------------------------------------------//
  NumericMatrix yhat_hurdle_post(nd, n_b);
  NumericMatrix m_hurdle_post(nd, n_b);
  NumericMatrix b_hurdle_post(nd, n_b);
  NumericVector eta_con_hurdle_post(nd);
  NumericVector eta_mod_hurdle_post(nd);

  NumericMatrix yhat_continuous_post(nd, n_c);
  NumericMatrix m_continuous_post(nd, n_c);
  NumericMatrix b_continuous_post(nd, n_c);
  NumericVector eta_con_continuous_post(nd);
  NumericVector eta_mod_continuous_post(nd);
  NumericVector sigma_continuous_post(nd);

  // Full-sample continuous estimation outputs (n_b rows)
  NumericMatrix m_continuous_est_post(nd, n_b);
  NumericMatrix b_continuous_est_post(nd, n_b);

  //-------------------------------------------------------------------------//
  // 4. MCMC LOOP                                                            //
  //-------------------------------------------------------------------------//
  Rcout << "\nBeginning Joint Hurdle-Continuous BCF (Path A) MCMC:\n";
  time_t tp;
  int time1 = time(&tp);

  size_t save_ctr = 0;
  for(size_t i = 0; i < (nd*thin + burn); i++) {

    Rcpp::checkUserInterrupt();
    if(i % status_interval == 0) {
      Rcout << "iteration: " << i << " continuous sigma: " << sigma_continuous << endl;
    }

    //-----------------------------------------------------------------------//
    //  HURDLE PART (Probit Link BCF)                                        //
    //-----------------------------------------------------------------------//
    // Update latent variables
    for(size_t k=0; k<n_b; k++) {
      if(y_hurdle[k] == 1) {
        y_lat_hurdle[k] = rtnormlo(allfit_hurdle[k], 1.0, 0.0);
      } else {
        y_lat_hurdle[k] = -rtnormlo(-allfit_hurdle[k], 1.0, 0.0);
      }
    }

    // Prognostic trees
    for(size_t j=0; j<ntree_con_hurdle; j++) {
      fit_basis(t_con_hurdle[j], xi_con_hurdle, di_con_hurdle, ftemp_hurdle.data(), node_pointers_con_hurdle[j], false, vanilla_hurdle);
      for(size_t k=0; k<n_b; k++) {
        allfit_hurdle[k] -= pi_con_hurdle.eta*ftemp_hurdle[k];
        allfit_con_hurdle[k] -= pi_con_hurdle.eta*ftemp_hurdle[k];
        r_con_hurdle[k] = (y_lat_hurdle[k] - allfit_hurdle[k])/pi_con_hurdle.eta;
      }
      double aa = bd_basis(t_con_hurdle[j], xi_con_hurdle, di_con_hurdle, pi_con_hurdle, gen, node_pointers_con_hurdle[j]);
      drmu_basis(t_con_hurdle[j], xi_con_hurdle, di_con_hurdle, pi_con_hurdle, gen);
      fit_basis(t_con_hurdle[j], xi_con_hurdle, di_con_hurdle, ftemp_hurdle.data(), node_pointers_con_hurdle[j], false, vanilla_hurdle);
      for(size_t k=0; k<n_b; k++) {
        allfit_hurdle[k] += pi_con_hurdle.eta*ftemp_hurdle[k];
        allfit_con_hurdle[k] += pi_con_hurdle.eta*ftemp_hurdle[k];
      }
    }

    // Moderating trees
    for(size_t j=0; j<ntree_mod_hurdle; j++) {
      fit_basis(t_mod_hurdle[j], xi_mod_hurdle, di_mod_hurdle, ftemp_hurdle.data(), node_pointers_mod_hurdle[j], false, false);
      for(size_t k=0; k<n_b; k++) {
        allfit_hurdle[k] -= pi_mod_hurdle.eta*ftemp_hurdle[k];
        allfit_mod_hurdle[k] -= pi_mod_hurdle.eta*ftemp_hurdle[k];
        r_mod_hurdle[k] = (y_lat_hurdle[k] - allfit_hurdle[k])/pi_mod_hurdle.eta;
      }
      double aa = bd_basis(t_mod_hurdle[j], xi_mod_hurdle, di_mod_hurdle, pi_mod_hurdle, gen, node_pointers_mod_hurdle[j]);
      drmu_basis(t_mod_hurdle[j], xi_mod_hurdle, di_mod_hurdle, pi_mod_hurdle, gen);
      fit_basis(t_mod_hurdle[j], xi_mod_hurdle, di_mod_hurdle, ftemp_hurdle.data(), node_pointers_mod_hurdle[j], false, false);
      for(size_t k=0; k<n_b; k++) {
        allfit_hurdle[k] += pi_mod_hurdle.eta*ftemp_hurdle[k];
        allfit_mod_hurdle[k] += pi_mod_hurdle.eta*ftemp_hurdle[k];
      }
    }

    // Scale updates for Hurdle part
    if(use_con_scale_hurdle) {
      for(size_t k=0; k<n_b; k++) ftemp_hurdle[k] = y_lat_hurdle[k] - allfit_mod_hurdle[k];
      double eta_old = pi_con_hurdle.eta;
      update_scale(ftemp_hurdle.data(), allfit_con_hurdle.data(), n_b, 1.0, pi_con_hurdle, gen);
      for(size_t k=0; k<n_b; k++) {
        allfit_hurdle[k] -= allfit_con_hurdle[k];
        allfit_con_hurdle[k] = allfit_con_hurdle[k] * pi_con_hurdle.eta / eta_old;
        allfit_hurdle[k] += allfit_con_hurdle[k];
      }
      pi_con_hurdle.sigma = 1.0/fabs(pi_con_hurdle.eta);
    }

    if(use_mod_scale_hurdle) {
      for(size_t k=0; k<n_b; k++) ftemp_hurdle[k] = y_lat_hurdle[k] - allfit_con_hurdle[k];
      double eta_old = pi_mod_hurdle.eta;
      update_scale(ftemp_hurdle.data(), allfit_mod_hurdle.data(), n_b, 1.0, pi_mod_hurdle, gen);
      for(size_t k=0; k<n_b; k++) {
        allfit_hurdle[k] -= allfit_mod_hurdle[k];
        allfit_mod_hurdle[k] = allfit_mod_hurdle[k] * pi_mod_hurdle.eta / eta_old;
        allfit_hurdle[k] += allfit_mod_hurdle[k];
      }
      pi_mod_hurdle.sigma = 1.0/fabs(pi_mod_hurdle.eta);
    }

    //-----------------------------------------------------------------------//
    //  CONTINUOUS PART (Log-Normal BCF)                                     //
    //-----------------------------------------------------------------------//
    // Prognostic trees
    for(size_t j=0; j<ntree_con_continuous; j++) {
      fit_basis(t_con_continuous[j], xi_con_continuous, di_con_continuous, ftemp_continuous.data(), node_pointers_con_continuous[j], false, vanilla_continuous);
      for(size_t k=0; k<n_c; k++) {
        allfit_continuous[k] -= pi_con_continuous.eta*ftemp_continuous[k];
        allfit_con_continuous[k] -= pi_con_continuous.eta*ftemp_continuous[k];
        r_con_continuous[k] = (y_continuous[k] - allfit_continuous[k])/pi_con_continuous.eta;
      }
      double aa = bd_basis(t_con_continuous[j], xi_con_continuous, di_con_continuous, pi_con_continuous, gen, node_pointers_con_continuous[j]);
      drmu_basis(t_con_continuous[j], xi_con_continuous, di_con_continuous, pi_con_continuous, gen);
      fit_basis(t_con_continuous[j], xi_con_continuous, di_con_continuous, ftemp_continuous.data(), node_pointers_con_continuous[j], false, vanilla_continuous);
      for(size_t k=0; k<n_c; k++) {
        allfit_continuous[k] += pi_con_continuous.eta*ftemp_continuous[k];
        allfit_con_continuous[k] += pi_con_continuous.eta*ftemp_continuous[k];
      }
    }

    // Moderating trees
    for(size_t j=0; j<ntree_mod_continuous; j++) {
      fit_basis(t_mod_continuous[j], xi_mod_continuous, di_mod_continuous, ftemp_continuous.data(), node_pointers_mod_continuous[j], false, false);
      for(size_t k=0; k<n_c; k++) {
        allfit_continuous[k] -= pi_mod_continuous.eta*ftemp_continuous[k];
        allfit_mod_continuous[k] -= pi_mod_continuous.eta*ftemp_continuous[k];
        r_mod_continuous[k] = (y_continuous[k] - allfit_continuous[k])/pi_mod_continuous.eta;
      }
      double aa = bd_basis(t_mod_continuous[j], xi_mod_continuous, di_mod_continuous, pi_mod_continuous, gen, node_pointers_mod_continuous[j]);
      drmu_basis(t_mod_continuous[j], xi_mod_continuous, di_mod_continuous, pi_mod_continuous, gen);
      fit_basis(t_mod_continuous[j], xi_mod_continuous, di_mod_continuous, ftemp_continuous.data(), node_pointers_mod_continuous[j], false, false);
      for(size_t k=0; k<n_c; k++) {
        allfit_continuous[k] += pi_mod_continuous.eta*ftemp_continuous[k];
        allfit_mod_continuous[k] += pi_mod_continuous.eta*ftemp_continuous[k];
      }
    }

    // Scale updates for Continuous part
    if(use_con_scale_continuous) {
      for(size_t k=0; k<n_c; k++) ftemp_continuous[k] = y_continuous[k] - allfit_mod_continuous[k];
      double eta_old = pi_con_continuous.eta;
      update_scale(ftemp_continuous.data(), allfit_con_continuous.data(), n_c, sigma_continuous, pi_con_continuous, gen);
      for(size_t k=0; k<n_c; k++) {
        allfit_continuous[k] -= allfit_con_continuous[k];
        allfit_con_continuous[k] = allfit_con_continuous[k] * pi_con_continuous.eta / eta_old;
        allfit_continuous[k] += allfit_con_continuous[k];
      }
      pi_con_continuous.sigma = sigma_continuous/fabs(pi_con_continuous.eta);
    }

    if(use_mod_scale_continuous) {
      for(size_t k=0; k<n_c; k++) ftemp_continuous[k] = y_continuous[k] - allfit_con_continuous[k];
      double eta_old = pi_mod_continuous.eta;
      update_scale(ftemp_continuous.data(), allfit_mod_continuous.data(), n_c, sigma_continuous, pi_mod_continuous, gen);
      for(size_t k=0; k<n_c; k++) {
        allfit_continuous[k] -= allfit_mod_continuous[k];
        allfit_mod_continuous[k] = allfit_mod_continuous[k] * pi_mod_continuous.eta / eta_old;
        allfit_continuous[k] += allfit_mod_continuous[k];
      }
      pi_mod_continuous.sigma = sigma_continuous/fabs(pi_mod_continuous.eta);
    }

    // Update continuous sigma
    double rss_continuous = 0.0;
    double restemp = 0.0;
    for(size_t k=0; k<n_c; k++) {
      restemp = y_continuous[k] - allfit_continuous[k];
      rss_continuous += restemp*restemp;
    }
    sigma_continuous = sqrt((nu*lambda + rss_continuous)/gen.chi_square(nu+n_c));
    pi_con_continuous.sigma = sigma_continuous/fabs(pi_con_continuous.eta);
    pi_mod_continuous.sigma = sigma_continuous/fabs(pi_mod_continuous.eta);

    //-----------------------------------------------------------------------//
    //  SAVE RESULTS                                                         //
    //-----------------------------------------------------------------------//
    if((i >= burn) && (i % thin == 0)) {

      // Hurdle draws
      eta_con_hurdle_post(save_ctr) = pi_con_hurdle.eta;
      eta_mod_hurdle_post(save_ctr) = pi_mod_hurdle.eta;
      for(size_t k=0; k<n_b; k++) {
        m_hurdle_post(save_ctr, k) = allfit_con_hurdle[k];
        b_hurdle_post(save_ctr, k) = pi_mod_hurdle.eta * fit_i_basis(k, t_mod_hurdle, xi_mod_hurdle, di_mod_hurdle, true);
        yhat_hurdle_post(save_ctr, k) = allfit_hurdle[k];
      }

      // Continuous draws (active subset)
      eta_con_continuous_post(save_ctr) = pi_con_continuous.eta;
      eta_mod_continuous_post(save_ctr) = pi_mod_continuous.eta;
      sigma_continuous_post(save_ctr) = sigma_continuous;
      for(size_t k=0; k<n_c; k++) {
        m_continuous_post(save_ctr, k) = allfit_con_continuous[k];
        b_continuous_post(save_ctr, k) = pi_mod_continuous.eta * fit_i_basis(k, t_mod_continuous, xi_mod_continuous, di_mod_continuous, true);
        yhat_continuous_post(save_ctr, k) = allfit_continuous[k];
      }

      // Out-of-sample prediction of continuous part for ALL observations (n_b rows)
      for(size_t k=0; k<n_b; k++) {
        m_continuous_est_post(save_ctr, k) = pi_con_continuous.eta * fit_i_basis(k, t_con_continuous, xi_con_continuous, di_con_continuous_est, vanilla_continuous);
        b_continuous_est_post(save_ctr, k) = pi_mod_continuous.eta * fit_i_basis(k, t_mod_continuous, xi_mod_continuous, di_mod_continuous_est, true);
      }

      save_ctr += 1;
    }
  }

  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;

  // std::vector handles dynamic array memory cleanup automatically and safely

  return(List::create(
      _["yhat_hurdle_post"] = yhat_hurdle_post,
      _["m_hurdle_post"] = m_hurdle_post,
      _["b_hurdle_post"] = b_hurdle_post,
      _["eta_con_hurdle"] = eta_con_hurdle_post,
      _["eta_mod_hurdle"] = eta_mod_hurdle_post,

      _["yhat_continuous_post"] = yhat_continuous_post,
      _["m_continuous_post"] = m_continuous_post,
      _["b_continuous_post"] = b_continuous_post,
      _["eta_con_continuous"] = eta_con_continuous_post,
      _["eta_mod_continuous"] = eta_mod_continuous_post,
      _["sigma_continuous"] = sigma_continuous_post,

      _["m_continuous_est_post"] = m_continuous_est_post,
      _["b_continuous_est_post"] = b_continuous_est_post
  ));
}
