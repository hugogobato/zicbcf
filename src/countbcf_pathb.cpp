#include "arma_config.h"
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "tree_samples.h"

using namespace Rcpp;

// Latent variable sampler for compound Poisson-Gamma (Tweedie with p=1.5)
// y_i: observed continuous semicontinuous outcome
// mu_i: current mean prediction (unused, kept for signature compatibility)
// phi: current dispersion parameter
// gen: random number generator
int draw_latent_N(double y_i, double mu_i, double phi, RNG& gen) {
  if (y_i == 0) return 0;
  double d = 4.0 * y_i / (phi * phi);
  if (d > 500000.0) d = 500000.0; // safety guard to prevent excessively large loops
  if (d <= 0 || std::isnan(d) || std::isinf(d)) return 1; // fall back
  
  double log_d = log(d);
  std::vector<double> log_probs;
  double max_log_prob = log_d;
  log_probs.push_back(log_d);
  
  double current_log_prob = log_d;
  int n = 1;
  while (true) {
    n++;
    current_log_prob += log_d - log(n) - log(n - 1);
    log_probs.push_back(current_log_prob);
    if (current_log_prob > max_log_prob) {
      max_log_prob = current_log_prob;
    }
    // Stop condition: if log_prob is much smaller than max_log_prob
    if (current_log_prob < max_log_prob - 25.0) {
      break;
    }
    if (n > 2000) break; // safety guard
  }
  
  // Exponentiate and normalize
  std::vector<double> probs;
  double sum_probs = 0.0;
  for (size_t i = 0; i < log_probs.size(); ++i) {
    double p = exp(log_probs[i] - max_log_prob);
    probs.push_back(p);
    sum_probs += p;
  }
  
  double u = gen.uniform() * sum_probs;
  double running_sum = 0.0;
  for (size_t i = 0; i < probs.size(); ++i) {
    running_sum += probs[i];
    if (u <= running_sum) {
      return i + 1;
    }
  }
  return probs.size();
}

// Zero-Inflated & Tweedie Bayesian Causal Forest (countbcf_pathb).
//
// [[Rcpp::export]]
List countbcf_pathb(arma::vec y_,
                    arma::vec offset_,
                    List bart_specs,
                    List bart_designs,
                    arma::mat random_des,
                    arma::mat random_var, arma::mat random_var_ix,
                    double random_var_df, arma::vec randeff_scales,
                    int burn, int nd, int thin,
                    int count_model,
                    double lambda, double nu,
                    double kappa_a, double kappa_b,
                    double leaf_c, double leaf_d,
                    double z_c, double z_d,
                    double kappa_prop_sd = 0.2,
                    bool return_trees = true,
                    bool save_trees = false,
                    bool est_mod_fits = false, bool est_con_fits = false,
                    bool prior_sample = false,
                    int status_interval = 100,
                    NumericVector lower_bd = NumericVector::create(0.0),
                    NumericVector upper_bd = NumericVector::create(0.0),
                    bool probit = false,
                    bool text_trace = true,
                    bool R_trace = false)
{

  bool randeff = true;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  }
  if(randeff) Rcout << "Using random effects." << std::endl;

  RNGScope scope;
  RNG gen;

  //-------------------------------------------------------------------------//
  // Read / format 'y'                                                       //
  //-------------------------------------------------------------------------//
  Rcout << "Reading in y\n\n";

  std::vector<double> y;
  std::vector<double> offset;

  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;

  int n_y0 = -1;

  for(NumericVector::iterator it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it;
    allys.sy2 += (*it)*(*it);
    if (*it == 0) n_y0 += 1;
  }

  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n;
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1));
  if(probit) shat = 1.0;

  double sigma = shat;

  for(NumericVector::iterator it = offset_.begin(); it != offset_.end(); ++it){
    offset.push_back(*it);
  }

  double kappa = 1;
  double kappa_acpt_rate = 0;

  std::vector<int> z;
  std::vector<double> log_w;
  std::vector<double> log_w_denom;
  std::vector<double> log_phi;
  std::vector<double> log_xi;
  std::vector<double> u_vec = y;

  // Latent variables for Tweedie model (count_model == 5)
  std::vector<int> N_latent(n, 0);

  for(std::size_t it = 0; it < y.size(); ++it){
    z.push_back(1);
    log_w.push_back(log(0.5));
    log_w_denom.push_back(log(2));
    log_phi.push_back(0);
    log_xi.push_back(0);
  }

  //-------------------------------------------------------------------------//
  // Read, format design info                                                //
  //-------------------------------------------------------------------------//

  Rcout << "Setting up designs\n\n";

  size_t num_designs = bart_designs.size();

  std::vector<std::vector<double> > x(num_designs);
  std::vector<std::vector<int> > groups(num_designs);
  std::vector<bool> group(num_designs);
  std::vector<xinfo> x_info(num_designs);
  std::vector<arma::mat> Omega(num_designs);
  std::vector<size_t> covariate_dim(num_designs);

  for(size_t i=0; i<num_designs; i++) {

    Rcout << "design " << i << endl;

    List dtemp = bart_designs[i];

    bool g = dtemp["group"];
    group[i] = g;

    Rcout << "group setting: " << g << endl;

    IntegerVector gt_ = dtemp["groups"];
    for(IntegerVector::iterator it=gt_.begin(); it!= gt_.end(); ++it) {
      groups[i].push_back(*it);
    }

    NumericVector xt_ = dtemp["X"];
    for(NumericVector::iterator it=xt_.begin(); it!= xt_.end(); ++it) {
      x[i].push_back(*it);
    }
    size_t p = x[i].size()/n;
    covariate_dim[i] = p;

    Rcout << "Instantiated covariate matrix " << i+1 << " with " << p << " columns" << endl;

    xinfo xi;
    xi.resize(p);
    List x_info_list = dtemp["info"];
    for(int j=0; j<p; ++j) {
      NumericVector tmp = x_info_list[j];
      std::vector<double> tmp2;
      for(size_t s=0; s<tmp.size(); ++s) {
        tmp2.push_back(tmp[s]);
      }
      xi[j] = tmp2;
    }
    x_info[i] = xi;

    Omega[i] = as<arma::mat>(dtemp["Omega"]);
  }

  //-------------------------------------------------------------------------//
  // Set up forests                                                          //
  //-------------------------------------------------------------------------//

  Rcout << "Setting up forests\n\n";

  size_t num_forests = bart_specs.size();
  std::vector<std::vector<tree> > trees(num_forests);
  std::vector<pinfo> prior_info(num_forests);
  std::vector<std::vector<double> > allfits(num_forests);
  std::vector<int> function_group(num_forests);

  std::vector<double> r_tree(n);
  std::fill(r_tree.begin(), r_tree.end(), 0.0);

  std::vector<double> log_mean(n);
  std::fill(log_mean.begin(), log_mean.end(), 0.0);

  Rcout << "Number of forests: " << num_forests <<  "\n\n" << endl;

  std::vector<dinfo> di(num_forests);
  std::vector<std::vector<std::vector<tree::tree_cp> > > node_pointers(num_forests);
  std::vector<double> sample_eta(num_forests);

  std::vector<tree_samples> final_tree_trace(num_forests);

  double* ftemp  = new double[n];

  size_t ntree_control = 250;
  if (bart_specs.size() > 0) {
    List spec0 = bart_specs[0];
    ntree_control = spec0["ntree"];
  }

  for(size_t i=0; i<num_forests; ++i) {

    List spec = bart_specs[i];

    double es = spec["sample_eta"];
    sample_eta[i] = es;

    int desi = spec["design_index"];
    size_t ntree = spec["ntree"];
    trees[i].resize(ntree);
    prior_info[i].vanilla = spec["vanilla"];

    int fg = spec["function_group"];
    function_group[i] = fg;

    Rcout << "forest " << i << " function_group=" << fg << " vanilla=" << prior_info[i].vanilla << endl;

    for(size_t j=0; j<ntree; ++j) trees[i][j].setm(zeros(Omega[desi].n_rows));

    prior_info[i].pbd = 1.0;
    prior_info[i].pb = .5;

    if ((count_model == 3) || (count_model == 4)){
      if (fg == 0){
        prior_info[i].c = leaf_c;
        prior_info[i].d = leaf_d;
      } else {
        prior_info[i].c = z_c;
        prior_info[i].d = z_d;
      }
    } else {
      prior_info[i].c = leaf_c;
      prior_info[i].d = leaf_d;
    }

    prior_info[i].alpha = spec["alpha"];
    prior_info[i].beta  = spec["beta"];
    prior_info[i].sigma = shat;

    prior_info[i].mu0 = as<arma::vec>(spec["mu0"]);

    prior_info[i].Sigma0 = as<arma::mat>(spec["Sigma0"]);
    prior_info[i].Prec0 = prior_info[i].Sigma0.i();
    prior_info[i].logdetSigma0 = log(det(prior_info[i].Sigma0));
    prior_info[i].eta = 1.0 / ntree_control;
    prior_info[i].gamma = 1;
    prior_info[i].scale_df = spec["scale_df"];

    dinfo dtemp;
    dtemp.n=n;
    dtemp.p = covariate_dim[desi];
    dtemp.x = &(x[desi])[0];
    dtemp.y = &r_tree[0];
    dtemp.u_i = &u_vec[0];
    dtemp.offset = &offset[0];
    dtemp.basis_dim = Omega[desi].n_rows;
    dtemp.omega = &(Omega[desi])[0];
    dtemp.groups = &(groups[desi])[0];
    dtemp.group = group[desi];

    node_pointers[i].resize(ntree);
    allfits[i].resize(n);
    std::fill(allfits[i].begin(), allfits[i].end(), 0.0);

    for(size_t j=0; j<ntree; ++j) {
      node_pointers[i][j].resize(n);
      fit_loglinear(trees[i][j], x_info[desi], dtemp, ftemp, node_pointers[i][j], true, prior_info[i].vanilla);
      for(size_t k=0; k<n; ++k) allfits[i][k] += ftemp[k];
    }

    prior_info[i].dart = spec["dart"];
    if(prior_info[i].dart) prior_info[i].dart_alpha = 1.0;
    std::vector<double> vp_mod(covariate_dim[desi], 1.0/covariate_dim[desi]);
    prior_info[i].var_probs = vp_mod;

    di[i] = dtemp;

    if (return_trees) {
      tree_samples ts(ntree, di[i].p, nd, di[i].basis_dim, x_info[desi]);
      final_tree_trace[i] = ts;
    }
  }

  Rcout << "Done." << endl;

  //-------------------------------------------------------------------------//
  // Build per-group log f^{(g)} and log_mean                                //
  //-------------------------------------------------------------------------//

  std::vector<std::vector<double> > log_f_per_group(3);
  for(int g=0; g<3; ++g) {
    log_f_per_group[g].resize(n);
    std::fill(log_f_per_group[g].begin(), log_f_per_group[g].end(), 0.0);
  }

  for(size_t s=0; s<num_forests; ++s) {
    int g = function_group[s];
    for(size_t k=0; k<n; ++k) {
      log_f_per_group[g][k] += allfits[s][k];
    }
  }

  for(size_t k=0; k<n; ++k) {
    if (count_model == 5) {
      log_mean[k] = offset[k] + 2.0 * log_f_per_group[0][k];
    } else {
      log_mean[k] = offset[k] + log_f_per_group[0][k];
    }
  }

  //-------------------------------------------------------------------------//
  // setup random effects (not used in count models)                         //
  //-------------------------------------------------------------------------//
  size_t random_dim = random_des.n_cols;
  int nr = 1;
  if(randeff) nr = n;

  arma::vec r(nr);
  arma::vec Wtr(random_dim);

  arma::mat WtW = random_des.t()*random_des;
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

  arma::vec eta(random_var_ix.n_cols);
  eta.fill(1.0);

  for(size_t k=0; k<nr; ++k) {
    r(k) = y[k];
    for(size_t j=0; j<num_forests; ++j) {
      r(k) -= allfits[j][k];
    }
  }

  Wtr = random_des.t()*r;
  arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random, Wtr/(sigma*sigma));
  arma::vec allfit_random = random_des*gamma;
  if(!randeff) allfit_random.fill(0);

  //-------------------------------------------------------------------------//
  // initialize tree fits                                                    //
  //-------------------------------------------------------------------------//
  double* allfit = new double[n];

  for(size_t i=0;i<n;i++) {
    allfit[i] = 0;
    for(size_t j=0; j<num_forests; ++j) {
      allfit[i] += allfits[j][i];
    }
    if(randeff) allfit[i] += allfit_random[i];
  }

  // output storage
  NumericVector sigma_post(nd);
  NumericVector kappa_post(nd);

  std::vector<std::vector<List> > forest_trace_R(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    forest_trace_R[j].resize(nd);
    for(size_t i=0; i<nd; ++i) {
      List init_list(trees[j].size());
      forest_trace_R[j][i] = init_list;
    }
  }

  std::vector<NumericMatrix> forest_fits_post(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    NumericMatrix postfits(nd, n);
    forest_fits_post[j] = postfits;
  }

  NumericMatrix yhat_post(nd, n);
  NumericMatrix etas_post(nd, num_forests);

  arma::mat gamma_post(nd, gamma.n_elem);
  arma::mat random_sd_post(nd, random_var.n_elem);

  std::vector<arma::cube> post_coefs(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    arma::cube tt(di[j].basis_dim, di[j].n, nd);
    tt.fill(0);
    post_coefs[j] = tt;
  }

  //-------------------------------------------------------------------------//
  // MCMC                                                                    //
  //-------------------------------------------------------------------------//

  Rcout << "\nBeginning MCMC:\n";
  time_t tp;
  int time1 = time(&tp);
  size_t save_ctr = 0;

  for(size_t i = 0; i < (nd*thin+burn); i++) {

    if(prior_sample) {
      for(int k=0; k<n; k++) y[k] = gen.normal(allfit[k], sigma);
    }

    if(lower_bd.size()>1) {
      for(int k=0; k<n; k++) {
        if(lower_bd[k]!=-INFINITY) y[k] = rtnormlo(allfit[k], sigma, lower_bd[k]);
      }
    }

    if(upper_bd.size()>1) {
      for(int k=0; k<n; k++) {
        if(upper_bd[k]!=INFINITY) y[k] = -rtnormlo(-allfit[k], sigma, -upper_bd[k]);
      }
    }

    Rcpp::checkUserInterrupt();
    if(i%status_interval == 0) {
      Rcout << "iteration: " << i << endl;
      if ((count_model == 2) || (count_model == 4)){
        Rcout << "dispersion parameter acceptance rate: " << kappa_acpt_rate << endl;
      }
      if (count_model == 5) {
        double sum_log_mean = 0.0;
        for (size_t k = 0; k < n; ++k) sum_log_mean += log_mean[k];
        Rcout << "  kappa=" << kappa << " avg_log_mean=" << sum_log_mean / n << endl;
      }
    }

    //-----------------------------------------------------------------------//
    //  update latent N_i (Tweedie compound Poisson, p=1.5)                  //
    //-----------------------------------------------------------------------//
    if (count_model == 5) {
      for (size_t k = 0; k < n; ++k) {
        double theta_k = log_mean[k] / 2.0;
        N_latent[k] = draw_latent_N(y[k], theta_k, kappa, gen);
      }
    }

    //-----------------------------------------------------------------------//
    //  update kappa (NB dispersion or Tweedie dispersion)                   //
    //-----------------------------------------------------------------------//
    if ((count_model == 2) || (count_model == 4)){

      double kappa_star = exp(gen.normal(log(kappa), kappa_prop_sd));

      double log_a_num = ll_loglinear(y, log_mean, kappa_star, count_model, true, log_w, n_y0)
        + R::dbeta(kappa_star / (1 + kappa_star), kappa_a, kappa_b, true) - 2 * log1p(kappa_star)
        + log(kappa_star);
      double log_a_denom = ll_loglinear(y, log_mean, kappa, count_model, true, log_w, n_y0)
        + R::dbeta(kappa / (1 + kappa), kappa_a, kappa_b, true) - 2 * log1p(kappa)
        + log(kappa);

      if(log(gen.uniform()) < log_a_num - log_a_denom){
        kappa = kappa_star;
        kappa_acpt_rate = ((kappa_acpt_rate * i) + 1) / (i + 1);
      } else {
        kappa_acpt_rate = ((kappa_acpt_rate * i) + 0) / (i + 1);
      }
    } else if (count_model == 5) {
      double sum_N = 0.0;
      double sum_V = 0.0;
      for (size_t k = 0; k < n; ++k) {
        if (y[k] > 0) {
          sum_N += N_latent[k];
        }
        double theta_k = log_mean[k] / 2.0;
        if (theta_k < -20.0) theta_k = -20.0;
        if (theta_k > 20.0) theta_k = 20.0;
        
        sum_V += 2.0 * exp(theta_k);
        if (y[k] > 0) {
          sum_V += 2.0 * y[k] * exp(-theta_k);
        }
      }
      double a_post = 1.0 + 2.0 * sum_N;
      double b_post = 1.0 + sum_V;
      
      double g_draw = gen.gamma(a_post, 1.0 / b_post);
      kappa = 1.0 / g_draw;
      if (kappa < 1e-4) kappa = 1e-4; // safety guard to prevent division by zero or explosion
    }

    //-----------------------------------------------------------------------//
    //  update Z_i (zero-inflated models)                                    //
    //-----------------------------------------------------------------------//
    if ((count_model == 3) || (count_model == 4)){
      drz_loglinear(z, log_w, kappa, log_mean, n_y0, count_model, gen);
    }

    //-----------------------------------------------------------------------//
    //  update xi (negative binomial)                                        //
    //-----------------------------------------------------------------------//
    if ((count_model == 2) || (count_model == 4)){
      drxi_loglinear(log_xi, kappa, log_mean, y, z, gen);
    }

    //-----------------------------------------------------------------------//
    //  update phi (zero-inflated)                                           //
    //-----------------------------------------------------------------------//
    if ((count_model == 3) || (count_model == 4)){
      drphi_loglinear(log_phi, log_w_denom, gen);
    }

    //-----------------------------------------------------------------------//
    // update trees                                                          //
    //-----------------------------------------------------------------------//
    for(size_t s = 0; s < num_forests; ++s) {
      int g = function_group[s];
      bool is_count_group = (g == 0);

      for(size_t j = 0; j < trees[s].size(); ++j) {

        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        for (size_t k = 0; k < n; k++){

          if (ftemp[k] != ftemp[k]){
            stop("nan in ftemp");
          }

          // remove current tree fit
          allfit[k]                -= prior_info[s].eta * ftemp[k];
          allfits[s][k]            -= prior_info[s].eta * ftemp[k];
          log_f_per_group[g][k]    -= prior_info[s].eta * ftemp[k];
          if (is_count_group) {
            if (count_model == 5) {
              log_mean[k]          -= 2.0 * prior_info[s].eta * ftemp[k];
            } else {
              log_mean[k]          -= prior_info[s].eta * ftemp[k];
            }
          }

          double omega_k = di[s].omega[k];

          // sufficient statistic vectors for the loglinear leaf update
          if (is_count_group) {
            // count component (group 0)
            if (count_model == 5) {
              double theta_k_current = log_mean[k] / 2.0;
              if (theta_k_current < -20.0) theta_k_current = -20.0;
              if (theta_k_current > 20.0) theta_k_current = 20.0;
              double exp_theta = exp(theta_k_current);
              double exp_neg_theta = exp(-theta_k_current);
              
              u_vec[k]  = (2.0 * N_latent[k] + (2.0 * y[k] * exp_neg_theta / kappa)) * omega_k;
              r_tree[k] = (2.0 * exp_theta / kappa) * omega_k;
            } else {
              u_vec[k]  = z[k] * y[k] * omega_k;
              r_tree[k] = z[k] * exp(log_xi[k] + log_mean[k]) * omega_k;
            }
          } else if (g == 1) {
            // f0 (zero-inflation odds)
            u_vec[k]  = (1 - z[k]) * omega_k;
            r_tree[k] = exp(log_phi[k] + log_f_per_group[1][k]) * omega_k;
          } else {
            // g == 2: f1 (NOT zero-inflation odds)
            u_vec[k]  = z[k] * omega_k;
            r_tree[k] = exp(log_phi[k] + log_f_per_group[2][k]) * omega_k;
          }

          if (r_tree[k] != r_tree[k]){
            stop("NaN in resid");
          }
        }

        bd_loglinear(trees[s][j], x_info[s], di[s], prior_info[s], gen, node_pointers[s][j]);
        drmu_loglinear(trees[s][j], x_info[s], di[s], prior_info[s], gen);
        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        // add new tree fit back
        for(size_t k=0; k<n; k++) {
          allfit[k]                += prior_info[s].eta * ftemp[k];
          allfits[s][k]            += prior_info[s].eta * ftemp[k];
          log_f_per_group[g][k]    += prior_info[s].eta * ftemp[k];
          if (is_count_group) {
            if (count_model == 5) {
              log_mean[k]          += 2.0 * prior_info[s].eta * ftemp[k];
            } else {
              log_mean[k]          += prior_info[s].eta * ftemp[k];
            }
          }
        }
      }
    }

    //-----------------------------------------------------------------------//
    // update w(x_i) using log_f_per_group[1] and log_f_per_group[2]         //
    //-----------------------------------------------------------------------//
    if ((count_model == 3) || (count_model == 4)){
      for (size_t k = 0; k < n; ++k){
        log_w_denom[k] = logsumexp(log_f_per_group[1][k], log_f_per_group[2][k]);
        log_w[k]       = log_f_per_group[2][k] - log_w_denom[k];
      }
    }

    //-----------------------------------------------------------------------//
    // update random effects (not used for count models)                     //
    //-----------------------------------------------------------------------//
    if(randeff) {

      for(size_t k=0; k<n; ++k) {
        r(k) = y[k] - allfit[k] + allfit_random[k];
        allfit[k] -= allfit_random[k];
      }

      Wtr = random_des.t()*r;

      arma::mat adj = diagmat(random_var_ix*eta);
      arma::mat Phi = adj*WtW*adj/(sigma*sigma) + Sigma_inv_random;
      Phi = 0.5*(Phi + Phi.t());
      arma::vec m = adj*Wtr/(sigma*sigma);
      gamma = rmvnorm_post(m, Phi);

      arma::mat adj2 = diagmat(gamma)*random_var_ix;
      arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
      arma::vec m2 = adj2.t()*Wtr/(sigma*sigma);
      Phi2 = 0.5*(Phi2 + Phi2.t());
      eta = rmvnorm_post(m2, Phi2);

      arma::vec ssqs   = random_var_ix.t()*(gamma % gamma);
      arma::rowvec counts = sum(random_var_ix, 0);
      for(size_t ii=0; ii<random_var_ix.n_cols; ++ii) {
        random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*2.0/(random_var_df/randeff_scales(ii)*randeff_scales(ii) + ssqs(ii));
      }
      Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

      allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;

      for(size_t k=0; k<n; ++k) {
        allfit[k] = allfit_random(k);
        for(size_t s=0; s<num_forests; ++s) {
          allfit[k] += allfits[s][k];
        }
      }
    }

    //-----------------------------------------------------------------------//
    // save results                                                          //
    //-----------------------------------------------------------------------//
    if( ((i>=burn) & (i % thin==0))) {

      gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
      random_sd_post.row(save_ctr) = (sqrt( eta % eta % random_var)).t();

      sigma_post(save_ctr) = sigma;
      kappa_post(save_ctr) = kappa;

      for(size_t k=0; k<n; k++) {
        yhat_post(save_ctr, k) = log_mean[k];
      }

      for(size_t s=0; s<num_forests; ++s) {
        etas_post(save_ctr,s) = prior_info[s].eta;
        for(size_t k=0; k<n; k++) {
          forest_fits_post[s](save_ctr, k) = allfits[s][k];
        }
        for(size_t j=0; j< trees[s].size(); ++j) {
          post_coefs[s].slice(save_ctr) += prior_info[s].eta*coef_basis(trees[s][j], x_info[s], di[s]);
          if (return_trees) {
            final_tree_trace[s].t[save_ctr][j] = trees[s][j];
            final_tree_trace[s].t[save_ctr][j].compress();
            final_tree_trace[s].t[save_ctr][j].scale(prior_info[s].eta);
          }
        }
      }

      save_ctr += 1;
    }
  }

  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;

  delete[] allfit;
  delete[] ftemp;

  return(List::create(_["yhat_post"] = yhat_post,
                      _["forest_fits"] = forest_fits_post,
                      _["function_group"] = function_group,
                      _["coefs"] = post_coefs,
                      _["etas"] = etas_post,
                      _["sigma"] = sigma_post,
                      _["kappa"] = kappa_post,
                      _["kappa_acceptance"] = kappa_acpt_rate,
                      _["gamma"] = gamma_post,
                      _["random_sd_post"] = random_sd_post,
                      _["tree_trace"] = final_tree_trace,
                      _["y_last"] = y
  ));
}
