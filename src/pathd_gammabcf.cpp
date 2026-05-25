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

// Log-likelihood of the Gamma BCF model (where theta = -log lambda)
double ll_gamma_bcf(const std::vector<double>& y, const std::vector<double>& theta, double kappa_c) {
  double ll = 0.0;
  size_t n = y.size();
  ll += n * kappa_c * log(kappa_c) - n * R::lgammafn(kappa_c);
  for (size_t i = 0; i < n; ++i) {
    double exp_theta = exp(theta[i]);
    if (exp_theta > 1e12) exp_theta = 1e12;
    if (exp_theta < 1e-12) exp_theta = 1e-12;
    ll += kappa_c * theta[i] + (kappa_c - 1.0) * log(y[i]) - kappa_c * y[i] * exp_theta;
  }
  return ll;
}

// [[Rcpp::export]]
List pathd_gammabcf(arma::vec y_,
                    arma::vec offset_,
                    List bart_specs,
                    List bart_designs,
                    arma::mat random_des,
                    arma::mat random_var, arma::mat random_var_ix,
                    double random_var_df, arma::vec randeff_scales,
                    int burn, int nd, int thin,
                    double kappa_a, double kappa_b,
                    double leaf_c, double leaf_d,
                    double kappa_prop_sd = 0.2,
                    bool return_trees = true,
                    bool save_trees = false,
                    bool est_mod_fits = false, bool est_con_fits = false,
                    bool prior_sample = false,
                    int status_interval = 100,
                    bool text_trace = true)
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

  for(NumericVector::iterator it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it;
    allys.sy2 += (*it)*(*it);
  }

  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n;
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1));

  double sigma = shat;

  for(NumericVector::iterator it = offset_.begin(); it != offset_.end(); ++it){
    offset.push_back(*it);
  }

  double kappa = 1.0; // Shape parameter kappa_c
  if (n > 1) {
    double var_y = (allys.sy2 - n * ybar * ybar) / (n - 1);
    if (var_y > 1e-8) {
      kappa = (ybar * ybar) / var_y;
      if (kappa < 0.1) kappa = 0.1;
      if (kappa > 50.0) kappa = 10.0;
    }
  }
  double kappa_acpt_rate = 0;

  std::vector<double> u_vec = y;
  std::vector<double> r_tree(n);

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

  std::vector<double> log_mean(n);
  std::fill(log_mean.begin(), log_mean.end(), 0.0);

  Rcout << "Number of forests: " << num_forests <<  "\n\n" << endl;

  std::vector<dinfo> di(num_forests);
  std::vector<std::vector<std::vector<tree::tree_cp> > > node_pointers(num_forests);
  std::vector<double> sample_eta(num_forests);

  std::vector<tree_samples> final_tree_trace(num_forests);

  double* ftemp  = new double[n];

  for(size_t i=0; i<num_forests; ++i) {
    List spec = bart_specs[i];

    double es = spec["sample_eta"];
    sample_eta[i] = es;

    int desi = spec["design_index"];
    size_t ntree = spec["ntree"];
    trees[i].resize(ntree);
    prior_info[i].vanilla = spec["vanilla"];

    for(size_t j=0; j<ntree; ++j) trees[i][j].setm(zeros(Omega[desi].n_rows));

    prior_info[i].pbd = 1.0;
    prior_info[i].pb = .5;

    // Per-spec GIG leaf prior (mu vs tau forests need different concentration);
    // fall back to global leaf_c / leaf_d when the spec omits them.
    if (spec.containsElementNamed("leaf_c")) {
      prior_info[i].c = as<double>(spec["leaf_c"]);
    } else {
      prior_info[i].c = leaf_c;
    }
    if (spec.containsElementNamed("leaf_d")) {
      prior_info[i].d = as<double>(spec["leaf_d"]);
    } else {
      prior_info[i].d = leaf_d;
    }

    prior_info[i].alpha = spec["alpha"];
    prior_info[i].beta  = spec["beta"];
    prior_info[i].sigma = shat;

    prior_info[i].mu0 = as<arma::vec>(spec["mu0"]);

    prior_info[i].Sigma0 = as<arma::mat>(spec["Sigma0"]);
    prior_info[i].Prec0 = prior_info[i].Sigma0.i();
    prior_info[i].logdetSigma0 = log(det(prior_info[i].Sigma0));
    prior_info[i].eta = 1;
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
  // Build initial log_mean (theta = sum of trees + offset)                  //
  //-------------------------------------------------------------------------//
  for(size_t k=0; k<n; ++k) {
    double sum_fits = 0.0;
    for(size_t s=0; s<num_forests; ++s) {
      sum_fits += allfits[s][k];
    }
    log_mean[k] = offset[k] + sum_fits;
  }

  //-------------------------------------------------------------------------//
  // setup random effects (not used in standard path but initialized)        //
  //-------------------------------------------------------------------------//
  size_t random_dim = random_des.n_cols;
  int nr = 1;
  if(randeff) nr = n;

  arma::vec r(nr);
  arma::vec Wtr(random_dim);

  arma::mat WtW = random_des.t()*random_des;
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

  arma::vec eta_re(random_var_ix.n_cols);
  eta_re.fill(1.0);

  arma::vec gamma = zeros(random_dim);
  arma::vec allfit_random = zeros(n);

  double* allfit = new double[n];
  for(size_t i=0;i<n;i++) {
    allfit[i] = log_mean[i];
  }

  // output storage
  NumericVector kappa_post(nd);
  std::vector<NumericMatrix> forest_fits_post(num_forests);
  for(size_t j=0; j<num_forests; ++j) {
    NumericMatrix postfits(nd, n);
    forest_fits_post[j] = postfits;
  }

  NumericMatrix yhat_post(nd, n);
  NumericMatrix etas_post(nd, num_forests);

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
    Rcpp::checkUserInterrupt();
    if(i%status_interval == 0) {
      Rcout << "iteration: " << i << endl;
      Rcout << "kappa acceptance rate: " << kappa_acpt_rate << " | current kappa: " << kappa << endl;
    }

    //-----------------------------------------------------------------------//
    //  update shape parameter kappa_c (using Metropolis-Hastings)           //
    //-----------------------------------------------------------------------//
    double kappa_star = exp(gen.normal(log(kappa), kappa_prop_sd));

    double log_a_num = ll_gamma_bcf(y, log_mean, kappa_star)
      + R::dbeta(kappa_star / (1.0 + kappa_star), kappa_a, kappa_b, true) - 2.0 * log1p(kappa_star)
      + log(kappa_star);
    double log_a_denom = ll_gamma_bcf(y, log_mean, kappa)
      + R::dbeta(kappa / (1.0 + kappa), kappa_a, kappa_b, true) - 2.0 * log1p(kappa)
      + log(kappa);

    if(log(gen.uniform()) < log_a_num - log_a_denom){
      kappa = kappa_star;
      kappa_acpt_rate = ((kappa_acpt_rate * i) + 1.0) / (i + 1.0);
    } else {
      kappa_acpt_rate = ((kappa_acpt_rate * i) + 0.0) / (i + 1.0);
    }

    //-----------------------------------------------------------------------//
    // update trees                                                          //
    //-----------------------------------------------------------------------//
    for(size_t s = 0; s < num_forests; ++s) {
      for(size_t j = 0; j < trees[s].size(); ++j) {

        // current tree fit
        fit_loglinear(trees[s][j], x_info[s], di[s], ftemp, node_pointers[s][j], false, prior_info[s].vanilla);

        for (size_t k = 0; k < n; k++){
          if (ftemp[k] != ftemp[k]){
            stop("nan in ftemp");
          }

          // remove current tree fit
          allfit[k]                -= prior_info[s].eta * ftemp[k];
          allfits[s][k]            -= prior_info[s].eta * ftemp[k];
          log_mean[k]              -= prior_info[s].eta * ftemp[k];

          double omega_k = di[s].omega[k];

          // GIG conjugate sufficient statistics:
          // u_vec = kappa * omega_k
          // r_tree = kappa * y * exp(log_mean) * omega_k
          u_vec[k]  = kappa * omega_k;
          double exp_log_mean = exp(log_mean[k]);
          if (exp_log_mean > 1e12) exp_log_mean = 1e12;
          if (exp_log_mean < 1e-12) exp_log_mean = 1e-12;
          r_tree[k] = kappa * y[k] * exp_log_mean * omega_k;

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
          log_mean[k]              += prior_info[s].eta * ftemp[k];
        }
      }
    }

    //-----------------------------------------------------------------------//
    // save results                                                          //
    //-----------------------------------------------------------------------//
    if( ((i>=burn) & (i % thin==0))) {
      kappa_post(save_ctr) = kappa;

      for(size_t k=0; k<n; k++) {
        yhat_post(save_ctr, k) = log_mean[k]; // Save log_mean = theta = -log lambda
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
                      _["coefs"] = post_coefs,
                      _["etas"] = etas_post,
                      _["kappa"] = kappa_post,
                      _["kappa_acceptance"] = kappa_acpt_rate,
                      _["tree_trace"] = final_tree_trace,
                      _["y_last"] = y
  ));
}
