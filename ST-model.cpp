#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * 3.14159265359);

// [[Rcpp::export]]
arma::mat tensor_rotate_X(arma::mat X,
                          arma::mat eigvecS,
                          arma::mat eigvecT,
                          arma::vec scales){
  
  int p = X.n_cols;
  int n = X.n_rows;
  int ns = eigvecS.n_cols;
  int nt = eigvecT.n_cols;
  arma::mat out(arma::size(X));
  
  for(int j = 0; j < p; j++){
    arma::mat X_j = arma::reshape(X.col(j), nt, ns);
    out.col(j) = arma::vectorise(eigvecT.t()*X_j*eigvecS);
  }
  
  for(int i = 0; i < n; i++){
    out.row(i) /= scales[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
double dmvnorm_kronecker_nugget(arma::colvec x, 
                                arma::colvec mu, 
                                arma::mat Sigma_S, //spatial
                                arma::mat Sigma_T, //temporal
                                double tau2, //variance
                                double sigma2, //nugget
                                bool logd = true){
  
  // x is vector of space time, filled by time first.
  // so that reshape(x, nt, ns) has a time series in each col.
  
  int ns = Sigma_S.n_cols;
  int nt = Sigma_T.n_cols;
  int k = ns*nt; //num of observations
  double out;
  
  //get eigen decomposition for both S T covs
  arma::vec eigval_S; arma::vec eigval_T;
  arma::mat eigvec_S; arma::mat eigvec_T;
  
  arma::eig_sym(eigval_S, eigvec_S, tau2*Sigma_S);
  arma::eig_sym(eigval_T, eigvec_T, Sigma_T);  
  
  //reshape and rotate response
  arma::mat y_mat = arma::reshape(x-mu, nt, ns);
  arma::mat y_new = eigvec_T.t()*y_mat*eigvec_S;
  arma::vec y_vec = arma::vectorise(y_new);
  
  //reshaped has diagonal covariance matrix.
  arma::vec eig_kron = arma::kron(eigval_S, eigval_T);
  arma::vec var_vec = eig_kron + sigma2;
  double quads = sum(pow(y_vec,2)/var_vec);
  
  out = -(k/2.0)*log2pi - 0.5 * quads - sum(log(var_vec))/2;     
  
  if (logd == false) {out = exp(out);}
  return(out);
}

// [[Rcpp::export]]
double dmvnorm_arma(arma::colvec x, 
                    arma::colvec mu, 
                    arma::mat Sigma, //spatial
                    bool logd = true){
  
  int k = x.n_rows; //num of observations
  double out;
  
  arma::mat Chol_Sigma = arma::trimatu(arma::chol(Sigma));
  arma::vec rooti = arma::solve(Chol_Sigma.t(), x-mu);
  double quads = arma::dot(rooti,rooti);
  
  out = -(k/2.0)*log2pi - 0.5*quads - sum(log(Chol_Sigma.diag()));
  
  if (logd == false) {out = exp(out);}
  return(out);
}

// [[Rcpp::export]]
double dmvnorm_arma_inv(arma::colvec x, 
                        arma::colvec mu, 
                        arma::mat Sigma_inv, //spatial
                        bool logd = true){
  
  int k = x.n_rows; //num of observations
  double out;
  
  arma::mat Chol_Sigma_inv = arma::trimatu(arma::chol(Sigma_inv));
  arma::vec rooti = Chol_Sigma_inv*(x-mu);
  double quads = arma::dot(rooti,rooti);
  
  out = -(k/2.0)*log2pi - 0.5*quads + sum(log(Chol_Sigma_inv.diag()));
  
  if (logd == false) {out = exp(out);}
  return(out);
}

// [[Rcpp::export]]
List ST_MCMC(arma::mat X_lin, //linear parts (just mean and covars, probably)
                  arma::colvec Y,
                  arma::mat dist_S,
                  arma::mat dist_T,
                  int n_samples,
                  int print_step,
                  bool tune,
                  int tune_every,
                  int tune_for,
                  int tune_after){
  
  //set up
  bool tune_sig2_bool = tune;
  bool tune_tau2_bool = tune;

  LogicalVector tune_phi_bool(2); for(int i = 0; i<2; i++){tune_phi_bool[i]=tune;}
  
  arma::mat par_covar_chol; //for empirical covariance of MCMC parameters
  
  //hyperparameters
  double sig2b = 10e6;
  
  //tuning parameters
  double tune_sig2 = .1;
  double tune_tau2 = 1;
  arma::vec tune_phi(2); tune_phi.ones();
  arma::mat did_accept_phi(2,n_samples); did_accept_phi.zeros();
  NumericVector did_accept_sig2(n_samples); did_accept_sig2.fill(0.0);
  NumericVector did_accept_tau2(n_samples); did_accept_tau2.fill(0.0);
  NumericVector did_accept_all(n_samples); did_accept_all.fill(0.0); //for simultaneous proposal
  
  //initial parameters
  arma::vec b_lin(X_lin.n_cols); b_lin.zeros(); //coefs
  arma::vec phi(2); phi.ones(); //ranges
  double sig2=1; double tau2=1; //scale and nugget
  
  arma::vec mn_lin = X_lin*b_lin;

  arma::vec b_lin_prior_mean(size(b_lin)); b_lin_prior_mean.zeros();

  //initial variance-covariance matrix
  arma::mat Sigma_S = exp(-phi[0]*dist_S);
  arma::mat Sigma_T = exp(-phi[1]*dist_T);
  
  //where to keep samples
  arma::mat b_lin_samples(X_lin.n_cols, n_samples);
  arma::mat phi_samples(2, n_samples);
  NumericVector sig2_samples(n_samples);
  NumericVector tau2_samples(n_samples);

  //some computation stuff
  arma::mat X_lin_t = X_lin.t()*X_lin;

  int nt = Sigma_T.n_cols;
  int ns = Sigma_S.n_cols;
  
  arma::mat I_lin; I_lin.eye(size(X_lin_t));
  int burn_in = tune_after + tune_for;
  
  //begin sampling
  for(int i = 0; i < n_samples; i++){
    
    // sample lin
    //get eigen decomposition for both S T covs
    arma::vec eigval_S; arma::vec eigval_T;
    arma::mat eigvec_S; arma::mat eigvec_T;
    
    arma::eig_sym(eigval_S, eigvec_S, tau2*Sigma_S);
    arma::eig_sym(eigval_T, eigvec_T, Sigma_T);  
    
    //reshape and rotate and scale "residuals"
    arma::vec new_vars = arma::kron(eigval_S, eigval_T) + sig2; //new ind. variances
    arma::mat y_mat = arma::reshape(Y, nt, ns); //reshape resids
    arma::mat y_new = eigvec_T.t()*y_mat*eigvec_S; //rotate resids
    arma::vec y_vec = arma::vectorise(y_new); //stretch back out
    y_vec = y_vec/sqrt(new_vars); //scale by sds
    
    //reshape design
    arma::mat X_lin_star = tensor_rotate_X(X_lin, eigvec_S, eigvec_T, sqrt(new_vars));
    arma::mat X_l_s_t = X_lin_star.t()*X_lin_star;
    arma::mat chol_lin = chol(X_l_s_t + I_lin/sig2b); //inverse cov
    arma::mat root_lin = trans(inv(trimatu(chol_lin))); //chol cov
    arma::mat Sig_lin = root_lin.t()*root_lin; //lin coef var

    arma::colvec mu_lin = Sig_lin*X_lin_star.t()*y_vec;
    b_lin = ((as<arma::rowvec>(rnorm(X_lin.n_cols)))*(root_lin) + mu_lin.t()).t();

    mn_lin = X_lin*b_lin;

    arma::vec mn_now = mn_lin;
    
    // if(i < burn_in){
    // sample sig2
    double sig2curr = sig2;
    double log_lk_curr = dmvnorm_kronecker_nugget(Y, 
                                                  mn_now, 
                                                  Sigma_S, //spatial
                                                  Sigma_T, //temporal
                                                  tau2, //variance
                                                  sig2curr); //nugget
    
    double sig2prop = exp(tune_sig2*rnorm(1)[0] + log(sig2curr));
    double log_lk_prop = dmvnorm_kronecker_nugget(Y, 
                                                  mn_now, 
                                                  Sigma_S, //spatial
                                                  Sigma_T, //temporal
                                                  tau2, //variance
                                                  sig2prop); //nugget
    
    double log_q_prop = log_lk_prop + R::dgamma(sig2prop, 0.01, 100, 1) + log(sig2prop);
    double log_q_curr = log_lk_curr + R::dgamma(sig2curr, 0.01, 100, 1) + log(sig2curr);
    
    double draw_r = runif(1)[0];
    if(log(draw_r) < log_q_prop - log_q_curr){
      sig2 = sig2prop;
      did_accept_sig2[i] = 1;
    } else {
      sig2 = sig2curr;
    }
    
    // sample tau2
    double tau2curr = tau2;
    log_lk_curr = dmvnorm_kronecker_nugget(Y, 
                                           mn_now, 
                                           Sigma_S, //spatial
                                           Sigma_T, //temporal
                                           tau2curr, //variance
                                           sig2); //nugget
    
    double tau2prop = exp(tune_tau2*rnorm(1)[0] + log(tau2curr));
    log_lk_prop = dmvnorm_kronecker_nugget(Y, 
                                           mn_now, 
                                           Sigma_S, //spatial
                                           Sigma_T, //temporal
                                           tau2prop, //variance
                                           sig2); //nugget
    
    log_q_prop = log_lk_prop + R::dgamma(tau2prop, 0.01, 100, 1) + log(tau2prop);
    log_q_curr = log_lk_curr + R::dgamma(tau2curr, 0.01, 100, 1) + log(tau2curr);
    
    draw_r = runif(1)[0];
    if(log(draw_r) < log_q_prop - log_q_curr){
      tau2 = tau2prop;
      did_accept_tau2[i] = 1;
    } else {
      tau2 = tau2curr;
    }

    // sample phi0
    double phi0curr = phi[0];
    log_lk_curr = dmvnorm_kronecker_nugget(Y, 
                                           mn_now, 
                                           Sigma_S, //spatial
                                           Sigma_T, //temporal
                                           tau2, //variance
                                           sig2); //nugget
    
    double phi0prop = exp(tune_phi[0]*rnorm(1)[0] + log(phi0curr));
    arma::mat Sigma_S_prop = exp(-phi0prop*dist_S);
    
    log_lk_prop = dmvnorm_kronecker_nugget(Y, 
                                           mn_now, 
                                           Sigma_S_prop, //spatial
                                           Sigma_T, //temporal
                                           tau2, //variance
                                           sig2); //nugget
    
    log_q_prop = log_lk_prop + R::dgamma(phi0prop, 0.01, 100, 1) + log(phi0prop);
    log_q_curr = log_lk_curr + R::dgamma(phi0curr, 0.01, 100, 1) + log(phi0curr);
    
    draw_r = runif(1)[0];
    if(log(draw_r) < log_q_prop - log_q_curr){
      phi[0] = phi0prop;
      Sigma_S = Sigma_S_prop;
      did_accept_phi(0,i) = 1;
    } else {
      phi[0] = phi0curr;
    }
    
    // sample phi1
    double phi1curr = phi[1];
    log_lk_curr = dmvnorm_kronecker_nugget(Y, 
                                           mn_now, 
                                           Sigma_S, //spatial
                                           Sigma_T, //temporal
                                           tau2, //variance
                                           sig2); //nugget
    
    double phi1prop = exp(tune_phi[1]*rnorm(1)[0] + log(phi1curr));
    arma::mat Sigma_T_prop = exp(-phi1prop*dist_T);
    
    log_lk_prop = dmvnorm_kronecker_nugget(Y, 
                                           mn_now, 
                                           Sigma_S, //spatial
                                           Sigma_T_prop, //temporal
                                           tau2, //variance
                                           sig2); //nugget
    
    log_q_prop = log_lk_prop + R::dgamma(phi1prop, 0.01, 100, 1) + log(phi1prop);
    log_q_curr = log_lk_curr + R::dgamma(phi1curr, 0.01, 100, 1) + log(phi1curr);
    
    draw_r = runif(1)[0];
    if(log(draw_r) < log_q_prop - log_q_curr){
      phi[1] = phi1prop;
      Sigma_T = Sigma_T_prop;
      did_accept_phi(1,i) = 1;
    } else {
      phi[1] = phi1curr;
    }
    
    // tune sig2 sampler every tune_every steps for tune_for iterations
    if(tune_sig2_bool){
      if((i < tune_for + tune_after) & (i > tune_after + 1)){
        if((i%tune_every == 0) & (i != 0)){
          Range which_tune(i-tune_every,i);
          double mean_accept = mean(did_accept_sig2[which_tune]);
          if(mean_accept > 0.5){
            tune_sig2 = tune_sig2*(1+0.1);
          } else if(mean_accept < 0.4){
            tune_sig2 = tune_sig2*(1-0.1);
          } else {
            tune_sig2_bool = FALSE;
            Rcpp::Rcout << "adaption complete for sig2 at iteration " << i << "." << std::endl;
          }
        }
      }
    }
    
    if(tune_tau2_bool){
      if((i < tune_for + tune_after) & (i > tune_after + 1)){
        if((i%tune_every == 0) & (i != 0)){
          Range which_tune(i-tune_every,i);
          double mean_accept = mean(did_accept_tau2[which_tune]);
          if(mean_accept > 0.5){
            tune_tau2 = tune_tau2*(1+0.1);
          } else if(mean_accept < 0.4){
            tune_tau2 = tune_tau2*(1-0.1);
          } else {
            tune_tau2_bool = FALSE;
            Rcpp::Rcout << "adaption complete for tau2 at iteration " << i << "." << std::endl;
          }
        }
      }
    }
    
    for(int j = 0; j<phi.n_elem; j++){
      if(tune_phi_bool[j]){
        if((i < tune_for + tune_after) & (i > tune_after + 1)){
          if((i%tune_every == 0) & (i != 0)){
            double mean_accept = sum(did_accept_phi(j,arma::span(i-tune_every,i)))/tune_every;
            if(mean_accept > 0.5){
              tune_phi[j] = tune_phi[j]*(1+0.1);
            } else if(mean_accept < 0.4){
              tune_phi[j] = tune_phi[j]*(1-0.1);
            } else {
              tune_phi_bool[j] = FALSE;
              Rcpp::Rcout << "adaption complete for phi " << j << " at iteration " << i << "." << std::endl;
            }
          }
        }
      }
    }
    
    // collect samples
    // lambda_samples[i] = lambda;
    b_lin_samples.col(i) = b_lin;
    phi_samples.col(i) = phi;
    sig2_samples[i] = sig2;
    tau2_samples[i] = tau2;
    
    if((i%print_step == 0) & (i != 0)){
      Rcpp::Rcout << "iteration: " << i << std::endl;
    }
  }
  
  List out; //returns samples
  out["tau2.samples"] = tau2_samples;
  out["phi.samples"] = phi_samples;
  out["sig2.samples"] = sig2_samples;
  out["b.lin.samples"] = b_lin_samples;
  out["did.accept"] = did_accept_all; 
  return out;
}

// [[Rcpp::export]]
arma::vec predict_new(arma::vec mn_new,
                      arma::vec mn_old,
                      arma::vec Y_old,
                      arma::mat Sigma_S_old, //spatial
                      arma::mat Sigma_T_old, //temporal
                      arma::mat Sigma_S_new, //spatial
                      arma::mat Sigma_T_new, //temporal
                      arma::mat Sigma_S_on, //spatial
                      arma::mat Sigma_T_on, //temporal
                      double tau2, //variance
                      double sigma2){ //nugget
  
  int ns = Sigma_S_old.n_cols;
  int nt = Sigma_T_old.n_cols;
  int nsn = Sigma_S_new.n_cols;
  int ntn = Sigma_T_new.n_cols;
  int k = nsn*ntn; //num of observations
  arma::vec out;
  
  //get eigen decomposition for both S T covs
  arma::vec eigval_S; arma::vec eigval_T;
  arma::mat eigvec_S; arma::mat eigvec_T;
  
  arma::eig_sym(eigval_S, eigvec_S, tau2*Sigma_S_old);
  arma::eig_sym(eigval_T, eigvec_T, Sigma_T_old);  
  
  //reshape and rotate response
  arma::mat y_mat = arma::reshape(Y_old - mn_old, nt, ns);
  arma::mat y_new = eigvec_T.t()*y_mat*eigvec_S;
  arma::vec y_vec = arma::vectorise(y_new);
  
  //reshaped has diagonal covariance matrix.
  arma::vec eig_kron = arma::kron(eigval_S, eigval_T);
  arma::vec var_vec = eig_kron + sigma2;
  arma::vec y_star = y_vec/var_vec;
  
  //the mean for prediction:
  arma::vec new_mn = mn_new + tau2*arma::vectorise(Sigma_T_on.t()*eigvec_T*arma::reshape(y_star, nt, ns)*eigvec_S.t()*Sigma_S_on);
  
  //the subtracted variance:
  arma::mat Sig_star = kron(eigvec_S.t()*Sigma_S_on, eigvec_T.t()*Sigma_T_on);
  for(int i = 0; i < Sig_star.n_rows; i++){
    Sig_star.row(i) /= sqrt(var_vec[i]);
  }
  
  arma::mat Chol_Sigma = arma::chol(tau2*kron(Sigma_S_new, Sigma_T_new) + sigma2*arma::eye(k,k) - pow(tau2,2)*Sig_star.t()*Sig_star);
  out = Chol_Sigma.t()*(as<arma::vec>(rnorm(k))) + new_mn;
  return(out);
}
