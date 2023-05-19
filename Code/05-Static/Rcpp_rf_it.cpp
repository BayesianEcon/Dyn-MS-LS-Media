#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <list>


// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include <RcppDist.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends( RcppDist , RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;




Function f("rnorm"); 


// // [[Rcpp::export]]
// arma::mat Th0Sampler (arma::mat alpha, arma::mat beta ) {
//   arma::mat theta_bar;
//   theta_bar = alpha + beta;
//   return theta_bar;
// }



// [[Rcpp::export]]
double mean_mat(NumericMatrix A){ 
  
  double n_i =  A.nrow() ;
  double n_j =  A.ncol();
  
  double result = 0;
  
  for (int i = 0; i < n_i; i++) {
    for (int j= 0; j < n_j; j++ ) {
      
      double aux = A(i,j);
      
      result = result + aux;
      
    }}
  
  return(result/(n_i*n_j));
}


// [[Rcpp::export]]
double sd_mat(NumericMatrix A){ 
  
  double n_i =  A.nrow() ;
  double n_j =  A.ncol();
  arma::vec val((n_i*n_j  - n_i)/2);
  
  double count = 0;
  
  for (int i = 0; i < n_i; i++) {
    for (int j= 0; j < i; j++ ) {
      
      double aux = A(i,j);
      val[count] = aux;
      count = count + 1;
      
    }}
  
  double sd_res = arma::stddev(val);
  return(sd_res);
}


// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// [[Rcpp::export]]
List Th0Sampler(arma::vec alpha, arma::vec beta, arma::vec m0, double ki, double a,
                arma:: mat S, double M) {
  
  arma::vec theta_bar = {0,0};
  theta_bar[0] = mean(alpha);
  theta_bar[1]  = mean(beta);
  
  arma::vec theta_ij_diff_1 =  alpha - mean(alpha);
  arma::vec theta_ij_diff_2 =  beta - mean(beta);
  
  auto theta_ij_diff =  arma::join_rows( theta_ij_diff_1, theta_ij_diff_2 );
  
  arma::mat theta_ij_diff_t = theta_ij_diff.t(); 
  
  arma::mat Psi(2,2,fill::zeros	);
  
  
  for (int i=0; i< M; i++) {
    
    arma::vec z = theta_ij_diff_t.col(i);
    arma::rowvec z_trans = z.t();
    arma::mat prod = z*z_trans;
    
    Psi = Psi + prod;
  }
  
  arma::vec th_bar_diff_e1 = (theta_bar - m0);
  arma::rowvec th_bar_diff_e2 =th_bar_diff_e1.t();
  arma::mat  th_bar_diff =  th_bar_diff_e1*th_bar_diff_e2;
  
  arma::vec mu0 = (ki*m0 + M*theta_bar)/(M+ki);
  double lambda = ki+M;
  arma::mat SS =  S + Psi + ((ki*M)/(M+ki))*th_bar_diff;
  double nu = a + M; 
  
  arma::mat invSS =  riwish( nu, SS);
  arma::mat theta0 = mvrnormArma(1, mu0, (1/lambda)*invSS);
  
  List L = List::create( Named("theta0")= theta0, Named("Ups") = invSS); 
  return L;
}


// [[Rcpp::export]]
double clog_sum_exp(arma::vec x) {
  
  auto b = max(x); 
  double result = (b + log(sum(exp(x-b))));
  return result;
  
}

// [[Rcpp::export]]

NumericVector csoftmax(NumericVector x) {
  
  NumericVector result = exp(x - clog_sum_exp(x));
  return result;
  
}


NumericVector logistic(NumericVector x) {
  NumericVector result = 1/(1+exp(-x));
  return result;
  
}

// [[Rcpp::export]]
double logistic_d(double x) {
  double result = 1/(1+exp(-x));
  return result;
}


// [[Rcpp::export]]
vec arma_sort(vec x, vec y, vec z) {
  // Order the elements of x by sorting y and z;
  // we order by y unless there's a tie, then order by z.
  // First create a vector of indices
  uvec idx = regspace<uvec>(0, x.size() - 1);
  // Then sort that vector by the values of y and z
  std::sort(idx.begin(), idx.end(), [&](int i, int j){
    if ( y[i] == y[j] ) {
      return z[i] < z[j];
    }
    return y[i] < y[j];
  });
  // And return x in that order
  return x(idx);
}



// [[Rcpp::export]]
double logbivnorm(arma::vec x,    arma::vec  mean,  arma::mat sigma) { 
  
  arma::vec diff = x- mean;
  arma::rowvec diff_t = diff.t();
  arma::mat inv_sigma = sigma;
  inv_sigma = inv(inv_sigma);
  
  double result = (-0.5*diff_t*inv_sigma*diff).eval()(0,0);
  return result;
};




// Foward Filtering Backward Sampling

// [[Rcpp::export]]
List FFBS(NumericVector Acc, NumericVector Bcc, NumericMatrix alef,NumericMatrix bet, NumericMatrix Eta, NumericMatrix P, double T) {
  
  for(int t = 1; t< T; t++) {
    for(int j = 0; j< 2; j++) {
      for(int i = 0; i < 2; i++){
        
        Acc[i]  = alef(t-1, i) + log(P(i, j)) + Eta(t,j);
        
        
      }
      
      alef(t,j) = clog_sum_exp(Acc);
    }
  }
  
  
  for(int t = 0; t< T; t++) {
    alef.row(t) = csoftmax(alef.row(t));
  }
  
  for(int t = (T-1); t > 0; t--) {
    for(int j = 0; j< 2; j++) {
      for(int i = 0; i < 2; i++){
        
        Bcc[i]  = bet(t, i) + log(P(j, i)) +  Eta(t,i);
        
      }
      
      bet(t-1, j) = clog_sum_exp(Bcc);
      
    }}
  
  // 
  for(int t = 0; t< T; t++) {
    
    
    bet(t,_) = csoftmax(bet(t,_));
    
  }
  
  
  List L = List::create( Named("alef") = alef,  Named("bet") = bet);
  return L;
  
}



// [[Rcpp::export]]

double tauSamplingEn(NumericVector ncomments, NumericVector  zi_2,  double a_tau, double b_tau) {
  
  double a_bar_tau =  a_tau +  sum(ncomments);
  double b_bar_tau =  b_tau +  sum(exp(zi_2));
  
  double tau = R::rgamma(a_bar_tau, 1/b_bar_tau);
  
  return tau;
}



// [[Rcpp::export]]

NumericVector  gammaSamplingEn(NumericVector leaning, NumericVector  zi_1, double phi, double gamma_0, double gamma_1, 
                               double a_gamma_0, double a_gamma_1, double b_gamma_0, double b_gamma_1, double prop_sigma) {
  
  Function f("dbeta");
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function g("dmvnorm");
  // 
  NumericVector  gamma = {gamma_0, gamma_1 };
  //
  arma::mat propdiag(2,2,fill::zeros);
  
  propdiag(0,0) = prop_sigma;
  propdiag(1,1) = prop_sigma;
  
  double gamma_prop_0 = R::rnorm(gamma_0,prop_sigma );
  
  //
  // double gamma_prop_1 = -1;
  // while(gamma_prop_1 < 0){
   double gamma_prop_1 = R::rnorm(gamma_1, prop_sigma );
    
  // }
  // 
  // arma::mat  gamma_prop_arma = mvrnormArma(1, gamma, propdiag );
  // 
  // NumericVector  gamma_prop = as<NumericVector>(wrap(gamma_prop_arma));
  
  NumericVector  gamma_prop(2) ;
  gamma_prop[0] = gamma_prop_0;
  gamma_prop[1] = gamma_prop_1;
  
  //gamma_prop[1] = 0;
  // 
  // if (gamma_prop[0] > gamma_prop[1] ) { double aux =  gamma_prop[1] ;   gamma_prop[0] =    gamma_prop[1];  gamma_prop[1] = aux ; }
  // 
  
  
  NumericVector  a_beta = phi*logistic(gamma_0 + abs(gamma_1)*zi_1);
  NumericVector  b_beta = phi*(1-logistic(gamma_0 + abs(gamma_1)*zi_1));
  // 
  NumericVector  a_beta_prop =  phi*logistic(gamma_prop[0] + abs(gamma_prop[1])*zi_1);
  NumericVector  b_beta_prop = phi*(1-logistic(gamma_prop[0] + abs(gamma_prop[1])*zi_1));
  
  double size = leaning.length();
  NumericVector  Beta(size);
  NumericVector Beta_prop(size);
  
  for(int t = 0; t< (size-1); t++) {
    Beta[t] = R::dbeta( leaning[t],  a_beta[t], b_beta[t],  true);
    Beta_prop[t] = R::dbeta( leaning[t], a_beta_prop[t],  b_beta_prop[t], true);}
  
  double prior = R::dnorm( gamma_0, a_gamma_0, b_gamma_0, true ) + R::dnorm( gamma_1, a_gamma_1, b_gamma_1, true );
  double prior_prop = R::dnorm( gamma_prop[0], a_gamma_0, b_gamma_0, true ) + R::dnorm( gamma_prop[1], a_gamma_1, b_gamma_1, true );
  
  // double norm_c_g1 = R::pnorm(gamma_1/sqrt(prop_sigma[1]), 0, 1 ,  true,  true) ;
  // double norm_c_g1_prop = R::pnorm(gamma_prop[1]/sqrt(prop_sigma[1]), 0, 1 ,  true,  true) ;
  
  // double norm_c_g0 = R::pnorm(gamma_0/sqrt(prop_sigma[0]), 0, 1 ,  true,  true) ;
  // double norm_c_g0_prop = R::pnorm(gamma_prop[0]/sqrt(prop_sigma[0]), 0, 1 ,  true,  true) ;
  
  double test =  sum(Beta_prop) - sum(Beta) ; //prior_prop - prior +
  
  // cout << test;
  // +  norm_c - norm_c_prop ;
  double rand = R::runif(0,1);
  
  if(test >  log(rand)){  return gamma_prop;  }else{ return gamma; }  ;
  
}



// [[Rcpp::export]]

double  phiSamplingEn(NumericVector leaning, NumericVector  zi_1, double phi, double gamma_0, double gamma_1, 
                      double a_phi, double b_phi, double prop_sigma) {
  
  double phi_prop = -1;
  
  while (phi_prop < 0) {
    
    phi_prop = R::rnorm(phi,prop_sigma );
  }
  
  // cout << phi_prop;
  
  NumericVector  a_beta = phi*logistic(gamma_0 + abs(gamma_1)*zi_1);
  NumericVector  b_beta = phi*(1-logistic(gamma_0 + abs(gamma_1)*zi_1));
  
  NumericVector  a_beta_prop =  phi_prop*logistic(gamma_0 + abs(gamma_1)*zi_1);
  NumericVector  b_beta_prop =  phi_prop*(1-logistic(gamma_0 + abs(gamma_1)*zi_1));
  
  double size = leaning.length();
  NumericVector  Beta(size);
  NumericVector Beta_prop(size);
  
  for(int t = 0; t< (size-1); t++) {
    Beta[t] = R::dbeta( leaning[t],  a_beta[t], b_beta[t],  true);
    Beta_prop[t] = R::dbeta( leaning[t], a_beta_prop[t],  b_beta_prop[t], true);}
  
  // gamma shape and scale
  double prior = R::dgamma(phi,a_phi,1/b_phi, true);
  double prior_prop = R::dgamma(phi_prop,a_phi,1/b_phi, true);
  
  double norm_c = R::pnorm(phi/(prop_sigma), 0, 1 ,  true,  true) ;
  double norm_c_prop = R::pnorm(phi_prop/(prop_sigma), 0, 1 ,  true,  true) ;
  
  
  double test =  sum(Beta_prop) - sum(Beta)  ;  //+  norm_c - norm_c_prop  ; //+ prior_prop - prior
  // cout << sum(Beta_prop);
  // cout <<   sum(Beta) ;
  
  
  double rand = R::runif(0,1);
  
  if(test > log(rand)){  return phi_prop;  }else{ return phi; }  ;
  
}


// [[Rcpp::export]]
NumericVector cseq(double first, double last, double by){
  
  int n = (last - first)/by + 1;
  NumericVector result(n);
  result[0] = first;
  
  
  for(int i = 1; i< n ; i++) {
    
    result[i] = result[i-1] + by;
    
  }
  
  return result;
}


static double const log2pi = std::log(2.0 * M_PI);



// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat const &x,  
                      arma::rowvec const &mean,  
                      arma::mat const &sigma, 
                      bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z      = (x.row(i) - mean) * rooti;    
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}


NumericMatrix  thetaijSamplingMALAEn(  NumericVector alpha_ij, NumericVector beta_ij   ,NumericVector w, NumericVector wl1, NumericVector ones,  NumericVector distance ,  NumericVector theta, arma::mat Ups, double M, double Time, NumericVector epsilon){
  
  NumericVector alpha = rep(alpha_ij, Time);
  NumericVector beta =  rep(beta_ij, Time);
  
  NumericVector nabla1 = w - exp(ones*alpha + wl1*beta - distance);
  NumericVector nabla2 = (w - exp(ones*alpha + wl1*beta - distance))*wl1;
  
  NumericMatrix sumGrad(M, 2);
  
  NumericVector nabla1_sub(100);
  NumericVector nabla2_sub(100);
  
  for(int i = 0; i < M; i++) {
    
    
    NumericVector idx = cseq(i,  (M*Time-1), M);
    nabla1_sub = nabla1[idx];
    nabla2_sub = nabla2[idx];
    
    
    sumGrad(i, 0)= sum(nabla1_sub)  ;
    sumGrad(i, 1)= sum(nabla2_sub)   ;
    
  };
  
  arma::mat invUps =  arma::inv(Ups);
  
  NumericMatrix theta_ij(M,2);
  theta_ij(_, 0) = alpha_ij;
  theta_ij(_, 1) = beta_ij;
  
  
  NumericMatrix gradPrior(M, 2);
  for(int i = 0; i < M; i++) {
    
    arma::vec diff =  (theta_ij(i,_) - theta) ;
    arma::vec res = -1 * invUps * diff;
    NumericVector result = as<NumericVector>(wrap(res));
    
    gradPrior(i, _) = result;
    
  };
  // 
  NumericMatrix theta_ij_aug(M,2);
  
  
  theta_ij_aug(_,0) = theta_ij(_,0) + (0.5*epsilon[0]*epsilon[0])*(gradPrior(_,0) + sumGrad(_, 0));
  theta_ij_aug(_,1) = theta_ij(_,1) + (0.5*epsilon[1]*epsilon[1])*(gradPrior(_,1) + sumGrad(_, 1));
  
  NumericMatrix theta_prop(M,2) ;
  
  for(int i = 0; i < M; i++) {
    
    arma::vec mu = theta_ij_aug(i,_);
    arma::mat sigma(2,2);
    sigma(0,0) = epsilon[0];
    sigma(1,1) = epsilon[1];
    
    NumericVector aux =  as<NumericVector>(wrap(mvrnormArma(1, mu,sigma) )) ;
    theta_prop(i,_ ) =  aux;
  }
  
  NumericVector alpha_prop = rep(theta_prop(_,0), Time);
  NumericVector beta_prop =  rep(theta_prop(_,1), Time);
  
  
  NumericVector pois(Time*M);
  NumericVector pois_prop(Time*M) ;
  
  for(int i = 0; i < Time*M; i++) {
    
    pois[i] =  R::dpois(w[i], exp(ones[i]*alpha[i] + wl1[i]*beta[i] - distance[i] ), true);
    pois_prop[i] = R::dpois(w[i], exp(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i] - distance[i] ), true);
    
  }
  
  NumericMatrix sumEL(M, 2);
  
  for(int i = 0; i < M; i++) {
    
    
    NumericVector idx = cseq(i,  (M*Time-1), M);
    NumericVector pois_sub = pois[idx];
    NumericVector pois_prop_sub = pois_prop[idx];
    
    sumEL(i, 0)= sum(pois_sub)  ;
    sumEL(i, 1)= sum(pois_prop_sub)   ;
  };
  
  NumericMatrix result(M,2);
  
  for(int i = 0; i < M; i++) {
    
    arma::vec thprop  = theta_prop(i,_);
    arma::vec thij  = theta_ij(i,_);
    
    
    double res_norm = logbivnorm(thprop, theta, Ups) -   logbivnorm(thij, theta, Ups);
    
    double res_beta = sumEL(i, 1)- sumEL(i, 0);
    double logr =  res_beta + res_norm ;
    
    double logunif = log(R::runif(0,1));
    
    if(logr > logunif){ result(i,_)  =  theta_prop(i,_);   }else{  result(i,_)  =  theta_ij(i,_);};
  };
  return result;
  
};


// [[Rcpp::export]]
NumericMatrix thetaijSamplingEn(  NumericVector alpha_ij, NumericVector beta_ij   ,NumericVector w, NumericVector wl1, NumericVector ones,  NumericVector distance ,  NumericVector theta, arma::mat Ups, double M, double Time, NumericVector epsilon){
  
  NumericVector alpha = rep(alpha_ij, Time);
  NumericVector beta =  rep(beta_ij, Time);
  
  NumericMatrix theta_ij(M,2);
  theta_ij(_, 0) = alpha_ij;
  theta_ij(_, 1) = beta_ij;
  
  
  NumericMatrix theta_prop(M,2) ;
  
  for(int i = 0; i < M; i++) {
    
    
    NumericVector mu =  {0,0};
    mu[0]  = theta_ij(i,0);
    mu[1]  = theta_ij(i,1);
    // 
    // mu[0]  = theta[0];
    // mu[1]  = theta[1];
    
    arma::mat sigma(2,2,fill::zeros);
    sigma(0,0) = epsilon[0];
    sigma(1,1) = epsilon[1];
    
    arma::mat aux =  mvrnormArma(1, mu, sigma)  ;
    theta_prop(i, 0 ) =  aux(0,0);
    theta_prop(i, 1 ) =  aux(0,1);
  };
  
  NumericVector alpha_prop = rep(theta_prop(_,0), Time);
  NumericVector beta_prop =  rep(theta_prop(_,1), Time);
  
  
  NumericVector pois(Time*M);
  NumericVector pois_prop(Time*M) ;
  
  
  for(int i = 0; i < Time*M; i++) {
    
    pois[i] =  R::dpois(w[i], exp(ones[i]*alpha[i] + wl1[i]*beta[i] - distance[i] ), true);
    pois_prop[i] = R::dpois(w[i], exp(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i] - distance[i] ), true);
    
    // pois[i] =  w[i]*(ones[i]*alpha[i] + wl1[i]*beta[i] ) - exp(ones[i]*alpha[i] + wl1[i]*beta[i] - distance[i]);
    // pois_prop[i] = w[i]*(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i]) - exp(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i] - distance[i]);
    // 
  };
  
  
  NumericMatrix sumEL(M, 2);
  
  for(int i = 0; i < M; i++) {
    
    
    NumericVector idx = cseq(i,  (M*Time-1), M);
    NumericVector pois_sub = pois[idx];
    NumericVector pois_prop_sub = pois_prop[idx];
    
    
    sumEL(i, 0)= sum(pois_sub)  ;
    sumEL(i, 1)= sum(pois_prop_sub)   ;
    
  };
  
  NumericMatrix result(M,2);
  
  for(int i = 0; i < M; i++) {
    
    arma::vec thprop  = theta_prop(i,_);
    arma::vec thij  = theta_ij(i,_);
    
    
    double res_norm = logbivnorm(thprop, theta, Ups) -   logbivnorm(thij, theta, Ups);
    double res_beta = sumEL(i, 1) - sumEL(i, 0);
    
    double logr =  res_beta + res_norm ;
    
    double logunif = log(R::runif(0,1));
    
    if(logr > logunif){
      double a_prop  =  theta_prop(i,0);
      double b_prop  =  theta_prop(i,1);
      
      result(i,0)  = a_prop;
      result(i,1)  = b_prop;  }else{
        
        double a  =  theta_ij(i,0); 
        double b  =  theta_ij(i,1); 
        
        result(i,0)  =  a;
        result(i,1)  =  b ;  };
      
  };
  
  return result;
  
  
}

// [[Rcpp::export]]
NumericMatrix cdiff(NumericVector A){
  
  int sizeA = A.length();
  
  NumericMatrix result(sizeA, sizeA);
  
  for(int i = 0; i < sizeA; i++) {
    double Ai = A(i);
    for(int j = 0; j < sizeA; j++) {
      double Aj = A(j);
      double aux =  (Ai-Aj);
      result(i,j) = aux;
      
    }}
  
  return result;
}




// [[Rcpp::export]]
List thetaijSamplingAdEn(  NumericVector alpha_ij, NumericVector beta_ij   ,NumericVector w, NumericVector wl1, NumericVector ones,  NumericVector distance ,  NumericVector theta, arma::mat Ups, NumericVector lam_ad, List Sigma_ad, double M, double Time){
  
  NumericVector alpha = rep(alpha_ij, Time);
  NumericVector beta =  rep(beta_ij, Time);
  
  NumericMatrix theta_ij(M,2);
  theta_ij(_, 0) = alpha_ij;
  theta_ij(_, 1) = beta_ij;
  
  
  NumericMatrix theta_prop(M,2) ;
  
  for(int i = 0; i < M; i++) {
    
    
    NumericVector mu =  {0,0};
    mu[0]  = theta_ij(i,0);
    mu[1]  = theta_ij(i,1);
    // 
    // mu[0]  = theta[0];
    // mu[1]  = theta[1];
    
    arma::mat sigma = Sigma_ad[i];
    double lm = lam_ad[i];
    
    arma::mat aux =  mvrnormArma(1, mu, lm*sigma)  ;
    theta_prop(i, 0 ) =  aux(0,0);
    theta_prop(i, 1 ) =  aux(0,1);
  };
  
  NumericVector alpha_prop = rep(theta_prop(_,0), Time);
  NumericVector beta_prop =  rep(theta_prop(_,1), Time);
  
  
  NumericVector pois(Time*M);
  NumericVector pois_prop(Time*M) ;
  
  
  for(int i = 0; i < Time*M; i++) {
    
    pois[i] =  R::dpois(w[i], exp(ones[i]*alpha[i] + wl1[i]*beta[i] - distance[i]*distance[i] ), true);
    pois_prop[i] = R::dpois(w[i], exp(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i] - distance[i]*distance[i] ), true);
    
    // pois[i] =  w[i]*(ones[i]*alpha[i] + wl1[i]*beta[i] ) - exp(ones[i]*alpha[i] + wl1[i]*beta[i] - distance[i]);
    // pois_prop[i] = w[i]*(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i]) - exp(ones[i]*alpha_prop[i] + wl1[i]*beta_prop[i] - distance[i]);
    // 
  };
  
  
  NumericMatrix sumEL(M, 2);
  
  for(int i = 0; i < M; i++) {
    
    
    NumericVector idx = cseq(i,  (M*Time-1), M);
    NumericVector pois_sub = pois[idx];
    NumericVector pois_prop_sub = pois_prop[idx];
    
    
    sumEL(i, 0)= sum(pois_sub)  ;
    sumEL(i, 1)= sum(pois_prop_sub)   ;
    
  };
  
  NumericMatrix result(M,2);
  NumericVector acceptance(M);
  
  for(int i = 0; i < M; i++) {
    
    arma::vec thprop  = theta_prop(i,_);
    arma::vec thij  = theta_ij(i,_);
    
    
    double res_norm = logbivnorm(thprop, theta, Ups) -   logbivnorm(thij, theta, Ups);
    double res_beta = sumEL(i, 1) - sumEL(i, 0);
    
    double logr =  res_beta + res_norm ;
    
    if(exp(logr) < 1 ){ acceptance[i] = exp(logr); }else{acceptance[i] = 1;};
    
    double logunif = log(R::runif(0,1));
    
    if(logr > logunif){
      double a_prop  =  theta_prop(i,0);
      double b_prop  =  theta_prop(i,1);
      
      result(i,0)  = a_prop;
      result(i,1)  = b_prop;  }else{
        
        double a  =  theta_ij(i,0); 
        double b  =  theta_ij(i,1); 
        
        result(i,0)  =  a;
        result(i,1)  =  b ;  };
      
  };
  
  
  List resutl_List = List::create( Named("theta_res")= result,  Named("Acceptance")=  acceptance ); 
  return resutl_List;
  
  
}

// [[Rcpp::export]]
List AdaptiveUpdate( NumericMatrix theta_ij, NumericVector lam_ad, List Sigma_ad, NumericMatrix mu_mat, NumericVector acceptance  ,  double acc_star, double M, double ite){
  
  for(int i = 0; i < M; i++) {
    
    double gamma = 1/pow(ite + 1, 0.55 );
    // double gamma = 1/pow(ite + 1, 1/(1+ lam_ad[i] ));
    double log_lam  = log(lam_ad[i]) + gamma*(acceptance[i] -  acc_star);
    
    NumericVector theta_ij_i = theta_ij(i,_);
    NumericVector mu_mat_i = mu_mat(i,_);
    
    
    NumericVector diff =  (theta_ij_i -  mu_mat_i );
    NumericVector mu = mu_mat(i,_) + gamma*diff;
    
    mu_mat(i,_) = mu;
    
    arma::vec dif_vec = diff;
    arma::rowvec dif_vec_t = dif_vec.t();
    arma::mat Sigma = Sigma_ad[i];
    
    NumericMatrix Sigma_up =  wrap(Sigma + gamma*(dif_vec*dif_vec_t - Sigma ) );
    
    Sigma_ad[i] = Sigma_up;
    
    lam_ad[i] = exp(log_lam);
    
    
  };
  
  
  List resutl_List = List::create( Named("lam_ad")= lam_ad, Named("mu_mat")=  mu_mat ,  Named("Sigma_ad")=  Sigma_ad );
  return resutl_List;
  
};


// [[Rcpp::export]]
NumericVector ecldist_point(NumericVector A, NumericVector B, int i){
  
  int sizeA = A.length();
  
  NumericVector result(sizeA-1);
  
  double Ai = A(i);
  double Bi = B(i);
  int idx = 0;
  
  for(int j = 0; j < sizeA; j++) {
    
    if( i != j){
      
      
      
      double Aj = A(j);
      double Bj = B(j);
      
      double aux = sqrt(pow(Ai-Aj,2) +  pow(Bi-Bj,2)) ;
      result[idx] = aux;
      
      idx = idx +1;
    }
  }
  
  return result;
};


// [[Rcpp::export]]
NumericMatrix ecldist(NumericVector A, NumericVector B){
  
  int sizeA = A.length();
  
  NumericMatrix result(sizeA, sizeA);
  
  for(int i = 0; i < sizeA; i++) {
    double Ai = A(i);
    double Bi = B(i);
    
    for(int j = 0; j < sizeA; j++) {
      double Aj = A(j);
      double Bj = B(j);
      
      double aux = sqrt(pow(Ai-Aj,2) +  pow(Bi-Bj,2)) ;
      result(i,j) = aux;
      
    }}
  
  return result;
};

// [[Rcpp::export]]
NumericMatrix ecldistp2(NumericVector A, NumericVector B){
  
  int sizeA = A.length();
  
  NumericMatrix result(sizeA, sizeA);
  
  for(int i = 0; i < sizeA; i++) {
    double Ai = A(i);
    double Bi = B(i);
    
    for(int j = 0; j < sizeA; j++) {
      double Aj = A(j);
      double Bj = B(j);
      
      double aux = pow(Ai-Aj,2) +  pow(Bi-Bj,2) ;
      result(i,j) = aux;
      
    }}
  
  return result;
};


// [[Rcpp::export]]

NumericMatrix gath(NumericMatrix A){
  
  int n_col = A.ncol();
  int n_row = A.nrow();
  int n_flat = n_col*n_row - n_col;
  
  NumericVector ix(n_flat);
  NumericVector jx(n_flat);
  NumericVector result(n_flat);
  int m = 0;
  
  for(int i = 0; i < n_row; i++) {
    for(int j = 0; j < n_col; j++) {
      
      if(i != j) {
        double aux = A(i, j);
        result(m) = aux;
        ix(m) = i+1;
        jx(m) = j+1;
        
        
        m = m+1;}
      
    }}
  
  NumericMatrix final(n_flat,3);
  final(_,0) = ix;
  final(_,1) = jx;
  final(_,2) = result;
  
  return final;
};

// [[Rcpp::export]]
NumericMatrix gath_tri(NumericMatrix A){
  
  int n_col = A.ncol();
  int n_row = A.nrow();
  int n_flat = (n_col*n_row - n_col)/2;
  
  NumericVector ix(n_flat);
  NumericVector jx(n_flat);
  NumericVector result(n_flat);
  int m = 0;
  
  for(int i = 0; i < n_row; i++) {
    for(int j = i; j < n_col; j++) {
      
      if(i != j){
        double aux = A(i, j);
        result(m) = aux;
        ix(m) = i+1;
        jx(m) = j+1;
        
        
        m = m+1;}
      
    }}
  
  NumericMatrix final(n_flat,3);
  final(_,0) = ix;
  final(_,1) = jx;
  final(_,2) = result;
  
  return final;
};






// [[Rcpp::export]]

DataFrame caggH(NumericVector ab_logl , NumericVector ab_logl_prop ,  NumericVector pois,  NumericVector pois_prop,  int N ){
  
  int len = ab_logl.length();
  int Time = len/N;
  
  NumericMatrix result(N,5);
  
  for(int i = 0; i < N; i++) {
    
    
    NumericVector idx = cseq(i,  (N*Time-1), N);
    
    NumericVector ab_logl_sub = ab_logl[idx];
    NumericVector  ab_logl_prop_sub = ab_logl_prop[idx];
    NumericVector pois_sub = pois[idx];
    NumericVector pois_prop_sub = pois_prop[idx];
    
    result(i, 0)=  i + 1  ;
    result(i, 1)= sum(ab_logl_sub)  ;
    result(i, 2)= sum(ab_logl_prop_sub)   ;
    result(i, 3)= sum(pois_sub)  ;
    result(i, 4)= sum(pois_prop_sub)   ;
    
  };
  
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["Sum_ab_logl"] = result(_, 1), _["Sum_ab_logl_prop"] = result(_, 2), _["Sum_pois"] = result(_, 3),   _["Sum_pois_prop"] = result(_, 4));
  return df;
  
}



// [[Rcpp::export]]

DataFrame caggPlane(NumericVector beta_grad , NumericVector c_grad ,  int N ){
  
  int len = beta_grad.length();
  int Time = len/N;
  
  NumericMatrix result(N,3);
  
  for(int i = 0; i < N; i++) {
    
    
    NumericVector idx = cseq(i,  (N*Time-1), N);
    
    NumericVector beta_grad_sub = beta_grad[idx];
    NumericVector  c_grad_sub = c_grad[idx];
    
    
    result(i, 0)=  i + 1  ;
    result(i, 1)= sum(beta_grad_sub)  ;
    result(i, 2)= sum(c_grad_sub)   ;
    
    
  };
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["Sum_beta_grad"] = result(_, 1), _["Sum_c_grad"] = result(_, 2));
  return df;
  
}


// [[Rcpp::export]]

DataFrame caggEdge(NumericVector lambda_grad_1 , NumericVector lambda_grad_2 ,  int N ){
  
  int len = lambda_grad_1.length();
  int Time = len/(N*N - N);
  
  NumericMatrix result(N,3);
  
  
  int k = 0;
  
  for(int t = 0; t < Time; t++) {
    
    int count = 0;
    
    for(int i = 0; i < N; i++) {
      
      
      NumericVector idx =cseq(0,  (N-2), 1) + count + k ;
      
      NumericVector lambda_grad_1_sub = lambda_grad_1[idx];
      NumericVector lambda_grad_2_sub = lambda_grad_2[idx];
      
      
      result(i, 0)=  i + 1  ;
      result(i, 1)= result(i, 1) + sum(lambda_grad_1_sub)  ;
      result(i, 2)=  result(i, 2) + sum(lambda_grad_2_sub)   ;
      count = count + N -1;
      
    };
    
    k =   k + (N*N -N);
    
  }
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["Sum_lambda_grad_1"] = result(_, 1), _["Sum_lambda_grad_2"] = result(_, 2));
  return df;
  
}



// [[Rcpp::export]]

DataFrame caggG(NumericVector pois , NumericVector pois_prop ,  int N ){
  
  
  int len = pois.length();
  int Time = len/(N*N - N);
  
  NumericMatrix result(N,3);
  
  int k = 0;
  
  for(int t = 0; t < Time; t++) {
    
    int count = 0;
    
    for(int i = 0; i < N; i++) {
      
      
      NumericVector idx = cseq(0,  (N-2), 1) + count + k ;
      
      NumericVector pois_sub = pois[idx];
      NumericVector pois_prop_sub = pois_prop[idx];
      
      
      result(i, 0)=  i + 1  ;
      result(i, 1)=  result(i, 1) + sum(pois_sub)  ;
      result(i, 2)=  result(i, 2) + sum(pois_prop_sub)   ;
      count = count + N -1;
      
    };
    
    k =   k + (N*N - N);
  }
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["Sum_pois"] = result(_, 1), _["Sum_pois_prop"] = result(_, 2));
  return df;
  
};


// [[Rcpp::export]]
DataFrame caggProdG(NumericVector g_a , NumericVector g_b ,  int M , int Times){
  
  
  NumericMatrix result(Times,3);
  
  int count = 0;
  
  for(int t = 0; t < Times; t++) {
    
    NumericVector idx = cseq(0,  (M-1), 1) + count ;
    
    NumericVector g_a_sub = g_a[idx];
    NumericVector g_b_sub = g_b[idx];
    
    
    result(t, 0)=  t + 1  ;
    result(t, 1)=   sum(g_a_sub)  ;
    result(t, 2)=   sum(g_b_sub)   ;
    count = count + M;
    
  }
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["prodG_a"] = result(_, 1), _["prodG_b"] = result(_, 2));
  return df;
  
};



// [[Rcpp::export]]
DataFrame caggProdH(NumericVector B_a , NumericVector B_b, NumericVector K_a , NumericVector K_b,   int N , int Times){
  
  
  NumericMatrix result(Times,5);
  
  int count = 0;
  
  for(int t = 0; t < Times; t++) {
    
    NumericVector idx = cseq(0,  (N-1), 1) + count ;
    
    NumericVector B_a_sub = B_a[idx];
    NumericVector B_b_sub = B_b[idx];
    NumericVector K_a_sub = K_a[idx];
    NumericVector K_b_sub = K_b[idx];
    
    
    result(t, 0)=  t + 1  ;
    result(t, 1)=   sum(B_a_sub)  ;
    result(t, 2)=   sum(B_b_sub)   ;
    result(t, 3)=   sum(K_a_sub)  ;
    result(t, 4)=   sum(K_b_sub)   ;
    
    count = count + N;
    
  }
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["prodB_a"] = result(_, 1), _["prodB_b"] = result(_, 2), _["prodK_a"] = result(_, 3), _["prodK_b"] = result(_, 4) );
  return df;
  
};




// [[Rcpp::export]]
NumericVector lgradF1(NumericVector dist1 , NumericVector distance,  NumericVector w , NumericVector wl1 ,  NumericVector ones , NumericVector  alpha, NumericVector beta, double Time){
  
  NumericVector alpha_ij =  rep(alpha, Time);
  NumericVector beta_ij =  rep(beta, Time);
  
  NumericVector result = dist1*(-1*w +  exp(ones*alpha_ij + beta_ij*wl1 - distance))/distance;
  
  return result;
}



// [[Rcpp::export]]
NumericVector lgradF2(NumericVector dist2 , NumericVector distance,  NumericVector w , NumericVector wl1 ,  NumericVector ones , NumericVector  alpha, NumericVector beta, int Time){
  
  NumericVector alpha_ij =  rep(alpha, Time);
  NumericVector beta_ij =  rep(beta, Time);
  
  NumericVector result = dist2*(-1*w +  exp(ones*alpha_ij + beta_ij*wl1 - distance))/distance;
  
  return result;
}



// [[Rcpp::export]]
NumericMatrix ziPropSampling(NumericMatrix zi,  double lambda ,  int N){
  
  NumericMatrix zi_prop(N,2) ;
  
  for(int i = 0; i < N; i++) {
    
    
    NumericVector mu =  {0,0};
    mu[0]  = zi(i,0);
    mu[1]  = zi(i,1);
    
    arma::mat sigma = arma::eye(2,2);
    double lm = lambda;
    
    arma::mat aux =  mvrnormArma(1, mu, lm*sigma)  ;
    zi_prop(i, 0 ) =  aux(0,0);
    zi_prop(i, 1 ) =  aux(0,1);
  };
  
  
  return zi_prop;
  
};




// [[Rcpp::export]]
NumericMatrix ziPropSamplingAd(NumericMatrix zi,  NumericVector lambda , List Sigma,   int N){
  
  NumericMatrix zi_prop(N,2) ;
  
  for(int i = 0; i < N; i++) {
    
    
    NumericVector mu =  {0,0};
    mu[0]  = zi(i,0);
    mu[1]  = zi(i,1);
    // 
    // mu[0]  = theta[0];
    // mu[1]  = theta[1];
    
    arma::mat sigma = Sigma[i];
    double lm = lambda[i];
    
    arma::mat aux =  mvrnormArma(1, mu, lm*sigma)  ;
    zi_prop(i, 0 ) =  aux(0,0);
    zi_prop(i, 1 ) =  aux(0,1);
  };
  
  
  return zi_prop;
  
};


// [[Rcpp::export]]


NumericMatrix procrustes_cpp(NumericMatrix X, NumericMatrix Y){
  
  
  arma::mat X_aux = as<arma::mat>(X);
  arma::mat Y_aux = as<arma::mat>(Y);
  arma::mat XY = X_aux.t()*Y_aux;
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  arma::svd(U,s,V,XY);
  
  arma::mat A = V*U.t();
  arma::mat result = Y_aux*A;
  
  
  return wrap(result);
}




#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
};





// [[Rcpp::export]]

NumericMatrix ZetaSamplingAd(NumericVector zi_s_1 , NumericVector  zi_s_2,  NumericVector lam_ad_zs,  NumericMatrix mu_mat_zs, List Sigma_ad_zs, NumericVector alpha_ij_doub, NumericVector beta_ij_doub,  NumericVector xi_state_s , NumericVector x_w,  NumericVector x_wl1, NumericVector leaning, NumericVector ncomments,  NumericVector x_ones , NumericVector state_i,NumericVector mu_s, double sigma_s, double gamma_0, double gamma_1, double phi, double tau, int M, int N){
  
  int Time_s = sum(xi_state_s);
  
  NumericVector state_x =  rep_each(xi_state_s, 2*M) ;
  NumericVector x_w_s =  x_w[ state_x > 0];
  NumericVector x_wl1_s =  x_wl1[ state_x > 0];
  NumericVector x_ones_s =  x_ones[ state_x > 0];
  NumericVector alpha_ij_doub_s = rep(alpha_ij_doub, Time_s);
  NumericVector beta_ij_doub_s = rep(beta_ij_doub, Time_s);
  NumericVector leaning_s = leaning[state_i > 0];
  NumericVector ncomments_s = ncomments[state_i > 0];
  
  NumericMatrix eucl_dist = ecldist(zi_s_1, zi_s_2 );
  NumericMatrix eucl =  gath(eucl_dist);
  NumericVector dist_s = eucl(_,2);
  NumericVector distan_s = rep(dist_s, Time_s) ;
  
  NumericMatrix zi_s(N, 3);
  zi_s(_ , 0) = zi_s_1;
  zi_s(_ , 1) = zi_s_2;
  
  // double eps = 0.002;
  // NumericMatrix zi_s_prop  =  ziPropSampling( zi_s ,  eps ,  N);
  NumericMatrix zi_s_prop  =  ziPropSamplingAd( zi_s ,  lam_ad_zs,  Sigma_ad_zs,  N);
  NumericVector zi_s_prop_1 =  zi_s_prop(_,0) ;
  NumericVector zi_s_prop_2 =  zi_s_prop(_,1) ;
  
  NumericMatrix eucl_prop_dist_s = ecldist(zi_s_prop_1, zi_s_prop_2 );
  NumericMatrix  eucl_prop_s =  gath(eucl_prop_dist_s);
  NumericVector dist_prop_s = eucl_prop_s(_ ,2 );
  NumericVector distan_prop_s = rep(dist_prop_s, Time_s);
  
  NumericVector lam_s = exp(x_ones_s*alpha_ij_doub_s  + x_wl1_s*beta_ij_doub_s - distan_s );
  NumericVector lam_prop_s = exp(x_ones_s*alpha_ij_doub_s  + x_wl1_s*beta_ij_doub_s - distan_prop_s );
  
  NumericVector zi_s_1_i_s =     rep(zi_s_1, Time_s);
  NumericVector zi_s_2_i_s =     rep(zi_s_2, Time_s);
  
  NumericVector zi_s_1_i_s_prop =     rep(zi_s_prop_1, Time_s);
  NumericVector zi_s_2_i_s_prop =     rep(zi_s_prop_2, Time_s);
  
  NumericVector DB_plane_s_a =  logistic(gamma_0+ gamma_1*zi_s_1_i_s)*phi;
  NumericVector DB_plane_s_b = (1-logistic(gamma_0+ gamma_1*zi_s_1_i_s))*phi;
  
  NumericVector DB_plane_s_a_prop =  logistic(gamma_0+ gamma_1*zi_s_1_i_s_prop)*phi;
  NumericVector DB_plane_s_b_prop = (1-logistic(gamma_0+ gamma_1*zi_s_1_i_s_prop))*phi;
  
  NumericVector DB_plane_s_lambda_like = tau*exp(zi_s_2_i_s);
  NumericVector DB_plane_s_lambda_like_prop = tau*exp(zi_s_2_i_s_prop);
  
  NumericVector DB_plane_s_ab_logl(Time_s*N) ;
  NumericVector  DBplane_s_ab_logl_prop(Time_s*N) ;
  NumericVector DB_plane_s_pois(Time_s*N) ;
  NumericVector  DBplane_s_pois_prop(Time_s*N) ;
  
  
  for(int i = 0; i < Time_s*N; i++) {
    
    
    DB_plane_s_ab_logl[i] =  R::dbeta(leaning_s[i], DB_plane_s_a[i], DB_plane_s_b[i]  , true);
    DBplane_s_ab_logl_prop[i] = R::dbeta(leaning_s[i], DB_plane_s_a_prop[i], DB_plane_s_b_prop[i]  , true);
    
    DB_plane_s_pois[i] =  R::dpois(ncomments_s[i], DB_plane_s_lambda_like[i] , true);
    DBplane_s_pois_prop[i] = R::dpois(ncomments_s[i], DB_plane_s_lambda_like_prop[i] , true);
    
  };
  
  DataFrame  sumH_s = caggH( DB_plane_s_ab_logl , DBplane_s_ab_logl_prop ,  DB_plane_s_pois , DBplane_s_pois_prop , N);
  
  NumericVector Sum_ab_logl_s =  sumH_s["Sum_ab_logl"];
  NumericVector Sum_ab_logl_prop_s =  sumH_s["Sum_ab_logl_prop"];
  NumericVector Sum_pois_s =  sumH_s["Sum_pois"];
  NumericVector Sum_pois_prop_s =  sumH_s["Sum_pois_prop"];
  
  NumericVector x_s_pois(Time_s*M*2) ;
  NumericVector x_s_pois_prop(Time_s*M*2) ;
  
  
  for(int i = 0; i < Time_s*M*2; i++) {
    
    
    x_s_pois[i] =  R::dpois(x_w_s[i], lam_s[i] , true);
    x_s_pois_prop[i] = R::dpois(x_w_s[i], lam_prop_s[i] , true);
    
  };
  
  
  DataFrame sumG_s = caggG(x_s_pois, x_s_pois_prop,  N);
  
  NumericVector SumG_pois_s =  sumG_s["Sum_pois"];
  NumericVector SumG_pois_prop_s =  sumG_s["Sum_pois_prop"];
  
  NumericVector aceptance_zs(N);
  
  for(int i = 0; i < N; i++) {
    
    
    double prior  =  R::dnorm(zi_s_prop_1[i] , mu_s[0] , sigma_s, true ) - R::dnorm(zi_s_1[i] , mu_s[0] , sigma_s, true ) +   R::dnorm(zi_s_prop_2[i] , mu_s[1] , sigma_s, true ) - R::dnorm(zi_s_2[i] , mu_s[1] , sigma_s, true );
    double logr_s =  prior  + Sum_ab_logl_prop_s[i] - Sum_ab_logl_s[i] + Sum_pois_prop_s[i] - Sum_pois_s[i] +  SumG_pois_prop_s[i] - SumG_pois_s[i];
    
    double runif_s =    log(R::runif(0,1))  ;
    
    if(exp(logr_s) < 1 ){ aceptance_zs[i] = exp(logr_s); }else{aceptance_zs[i] = 1;};
    
    
    if(logr_s > runif_s){
      double zi_s_aux_1  =  zi_s_prop_1[i];
      double zi_s_aux_2  =   zi_s_prop_2[i];
      
      zi_s_1[i]  = zi_s_aux_1;
      zi_s_2[i]  = zi_s_aux_2; }
    
  };
  
  
  
  zi_s(_ , 0) = zi_s_1;
  zi_s(_ , 1) = zi_s_2;
  zi_s(_,2) =  aceptance_zs; 
  
  return zi_s;
  
  
}


// [[Rcpp::export]]
arma::mat rdirichlet_cpp(int num_samples,
                         arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::mat distribution = arma::zeros(num_samples, distribution_size);
  
  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
      double cur = R::rgamma(alpha_m[j],1.0);
      distribution(i,j) = cur;
      sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
      distribution(i,j) = distribution(i,j)/sum_term;
    }
  }
  return(distribution);
};


// [[Rcpp::export]]
NumericMatrix sampleTransP(NumericVector xi_state1, NumericVector xi_state2, double omega_lower_a, double omega_lower_b, NumericMatrix P, double Time){
  
  NumericVector xi_a_a(Time-1) ;
  NumericVector xi_b_b(Time-1) ;
  NumericVector xi_a_b(Time-1) ;
  NumericVector xi_b_a(Time-1) ;  
  
  for (int i = 0; i < (Time -1); ++i) {
    
    xi_a_a[i] = xi_state1[i]* xi_state1[i+1];  
    xi_b_b[i] = xi_state2[i]* xi_state2[i+1];
    
    xi_a_b[i] = xi_state2[i]*xi_state1[i+1];
    xi_b_a[i]= xi_state1[i]*xi_state2[i+1];
    
  }
  
  double omega_upper_a_a = omega_lower_a -1 + sum(xi_a_a);
  double omega_upper_b_a = omega_lower_a -1 + sum(xi_b_a);
  
  double omega_upper_b_b = omega_lower_b -1 + sum(xi_b_b);
  double omega_upper_a_b = omega_lower_a -1 + sum(xi_a_b);
  
  
  arma::mat p1 = rdirichlet_cpp(1, {omega_upper_a_a, omega_upper_b_a});
  arma::mat p2 = rdirichlet_cpp(1, {omega_upper_a_b, omega_upper_b_b});
  
  
  P(0,0) = p1(0,0);
  P(1,0) = p1(0,1);
  
  P(0,1) = p2(0,0);
  P(1,1) = p2(0,1);
  
  return P;
  
}


// [[Rcpp::export]]

List ZetaSamplingAdC(NumericVector zi_s_1 , NumericVector  zi_s_2,  NumericVector lam_ad_zs,  NumericMatrix mu_mat_zs, List Sigma_ad_zs, NumericVector beta,  NumericVector xi_state_s , NumericVector x_i, NumericVector DBplane_i ,  NumericVector x_w, NumericVector leaning, NumericVector ncomments,  NumericVector x_ones ,NumericVector mu_s, double sigma_s, double gamma_0, double gamma_1, double phi, int M, int N, double acc_star, double ite){
  
  int Time_s = sum(xi_state_s);
   zi_s_2 = rep(0, N);
  
  
  // NumericVector state_x =  rep_each(xi_state_s, 2*M) ; // len: 380 x 100 = 380000
  // NumericVector state_i =  rep_each(xi_state_s, N); //len: 100x20
  // 
  NumericVector x_i_s =  clone(x_i); // len: 2M x Ts
  NumericVector x_w_s =  clone(x_w); // len: 2M x Ts
  NumericVector x_ones_s =  clone(x_ones); // len: 2M x Ts
  // 
  
  NumericVector DBplane_i_s = clone(DBplane_i); // 20 x Ts
  NumericVector leaning_s = clone(leaning);  // 20 x Ts
  NumericVector ncomments_s = clone(ncomments);  // 20 x Ts
  
  
  NumericMatrix zi_s(N, 3);
  
  // zi_s_1[0] = 0;
  // zi_s_2[0] = 0;
  // 
  // zi_s_2[8] = -2.5 ;
  // 
  // 
  
  zi_s(_ , 0) = zi_s_1;
  zi_s(_ , 1) = zi_s_2;
  
  
  
  NumericVector zi_s_prop_1(N);
  NumericVector zi_s_prop_2(N);
  
  NumericVector range = cseq(0, N-1, 1) ; //19
  NumericVector sample_v  = sample( range,  N, false);
  
  for(int k = 0; k < N; k++) {
    
    //double i = sample_v[k];
    double i = k;
    
    double beta_i = beta[i];
    NumericVector beta_ni = clone(beta) ;
    beta_ni.erase(i);
    
    NumericVector dist_s = ecldist_point(zi_s_1, zi_s_2, i ); // 19
    NumericVector distan_s = rep(dist_s, Time_s ); // 19 x T_s
    
    NumericVector mu =  {0,0};
    double mu_xx  = zi_s_1[i];
    double mu_yy  = zi_s_2[i];
    
    mu[0]  = mu_xx;
    mu[1]  = mu_yy;
    
    arma::mat sigma = Sigma_ad_zs[i];
    double lm = lam_ad_zs[i];
    
    // 
    double zi_prop_x = 0;
    double zi_prop_y = 0;
    //   
    // arma::mat EYE = arma::eye(2,2);
    // 
    // // int ctrl = 1;
    // // do{   
    // 
    //  if(ite < 5000){
    // //   
    //  zi_prop_x = R::runif( -1,1);
    //  zi_prop_y = R::runif( -1,1);
    // //   
    // }else{
    
    //arma::mat aux =  mvrnormArma(1, mu, 0.0001*EYE );
    
    arma::mat  aux =  mvrnormArma(1, mu, sigma*lm )  ; 
    zi_prop_x = aux(0,0);
    zi_prop_y = aux(0,1);
    zi_prop_y = 0;
    
    
    // 
    // }
    
    
    //
    // 
    
    //cout << aux;
    
    // if( abs(zi_prop_x) < 1 && abs(zi_prop_y) < 1){ ctrl = 1;}else{  ctrl = 0 ;};  } while (ctrl < 1);
    
    NumericVector  zi_s_prop_1 = clone(zi_s_1);
    NumericVector zi_s_prop_2  =  clone(zi_s_2);
     zi_s_prop_2  =  rep(0, N);
    
    
    // if(i == 0){zi_prop_x = 0;};
    
    zi_s_prop_1[i] = zi_prop_x;
    zi_s_prop_2[i] = zi_prop_y;
    
    NumericVector dist_s_prop = ecldist_point(zi_s_prop_1, zi_s_prop_2, i ); //19
    NumericVector distan_s_prop = dist_s_prop; // 19xT_s
    
    double idx = i +1;
    NumericVector x_ones_s_i = x_ones_s[x_i_s == idx]; //19x72
    NumericVector x_w_s_i =  x_w_s[x_i_s == idx]; //19x72
    
    NumericVector x_s_pois(Time_s*(N-1)) ;
    NumericVector x_s_pois_prop(Time_s*(N-1)) ;
    
    for(int j = 0; j < Time_s*(N-1); j++) {
      
      double lam_s_aux = exp(beta_i+ beta_ni[j] - pow(distan_s[j],2) );
      double lam_s_prop_aux = exp(beta_i + beta_ni[j] - pow(distan_s_prop[j],2));
      
      double x_s_pois_aux =   R::dpois(x_w_s_i[j], lam_s_aux , true);
      double x_s_pois_prop_aux = R::dpois(x_w_s_i[j], lam_s_prop_aux , true);
      
      x_s_pois[j] =  x_s_pois_aux;
      x_s_pois_prop[j] = x_s_pois_prop_aux;
      
    }
    
    double SumG_pois_s = sum(x_s_pois);
    double SumG_pois_prop_s =  sum(x_s_pois_prop);
    
    double DB_plane_s_a =  logistic_d(gamma_0+ abs(gamma_1)*mu[0])*phi;
    double DB_plane_s_b = (1-logistic_d(gamma_0+ abs(gamma_1)*mu[0]))*phi;
    
    double DB_plane_s_a_prop =  logistic_d(gamma_0+ abs(gamma_1)*zi_prop_x)*phi;
    double DB_plane_s_b_prop = (1- logistic_d(gamma_0+ abs(gamma_1)*zi_prop_x))*phi;
    
    
    NumericVector DB_plane_s_ab_logl(Time_s) ;
    NumericVector  DBplane_s_ab_logl_prop(Time_s) ;
    
    NumericVector leaning_s_i = leaning_s[DBplane_i_s == idx ]; //T_s
    NumericVector ncomments_s_i = ncomments_s[DBplane_i_s == idx]; //T_s
    
    
    for(int j = 0; j < Time_s; j++) {
      
      DB_plane_s_ab_logl[j] =  R::dbeta(leaning_s_i[j], DB_plane_s_a, DB_plane_s_b  , true);
      DBplane_s_ab_logl_prop[j] = R::dbeta(leaning_s_i[j], DB_plane_s_a_prop, DB_plane_s_b_prop  , true);
      
    };
    
    double Sum_ab_logl_s = sum(DB_plane_s_ab_logl);
    double Sum_ab_logl_prop_s = sum(DBplane_s_ab_logl_prop);
    
    double prior  =  R::dnorm(zi_s_prop_1[i] , mu_s[0] , sigma_s, true ) - R::dnorm(mu[0] , mu_s[0] , sigma_s, true ); //+   R::dnorm(zi_s_prop_2[i] , mu_s[1] , sigma_s, true ) - R::dnorm(mu[1] , mu_s[1] , sigma_s, true );
    double logr_s =  prior  + (Sum_ab_logl_prop_s - Sum_ab_logl_s) +  (SumG_pois_prop_s - SumG_pois_s);
    double runif_s =    log(R::runif(0,1))  ;
    
    double aceptance_zs;
    
    if( logr_s < 0 ){ 
      // if(Rcpp::traits::is_infinite<REALSXP>(logr_s)){
      //   logr_s = -1000;}
      aceptance_zs = exp(logr_s); }else{
        aceptance_zs = 1; };
      
      
      if(logr_s > runif_s){
        zi_s_1[i]  = zi_prop_x;
        zi_s_2[i]  = zi_prop_y;
      }else{
        
        zi_s_1[i]  = mu_xx;
        zi_s_2[i]  = mu_yy;
        
      };
      // 
      // zi_s_1[0] = 0;
      // zi_s_2[0] = 0;
      // 
      // zi_s_2[8] = -2.5 ;
      
      
      zi_s(_ , 0) = zi_s_1;
      zi_s(_ , 1) = zi_s_2;
      zi_s(i,2) =  aceptance_zs;
      
      NumericMatrix zi_s_sub(N, 2);
      zi_s_sub(_ , 0) = zi_s_1;
      zi_s_sub(_ , 1) = zi_s_2;
      
      // 
      double gamm = 1/pow(ite + 2, 0.55 );
      
      double log_lam  = log(lm) + gamm*(aceptance_zs -  acc_star);
      
      NumericVector theta_ij_i = zi_s_sub(i,_);
      
      NumericVector mu_mat_i = mu_mat_zs(i,_);
      
      
      NumericVector diff =  (theta_ij_i -  mu_mat_i );
      NumericVector mu_x = mu_mat_i + gamm*diff;
      
      mu_mat_zs(i,_) = mu_x;
      
      arma::vec dif_vec = diff;
      arma::rowvec dif_vec_t = dif_vec.t();
      arma::mat Sigma = sigma;
      
      NumericMatrix Sigma_up =  wrap(Sigma + gamm*(dif_vec*dif_vec_t - Sigma ) );
      
      Sigma_ad_zs[i] = Sigma_up;
      
      lam_ad_zs[i] = exp(log_lam);
      
      // 
  };
  
  List result = List::create( Named("zi_s")= zi_s, Named("lam_ad_zs") = lam_ad_zs,  Named("mu_mat_zs") = mu_mat_zs, Named("Sigma_ad_zs") = Sigma_ad_zs);
  return result;
  // 
  
}



// [[Rcpp::export]]

NumericVector vec_nbet(NumericVector beta, double N){
  
  NumericVector res;
  
  for( int i = 0; i < N; i++){
    
    NumericVector  bb = clone(beta);
    bb.erase(i);
    
    for( int j = 0; j < (N-1); j++){
    
    double tt = bb[j];
    res.push_back(tt);}
    
  }
  
  return(res);
  
}



// [[Rcpp::export]]

NumericVector vec_nbet_tri(NumericVector beta, double N){
  
  NumericVector res;
  
  for( int i = 0; i < (N-1); i++){
    
    NumericVector  bb = clone(beta);
    
    for( int j = i + 1; j < (N); j++){
      
      double tt = bb[j];
      res.push_back(tt);
      
    }
    
  }
  
  return(res);
  
}



// [[Rcpp::export]]

NumericVector vec_nbet_tri_a(NumericVector beta, double N){
  
  NumericVector res;
  
  for( int i = 0; i < (N-1); i++){
    
    NumericVector  bb = rep(beta(i), N-i-1);
    double len = bb.length();
    
    for( int j = 0 ; j < len; j++){
      
      double tt = bb(j);
      res.push_back(tt);
      
    }
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]

NumericVector  betaSamplingEn(NumericVector vec_y,  NumericVector  beta,  NumericVector vec_dist, double beta_a,  double beta_b,  double N, double prop_sigma) {
  
  for(int i = 0; i < N; i++) {
    
    int ind = i;
    
    NumericVector beta_prop = clone(beta);
    double beta_x = beta[ind];
    double beta_prop_x = R::rnorm(beta_x,prop_sigma );
    beta_prop[ind] = beta_prop_x;
    
    int nn = vec_y.length();
    
    NumericVector vec_beta = rep_each(beta, N-1);
    NumericVector vec_beta_prop = rep_each(beta_prop, N-1);
    
    NumericVector vec_nbeta =  vec_nbet( beta,  N);
    
    NumericVector pois(nn);
    NumericVector pois_prop(nn);
    
    
    for(int k = 0; k < nn; k++) {
      
      double lambda_aux = exp(vec_beta[k] + vec_nbeta[k]  - pow(vec_dist[k],2));
      double lambda_aux_prop = exp(vec_beta_prop[k]  + vec_nbeta[k] - pow(vec_dist[k],2));
      
      double pois_aux = R::dpois(vec_y[k], lambda_aux, true);
      double pois_prop_aux =  R::dpois(vec_y[k], lambda_aux_prop, true);
      
      pois[k]  = pois_aux;  
      pois_prop[k]  = pois_prop_aux;
      
    }
    
    double prior = R::dnorm( beta_x, beta_a, beta_b, true );
    double prior_prop = R::dnorm( beta_prop_x, beta_a, beta_b, true );
    
    double test = prior_prop - prior + sum(pois_prop) - sum(pois) ;
    double rand = R::runif(0,1);
    if(test >  log(rand)){  beta =  beta_prop;  } ;
    
  }
  
  return(beta);
};




// [[Rcpp::export]]
double rinvgamma(double shape, double scale) {

    //GetRNGstate();
    double x  = R::rgamma( shape, 1/scale);
    //PutRNGstate();
    
    return(1/x);
}


// [[Rcpp::export]]
double invSamplingEn(NumericVector zi_s_1 ,  double a_prior_ig , double b_prior_ig, int N){
  
 
 double a_star = a_prior_ig + N/2;
 double b_star = b_prior_ig + sum(pow(zi_s_1,2))/2;
 
  double sigma = rinvgamma(a_star, b_star) ;

  return(sqrt(sigma));
  
  
  };






// [[Rcpp::export]]
List MCMC(NumericVector princ_w,  NumericVector princ_ones, NumericVector x_w,  NumericVector x_ones, NumericVector mu_ad,  double lm, arma::mat sigma_ad, NumericVector leaning, NumericVector ncomments ,  NumericVector lam_ad_za, NumericMatrix  mu_mat_za,  List Sigma_ad_za,  NumericVector lam_ad_zb, NumericMatrix  mu_mat_zb,  List Sigma_ad_zb,  double alpha, NumericVector beta, NumericVector xi_state1, NumericVector xi_state2  , NumericVector zi_a_1, NumericVector zi_a_2,  NumericVector zi_b_1, NumericVector zi_b_2 ,  double mu_alpha, double mu_beta, double sigma_alpha, double sigma_beta , NumericVector mu_a ,  double sigma_a,  NumericVector mu_b ,  double sigma_b, double phi, double gamma_0,  double gamma_1,  double a_phi, double b_phi, double a_gamma_0, double b_gamma_0 , double a_gamma_1, double b_gamma_1,   double omega_lower_a, double omega_lower_b, NumericMatrix P,  int M, int N, int Time  ,
          NumericVector x_i, NumericVector DBplane_i, int Iterations){
  
  NumericMatrix th0_ite(Iterations, 3);
  NumericMatrix phi_gamma_tau_ite(Iterations, 5);
  
  NumericMatrix zi_a_ite(Iterations*N, 4);
  NumericMatrix zi_b_ite(Iterations*N, 4);
  
  NumericMatrix P_ite(Iterations, 5);
  NumericMatrix xi_ite(Iterations*Time, 4);
  NumericMatrix HETA(Iterations, 1);
  NumericMatrix beta_ite(Iterations, N);
  NumericMatrix sigma_ite(Iterations, 2);

  
  zi_a_2 = rep(0, N);
  
  int index_z = 0;
  //int index_x = 0;
  
  
  // double ite = 0;
  //
  //   // percentage charging
  
  for(double ite = 0; ite < Iterations; ite++) {
    double prct = (ite+1)/Iterations;
    printProgress(prct);
    
    try {
      
      
      // start sampling
      
      NumericMatrix eucl_a_mat =  ecldist(zi_a_1, zi_a_2 );
      NumericMatrix eucl_a = gath(eucl_a_mat);
      NumericVector distance_a = rep(eucl_a(_,2) , Time);
      
      beta =  betaSamplingEn(x_w,    beta,   distance_a, mu_beta,  sigma_beta,   N,  0.02);
      beta_ite(ite, _ ) = beta;
        
        //phi = 50;
        
        phi = phiSamplingEn(leaning, zi_a_1, phi, gamma_0, gamma_1, a_phi, b_phi,  3 );//3
      
      // phi = 50;
      //
      NumericVector gamma = gammaSamplingEn(leaning, zi_a_1, phi, gamma_0, gamma_1,  a_gamma_0, b_gamma_0 ,  a_gamma_1, b_gamma_1,  0.15);
      
      gamma_0 = gamma[0];
      gamma_1 = gamma[1];
      
      NumericVector phi_gamma_tau_aux = {phi, gamma_0, gamma_1, 0, ite +1};
      
      phi_gamma_tau_ite(ite,_ ) = phi_gamma_tau_aux;
      
      // 
      List ZAobj = ZetaSamplingAdC(zi_a_1, zi_a_2,  lam_ad_za,   mu_mat_za, Sigma_ad_za,  beta,  xi_state1 , x_i, DBplane_i , x_w ,   leaning ,  ncomments , x_ones ,mu_a, sigma_a, gamma_0, gamma_1, phi, M,  N, 0.20, ite);
      
      NumericMatrix ziRes_a =   ZAobj["zi_s"];
      
      zi_a_1 = ziRes_a(_, 0);
      zi_a_2 = ziRes_a(_, 1);
      
      NumericVector acceptance_za = ziRes_a(_, 2);
      
      NumericVector aux_vec_za =  ZAobj["lam_ad_zs"];
      lam_ad_za =  aux_vec_za;
      NumericMatrix aux_mat_za = ZAobj["mu_mat_zs"];
      mu_mat_za =   aux_mat_za;
      List aux_list_za = ZAobj["Sigma_ad_zs"];
      Sigma_ad_za = aux_list_za;
      
      
      zi_a_1 = zi_a_1 - mean(zi_a_1);
      zi_a_2 = zi_a_2 - mean(zi_a_2);
      
      //
      double zi_a_1x = zi_a_1[ 36 ];

      if(zi_a_1x < 0){
        
        double mean_za_1 = mean(zi_a_1);
        
        NumericVector zi_a_1_c_rf = -1*(zi_a_1 - mean_za_1);
        
        zi_a_1 = zi_a_1_c_rf ;
        
      }
      
      
      for(int i = 0; i < N; i++) {
        
        double it = i;
        double zi_a_1_ax =  zi_a_1[i];
        double zi_a_2_ax =  zi_a_2[i];
        
        double zi_b_1_ax = zi_b_1[i];
        double zi_b_2_ax = zi_b_2[i];
        
        NumericVector aux_vec_za = { zi_a_1_ax, zi_a_2_ax, it + 1,  ite +1  };
        NumericVector aux_vec_zb = { zi_b_1_ax, zi_b_2_ax, it + 1 ,  ite +1  };
        
        zi_a_ite(i + index_z, _ ) =  aux_vec_za;
        zi_b_ite(i + index_z, _ ) =  aux_vec_zb;
        
      }
      
      index_z = index_z + N;
      
      double sigma_a =  invSamplingEn(zi_a_1 ,   s_a ,  s_b,  N);
            sigma_ite(ite, 0) = sigma_a;

     eucl_a_mat =  ecldist(zi_a_1, zi_a_2 );
     eucl_a = gath_tri(eucl_a_mat);
   
     distance_a = rep(eucl_a(_,2) , Time);
   
 
    NumericVector vec_beta_t = vec_nbet_tri_a(beta, N);
    NumericVector vec_nbeta_t =  vec_nbet_tri( beta,  N);

     
     NumericVector lambda_a = exp(vec_beta_t + vec_nbeta_t - distance_a*distance_a);
   
     NumericVector g_a(M);
   
     for(int i = 0; i < M; i++) {
       g_a(i) = R::dpois(princ_w(i), lambda_a(i), true );   
     }
   
    // DataFrame ProdG =  caggProdG( g_a ,  g_a ,   M , Time);
   
   
    NumericVector zi_a_1_i =     rep(zi_a_1, Time);
   
     NumericVector DB_plane_a_a =  logistic(gamma_0+ abs(gamma_1)*zi_a_1_i)*phi;
     NumericVector DB_plane_b_a = (1-logistic(gamma_0+ abs(gamma_1)*zi_a_1_i))*phi;
   
   
     NumericVector B_a(N);
   
      
     for(int i = 0; i < N; i++) {
       B_a(i) = R::dbeta(leaning(i), DB_plane_a_a(i), DB_plane_b_a(i)  , true);
     }
   
   
    // DataFrame ProdH = caggProdH( B_a ,  B_a,  K_a ,  K_a,    N,   Time);
      
   
     double likelihood_a =  sum(g_a)  + sum(B_a);
   
    HETA(ite,0) = sum(likelihood_a);

      
    }catch(std::runtime_error & e) {std::cout << "Exception caught!" << std::endl;
      forward_exception_to_r(e); break;};
    
  }
  
  
  
  List RESULT = List::create( Named("beta_it")= beta_ite,   Named("phi_gamma_tau_it") = phi_gamma_tau_ite, Named("zi_a_it") = zi_a_ite, Named("lLik") = HETA,   Named("Sigma_z_ite") =  sigma_ite);
  return RESULT;
}


////////////


// 


// [[Rcpp::export]]
NumericVector shuffling(NumericVector x, int size){ NumericVector aux= Rcpp::sample( x,  size, false);
  return(aux);
}


// [[Rcpp::export]]

NumericMatrix procrustes_preprocessing(NumericMatrix ref_position_c, NumericVector ref_position_c_mean , NumericVector result_x ,  NumericVector result_y ,   int N, int iterations ){
  
  NumericMatrix RESULT(N*iterations, 2);
  NumericVector v_x(N*iterations);
  NumericVector v_y(N*iterations);
  int idx = 0;
  
  for(int ite = 0; ite < iterations; ite++) {
    
    
    // double prct = (ite+1)/iterations;
    // printProgress(prct);
    
    NumericVector index =  cseq(0 + idx, (N-1) +idx ,  1);
    
    NumericVector position_x = result_x[index];
    NumericVector position_y = result_y[index];
    
    position_x =  position_x - mean(position_x);
    position_y =  position_y - mean(position_y);
    
    NumericMatrix position(N,2) ;
    
    position(_,0) = position_x;
    position(_,1) = position_y;
    
    position = procrustes_cpp(ref_position_c , position);
    position_x =  position(_,0);
    position_y = position(_,1);
    
    position_x = position_x + ref_position_c_mean[0];
    position_y = position_y + ref_position_c_mean[1];
    
    v_x[index] = position_x;
    v_y[index] = position_y;
    
    idx = idx + N;
    Rcout << "The iteration is : " << ite << "\n";
    
  }
  
  
  RESULT(_,0) = v_x;
  RESULT(_,1) = v_y;
  
  return RESULT;
  
}
