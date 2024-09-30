
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <list>
#include <RcppDist.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends( RcppDist , RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;



Function f("rnorm"); 
Environment pkg = Environment::namespace_env("stats");
Function rmulti = pkg["rmultinom"];


// [[Rcpp::export]]
NumericMatrix dist_fcpp(NumericMatrix A, NumericMatrix B)
{
  int n = A.nrow();
  int m = B.nrow();
  int K = A.ncol();
  NumericMatrix result(n, m);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {

        double aux = 0;
        for (int k = 0; k < K; k++)
        {
          double Ai = A(i, k);
          double Bj = B(j, k);
          aux = aux + pow(Ai - Bj, 2);
        }
        result(i, j) = sqrt(aux);
    }
  }
  return result;
};

// [[Rcpp::export]]
NumericMatrix dist_lcpp(NumericMatrix A, NumericMatrix B)
{
  int n = A.nrow();
  int m = B.nrow();
  int K = A.ncol();
  NumericMatrix result(n, m);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      if (j > i)
      {
        double aux = 0;
        for (int k = 0; k < K; k++)
        {
          double Ai = A(i, k);
          double Bj = B(j, k);
          aux = aux + pow(Ai - Bj, 2);
        }
        result(i, j) = sqrt(aux);
      }
    }
  }
  return result;
};

// [[Rcpp::export]]
double avg_dist_lcpp(NumericMatrix A, NumericMatrix B)
{
  int n = A.nrow();
  int m = B.nrow();
  int K = A.ncol();
  
  double sum_result;
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      if (j > i)
      {
        double aux = 0;
        for (int k = 0; k < K; k++)
        {
          double Ai = A(i, k);
          double Bj = B(j, k);
          aux = aux + pow(Ai - Bj, 2);
        }
         sum_result = sum_result + sqrt(aux);
      }
    }
  }
  
  double result = sum_result/((n/2.00)*(n-1));
  return result;
};


// [[Rcpp::export]]
NumericVector dist_i_lcpp(NumericMatrix A, NumericMatrix B, int i)
{
  int n = A.nrow();
  int m = B.nrow();
  int K = A.ncol();
  NumericVector result(n-1);
  
  int idx = 0;

    for (int j = 0; j < m; j++)
    {
      if (j != i)
      {
        double aux = 0;
        for (int k = 0; k < K; k++)
        {
          double Ai = A(i, k);
          double Bj = B(j, k);
          aux = aux + pow(Ai - Bj, 2);
        }
        result(idx) = sqrt(aux);
        idx = idx + 1;
      }
    }
  
  return result;
};


// [[Rcpp::export]]
NumericMatrix dist2_lcpp(NumericMatrix A, NumericMatrix B)
{
  int n = A.nrow();
  int m = B.nrow();
  int K = A.ncol();
  NumericMatrix result(n, m);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      if (j > i)
      {
        double aux = 0;
        for (int k = 0; k < K; k++)
        {
          double Ai = A(i, k);
          double Bj = B(j, k);
          aux = aux + pow(Ai - Bj, 2);
        }
        result(i, j) = (aux);
      }
    }
  }
  return result;
};



// [[Rcpp::export]]
NumericVector dist(NumericVector x) {
  int n = x.size();
  NumericVector out(n * (n - 1) / 2);
  int k = 0;
  
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      double d = 0;
        d += pow(x[i] - x[j], 2);
      out[k] = sqrt(d);
      k++;
    }
  }
  
  return unique(out);
}

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
double median_mat(NumericMatrix A){ 
  
  double n_i =  A.nrow() ;
  double n_j =  A.ncol();
  double m = (n_i/2)*(n_i-1);
  
  NumericVector res(m);
  double idx = 0;
  
  for (int i = 0; i < n_i; i++) {
    for (int j= i+1; j < n_j; j++ ) {
      
      double aux = A(i,j);
      
      res(idx) = aux;
      idx = idx +1;
      
    }}
  
  double med = median(res);
  return(med);
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


// [[Rcpp::export]]

NumericVector  gammaSamplingEn(NumericVector leaning, NumericVector  zi_1, double phi, double gamma_0, double gamma_1, 
                               double a_gamma_0, double a_gamma_1, double b_gamma_0, double b_gamma_1,  double interp_eq, NumericVector prop_sigma) {
  

  NumericVector  gamma(2);
  
  gamma[0] = gamma_0;
  gamma[1] = gamma_1*interp_eq;
  
  
  arma::mat propdiag(2,2,fill::zeros);

  propdiag(0,0) = prop_sigma[0];
  propdiag(1,1) = prop_sigma[1];
  
  NumericVector gamma_prop(2);
  GetRNGstate();
  double gamma_prop_0 = R::rnorm(gamma_0, prop_sigma[0]);
  PutRNGstate();
  
  GetRNGstate();
  double gamma_prop_1 = R::rnorm(gamma_1, prop_sigma[1]);
  PutRNGstate();

   gamma_prop(0) = gamma_prop_0;
   gamma_prop(1) = gamma_prop_1*interp_eq;
  
  
  // double gamma_prop_1 = 0;
  // 
  // int ctrl = 1;
  // do{  gamma_prop_1 = R::rnorm(gamma_1, prop_sigma[1] ); if(gamma_prop_1 <= 0){ ctrl = 0;}else{  ctrl = 1 ;};  } while (ctrl < 1);
  
  
   //arma::mat  gamma_prop_arma = mvrnormArma(1, gamma, propdiag );
   //NumericVector  gamma_prop = as<NumericVector>(wrap(gamma_prop_arma)); 
  
  // 
  // if (gamma_prop[0] > gamma_prop[1] ) { double aux =  gamma_prop[1] ;   gamma_prop[0] =    gamma_prop[1];  gamma_prop[1] = aux ; }
  // 
  
  NumericVector  a_beta = phi*logistic(gamma_0 +interp_eq*(gamma_1)*zi_1);
  NumericVector  b_beta = phi*(1-logistic(gamma_0 +interp_eq*(gamma_1)*zi_1));
  // 
  NumericVector  a_beta_prop =  phi*logistic(gamma_prop[0] + interp_eq*(gamma_prop[1])*zi_1);
  NumericVector  b_beta_prop =  phi*(1-logistic(gamma_prop[0] +  interp_eq*(gamma_prop[1])*zi_1));
  
  double size = b_beta_prop.length();
  NumericVector  Beta(size);
  NumericVector Beta_prop(size);
  
  for(int t = 0; t< (size); t++) {
    Beta[t] = R::dbeta( leaning[t],  a_beta[t], b_beta[t],  true);
    Beta_prop[t] = R::dbeta( leaning[t], a_beta_prop[t],  b_beta_prop[t], true);}
  
  double prior = R::dnorm( gamma_0, a_gamma_0, b_gamma_0, true ) + (R::dnorm( gamma_1, a_gamma_1, b_gamma_1, true ))*interp_eq;
  double prior_prop = R::dnorm( gamma_prop[0], a_gamma_0, b_gamma_0, true ) + (R::dnorm( gamma_prop[1], a_gamma_1, b_gamma_1, true ))*interp_eq;
  
  // double norm_c_g1 = R::pnorm(gamma_1/sqrt(prop_sigma[1]), 0, 1 ,  true,  true) ;
  // double norm_c_g1_prop = R::pnorm(gamma_prop[1]/sqrt(prop_sigma[1]), 0, 1 ,  true,  true) ;
  
  // double norm_c_g0 = R::pnorm(gamma_0/sqrt(prop_sigma[0]), 0, 1 ,  true,  true) ;
  // double norm_c_g0_prop = R::pnorm(gamma_prop[0]/sqrt(prop_sigma[0]), 0, 1 ,  true,  true) ;
  
  double test = sum(Beta_prop) - sum(Beta)  + prior_prop - prior  ;
  //cout << prior << "  " << gamma_0 <<  "  " << a_gamma_0 << "  " << b_gamma_0 << "  " << gamma_1 << "  " << a_gamma_1 << "  " << b_gamma_1;
  
  // +  norm_c - norm_c_prop ;
  GetRNGstate();
  double rand = R::runif(0,1);
  PutRNGstate();
  
  if(test >  log(rand)){  return gamma_prop;  }else{ return gamma; }  ;
}
 


// [[Rcpp::export]]

double  phiSamplingEn(NumericVector leaning, NumericVector  zi_1, double phi, double gamma_0, double gamma_1, 
                      double a_phi, double b_phi, double interp_eq , double prop_sigma) {
  
  double phi_prop = -1;
  while(phi_prop <= 0){
    GetRNGstate();
    phi_prop = R::rnorm(phi,prop_sigma );
    PutRNGstate();
  };
  
  NumericVector  a_beta = phi*logistic(gamma_0 + interp_eq*(gamma_1)*zi_1);
  NumericVector  b_beta = phi*(1-logistic(gamma_0 + interp_eq*(gamma_1)*zi_1));
  
  NumericVector  a_beta_prop =  phi_prop*logistic(gamma_0 +interp_eq*(gamma_1)*zi_1);
  NumericVector  b_beta_prop =  phi_prop*(1-logistic(gamma_0 +interp_eq*(gamma_1)*zi_1));
  
  double size = b_beta_prop.length();
  NumericVector  Beta(size);
  NumericVector Beta_prop(size);
  
  for(int t = 0; t< (size); t++) {
    Beta[t] = R::dbeta( leaning[t],  a_beta[t], b_beta[t],  true);
    Beta_prop[t] = R::dbeta( leaning[t], a_beta_prop[t],  b_beta_prop[t], true);}
  
  // gamma shape and scale
  double prior = R::dgamma(phi,a_phi,1/b_phi, true);
  double prior_prop = R::dgamma(phi_prop,a_phi,1/b_phi, true);
  
  double norm_c = R::pnorm(phi/prop_sigma, 0, 1 ,  true,  true) ;
  double norm_c_prop = R::pnorm(phi_prop/prop_sigma, 0, 1 ,  true,  true) ;
  
  
  double test = sum(Beta_prop) - sum(Beta) + prior_prop - prior +  norm_c - norm_c_prop ;
  GetRNGstate();
  double rand = R::runif(0,1);
  PutRNGstate();
  
  if(test > log(rand)){  return phi_prop;  }else{ return phi; }  ;
  
}


// [[Rcpp::export]]
NumericVector cseq(double first, double last, double by){
  
  int n = (last - first)/by + 1;
  NumericVector result(n);
  result(0) = first;
  
  
  for(int i = 1; i< n ; i++) {
    
    result(i) = result(i-1) + by;
    
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
NumericMatrix ecldist(NumericVector A){
  
  int sizeA = A.length();
  
  NumericMatrix result(sizeA, sizeA);
  
  for(int i = 0; i < sizeA; i++) {
    double Ai = A(i);
    
    for(int j = 0; j < sizeA; j++) {
      double Aj = A(j);
      
      double aux = sqrt(pow(Ai-Aj,2) ) ;
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
DataFrame cagglogG(NumericVector g_a , NumericVector g_b ,  int M , int Times){
  
  
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
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["logG_a"] = result(_, 1), _["logG_b"] = result(_, 2));
  return df;
  
};



// [[Rcpp::export]]
DataFrame cagglogH(NumericVector B_a , NumericVector B_b,   int N , int Times){
  
  
  NumericMatrix result(Times,5);
  
  int count = 0;
  
  for(int t = 0; t < Times; t++) {
    
    NumericVector idx = cseq(0,  (N-1), 1) + count ;
    
    NumericVector B_a_sub = B_a[idx];
    NumericVector B_b_sub = B_b[idx];
    
    
    
    result(t, 0)=  t + 1  ;
    result(t, 1)=   sum(B_a_sub)  ;
    result(t, 2)=   sum(B_b_sub)   ;
    
    
    count = count + N;
    
  }
  
  DataFrame df = DataFrame::create( Named("i") = result(_, 0) , _["logB_a"] = result(_, 1), _["logB_b"] = result(_, 2) );
  return df;
  
};




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
arma::mat rdirichlet_cpp(int num_samples,
                         arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::mat distribution = arma::zeros(num_samples, distribution_size);
  
  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
      GetRNGstate() ;
      double cur = R::rgamma(alpha_m(j),1.0);
      PutRNGstate();
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

NumericVector vec_nbet(NumericVector x, double N){
  
  NumericVector res;
  NumericVector x_aux = clone(x);
  
  for( int i = 0; i < (N); i++){
    
    NumericVector bb = clone(x_aux);
    bb.erase(i);
    NumericVector bb_aux  = bb;
    
    
    for( int j = 0; j < (N-1); j++){
      
      double tt = bb_aux(j);
      res.push_back(tt);}
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]

List ZetaSamplingAdC(NumericMatrix zi_s, NumericVector lam_ad_zs,  NumericMatrix mu_mat_zs, List Sigma_ad_zs, NumericVector beta,  NumericMatrix xi_state , NumericVector x_i, NumericVector DBplane_i ,  NumericVector x_w,  NumericVector leaning, NumericVector x_ones , NumericVector mu_s, NumericVector sigma_state, double gamma_0, double gamma_1, double phi, int M, int N, int D, int K, int k_state, double acc_star, double interp_eq, double ite){
  
  double sigma_s = sigma_state(k_state);
  NumericVector xi_state_s = xi_state(_, k_state);
  int Time_s = sum(xi_state_s);
  NumericVector state_x =  rep_each(xi_state_s, 2*M) ; // len: 380 x 100 = 380000
  NumericVector state_i =  rep_each(xi_state_s, N); //len: 100x20
  
  NumericVector vec_beta_a = rep_each(beta, N-1 );
  NumericVector vec_beta_a_s = rep(vec_beta_a, Time_s);
  
  NumericVector vec_beta_b = vec_nbet(beta, N);
  NumericVector vec_beta_b_s = rep(vec_beta_b, Time_s);
  
  NumericVector x_i_s =  x_i[ state_x > 0]; // len: 2M x Ts
  NumericVector x_w_s =  x_w[ state_x > 0]; // len: 2M x Ts
  NumericVector x_ones_s =  x_ones[ state_x > 0]; // len: 2M x Ts
  
  
  NumericVector DBplane_i_s = DBplane_i[state_i > 0]; // 20 x Ts
  NumericVector leaning_s = leaning[state_i > 0];  // 20 x Ts
  
  
  NumericMatrix zi_s_prop(N, D);
  
  NumericVector range = cseq(0, N-1, 1) ; //19
  NumericVector sample_v  = sample( range,  N, false);
  // int k = 0;
  
  for(int k = 0; k < N; k++) {
    
    double i = sample_v[k];
    //double i = k;
    
    NumericVector dist_s = dist_i_lcpp(zi_s, zi_s, i ); // 19
    NumericVector distan_s = rep(dist_s, Time_s ); // 19 x T_s
    
    NumericMatrix zi_s_prop_aux = zi_s;
    NumericMatrix zi_s_prop = clone(zi_s_prop_aux);
    NumericVector mean_zs = zi_s_prop(i,_);

    NumericMatrix sigma = as<NumericMatrix>(Sigma_ad_zs[i]);
    double lm = lam_ad_zs[i];
    
    NumericVector z_i_x(D);
    
    for (int d = 0; d < D; d++)
    {
      
      double sigma_d = sigma(d, d);
      GetRNGstate();
      double z_i_prop_x = R::rnorm(mean_zs(d), sigma_d * lm);
      PutRNGstate();
      z_i_x(d) = z_i_prop_x;
      
    }
    
    zi_s_prop(i,_) = clone(z_i_x);

    NumericVector dist_s_prop = dist_i_lcpp(zi_s_prop, zi_s_prop, i ); //19
    NumericVector distan_s_prop = rep(dist_s_prop, Time_s); // 19xT_s
    
    double idx = i +1;
    NumericVector x_ones_s_i = x_ones_s[x_i_s == idx]; //19x72
    NumericVector x_w_s_i =  x_w_s[x_i_s == idx]; //19x72
    
    NumericVector x_s_pois(Time_s*(N-1)) ;
    NumericVector x_s_pois_prop(Time_s*(N-1)) ;
    
    NumericVector vec_beta_a_s_i =  vec_beta_a_s[x_i_s == idx];
    NumericVector vec_beta_b_s_i =  vec_beta_b_s[x_i_s == idx];
    
    for(int j = 0; j < Time_s*(N-1); j++) {
      
      double beta_a = vec_beta_a_s_i[j];
      double beta_b = vec_beta_b_s_i[j];
      double n_distan_s = distan_s[j];
      double n_distan_s_prop = distan_s_prop[j];
      
      double lam_s_aux = exp(beta_a + beta_b -  pow(n_distan_s,2) );
      double lam_s_prop_aux = exp(beta_a + beta_b    - pow(n_distan_s_prop,2) );
      
      double x_s_pois_aux =      R::dpois(x_w_s_i[j], lam_s_aux , true);
      double x_s_pois_prop_aux = R::dpois(x_w_s_i[j], lam_s_prop_aux , true);
      
      x_s_pois[j] =  x_s_pois_aux;
      x_s_pois_prop[j] = x_s_pois_prop_aux;
      
    }
    
    double SumG_pois_s = sum(x_s_pois);
    double SumG_pois_prop_s =  sum(x_s_pois_prop);
    

    double DB_plane_s_a =  logistic_d(gamma_0+interp_eq*(gamma_1)*mean_zs(0))*phi;
    double DB_plane_s_b = (1-logistic_d(gamma_0+interp_eq*(gamma_1)*mean_zs(0)))*phi;
    
    double DB_plane_s_a_prop =  logistic_d(gamma_0+interp_eq*(gamma_1)*z_i_x(0))*phi;
    double DB_plane_s_b_prop = (1- logistic_d(gamma_0+interp_eq*(gamma_1)*z_i_x(0)))*phi;
    
    NumericVector DB_plane_s_ab_logl(Time_s) ;
    NumericVector  DBplane_s_ab_logl_prop(Time_s) ;
    
    NumericVector leaning_s_i = leaning_s[DBplane_i_s == idx ]; //T_s
    
    
    for(int j = 0; j < Time_s; j++) {
      
      DB_plane_s_ab_logl[j] =  R::dbeta(leaning_s_i[j], DB_plane_s_a, DB_plane_s_b  , true);
      DBplane_s_ab_logl_prop[j] = R::dbeta(leaning_s_i[j], DB_plane_s_a_prop, DB_plane_s_b_prop  , true);
      
    };
    
    double Sum_ab_logl_s = sum(DB_plane_s_ab_logl);
    double Sum_ab_logl_prop_s = sum(DBplane_s_ab_logl_prop); 
    
    double prior = 0;
    double prior_prop = 0;
    
    for (int d = 0; d < D; d++)
    {
      double prior_aux = R::dnorm(mean_zs(d), mu_s(d), sigma_s, true);
      double prior_prop_aux = R::dnorm(z_i_x(d), mu_s(d), sigma_s, true);
      prior = prior + prior_aux;
      prior_prop = prior_prop + prior_prop_aux;
    }    
    
    
    double logr_s =  (prior_prop - prior)  + (Sum_ab_logl_prop_s - Sum_ab_logl_s) +  (SumG_pois_prop_s - SumG_pois_s);
    GetRNGstate();
    double runif_s =    log(R::runif(0,1))  ;
    PutRNGstate();
    double aceptance_zs;
    if(exp(logr_s) < 1 ){ aceptance_zs = exp(logr_s); }else{aceptance_zs = 1;};
    
    if(logr_s > runif_s){
      zi_s(i,_ )  = clone(z_i_x);
    }
    

    double gamm = 1/pow(ite + 1, 0.52 );
    
    double log_lam  = log(lm) + gamm*(aceptance_zs -  acc_star);
    
    NumericVector theta_ij_i = zi_s(i, _);
    NumericVector mu_mat_i = mu_mat_zs(i, _);
    NumericVector diff = (theta_ij_i - mu_mat_i);
    
    NumericVector mu_x = mu_mat_i + gamm * diff;
    
    mu_mat_zs(i, _) = mu_x;
    
    arma::vec dif_vec = diff;
    arma::rowvec dif_vec_t = dif_vec.t();
    
    arma::mat Sigma = as<arma::mat>(sigma);
    
    NumericMatrix Sigma_up = wrap(Sigma + gamm * (dif_vec * dif_vec_t - Sigma));
    
    Sigma_ad_zs(i) = Sigma_up;
    lam_ad_zs(i) = exp(log_lam);

  };
  
  List result = List::create( Named("zi_s")= zi_s, Named("lam_ad_zs") = lam_ad_zs,  Named("mu_mat_zs") = mu_mat_zs, Named("Sigma_ad_zs") = Sigma_ad_zs);
  return result;

}




// [[Rcpp::export]]

NumericVector vec_nbet_tri(NumericVector beta, double N){
  
  NumericVector res;
  
  for( int i = 0; i < (N-1); i++){
    
    NumericVector  bb = rep(beta,1);
    
    for( int j = i + 1; j < (N); j++){
      
      double tt = bb(j);
      res.push_back(tt);
      
    }
    
  }
  
  return(res);
  
}


NumericVector my_pois(double lambda){
  // calling rnorm()
  Function f("rpois");   
  
  // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
  return f(1,lambda);
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

NumericVector  betaSamplingEn(NumericVector vec_y,  NumericVector  beta,  NumericVector vec_dist, double beta_a,  double beta_b,  double N, double Time, NumericVector x_i,double prop_sigma) {
  
  NumericVector range = cseq(0, N-1, 1) ; //19
  NumericVector sample_v  = sample( range,  N, false);
  
  for(int i = 0; i < N; i++) {
    
    int ind = sample_v[i];
    
    
    
    
    NumericVector beta_prop = clone(beta);
    double beta_x = beta[ind];
    GetRNGstate();
    double beta_prop_x = R::rnorm(beta_x,prop_sigma );
    PutRNGstate();
    beta_prop[ind] = beta_prop_x;
    
    int nn = (N-1)*Time;
    
    double idx = ind +1;
    
    
    
    NumericVector vec_beta = rep(beta_x, N-1);
    NumericVector vec_beta_prop = rep(beta_prop_x, N-1);
    
    NumericVector vec_beta_t_i = rep(vec_beta, Time);
    NumericVector vec_beta_prop_t_i = rep(vec_beta_prop,Time);
    
    NumericVector beta_aux = clone(beta);
    beta_aux.erase(ind);
    
    NumericVector vec_nbeta =  beta_aux ;
    NumericVector vec_nbeta_t_i =  rep( vec_nbeta,  Time);
    
    NumericVector vec_y_i =  vec_y[x_i == idx]; //19x72
    NumericVector vec_dist_i =  vec_dist[x_i == idx]; //19x72
    
    NumericVector pois(nn);
    NumericVector pois_prop(nn);
    
    
    for(int k = 0; k < nn; k++) {
      
      double lambda_aux = exp(vec_beta_t_i(k) + vec_nbeta_t_i(k)  - pow(vec_dist_i(k),2) );
      double lambda_aux_prop = exp(vec_beta_prop_t_i(k)  + vec_nbeta_t_i(k) - pow(vec_dist_i(k),2));
      
      double pois_aux = R::dpois(vec_y_i(k), lambda_aux, true);
      double pois_prop_aux =  R::dpois(vec_y_i(k), lambda_aux_prop, true);
      
      
      //if(arma::is_finite(pois_aux)){}else{pois_aux = -1*pow(10, 26);};
      //if(arma::is_finite(pois_aux)){}else{pois_prop_aux = -1*pow(10, 26);};
      
      pois(k)  = pois_aux;  
      pois_prop(k)  = pois_prop_aux;
      
    }
    
    
    double prior = R::dnorm( beta_x, beta_a, beta_b, true );
    double prior_prop = R::dnorm( beta_prop_x, beta_a, beta_b, true );
    
    double test = prior_prop - prior + sum(pois_prop) - sum(pois) ;
    
    GetRNGstate();
    double rand = R::runif(0,1);
    PutRNGstate();
    
    if(test >  log(rand)){  beta =  clone(beta_prop);  } ;
    
  }
  
  
  return(beta);
  

};



// [[Rcpp::export]]

List  betaSamplingEnAdC(NumericVector vec_y,  NumericVector lam_ad_beta,  NumericMatrix  mu_mat_beta, List Sigma_ad_beta,  NumericVector  beta,  NumericVector vec_dist, double beta_a,  double beta_b,  double N, double Time, NumericVector x_i,double acc_star, double ite) {
  
  NumericVector range = cseq(0, N-1, 1) ; //19
  NumericVector sample_v  = sample( range,  N, false);
  
  for(int i = 0; i < N; i++) {
    
    int ind = sample_v[i];
    
    arma::mat sigma = Sigma_ad_beta[ind];
    double lm = lam_ad_beta[ind];
    
    NumericVector beta_prop = clone(beta);
    double beta_x = beta[ind];
    double sig_b = sigma(0,0);
    GetRNGstate();
    double beta_prop_x = R::rnorm(beta_x, sig_b*lm );
    PutRNGstate();
    beta_prop[ind] = beta_prop_x;
    
    int nn = (N-1)*Time;
    
    double idx = ind +1;
    
    
    
    NumericVector vec_beta = rep(beta_x, N-1);
    NumericVector vec_beta_prop = rep(beta_prop_x, N-1);
    
    NumericVector vec_beta_t_i = rep(vec_beta, Time);
    NumericVector vec_beta_prop_t_i = rep(vec_beta_prop,Time);
    
    NumericVector beta_aux = clone(beta);
    beta_aux.erase(ind);
    
    NumericVector vec_nbeta =  beta_aux ;
    NumericVector vec_nbeta_t_i =  rep( vec_nbeta,  Time);
    
    NumericVector vec_y_i =  vec_y[x_i == idx]; //19x72
    NumericVector vec_dist_i =  vec_dist[x_i == idx]; //19x72
    
    NumericVector pois(nn);
    NumericVector pois_prop(nn);
    
    
    for(int k = 0; k < nn; k++) {
      
      double lambda_aux = exp(vec_beta_t_i(k) + vec_nbeta_t_i(k)  - pow(vec_dist_i(k),2) );
      double lambda_aux_prop = exp(vec_beta_prop_t_i(k)  + vec_nbeta_t_i(k) - pow(vec_dist_i(k),2));
      
      double pois_aux = R::dpois(vec_y_i(k), lambda_aux, true);
      double pois_prop_aux =  R::dpois(vec_y_i(k), lambda_aux_prop, true);
      
      
      //if(arma::is_finite(pois_aux)){}else{pois_aux = -1*pow(10, 26);};
      //if(arma::is_finite(pois_aux)){}else{pois_prop_aux = -1*pow(10, 26);};
      
      pois(k)  = pois_aux;  
      pois_prop(k)  = pois_prop_aux;
      
    }
    
    
    double prior = R::dnorm( beta_x, beta_a, beta_b, true );
    double prior_prop = R::dnorm( beta_prop_x, beta_a, beta_b, true );
    
    double test = prior_prop - prior + sum(pois_prop) - sum(pois) ;
    
    GetRNGstate();
    double rand = R::runif(0,1);
    PutRNGstate();
    
    if(test >  log(rand)){  beta =  clone(beta_prop);  } ;
    
    
   double aceptance_beta;
    if(exp(test) < 1 ){ aceptance_beta = exp(test); }else{aceptance_beta = 1;};
    
    
    double gamm = 1/pow(ite + 1, 0.52 );
    
    double log_lam  = log(lm) + gamm*(aceptance_beta -  acc_star);
    
    double theta_ij_i =  beta[ind];
    double mu_mat_i = mu_mat_beta(ind,0);
    double diff =  (theta_ij_i -  mu_mat_i );
    double mu_x = mu_mat_i + gamm*diff;
    mu_mat_beta(ind,0) = mu_x;
   
   double sigma_i = sigma(0,0);
   double Sigma_up_i =  sigma_i + gamm*(diff*diff - sigma_i);
   sigma(0,0) = Sigma_up_i;
   
   Sigma_ad_beta[ind] = sigma;
    
   lam_ad_beta[ind] = exp(log_lam);
    
    
    
  }
  
  
  
  List result = List::create( Named("beta")= beta, Named("lam_ad_beta") = lam_ad_beta,  Named("mu_mat_beta") = mu_mat_beta, Named("Sigma_ad_beta") = Sigma_ad_beta);
  return result;
  

};


// [[Rcpp::export]]
double rinvgamma(double shape, double scale) {

    GetRNGstate();
    double x  = R::rgamma( shape, 1/scale);
    PutRNGstate();
    
    return(1/x);
}


// [[Rcpp::export]]
double invSamplingEn(NumericMatrix zi_s ,  double a_prior_ig , double b_prior_ig, double N, double D){
  
 double a_star = a_prior_ig + (N/2)*D;
  
  // cout << a_star << endl;
  
  double sum_v = sum(pow(zi_s,2));
  
 double b_star = b_prior_ig + sum_v/2;
 
  double sigma = rinvgamma(a_star, b_star) ;

  return(sqrt(sigma));
  
  };


// [[Rcpp::export]]
NumericMatrix sampleTransP(NumericMatrix S_int, NumericVector omega ,  double Time, double K){
  
  NumericMatrix S = clone(S_int);
  NumericMatrix P(K,K);
  
  for (int i = 0; i < K; i++) {
    
    NumericVector Count(K);
    Count = rep(0, K);
    
    for (int t = 0; t < (Time-1); t++) {
      
      NumericVector St = S(t,_);
      double St_i = St(i);
      
      NumericVector Stp1 = S(t+1,_);
      
      for (int j = 0; j < K; j++) {
        
        double Stp1_j = Stp1(j);
        
        double aux_j = Count(j);
        
        aux_j = aux_j + St_i*Stp1_j;
        
        Count(j) = aux_j;
        
      }
    }
    
    NumericVector omega_star = omega + Count;
    arma::vec omega_star_arma = as<arma::vec>(omega_star);
    
    arma::mat p_i_arma = rdirichlet_cpp(1, omega_star);
    NumericVector p_i  = wrap(p_i_arma);
    
    P(_, i) = p_i;
    
  }
  
  
  return P;
  
}


// [[Rcpp::export]]
NumericMatrix FFBS_v1( NumericMatrix Eta, NumericMatrix P_in, double T, double K) {
  
  // Eta is a numeric matrix representing emission probabilities.
  // P is a numeric matrix representing transition probabilities.
  // T is a double representing the number of time steps.
  // K is a double representing the number of states.
  
  arma::mat Pt = as<arma::mat>(P_in);
  arma::mat P_arma = Pt.t();
  NumericMatrix P = wrap(P_arma);
  
  NumericVector Acc(K);
  NumericVector Bcc(K);
  NumericMatrix alef(T, K);
  
  NumericVector x10 = rep(1.00/K, K);
  NumericVector eta0 = Eta(0,_);
  alef(0,_) = log(x10) + eta0;
  
  NumericMatrix bet(T, K);
  NumericVector one_aux = rep(1.00, K);
  bet(T - 1, _) = one_aux;
  
  // It computes the forward probabilities alef using the forward recursion.
  
  for(int t = 1; t< T; t++) {
    for(int j = 0; j< K; j++) {
      for(int i = 0; i < K; i++){
        
        Acc(i)  = alef(t-1, i) + log(P(i, j)) + Eta(t,j);
        
        
      }
      
      alef(t,j) = clog_sum_exp(Acc);
    }
  }
  
  
  for(int t = 0; t< T; t++) {
    alef.row(t) = csoftmax(alef.row(t));
  }
  
  // It computes the backward probabilities bet using the backward recursion.
  
  for(int t = (T-1); t > 0; t--) {
    for(int j = 0; j< K; j++) {
      for(int i = 0; i < K; i++){
        
        Bcc(i)  = bet(t, i) + log(P(j, i)) +  Eta(t,i);
        
      }
      
      bet(t-1, j) = clog_sum_exp(Bcc);
      
    }}
  
  // 
  for(int t = 0; t< T; t++) {
    
    
    bet(t,_) = csoftmax(bet(t,_));
    
  }
  
  // It computes the joint probabilities ghimel.
  NumericMatrix ghimel(T,K);
  
  for(int k = 0; k < K; k++) {
    
    NumericVector aux1 = alef(_,k);
    NumericVector aux2 = bet(_,k);
    
    ghimel(_, k) = aux1*aux2;
    
  }
  
  //States Sampling
  NumericMatrix Xstate(T,K);
  
  // It samples states using the joint probabilities.
  for(int t = 0; t < T; t++) {
    
    NumericVector aux_i = ghimel(t,_);
    NumericVector aux = aux_i/sum(aux_i);
    
    GetRNGstate();
    NumericVector res = rmulti(1,1, aux);
    PutRNGstate();
    
    Xstate(t, _) = res;
    
  }
  
  return Xstate;
}

// [[Rcpp::export]]
NumericMatrix FFBS( NumericMatrix Eta, NumericMatrix P_in, double T, double K) {
  
  // Eta is a numeric matrix representing emission probabilities.
  // P is a numeric matrix representing transition probabilities.
  // T is a double representing the number of time steps.
  // K is a double representing the number of states.
  
  arma::mat Pt = as<arma::mat>(P_in);
  arma::mat P_arma = Pt.t();
  NumericMatrix P = wrap(P_arma);
  
  // NumericMatrix P = clone(P_in);
  NumericVector Acc(K);
  NumericVector Bcc(K);
  NumericMatrix alef(T, K);
  
  // \alpha_1(x_1 ) = p(x_1) * p(y_1 | x_1) [[ initialization ]]
  
  NumericVector p10 = rep(1/K, K);
  NumericVector eta0 = Eta(0,_);
  alef(0,_) = log(p10) + eta0;
  
  
  // It computes the forward probabilities alef using the forward recursion.
  
  // \alpha_t(x_t = j ) = p(y_t | x_t = j ) * \sum_{i=1}^N p(x_{t-1} = i, y_1, ..., y_{t-1}) * p(x_t  = j| x_{t-1} = i)
  //   = p(y_t | x_t = j) * \sum_{i=1}^N \alpha_t(x_t-1 = i) * p(x_t = j | x_{t-1} = i)    for t=2, .., T
  
  for(int t = 1; t< T; t++) {
    for(int j = 0; j< K; j++) {
      
      for(int i = 0; i < K; i++){
        
        Acc(i)  = alef(t-1, i) + log(P(i, j)) + Eta(t,j);
        
      }
      
      alef(t,j) = clog_sum_exp(Acc);
    }
  }
  
  
  for(int t = 0; t< T; t++) {
    alef.row(t) = csoftmax(alef.row(t));
  }
  
  // Backward Sampling
  
  NumericMatrix Xstate(T,K);
  
  NumericVector alef_T = alef(T-1,_);
  
  GetRNGstate();
  NumericVector res = rmulti(1,1, alef_T);
  PutRNGstate();
  
  Xstate(T-1, _) = res;
  
  for(int t = (T-2); t >= 0; t--) {
    
    NumericVector alefp1 = alef(t+1,_);
    NumericVector alef_current = alef(t,_);
    NumericVector SamplDistr(K);
    
    NumericVector state = Xstate(t+1, _);
    
    NumericVector range = cseq(0, K-1,1);

    NumericVector jth_v = range[state == 1];
    
    double jth = jth_v(0);
    
    // p(x_{t-1} = i | x_t = j, y_{1:T}) = p(x_{t-1} = i | y_{1:t-1}) * p(x_t = j | x_{t-1} = i) / Z_t
    // where Z_t = \sum_{i=1}^N p(x_{t-1} = i | y_{1:t-1}) * p(x_t = j | x_{t-1} = i)
    
    for(int i = 0; i < K; i++){
      double aux =  P(i, jth)*alef_current(i);
      SamplDistr(i) = aux;
    }
    
    SamplDistr = SamplDistr/sum(SamplDistr); //propto in Murphy's Book
    
    GetRNGstate();
    NumericVector res = rmulti(1,1, SamplDistr);
    PutRNGstate();
    
    Xstate(t, _) = clone(res);
    
  }
  
  return Xstate;
}


// [[Rcpp::export]]
NumericVector  PosteriorPredictive( NumericVector beta, NumericVector zi_a_1 , NumericVector zi_b_1,  NumericVector xi_state1 , double N, double Time){
  
  NumericVector EY(Time);
  NumericVector VarY(Time);
  NumericVector DiY(Time);
  arma::vec v_ones = arma::colvec(N,  fill::ones);
  
  NumericVector state = clone(xi_state1);
  
  for(int t = 0; t < Time; t++) {
    
    arma::mat Q = arma::mat(N,N,  fill::zeros);
    
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        
        if(j > i){
          double chk = state(t);
          
          if(chk == 1){
            
            double lambda_y  = exp(beta(i)  + beta(j) - pow(zi_a_1(i) - zi_a_1(j),2));
            
            GetRNGstate();
            double rp =  R::rpois(lambda_y);
            PutRNGstate();
            
            Q(i,j)  = rp;
            
          }else{
            
            double lambda_y  = exp(beta(i) + beta(j) - pow(zi_b_1(i) - zi_b_1(j),2));
            GetRNGstate();
            double rp =  R::rpois(lambda_y);
            PutRNGstate();
            Q(i,j)  = rp;
          }
          
        }
        
      }
    }
    
    Q  = Q+Q.t();
    
    arma::colvec vr = Q*v_ones;
    // cout << vr;
    Rcpp::NumericVector scr = Rcpp::NumericVector(Rcpp::wrap(vr));
    
    double ey = mean(scr);
    double vary = var(scr);
    double diy = vary/ey;
    
    EY(t) = ey;
    VarY(t) = vary;
    DiY(t)= diy;
  }
  //
  //
  double EY_tot = mean(EY);
  double VarY_tot = mean(VarY);
  double DiY_tot = mean(DiY);
  
  NumericVector nw_quantities = {EY_tot, VarY_tot, DiY_tot };
  
  return(nw_quantities);
  
}


NumericVector  PosteriorPredictive_alpha( double alpha, double N, double Time){
  
  NumericVector EY(Time);
  NumericVector VarY(Time);
  NumericVector DiY(Time);
  arma::vec v_ones = arma::colvec(N,  fill::ones);
  
  for(int t = 0; t < Time; t++) {
    
    arma::mat Q = arma::mat(N,N,  fill::zeros);
    
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        
        
        if( j > i){
          
          
          double lambda_y  = exp(alpha);
          
          GetRNGstate();
          double rp =  R::rpois(lambda_y);
          PutRNGstate();
          
          Q(i,j)  = rp;
          
        }
        
      }
      
    }
    
    Q  = Q+Q.t();
    
    arma::colvec vr = Q*v_ones;
    // cout << vr;
    Rcpp::NumericVector scr = Rcpp::NumericVector(Rcpp::wrap(vr));
    
    
    double ey = mean(scr);
    double vary = var(scr);
    double diy = vary/ey;
    
    EY(t) = ey;
    VarY(t) = vary;
    DiY(t)= diy;
  }
  //
  //
  double EY_tot = mean(EY);
  double VarY_tot = mean(VarY);
  double DiY_tot = mean(DiY);
  
  NumericVector nw_quantities = {EY_tot, VarY_tot, DiY_tot };
  
  return(nw_quantities);
  
}

// [[Rcpp::export]]
double  alphaSamplingEn(NumericVector vec_y,  double  alpha,  double alpha_a,  double alpha_b,  double N, double Time, double M, double prop_sigma) {
  
  
  GetRNGstate();
  double alpha_prop = R::rnorm(alpha , prop_sigma );
  PutRNGstate();
  
  NumericVector pois(M*Time);
  NumericVector pois_prop(M*Time);
  
  for(int i = 0; i < M*Time; i++) {
    
    double lambda_aux = exp(alpha);
    double lambda_aux_prop = exp(alpha_prop);
    
    double pois_aux = R::dpois(vec_y(i), lambda_aux, true);
    double pois_prop_aux =  R::dpois(vec_y(i), lambda_aux_prop, true);
    
    
    //if(arma::is_finite(pois_aux)){}else{pois_aux = -1*pow(10, 26);};
    //if(arma::is_finite(pois_aux)){}else{pois_prop_aux = -1*pow(10, 26);};
    
    pois(i)  = pois_aux;  
    pois_prop(i)  = pois_prop_aux;
    
  }
  
  
  double prior = R::dnorm( alpha, alpha_a, alpha_b, true );
  double prior_prop = R::dnorm( alpha_prop, alpha_a, alpha_b, true );
  
  double test = prior_prop - prior + sum(pois_prop) - sum(pois) ;
  
  GetRNGstate();
  double rand = R::runif(0,1);
  PutRNGstate();

  
  if(test >  log(rand)){  alpha =  alpha_prop;  } ;
  
  return(alpha);
  
  
} 


// [[Rcpp::export]]
NumericVector  Likelihood_k(NumericVector princ_w, NumericVector leaning,  NumericVector beta, NumericMatrix zi_k, double gamma_0, double gamma_1, double phi, int N, int Time, int interp_eq){
  
  NumericVector likelihood_k(Time);
  NumericMatrix dist = dist_lcpp(zi_k, zi_k);
  
  double index_1 = 0;
  double index_2 = 0;
  
  for(int t = 0; t < Time; t++) {
    
    double sum_g = 0;
    double sum_b = 0;
    
    
    for(int i = 0; i < N; i++) {
      
        for(int j = i+1; j < N; j++) {

        double lambda = exp(beta(i) + beta(j) - pow(dist(i, j),2));
        double p_lk = R::dpois(princ_w(index_1), lambda, true);
        
        sum_g = sum_g + p_lk;
        index_1 = index_1 + 1;
      }
      
      double zi_k_1 = zi_k(i, 0);
      double DB_plane_a =  logistic_d(gamma_0+ interp_eq*gamma_1*zi_k_1)*phi;
      double DB_plane_b =  (1-logistic_d(gamma_0+ interp_eq*gamma_1*zi_k_1))*phi;
      
      double b_lk = R::dbeta(leaning(index_2), DB_plane_a, DB_plane_b, true);
      sum_b = sum_b + b_lk;
      
      index_2 = index_2 +1;
    }
    
    likelihood_k(t) = sum_g + sum_b;
  }
  
  return likelihood_k;
}
  


//MCMC function

// [[Rcpp::export]]
List MCMC(NumericVector princ_w,  NumericVector princ_ones, NumericVector x_w,  NumericVector x_ones,  NumericVector leaning,  NumericVector lam_ad_beta, NumericMatrix  mu_mat_beta,  List Sigma_ad_beta, NumericMatrix lam_ad_z, List  mu_mat_z,  List Sigma_ad_z,  NumericVector beta, NumericMatrix xi_state , List zi_in,   double mu_beta,  double sigma_beta , NumericVector mu_z ,  NumericVector sigma_z,  double s_a, double s_b, double phi, double gamma_0,  double gamma_1,  double a_phi, double b_phi, double a_gamma_0, double b_gamma_0 , double a_gamma_1, double b_gamma_1,  NumericVector omega, NumericMatrix P,  int N, int D, int K, int Time  ,
          NumericVector x_i, NumericVector DBplane_i, double prop_sd_rg = 0.0025,  double prop_sd_gamma = 0.00015, double prop_sd_phi = 40, double acc_beta = 0.25 , double acc_zeta = 0.25 , double pivot = 14, double sign = -1 , double interp_eq = 1, double ms_eq = 1, double rg_eq = 0,  int Iterations =10000){
  
  
  List zi = clone(zi_in);
  NumericMatrix nw_quantities_list(Iterations, 3);
  
  double M = (N/2.00)*(N-1);
  
    ////////////////////////
    //FULL SAMPLER 
    ////////////////////////
    
    NumericMatrix beta_ite(Iterations, N);
    NumericMatrix phi_gamma_ite(Iterations, 5);
    NumericMatrix sigma_ite(Iterations, K);
    
    List zi_ite(Iterations);
    
    NumericMatrix P_ite(Iterations, K*K);
    List xi_ite(Iterations); 
    
    NumericMatrix Eta(Time, K);
    List HETA(Iterations);
    NumericMatrix dist_m(2*M*Time, K);
    
    [[maybe_unused]] int index_z = 0;
    int index_x = 0;
    
    
    for(double ite = 0; ite < Iterations; ite++) {
      double prct = (ite+1)/Iterations;
      printProgress(prct);

      // start sampling
      
      NumericVector distance(2*M*Time);
      
      for(double k = 0; k < K; k++) {
        
        NumericMatrix zi_k = zi(k);
        NumericMatrix eucl_mat =  dist_fcpp(zi_k, zi_k );
        NumericMatrix eucl = gath(eucl_mat);
        NumericVector ecl = eucl(_,2);
        NumericVector distance_v = rep(ecl , Time);
        
        NumericVector xi_state_k = xi_state(_,k);
        NumericVector state_k_ij =  rep_each(xi_state_k, 2*M) ;
        
        distance= distance + distance_v*state_k_ij;
         
      }
      
      List betaObj = betaSamplingEnAdC( x_w,   lam_ad_beta,    mu_mat_beta,  Sigma_ad_beta,    beta,   distance,  mu_beta,   sigma_beta,  N,  Time,  x_i, acc_beta, ite) ;
    
      NumericVector beta_aux =   betaObj("beta");
      beta = beta_aux ;

       NumericVector aux_lam =  betaObj("lam_ad_beta");
       lam_ad_beta =  aux_lam;
       NumericMatrix aux_mat_beta = betaObj("mu_mat_beta");
       mu_mat_beta =   aux_mat_beta;
       List aux_list_beta = betaObj("Sigma_ad_beta");
       Sigma_ad_beta = aux_list_beta;
      
       beta_ite(ite, _ ) = beta;
       
      NumericVector zi_1_i(Time*N);
      
      for(double k = 0; k < K; k++) {
        
        NumericMatrix zi_k = zi(k);
        NumericVector zi_k_1 = zi_k(_, 0);
        NumericVector zi_k_1_v =  rep(zi_k_1, Time);
        NumericVector xi_state_k = xi_state(_,k);
        NumericVector state_k_i =  rep_each(xi_state_k, N) ;
        
        zi_1_i= zi_1_i + zi_k_1_v*state_k_i;
        
      }
      
      // GAMMAS and PHI
      
       phi = phiSamplingEn(leaning, zi_1_i, phi, gamma_0, gamma_1, a_phi, b_phi, interp_eq, prop_sd_phi );

       NumericVector gamma = gammaSamplingEn(leaning, zi_1_i, phi, gamma_0, gamma_1,  a_gamma_0, a_gamma_1,  b_gamma_0 ,  b_gamma_1, interp_eq, {prop_sd_gamma, prop_sd_gamma});
      //
       gamma_0 = gamma(0);
       gamma_1 = gamma(1);
      
       NumericVector phi_gamma_aux = {phi, gamma_0, gamma_1, 0, ite +1};
      
       phi_gamma_ite(ite,_ ) = phi_gamma_aux;
      
      NumericVector avg_dist(K);
      
      for(double k = 0; k < K; k++) {

        //draw Z and update 
        NumericVector xi_state_k = xi_state(_,k);
        NumericMatrix zi_s_aux = zi(k);
        NumericMatrix zi_s = clone(zi_s_aux);
      
      bool test = (sum(xi_state_k) == 0);
      cout << test;
      
      if(sum(xi_state_k) == 0){
      
      NumericMatrix zi_null(N, D);
      
      for (int i = 0; i < N; i++){
        for (int d = 0; d < D; d++){
          
          GetRNGstate();
          double zprop = R::rnorm(0, 1);
          PutRNGstate();
          zi_null(i, d) = zprop;
          
         }
      }
      
      zi_s = clone(zi_null);
        
      }else{
        
        
        NumericVector lam_ad_zs  = lam_ad_z(_, k);
        NumericMatrix mu_mat_zs = mu_mat_z(k);
        List Sigma_ad_zs = Sigma_ad_z(k);
        
        List Zobj_s =  ZetaSamplingAdC(zi_s,  lam_ad_zs,   mu_mat_zs, Sigma_ad_zs,  beta,   xi_state ,  x_i,  DBplane_i ,   x_w,   leaning,  x_ones ,  mu_z,  sigma_z,  gamma_0,  gamma_1,  phi,  M,  N,  D,  K,  k, acc_zeta,  interp_eq,  ite);
        
        NumericMatrix ziRes =   Zobj_s("zi_s");
        zi_s = clone(ziRes);
        
        NumericVector aux_vec_zs =  Zobj_s("lam_ad_zs");
        lam_ad_zs =  aux_vec_zs;
        NumericMatrix aux_mat_zs = Zobj_s("mu_mat_zs");
        mu_mat_zs =   aux_mat_zs;
        List aux_list_zs = Zobj_s("Sigma_ad_zs");
        Sigma_ad_zs = aux_list_zs;
        
        lam_ad_z(_, k) = clone(lam_ad_zs);
        mu_mat_z(k) = clone(mu_mat_zs);
        Sigma_ad_z(k) = clone(Sigma_ad_zs);
      }
        
        
        for (int d = 0; d < D; d++)
        {
          NumericVector zeta_v = zi_s(_, d);
          NumericVector zeta_demean = zeta_v - mean(zeta_v);
          zi_s(_, d) = clone(zeta_demean);
        }
        
        
       // if(ite > 1000){ 
        // List zi_ref =   zi_ite(ite -1);
        // NumericMatrix  zi_ref_k = zi_ref(k);
        // NumericMatrix zi_s_aux = procrustes_cpp(zi_ref_k, zi_s);
        // zi_s = clone(zi_s_aux);
        // }
  
         NumericVector zi_s_1 = zi_s(_, 0);
         double zi_s_1xx = zi_s_1(pivot-1);
        
         if(sign*zi_s_1xx < 0){
          zi_s_1 = -1*(zi_s_1);
          zi_s(_, 0) = clone(zi_s_1) ;
         }
        
        zi(k) = clone(zi_s);
        
        NumericVector zi_s_1x = zi_s(_, 0);
        NumericMatrix dist_s = ecldist(zi_s_1x);
        NumericMatrix dist_s_m = gath_tri(dist_s);
        NumericVector dist_s_v = dist_s_m(_,2);
        
        avg_dist(k) = median(dist_s_v);
      }
      
      IntegerVector indx_v =  seq_along(avg_dist) -1;
      std::sort(indx_v.begin(), indx_v.end(), [&](int i, int j){return avg_dist[i] > avg_dist[j];});
      NumericVector indx_vv = wrap(indx_v);
            
      List zi_aux(K);
      
      for (int k = 0; k < K; k++) {
        NumericMatrix aux_zi_k = zi(k);
        double idx =indx_vv[k];
        zi_aux(idx) = clone(aux_zi_k);
      }
      
      zi = clone(zi_aux); 
      
      // NumericMatrix zi_k1 = zi(0);
      // NumericMatrix zi_k2 = zi(1);
      
      // cout << " dist "<< avg_dist << endl;

      zi_ite(ite) = clone(zi);
      
      for(double k = 0; k < K; k++) {
         NumericMatrix zi_s_aux = zi(k);
         NumericMatrix zi_s = clone(zi_s_aux);
         double sigma_k =  invSamplingEn(zi_s ,   s_a ,  s_b,  N, D);
         sigma_z(k) = sigma_k;
       }
      
        sigma_ite(ite, _) = clone(sigma_z);
     
        // //// TRANSITION PROBABILITIES
        // 
        NumericMatrix P_aux = sampleTransP(xi_state,  omega ,   Time,  K);
        P = clone(P_aux);

        // cout << P<<endl;
        
        arma::mat P_arma = as<arma::mat>(P);
        arma::vec P_arma_vec = vectorise(P_arma);
        NumericVector P_vec = wrap(P_arma_vec);

        P_ite(ite, _) = P_vec;

      //   //// LIKELIHOOD COMPUTATION
      
      for(double k = 0; k < K; k++) {
      NumericMatrix zi_k_aux =  zi(k);
      NumericMatrix zi_k = clone(zi_k_aux);
      NumericVector lik_k = Likelihood_k( princ_w,  leaning,   beta,  zi_k,  gamma_0,  gamma_1,  phi,  N,  Time, interp_eq);
      Eta(_,k) = clone(lik_k);
      }
      
      HETA(ite) = clone(Eta);
      // 
       // cout << Eta;
      // 
      // //// FFBS ALGORITHM
      // 
      // double min_state = -1;
      // while(min_state <= 0){
      NumericMatrix xi_state_aux = FFBS(Eta,  P,  Time ,  K);
      xi_state = clone(xi_state_aux);
      
      // 
      // arma::mat xi_state_arma = as<arma::mat>(xi_state);
      // min_state = min(arma::sum(xi_state_arma, 0));
      // }

      xi_ite(ite) = clone(xi_state);
      
      // cout << xi_state;
      // // 
    }
    // 
    
    List RESULT = List::create( Named("beta_it")= beta_ite, Named("phi_gamma_ite")= phi_gamma_ite,   Named("zi_ite")= zi_ite,  Named("P_ite")= P_ite,  Named("xi_ite")= xi_ite, Named("sigma_ite")= sigma_ite, Named("logLik")= HETA);
    return RESULT;
  
  }
  
  
  //post processing
  
  
  // [[Rcpp::export]]
List update_zi(NumericMatrix zi_a_it, NumericMatrix zi_b_it) {
  int ncol = zi_a_it.ncol();
  int nrow = zi_a_it.nrow();
  
  for (int i = 0; i < ncol; i++) {
      for (int it = 0; it < nrow; it++) {
        double a = zi_a_it(it, i);
        double b = zi_b_it(it, i);
        
        if (std::abs(a) < std::abs(b)) {
          zi_a_it(it, i) = b;
          zi_b_it(it, i) = a;
        }
      }
  }
  
  List result = List::create(zi_a_it, zi_b_it);
  return result;
}


// [[Rcpp::export]]
NumericMatrix update_x(NumericMatrix x_it, NumericMatrix zi_a_it, NumericMatrix zi_b_it) {
  int ncol = x_it.ncol();
  int nrow = x_it.nrow();
  
  int current_iteration = 0;
  int total_iterations = ncol * nrow;

  
  for (int t = 0; t < ncol; t++) {
    for (int it = 0; it < nrow; it++) {
      NumericVector dist_a = dist(zi_a_it(it, _));
      NumericVector dist_b = dist(zi_b_it(it, _));
      
      if (median(dist_a) < median(dist_b)) {
        x_it(it, t) = 1 - x_it(it, t);
      }
      
      current_iteration ++;
    if (current_iteration % 100 == 0) {
        Rcpp::Rcout << "\rProgress: " << current_iteration << "/" << total_iterations << " [" << std::round(current_iteration / (double)total_iterations * 100.0) << "%]";
      }
      
      
    }
  }
  
  return x_it;
  

  
}