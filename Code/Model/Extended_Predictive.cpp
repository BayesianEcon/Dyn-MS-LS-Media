
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

Environment pkg = Environment::namespace_env("stats");
Function rmulti = pkg["rmultinom"];

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

// [[Rcpp::export]]
double logistic_d(double x) {
  double result = 1/(1+exp(-x));
  return result;
}


// [[Rcpp::export]]

NumericVector Prediction_Leaning(NumericVector leaning , NumericVector time, NumericVector beta, double delta, double N, double T){
  
  NumericVector EY(T);
  NumericVector VarY(T);
  NumericVector DiY(T);
  arma::vec v_ones = arma::colvec(N,  arma::fill::ones);
  
  double count = 0;
  double M =  (N-1)*(N/2);
  
  for(int t = 0; t <  T; t++ ){
    arma::mat Y(N,N, fill::zeros);
    
    double idx = 0;
    
    NumericVector selector = cseq(count, count + M -1 , 1);
    count = count + M;
    
    NumericVector leaning_t = leaning[selector];
    //NumericVector leaning_t = leaning[time == t+1];
    
    for(int i = 0; i <  N; i++ ){
      for(int j = i+1; j <  N; j++ ){
        
          double lambda = exp(beta(i) + beta(j) + delta*leaning_t[idx]);
          GetRNGstate() ;
          double y = R::rpois(lambda);
          PutRNGstate() ;
          
          Y(i,j)  = y;
          
          idx = idx +1;
        
      }}
    
    Y  = Y+Y.t();

    arma::colvec vr = Y*v_ones;
    Rcpp::NumericVector scr = Rcpp::NumericVector(Rcpp::wrap(vr));
    
    double ey = mean(scr);
    double vary = var(scr);
    double diy = vary/ey;
    
    //return {ey, vary, diy};
    
    EY(t) = ey;
    VarY(t) = sqrt(vary);
    DiY(t)= diy;
    
  }
  
  double m_ey =mean(EY);
  double m_VarY =mean(VarY);
  double m_DiY =mean(DiY);
  
  NumericVector res = {m_ey, m_VarY, m_DiY};
  return res;
  
}


// [[Rcpp::export]]
NumericVector  Prediction_MS( NumericVector beta, List zi_in, NumericMatrix xi_state , double Npages, double Time, double K, double D){
  
  NumericVector EY(Time);
  NumericVector VarY(Time);
  NumericVector DiY(Time);
  arma::vec v_ones = arma::colvec(Npages,  fill::ones);
  NumericVector n_states = cseq(0, K-1, 1);
  NumericMatrix state = clone(xi_state);
  
  for(int t = 0; t < Time; t++) {
    
    arma::mat Q = arma::mat(Npages,Npages,  fill::zeros);
    
    NumericVector chk_time = state(t,_);

    NumericVector chk_v = n_states[chk_time == 1];
    double chk = chk_v(0);

    NumericMatrix zi_k = zi_in(chk);
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
            
            double dist2 = 0;
            for(int d = 0; d < D; d++) { 
              dist2 = dist2 +  pow(zi_k(i, d) - zi_k(j, d),2);
            }

            double lambda_y  = exp(beta(i)  + beta(j) - dist2);
            
            GetRNGstate();
            double rp =  R::rpois(lambda_y);
            PutRNGstate();
            
            Q(i,j)  = rp;
            
        
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
    VarY(t) = sqrt(vary);
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

NumericMatrix PredictiveIterationLeaning(NumericVector leaning , NumericVector time, NumericMatrix beta_m, NumericVector delta_v, double N, double T, double start, double end , double by){
  
  
  NumericVector all_it = cseq(start, end, by);
  double Iterations = all_it.length();
  NumericMatrix RES(Iterations ,3);
  
  for(int it = 0; it <  Iterations; it++ ){
    
  double prct = (it+1)/Iterations;
  printProgress(prct);
    
  double idx = all_it[it];
  NumericVector beta = beta_m(idx-1, _);
  double delta = delta_v(idx-1);
    
  NumericVector res = Prediction_Leaning( leaning ,  time,  beta,  delta,  N,  T);
  RES(it, _) = res;

  }
  
  
  return(RES);

}




// [[Rcpp::export]]

NumericMatrix PredictiveIterationMS(NumericMatrix beta_m, List zi_l, List xi_l , double Npages, double Time, double K, double D, double start, double end , double by){
  
  NumericVector all_it = cseq(start, end, by);
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(Iterations ,3);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    NumericVector beta = beta_m(idx-1, _);
    List zi_in = zi_l[idx-1];
    NumericMatrix xi_state_in = xi_l[idx-1];

    NumericVector res = Prediction_MS(  beta,  zi_in ,  xi_state_in ,  Npages,  Time, K, D);

    RES(it, _) = res;
    
  }
  
  
  return(RES);
  
}



// [[Rcpp::export]]

NumericVector PredictiveIterationScatterMS(NumericVector beta, NumericVector zi_a_1 , NumericVector zi_b_1,  NumericVector xi_state1_m , double Npages, double Time){
  
  
  double M = (Npages/2)*(Npages-1);
  NumericVector RES(M*Time);
  
  double idx = 0;
  
  
  for(int t = 0; t < Time; t++) {
    
    double chk = xi_state1_m(t);
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        
        if(chk == 1){
          
          double lambda_y  =  exp(beta(i)  + beta(j) - pow(zi_a_1(i) - zi_a_1(j),2));
          
          // GetRNGstate();
          // double rp =  R::rpois(lambda_y);
          // PutRNGstate();
          
          RES(idx) = lambda_y;          
        }else{
          
          double lambda_y  =  exp( beta(i) + beta(j) - pow(zi_b_1(i) - zi_b_1(j),2));
          // GetRNGstate();
          // double rp =  R::rpois(lambda_y);
          // PutRNGstate();
          RES(idx) = lambda_y;
        }
        
        idx = idx +1;
        
      }
    }
  }
   
  

  
  return(RES);
  
  
}


// [[Rcpp::export]]

NumericVector PredictiveIterationScatterLeaning(NumericVector beta, double delta, NumericVector leaning, double Npages, double Time){
  
  
  double M = (Npages/2)*(Npages-1);
  NumericVector RES(M*Time);
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        double lambda_y  =  exp( beta(i)  + beta(j)  + delta*leaning[idx]);
        
        // GetRNGstate();
        // double rp =  R::rpois(lambda_y);
        // PutRNGstate();
        
        RES(idx) = lambda_y;
        
        
        
        idx = idx +1;
        
      }
    }
  }
  
  
  
  
  return(RES);
  
}




// [[Rcpp::export]]

List Prediction_Alpha(double alpha, double N, double T){
  
  NumericVector EY(T);
  NumericVector VarY(T);
  NumericVector DiY(T);
  arma::vec v_ones = arma::colvec(N,  arma::fill::ones);
  NumericVector Yv(T*(N/2)*(N-1));
  
  double count = 0;
  double M =  (N-1)*(N/2);
  
  double idx2 = 0;
  
  for(int t = 0; t <  T; t++ ){
    arma::mat Y(N,N, fill::zeros);
    
    double idx = 0;
    
    NumericVector selector = cseq(count, count + M -1 , 1);
    count = count + M;
    

    for(int i = 0; i <  N; i++ ){
      for(int j = i+1; j <  N; j++ ){
        
        double lambda = exp(alpha);
        GetRNGstate() ;
        double y = R::rpois(lambda);
        PutRNGstate() ;
        
        Y(i,j)  = y;
        Yv(idx2) = y;
        
        
        idx = idx +1;
        idx2 = idx2 +1;
        
        
      }}
    
    Y  = Y+Y.t();
    
    arma::colvec vr = Y*v_ones;
    Rcpp::NumericVector scr = Rcpp::NumericVector(Rcpp::wrap(vr));
    
    double ey = mean(scr);
    double vary = var(scr);
    double diy = vary/ey;
    
    //return {ey, vary, diy};
    
    EY(t) = ey;
    VarY(t) = sqrt(vary);
    DiY(t)= diy;
    
  }
  
  double m_ey =mean(EY);
  double m_VarY =mean(VarY);
  double m_DiY =mean(DiY);
  
  NumericVector res = {m_ey, m_VarY, m_DiY};
  List Res  = List::create( Named("res")= res, Named("Yv")= Yv);
  return Res;
  
}





// [[Rcpp::export]]

List PredictiveIterationAlpha( NumericVector alpha_v, double N, double T, double start, double end , double by){
  
  
  NumericVector all_it = cseq(start, end, by);
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(Iterations ,3);
  NumericMatrix Ymat = no_init_matrix( T*(N/2)*(N-1), Iterations);
  
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    double alpha = alpha_v(idx-1);
    
    List resobj = Prediction_Alpha(alpha,  N,  T);
    NumericVector res = resobj("res");
    NumericVector Yvs = resobj("Yv");
    
    
    RES(it, _) = res;
    Ymat(_,it) = Yvs;
    
  }
  
  
  List RESULT = List::create( Named("RES")= RES, Named("Ymat")= Ymat);
  return(RESULT);
  
}


//////////////////////////

// [[Rcpp::export]]
NumericVector  PredictionScore_MS(NumericVector y_princ, NumericVector beta, List zi_in ,  NumericMatrix xi_state , double Npages, double Time, double K, double D, bool lg){
  
  
  double len = y_princ.length();
  NumericVector res(len);
  NumericVector n_states = cseq(0, K-1, 1);
  
  NumericMatrix state = clone(xi_state);
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    NumericVector chk_time = state(t,_);
    NumericVector chk_v = n_states[chk_time == 1];
    double chk = chk_v(0);
    
    NumericMatrix zi_k = zi_in(chk);
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        double dist2 = 0;
        for(int d = 0; d < D; d++) { 
          dist2 = dist2 +  pow(zi_k(i, d) - zi_k(j, d),2);
        }
        
        double yit = y_princ[idx];
       
            double lambda_y  = exp(beta(i)  + beta(j) - dist2);
            double lkpois = R::dpois(yit, lambda_y, lg);
            res[idx] = lkpois;
            
        idx = idx +1;
        
        
      }
    }
  }
  //
  //

  
  return(res);
  
}


// [[Rcpp::export]]
NumericVector  PredictionScore_AlphaEll(NumericVector y_princ, NumericVector ell, NumericVector beta, double delta , double Npages, double Time, bool lg){
  
  
  double len = y_princ.length();
  NumericVector res(len);
  

  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    

    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        double yit = y_princ[idx];
        
          
          double lambda_y  = exp(beta(i)  + beta(j) - delta*ell[idx]);
          double lkpois = R::dpois(yit, lambda_y, lg);
          res[idx] = lkpois;
          
       
        
        
        idx = idx +1;
        
        
      }
    }
  }
  //
  //
  
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector  PredictionScore_Alpha(NumericVector y_princ, double alpha, double Npages, double Time, bool lg){
  
  
  double len = y_princ.length();
  NumericVector res(len);
  
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
        
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        double yit = y_princ[idx];
                   
            double lambda_y  = exp(alpha);
            double lkpois = R::dpois(yit, lambda_y, lg);
            res[idx] = lkpois;
            
                
          idx = idx +1;
        
        
      }
    }
  }
  //
  //

  
  return(res);
  
}



// [[Rcpp::export]]
NumericVector  PredictionScore_Plane(NumericVector plane_ell, double gamma_0, double gamma_1,  double phi, List zi_in ,  NumericMatrix xi_state , double Npages, double Time, double K, bool lg){
  
  
  double len = plane_ell.length();
  NumericVector res(len);
  NumericVector n_states = cseq(0, K-1, 1);
  
  NumericMatrix state = clone(xi_state);
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    NumericVector chk_time = state(_,t);
    NumericVector chk_v = n_states[state == 1];
    double chk = chk_v(0);
    
    NumericMatrix zi_k = zi_in(chk);
    
    for(int i = 0; i < Npages; i++) {
     
        double plane_ell_it = plane_ell[idx];
            
            double DB_plane_a_a  = logistic_d(gamma_0+ (gamma_1)*zi_k(i, 0))*phi;
            double DB_plane_a_b  =  (1-logistic_d(gamma_0+(gamma_1)*zi_k(i, 0)))*phi;

            double lkbeta =  R::dbeta(plane_ell_it, DB_plane_a_a,DB_plane_a_b , lg);
            res[idx] = lkbeta;
            
           idx++;
              
    }
  }
  //
  //

  
  return(res);
  
}








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

NumericMatrix logPredictiveScoreMS(NumericVector  y_princ,  NumericMatrix beta_m, List zi_l ,  List xi_l , double Npages, double Time, double K, double D, double start, double end , double by, bool lg = true){
  
  
  NumericVector all_it = cseq(start-1, end-1, by);
  double len = y_princ.length();
  double Iterations = all_it.length();
  
  NumericMatrix RES(len , Iterations);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it(it);
    NumericVector beta = beta_m(idx-1, _);
    List zi_in = zi_l(idx-1);
    NumericMatrix xi_state_in = xi_l(idx-1);
    
    NumericVector res = PredictionScore_MS( y_princ, beta,  zi_in,   xi_state_in ,  Npages,  Time, K, D, lg);
    RES(_, it) = res;
    
  }
  
  
  return(RES);
  
}

// [[Rcpp::export]]

NumericMatrix logPredictiveScoreAlpha(NumericVector  y_princ,  NumericVector alpha_m,  double Npages, double Time, double start, double end , double by, bool lg = true){
  
  
  NumericVector all_it = cseq(start-1, end-1, by);
  double len = y_princ.length();
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(len , Iterations);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    double alpha = alpha_m(idx);

    
    NumericVector res = PredictionScore_Alpha( y_princ, alpha,  Npages,  Time, lg);
    RES(_, it) = res;
    
  }
  
  
  return(RES);
  
}


// [[Rcpp::export]]

NumericMatrix logPredictiveScoreAlphaEll(NumericVector  y_princ, NumericVector ell, NumericMatrix beta_m, NumericVector delta_m, double Npages, double Time, double start, double end , double by, bool lg = true){
  
  
  NumericVector all_it = cseq(start-1, end-1, by);
  double len = y_princ.length();
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(len , Iterations);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    NumericVector beta = beta_m(idx, _);
    double delta = delta_m(idx);
    
    NumericVector res = PredictionScore_AlphaEll( y_princ, ell,  beta,  delta ,  Npages,  Time, lg);
    RES(_, it) = res;
    
  }
  
  
  return(RES);
  
}

// [[Rcpp::export]]

NumericMatrix logPredictiveScorePlane(NumericVector plane_ell, NumericVector gamma_0_ite, NumericVector gamma_1_ite,  NumericVector phi_ite, List zi_l ,  List xi_l , double Npages, double Time, double K, double start, double end , double by, bool lg = true){
  
  
  NumericVector all_it = cseq(start-1, end-1, by);
  double len = plane_ell.length();
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(len , Iterations);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];

    double gamma_0 = gamma_0_ite(idx);
    double gamma_1 = gamma_1_ite(idx);
    double phi = phi_ite(idx);

    List zi_in = zi_l(idx-1);
    NumericMatrix xi_state_in = xi_l(idx-1);
    
    NumericVector res = PredictionScore_Plane( plane_ell,  gamma_0,  gamma_1,   phi, zi_in,   xi_state_in ,  Npages,  Time, K, lg);
    RES(_, it) = res;
    
  }
  
  return(RES);
  
}

// [[Rcpp::export]]

NumericVector vec_nbet_tri(NumericVector beta, double Npages){
  
  NumericVector res;
  
  for( int i = 0; i < (Npages-1); i++){
    
    NumericVector  bb = rep(beta,1);
    
    for( int j = i + 1; j < (Npages); j++){
      
      double tt = bb(j);
      res.push_back(tt);
      
    }
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
NumericVector vec_nbet_tri_a(NumericVector beta, double Npages){
  
  NumericVector res;
  
  for( int i = 0; i < (Npages-1); i++){
    
    NumericVector  bb = rep(beta(i), Npages-i-1);
    double len = bb.length();
    
    for( int j = 0 ; j < len; j++){
      
      double tt = bb(j);
      res.push_back(tt);
    }
  }
  
  return(res);
  
}

// [[Rcpp::export]]
NumericMatrix  PredictionScore_MS_metrics(NumericVector y_princ, NumericVector beta,  List zi_in ,  NumericMatrix xi_state , double Npages, double Time, double K, double D){
  
  arma::vec v_ones = arma::colvec(Npages,  arma::fill::ones);
  double len = y_princ.length();
  NumericMatrix res(Time, 3);
  NumericVector n_states = cseq(0, K-1, 1);
  NumericMatrix state = clone(xi_state);
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    arma::mat Y(Npages, Npages);
    
    NumericVector chk_time = state(t,_);
    NumericVector chk_v = n_states[chk_time == 1];
    double chk = chk_v(0);
    
    NumericMatrix zi_k = zi_in(chk);
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        double dist2 = 0;
        for(int d = 0; d < D; d++) { 
          dist2 = dist2 +  pow(zi_k(i, d) - zi_k(j, d),2);
        }
        
        double yit = y_princ[idx];
        
        double lambda_y  = exp(beta(i)  + beta(j) - dist2);
        double lkpois = R::rpois(lambda_y);
        
        Y(i, j) = lkpois;
        Y(j, i) = lkpois;
        
        idx = idx +1;

      }
    }
    
    arma::colvec vr = Y*v_ones;
    Rcpp::NumericVector scr = Rcpp::NumericVector(Rcpp::wrap(vr));
    
    double ey = mean(scr);
    double vary = var(scr);
    double diy = vary/ey;
    
    res(t,0) = ey;
    res(t,1) = vary;
    res(t,2) = diy;
  }
  
  return(res);
  
}


// [[Rcpp::export]]
List logPredictiveScoreMS_metrics(NumericVector  y_princ,  NumericMatrix beta_m, List zi_l ,  List xi_l , double Npages, double Time, double K, double D, double start, double end , double by){
  
  NumericVector all_it = cseq(start-1, end-1, by);
  double len = y_princ.length();
  double Iterations = all_it.length();
  NumericMatrix RES1 = no_init_matrix(Iterations , Time);
  NumericMatrix RES2 = no_init_matrix(Iterations , Time);
  NumericMatrix RES3 = no_init_matrix(Iterations , Time);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    NumericVector beta = beta_m(idx, _);
    List zi_in = zi_l(idx-1);
    NumericMatrix xi_state_in = xi_l(idx-1);
    
    NumericMatrix res = PredictionScore_MS_metrics(y_princ, beta,  zi_in,   xi_state_in ,  Npages,  Time, K, D);
    
    NumericVector res1_v = res(_,0);
    RES1(it, _) = res1_v;
    
    NumericVector res2_v = res(_,1);
    RES2(it, _) = res2_v;
    
    NumericVector res3_v = res(_,2);
    RES3(it, _) = res3_v;
  }
  
  List RES  = List::create( Named("EY")= RES1, Named("Var")= RES2,Named("DI")= RES3);
  return(RES);
  
}
// 
// 
// // [[Rcpp::export]]
// List OutofSample(NumericMatrix beta_m, List zi_l ,  List xi_l , NumericMatrix Plong, double Npages, double Time, double K, double D, double start, double end , double by){
//   
//   NumericVector all_it = cseq(start-1, end-1, by);
//   double Iterations = all_it.length();
//   List Ypred_list(Iterations);
//   
//   NumericVector n_states = cseq(0, K-1, 1);
//   
// 
//   for(int it = 0; it <  Iterations; it++ ){
//     
//     double prct = (it+1)/Iterations;
//     printProgress(prct);
//     
//     double idx = all_it[it];
//     NumericVector beta = beta_m(idx, _);
//     List zi_in = zi_l(idx-1);
//     NumericMatrix xi_state_in = xi_l(idx-1);
//     NumericVector P_vec = Plong(idx-1, _);
//     arma::vec P_vec_arma = as<arma::vec>(P_vec);
//     arma::mat P_arma(P_vec_arma.memptr(), Npages, Npages, false, false);
//       
//     NumericVector xi_state_T =  xi_state_in(Time, _);
//     arma::vec xi_state_T_arma = as<arma::vec>(xi_state_T);
//     arma::vec Prob = P_arma*xi_state_T_arma;
//     NumericVector Prob_num = wrap(Prob);
//     
//     GetRNGstate();
//     NumericVector next_state = rmulti(1,1, Prob_num);
//     PutRNGstate();
//     
//     arma::mat Ypred(Npages, Npages);
//     NumericVector chk_v = n_states[next_state == 1];
//     double chk = chk_v(0);
//     
//     NumericMatrix zi_k = zi_l(chk);
//     
//     for(int i = 0; i < Npages; i++) {
//       for(int j = i+1; j < Npages; j++) {
//         
//         double dist2 = 0;
//         for(int d = 0; d < D; d++) {
//           dist2 = dist2 +  pow(zi_k(i, d) - zi_k(j, d),2);
//         }
//         
//         double lambda_y  = exp(beta(i)  + beta(j) - dist2);
//         double yres = R::rpois(lambda_y);
//         
//         Ypred(i, j) = yres;
//         Ypred(j, i) = yres;
//         
//       }
//     }
//     
//     Ypred_list(it) = clone(Ypred);
// 
//   
//   }
//   
//   List RES  = List::create( Named("Ypred_list")= Ypred_list);
//   return(RES);
//   
// }
