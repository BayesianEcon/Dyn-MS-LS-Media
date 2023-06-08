
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
NumericVector  Prediction_MS( NumericVector beta, NumericVector zi_a_1 , NumericVector zi_b_1,  NumericVector xi_state1 , double Npages, double Time){
  
  NumericVector EY(Time);
  NumericVector VarY(Time);
  NumericVector DiY(Time);
  arma::vec v_ones = arma::colvec(Npages,  fill::ones);
  
  NumericVector state = clone(xi_state1);
  
  for(int t = 0; t < Time; t++) {
    
    arma::mat Q = arma::mat(Npages,Npages,  fill::zeros);
    
    double chk = state(t);
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
       
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

NumericMatrix PredictiveIterationMS(NumericMatrix beta_m, NumericMatrix zi_a_1_m , NumericMatrix zi_b_1_m,  NumericMatrix xi_state1_m , double Npages, double Time, double start, double end , double by){
  
  
  NumericVector all_it = cseq(start, end, by);
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(Iterations ,3);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    NumericVector beta = beta_m(idx-1, _);
    NumericVector zi_a_1 = zi_a_1_m(idx-1, _);
    NumericVector zi_b_1 = zi_b_1_m(idx-1, _);
    NumericVector xi_state1 = xi_state1_m(idx-1, _);
    
    NumericVector res = Prediction_MS(  beta,  zi_a_1 ,  zi_b_1,   xi_state1 ,  Npages,  Time);
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
NumericVector  PredictionScore_MS(NumericVector y_princ, NumericVector beta, NumericVector zi_a_1 , NumericVector zi_b_1,  NumericVector xi_state1 , double Npages, double Time, bool lg){
  
  
  double len = y_princ.length();
  NumericVector res(len);
  
  NumericVector state = clone(xi_state1);
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    
    double chk = state(t);
    
    for(int i = 0; i < Npages; i++) {
      for(int j = i+1; j < Npages; j++) {
        
        double yit = y_princ[idx];
       
          if(chk == 1){
            
            double lambda_y  = exp(beta(i)  + beta(j) - pow(zi_a_1(i) - zi_a_1(j),2));
            double lkpois = R::dpois(yit, lambda_y, lg);
            res[idx] = lkpois;
            
        }else{
            
            double lambda_y  = exp(beta(i) + beta(j) - pow(zi_b_1(i) - zi_b_1(j),2));
            double lkpois = R::dpois(yit, lambda_y, lg);
            res[idx] = lkpois;

          }
          
          
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
NumericVector  PredictionScore_Plane(NumericVector plane_ell, double gamma_0, double gamma_1,  double phi, NumericVector zi_a_1 , NumericVector zi_b_1,  NumericVector xi_state1 , double Npages, double Time, bool lg){
  
  
  double len = plane_ell.length();
  NumericVector res(len);
  
  NumericVector state = clone(xi_state1);
  
  double idx = 0;
  
  for(int t = 0; t < Time; t++) {
    
    
    double chk = state(t);
    
    for(int i = 0; i < Npages; i++) {
     
        
        double plane_ell_it = plane_ell[idx];
        
        
          if(chk == 1){
            
            double DB_plane_a_a  = logistic_d(gamma_0+ (gamma_1)*zi_a_1(i))*phi;
            double DB_plane_a_b  =  (1-logistic_d(gamma_0+(gamma_1)*zi_a_1(i)))*phi;

            
            double lkbeta =  R::dbeta(plane_ell_it, DB_plane_a_a,DB_plane_a_b , lg);
            res[idx] = lkbeta;
            
        }else{
            
            double DB_plane_b_a  = logistic_d(gamma_0+ (gamma_1)*zi_b_1(i))*phi;
            double DB_plane_b_b  =  (1-logistic_d(gamma_0+ (gamma_1)*zi_b_1(i)))*phi;

            
            double lkbeta =  R::dbeta(plane_ell_it, DB_plane_b_a,DB_plane_b_b , lg);
            res[idx] = lkbeta;

          }
          
          
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

NumericMatrix logPredictiveScoreMS(NumericVector  y_princ,  NumericMatrix beta_m, NumericMatrix zi_a_1_m , NumericMatrix zi_b_1_m,  NumericMatrix xi_state1_m , double Npages, double Time, double start, double end , double by, bool lg = true){
  
  
  NumericVector all_it = cseq(start-1, end-1, by);
  double len = y_princ.length();
  double Iterations = all_it.length();
  NumericMatrix RES = no_init_matrix(len , Iterations);
  
  for(int it = 0; it <  Iterations; it++ ){
    
    double prct = (it+1)/Iterations;
    printProgress(prct);
    
    double idx = all_it[it];
    NumericVector beta = beta_m(idx, _);
    NumericVector zi_a_1 = zi_a_1_m(idx, _);
    NumericVector zi_b_1 = zi_b_1_m(idx, _);
    NumericVector xi_state1 = xi_state1_m(idx, _);
    
    NumericVector res = PredictionScore_MS( y_princ, beta,  zi_a_1 ,  zi_b_1,   xi_state1 ,  Npages,  Time, lg);
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

NumericMatrix logPredictiveScorePlane(NumericVector plane_ell, NumericVector gamma_0_ite, NumericVector gamma_1_ite,  NumericVector phi_ite, NumericMatrix zi_a_1_m , NumericMatrix zi_b_1_m,  NumericMatrix xi_state1_m , double Npages, double Time, double start, double end , double by, bool lg = true){
  
  
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


    NumericVector zi_a_1 = zi_a_1_m(idx, _);
    NumericVector zi_b_1 = zi_b_1_m(idx, _);
    NumericVector xi_state1 = xi_state1_m(idx, _);;
    
    NumericVector res = PredictionScore_Plane( plane_ell,  gamma_0,  gamma_1,   phi,  zi_a_1 ,  zi_b_1,   xi_state1 ,  Npages,  Time, lg);
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






