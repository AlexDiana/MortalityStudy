
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


arma::vec mvrnormArma(arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return mu + arma::trans(Y * arma::chol(Sigma));
}

// [[Rcpp::export]]
arma::vec mrt2(arma::vec mean, arma::mat Sigma, 
               double df){
  arma::vec zeroVec = arma::zeros(mean.size());
  arma::vec y = mvrnormArma(zeroVec, Sigma);
  double u = R::rchisq(df);
  arma::vec x = sqrt(df / u) * y + mean;
  return x;
}

// [[Rcpp::export]]
double dmt_cpp(arma::vec x, double nu, arma::vec mu, arma::mat Sigma, 
               bool returnLog){
  int p = x.size();
  
  double logratio = R::lgammafn(( nu + p ) / 2) - R::lgammafn( nu / 2 );
  
  arma::vec product = (arma::trans(x - mu) * arma::inv(Sigma) * (x - mu));
  double lognum = (- ( nu + p ) / 2) * log(1 + (1 / nu) * product[0]);
  double logdetSigma = arma::log_det(Sigma).real();
  
  double logden = (p / 2.0) * log(M_PI * nu) + (0.5) * logdetSigma;
  // double logden = (p / 2.0) * log(M_PI * nu) + (0.5) * log(arma::det(Sigma));
  
  double loglikelihood = logratio + lognum - logden;
  
  if(returnLog){
    return loglikelihood;
  } else {
    return exp(loglikelihood);
  }
}

// [[Rcpp::export]]
arma::mat createTerm4_cpp(arma::vec gtx, int X, int Y){
  
  arma::mat term4 = arma::zeros(X, Y);
  
  for(int x = 0; x < X; x++){
    for(int y = 0; y < Y; y++){
      term4(x,y) = gtx[(y - x) - ( - (X-1))];
    }  
  }
  
  return(term4);
  
}

// [[Rcpp::export]]
arma::vec gr_loglik_term4_gtx_m1_fun_cpp(arma::vec param, 
                                         arma::mat d, 
                                         arma::mat E, 
                                         arma::mat cxt,
                                         arma::vec ik_m1){
    
    int X = d.n_rows;
    int Y = d.n_cols;
    
    arma::vec gtx = param;
    gtx.resize(param.n_elem + 1); 
    gtx(param.n_elem) = 0;
  
    arma::mat gtx_mat = createTerm4_cpp(gtx, X, Y);
    
    double sumParams = accu(gtx_mat);
    gtx_mat(0,Y-1) = -sumParams;
    
    arma::vec term1 = arma::zeros(X+Y-2);
    arma::vec term2 = arma::zeros(X+Y-2);
    for(int g = 0; g < (X+Y-2); g++){
      term1[g] += sum(diagvec(d, g - X + 1));
      term2[g] -= sum(diagvec(exp(cxt + gtx_mat + log(E)), g - X + 1));
    }

    arma::vec term3 = d(0, Y-1) * (- ik_m1);
    
    arma::vec term4 = - exp(cxt(0, Y-1) + gtx_mat(0, Y-1) + 
      log(E(0, Y-1))) * (-ik_m1);

    return(- (term1 + term2 + term3 + term4));
  
}


// [[Rcpp::export]]
arma::mat hess_loglik_term4_gtx_m1_fun_cpp(arma::vec param, 
                                         arma::mat d, 
                                         arma::mat E, 
                                         arma::mat cxt,
                                         arma::vec ik_m1){
    
    int X = d.n_rows;
    int Y = d.n_cols;
    
    arma::vec gtx = param;
    gtx.resize(param.n_elem + 1); 
    gtx(param.n_elem) = 0;
  
    arma::mat gtx_mat = createTerm4_cpp(gtx, X, Y);
    
    double sumParams = accu(gtx_mat);
    gtx_mat(0,Y-1) = -sumParams;
    
    arma::mat hess = arma::zeros(X+Y-2, X+Y-2);
    
    for (int g1 = 0; g1 < (X+Y-2); g1++) {
      
    // term2
      
      hess(g1,g1) -= 
        sum(diagvec(exp(cxt + gtx_mat + log(E)), g1 - X  + 1));
        
        for (int g2 = 0; g2 < (X+Y-2); g2++) {
          
    // term3
          
          hess(g1,g2) -=
            ik_m1[g1] * ik_m1[g2] * exp(cxt(0, Y-1) + gtx_mat(0, Y-1) + 
            log(E(0, Y-1)));
            
        }  
    }
    
    return(- hess);
  
}



// loglikelihood of 

// [[Rcpp::export]]
double loglik_LC(
    arma::vec a,
    arma::vec b,
    arma::vec k,
    arma::cube d,
    arma::cube E) {
  
  double loglik = 0;
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  for(int x = 0; x < X; x++){
    
    for(int t = 0; t < Y; t++){
      
      for(int p = 0; p < P; p++){
        
        double term1 = 
          d(x,t,p) * (a[x] + b[x] * k[t] + log(E(x,t,p)));
        double term2 = exp(a[x] + b[x] * k[t] + log(E(x,t,p)));
        loglik += term1 - term2;
        
      }
      
    }
    
  }
  
  return loglik;
}

// [[Rcpp::export]]
double loglik_LCp(
    arma::mat a,
    arma::mat b,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  double loglik = 0;
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  for(int x = 0; x < X; x++){
    
    for(int t = 0; t < Y; t++){
      
      for(int p = 0; p < P; p++){
        
        if(arma::is_finite(d(x,t,p))){
          
          double mxtp = a(x,p) + b(x,p) * k(t,p);
          
          double term1 = d(x,t,p) * (mxtp + log(E(x,t,p)));
          double term2 = exp(mxtp + log(E(x,t,p)));
          
          loglik += term1 - term2;
        }
        
        
      }
      
    }
    
  }
  
  return loglik;
}

// [[Rcpp::export]]
double loglik_LCp_x(
    int x,
    arma::mat a,
    arma::mat b,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  double loglik = 0;
  
  int Y = E.n_cols;
  int P = E.n_slices;
  
  for(int t = 0; t < Y; t++){
    
    for(int p = 0; p < P; p++){
      
      if(arma::is_finite(d(x,t,p))){
        double mxtp = a(x,p) + b(x,p) * k(t,p);
        
        double term1 = d(x,t,p) * (mxtp + log(E(x,t,p)));
        double term2 = exp(mxtp + log(E(x,t,p)));
        
        loglik += term1 - term2;
      }
      
    }
    
  }
  
  return loglik;
}

// [[Rcpp::export]]
double loglik_LCp_t(
    int t,
    arma::mat a,
    arma::mat b,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  double loglik = 0;
  
  int X = E.n_rows;
  int P = E.n_slices;
  
  for(int x = 0; x < X; x++){
    
    for(int p = 0; p < P; p++){
      
      if(arma::is_finite(d(x,t,p))){
        
        double mxtp = a(x,p) + b(x,p) * k(t,p);
        
        double term1 = d(x,t,p) * (mxtp + log(E(x,t,p)));
        double term2 = exp(mxtp + log(E(x,t,p)));
        
        loglik += term1 - term2;
        
      }
    }
    
  }
  
  return loglik;
}

// gradient of loglikelihood of ap for a given x given a

// [[Rcpp::export]]
arma::vec gr_loglik_ap_x(
    int x,
    arma::vec a,
    arma::mat ap,
    arma::mat bp,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(P);
  
  for(int t = 0; t < Y; t++){
    
    for(int p = 0; p < P; p++){
      
      double m_xtp = (a[x] + ap(x,p)) + bp(x,p) * k(t,p) + log(E(x, t, p));
      
      double term2 = exp(m_xtp);
      
      gr_loglik[p] += d(x, t, p) - term2;
      
    }
    
  }
  
  return gr_loglik;
}

// gradient of loglikelihood of bp for a given x given a

// [[Rcpp::export]]
arma::vec gr_loglik_bp_x(
    int x,
    arma::mat ap,
    arma::vec b,
    arma::mat bp,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(P);
  
  for(int t = 0; t < Y; t++){
    
    for(int p = 0; p < P; p++){
      
      double m_xtp = ap(x,p) + (b[x] + bp(x,p)) * k(t,p) + log(E(x, t, p));
      
      double term2 = exp(m_xtp);
      
      gr_loglik[p] += d(x, t, p) * k(t,p) - term2 * k(t,p);
      
    }
    
  }
  
  return gr_loglik;
}

// gradient of loglikelihood of kp for a given t given k

// [[Rcpp::export]]
arma::vec gr_loglik_kp_t(
    int t,
    arma::mat ap,
    arma::mat bp,
    arma::vec k,
    arma::mat kp,
    arma::cube d,
    arma::cube E) {
  
  int X = E.n_rows;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(P);
  
  for(int x = 0; x < X; x++){
    
    for(int p = 0; p < P; p++){
      
      double m_xtp = ap(x,p) + bp(x,p) * (k[t] + kp(t,p)) + log(E(x, t, p));
      
      double term2 = exp(m_xtp);
      
      gr_loglik[p] += d(x, t, p) * bp(x,p) - term2 * bp(x,p);
      
    }
    
  }
  
  return gr_loglik;
}

// gradient of loglikelihood of ab for a given x given abp

// [[Rcpp::export]]
arma::vec gr_loglik_ab_x(
    int x,
    arma::vec a,
    arma::vec b,
    arma::mat ap,
    arma::mat bp,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(2);
  
  for(int t = 0; t < Y; t++){
    
    for(int p = 0; p < P; p++){
      
      double m_xtp = (a[x] + ap(x,p)) + (b[x] + bp(x,p)) * k(t,p) + log(E(x, t, p));
      
      double term2 = exp(m_xtp);
      
      // a
      gr_loglik[0] += d(x, t, p) - term2;
      // b
      gr_loglik[1] += k(t,p) * d(x, t, p) - term2 * k(t, p);
      
      
    }
    
  }
  
  return gr_loglik;
}


// gradient of loglikelihood of k for a given t given kp

// [[Rcpp::export]]
double gr_loglik_k_t(
    int t,
    arma::mat ap,
    arma::mat bp,
    arma::vec k,
    arma::mat kp,
    arma::cube d,
    arma::cube E) {
  
  int X = E.n_rows;
  int P = E.n_slices;
  
  double gr_loglik = 0;
  
  for(int x = 0; x < X; x++){
    
    for(int p = 0; p < P; p++){
      
      double m_xtp = ap(x,p) + bp(x,p) * (k[t] + kp(t,p)) + log(E(x, t, p));
      
      double term2 = exp(m_xtp);
      
      gr_loglik += bp(x, p) * d(x, t, p) - term2 * bp(x, p);
      
    }
    
  }
  
  return gr_loglik;
}

// gradient of loglikelihood of k given kp

// [[Rcpp::export]]
arma::vec gr_loglik_k(
    arma::mat a,
    arma::mat b,
    arma::vec k0,
    arma::mat kp,
    arma::cube d,
    arma::cube E) {
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(Y);
  
  for(int t = 0; t < Y; t++){
    
    for(int x = 0; x < X; x++){
      
      for(int p = 0; p < P; p++){
        
        double m_xtp = a(x, p) + b(x, p) * (k0[t] + kp(t,p)) + log(E(x, t, p));
        
        double term2 = exp(m_xtp);
        
        gr_loglik[t] += b(x, p) * d(x, t, p) - term2 * b(x,p);
        
        
      }
      
    }
    
  }
  
  return gr_loglik;
}



// gradient of loglikelihood of ab given abp

// [[Rcpp::export]]
arma::vec gr_loglik_ab(
    arma::vec a,
    arma::vec b,
    arma::mat ap,
    arma::mat bp,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(2 * X);
  
  for(int x = 0; x < X; x++){
    
    for(int t = 0; t < Y; t++){
      
      for(int p = 0; p < P; p++){
        
        double m_xtp = (a[x] + ap(x,p)) + (b[x] + bp(x,p)) * k(t,p) + log(E(x, t, p));
        
        double term2 = exp(m_xtp);
        
        // a
        gr_loglik[x] += d(x, t, p) - term2;
        // b
        gr_loglik[X + x] += k(t,p) * d(x, t, p) - term2 * k(t, p);
        
        
      }
      
    }
    
  }
  
  return gr_loglik;
}

// loglikelihood for a specific x

// [[Rcpp::export]]
double loglik_abp(
    int x,
    arma::mat a,
    arma::mat b,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  double loglik = 0;
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  for(int t = 0; t < Y; t++){
    
    for(int p = 0; p < P; p++){
      
      double term1 = 
        d(x,t,p) * (a(x,p) + b(x,p) * k(t,p) + log(E(x,t,p)));
      double term2 = exp(a(x,p) + b(x,p) * k(t,p) + log(E(x,t,p)));
      loglik += term1 - term2;
      
    }
    
  }
  
  return loglik;
}

// gradient of axp for a specific x

// [[Rcpp::export]]
arma::vec gr_loglik_abp(
    int x,
    arma::vec a,
    arma::mat ap,
    arma::mat b,
    arma::mat k,
    arma::cube d,
    arma::cube E) {
  
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;

  arma::vec ap_grad = arma::zeros(P);
  
  for(int t = 0; t < Y; t++){
    
    for(int p = 0; p < P; p++){
      
      double term1 = 
        d(x,t,p) * (a[x] + ap(x,p) + b(x,p) * k(t,p) + log(E(x,t,p)));
      double term2 = exp(a[x] + ap(x,p) + b(x,p) * k(t,p) + log(E(x,t,p)));
      
      ap_grad[p] += d(x, t, p) - term2;
      
    }
    
  }
  
  return ap_grad;
}


