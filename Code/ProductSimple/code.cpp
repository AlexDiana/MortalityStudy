
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

// [[Rcpp::export]]
double loglik_cpm1_ax_cpp(
    arma::vec cp,
    arma::cube d,
    arma::cube E,
    arma::vec ax, 
    arma::cube m_terms){
  
  int X = d.n_rows;
  int Y = d.n_cols;
  int P = d.n_slices;
  
  double loglik = 0;
  
  for(int x = 0; x < X; x++){
    
    for(int t = 0; t < Y; t++){
      
      for(int p = 0; p < P; p++){
        
        if(arma::is_finite(d(x,t,p))){
          
          double mxtp = ax[x] * cp[p] + m_terms(x,t,p) + log(E(x,t,p));
          
          loglik += (d(x,t,p) * mxtp - exp(mxtp));  
          
        }
        
      }
      
    }
    
  }
  
  return(loglik);
  
}

// [[Rcpp::export]]
double loglik_ktm1_cp_cxp_cpp(
  arma::vec kt,
  arma::cube d, 
  arma::cube E, 
  arma::mat cxp, 
  arma::cube m_terms){
  
  int X = d.n_rows;
  int Y = d.n_cols;
  int P = d.n_slices;
  
  double loglik = 0;
  
  for(int x = 0; x < X; x++){
    
    for(int t = 0; t < Y; t++){
      
      for(int p = 0; p < P; p++){
        
        if(arma::is_finite(d(x,t,p))){
          
          double mxtp = kt[t] * cxp(t,p) + m_terms(x,t,p) + log(E(x,t,p));
          
          loglik += (d(x,t,p) * mxtp - exp(mxtp));  
          
        }
        
      }
      
    }
    
  }
  
 return loglik; 
}  

// arma::vec gr_loglik_ktm1_cp_cxp_cpp(
//     arma::vec kt,
//     arma::cube d, 
//     arma::cube E, 
//     arma::mat cxp, 
//     arma::cube m_terms){
//     
//     
//     int X = d.n_rows;
//   int Y = d.n_cols;
//   int P = d.n_slices;
//   
//   arma::vec grad = arma::zeros(Y - 1);
//   
//   for(int x = 0; x < X; x++){
//     
//     for(int t = 0; t < Y; t++){
//       
//       for(int p = 0; p < P; p++){
//         
//         if(t != (Y-1)){
//           
//           grad[t] += 
//           
//         }
//         
//       }
//       
//     }
//     
//   }
//     
//     kt <- c(ktm1, - sum(ktm1))
//     
//     ktcp <- matrix(kt, Y, P, byrow = F) * cxp
//     
//     ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))
//     
//     mxtp <- ktcp_mat + m_terms + log(E)
//     
//     term1 <- d * mxtp
//     term2 <- exp(mxtp)
//     
//     cxp_array <- aperm(array(cxp, dim = c(Y, P, X)), perm = c(3,1,2))
//     
//     term1_grad <- apply(d * cxp_array, 2, function(x){
//       sum(x, na.rm = T)
//     })
//     
//     term2_grad <- apply(exp(mxtp) * cxp_array, 2, function(x){
//       sum(x, na.rm = T)
//     })
//     
//     term1_1_grad <- term1_grad[1:(Y-1)]
//     term1_2_grad <- term1_grad[Y]
//     
//     term2_1_grad <- term2_grad[1:(Y-1)]
//     term2_2_grad <- term2_grad[Y]
//     
//     - (term1_1_grad - term2_1_grad - term1_2_grad + term2_2_grad)
//     
//     
//   }  
//   
// }

// gradient of loglikelihood of ap for a given x given a

// [[Rcpp::export]]
double loglik_gcm2_cpp(
    arma::vec gtx, 
    arma::cube& d, 
    arma::cube& E, 
    arma::cube& m_terms){
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  double loglik = 0;
  
  for(int x = 0; x < X; x++){
    
    for(int t = 0; t < Y; t++){
      
      for(int p = 0; p < P; p++){
        
        if(arma::is_finite(d(x,t,p))){
          
          double mxtp = gtx[(t - x) - ( - (X-1))] + m_terms(x,t,p) + log(E(x,t,p));
          
          loglik += -(d(x,t,p) * mxtp - exp(mxtp));  
          
        }
        
      }
      
    }
    
  }
  
  return loglik;
  
}

// [[Rcpp::export]]
arma::vec gr_loglik_gcm2_cpp(
    arma::vec gtx, 
    arma::cube& d, 
    arma::cube& E, 
    arma::cube& m_terms){
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec grad = arma::zeros(X + Y - 3);
  
  // int n_c = X + Y - 2;
  // 
  // for(int cc = 1; cc < n_c; cc++){
  //   
  //   for(int i = 0; i < ){
  //     
  //     for(int p = 0; p < P; p++){
  //       
  //       if(arma::is_finite(d(x,t,p))){
  //         
  //         double mxtp = gtx[cc] + m_terms(x,t,p) + log(E(x,t,p));
  //         
  //         if(cc < (n_c - 1)){
  //           
  //           grad[cc] += - (d(x,t,p) - exp(mxtp));
  //           
  //         } else {
  //           
  //           grad[cc] += (d(x,t,p) - exp(mxtp));
  //           
  //         }
  //         
  //       }
  //       
  //     }
  //     
  //   }
  //   
  // }
  
  for(int x = 0; x < X; x++){

    for(int t = 0; t < Y; t++){

      for(int p = 0; p < P; p++){

        if(arma::is_finite(d(x,t,p))){

          int idx_current = (t - x) - ( - (X-1));

          double mxtp = gtx[idx_current] + m_terms(x,t,p) + log(E(x,t,p));
          // double mxtp = gtx[(t - x) - ( - (X-1))] + m_terms(x,t,p) + log(E(x,t,p));

          // if(!(t == 0 & x == (X-1))){
          if(idx_current != 0){

            // if(!(t == (Y-1) & x == 0)){
            if(idx_current != (X+Y-2)){

              grad[idx_current - 1] += - (d(x,t,p) - exp(mxtp));
              // grad[(t - x) - ( - (X-1)) - 1] += - (d(x,t,p) - exp(mxtp));

            } else {

              for(int l = 0; l < grad.size(); l++){

                grad[l] += (d(x,t,p) - exp(mxtp));

              }

            }

          }

        }

      }

    }

  }
  
  return grad;
  
}

// [[Rcpp::export]]
arma::mat createTerm1(arma::vec gtx, 
                      arma::cube& d, 
                      arma::cube& E, 
                      arma::cube& m_terms){

  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;

  int n_c = X + Y - 1;

  arma::mat term12 = arma::zeros(n_c - 2, P);

  for(int x = 0; x < X; x++){

    for(int t = 0; t < Y; t++){

      for(int p = 0; p < P; p++){

        if(arma::is_finite(d(x,t,p))){

          int idx_current = (t - x) - ( - (X-1));

          double mxtp = gtx[idx_current] + m_terms(x,t,p) + log(E(x,t,p));
          // double mxtp = gtx[(t - x) - ( - (X-1))] + m_terms(x,t,p) + log(E(x,t,p));

          // if(!(t == 0 & x == (X-1))){
          if(idx_current != 0){

            // if(!(t == (Y-1) & x == 0)){
            if(idx_current != (X+Y-2)){

              term12(idx_current - 1, p) += 
                d(x,t,p) - exp(mxtp);

              // grad[idx_current - 1] += - (d(x,t,p) - exp(mxtp));
              // grad[(t - x) - ( - (X-1)) - 1] += - (d(x,t,p) - exp(mxtp));

            } 
            
          }

        }

      }

    }

  }

  return term12;
}

// // [[Rcpp::export]]
// arma::mat createTerm21(
//   arma::cube gtx_array,
//   arma::cube m_terms,
// ){
//   
//   
// }
// 
// term2_1 <-  sapply(1:P, function(p){
//   sapply(2:(X+Y-2), function(g){
//     - sum(sdiag(exp(m_terms[,,p] + gtx_array[,,p] + log(E[,,p])), k = g - X))
//   })
// })

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



// [[Rcpp::export]]
arma::vec gr_loglik_km1(
    arma::mat a,
    arma::mat b,
    arma::vec k0,
    arma::mat kp,
    arma::cube d,
    arma::cube E) {
  
  int X = E.n_rows;
  int Y = E.n_cols;
  int P = E.n_slices;
  
  arma::vec gr_loglik = arma::zeros(Y - 1);
  
  for(int x = 0; x < X; x++){
    
    for(int p = 0; p < P; p++){
      
      for(int t = 0; t < (Y-1); t++){
        
        if(arma::is_finite(d(x, t, p))){
          
          double m_xtp = a(x, p) + b(x, p) * (k0[t] + kp(t,p)) + log(E(x, t, p));
          
          double term2 = exp(m_xtp);
          
          gr_loglik[t] += b(x, p) * d(x, t, p) - term2 * b(x,p);
          
        }
        
        if(arma::is_finite(d(x, Y - 1, p))) {
          
          double m_xym1p = a(x, p) + b(x, p) *
            (k0[Y - 1] + kp(Y - 1, p)) + log(E(x, Y - 1, p));
          
          double term2m1 = exp(m_xym1p);
          
          gr_loglik[t] += - b(x, p) * d(x, Y - 1, p) - (-term2m1 * b(x,p));
          
          
        }
        
      }
      
    }
    
  }
  
  return (- gr_loglik);
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


