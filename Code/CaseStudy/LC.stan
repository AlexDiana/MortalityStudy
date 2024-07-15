data {
  // number of sites
  int<lower = 1> X;
  int<lower = 1> Y;
  
  // data
  int<lower = 0> d[X, Y];
  real<lower = 0> E[X, Y];
  
  // priors
  // real prior_beta_psi;
  // real<lower = 0> prior_beta_psi_sd;
  // real prior_beta_theta;
  // real<lower = 0> prior_beta_theta_sd;
  // real<lower = 0> prior_atheta0;
  // real<lower = 0> prior_btheta0;
  // real<lower = 0> prior_ap;
  // real<lower = 0> prior_bp;
  // real<lower = 0> prior_aq;
  // real<lower = 0> prior_bq;
  
}

parameters {
  
  vector[X] ax;
  vector[X - 1] bxm1;
  vector[Y - 1] ktm1;
  
}

transformed parameters {
  
  vector[X] bx;
  for(x in 1:(X-1)){
    bx[x] = bxm1[x];
  }
  bx[X] = 1 - sum(bxm1);
  
  vector[Y] kt;
  for(t in 1:(Y-1)){
    kt[t] = ktm1[t];
  }
  kt[Y] = - sum(ktm1);
  
}

model {
  
  for(x in 1:X){
    
    for(t in 1:Y){
      
      target += 
      d[x,t] * log(inv_logit(ax[x] + bx[x] * kt[t])) + 
      (E[x,t] - d[x,t]) * log(1 - inv_logit(ax[x] + bx[x] * kt[t]));
      
    }
    
  }
  
}
