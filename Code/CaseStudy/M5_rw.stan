data {
  // number of sites
  int<lower = 1> X;
  int<lower = 1> Y;
  
  // data
  int<lower = 0> d[X, Y];
  real<lower = 0> E[X, Y];
  
}

parameters {
  
  vector[X] ax;
  vector[Y - 1] ktm1;
  vector[X + Y - 3] gtxm1XY;
  
  real alpha;
  real<lower = 0> sigmak;
}

transformed parameters {
  
  vector[Y] kt;
  for(t in 1:(Y-1)){
    kt[t] = ktm1[t];
  }
  kt[Y] = - sum(ktm1);
  
  vector[X + Y - 1] gtx;
  gtx[1] = 0;
  for(g in 2:(X + Y - 2)){
    gtx[g] = gtxm1XY[g-1];
  }
  gtx[X + Y - 1] = - sum(gtxm1XY);
  
}

model {
  
  // prior
  
  for(t in 2:Y){
    
    target += normal_lpdf(kt[t] | kt[t - 1] + alpha, sigmak);
    
  }
  
  // likelihood
  
  for(x in 1:X){
    
    for(t in 1:Y){
      
      target += 
        d[x,t] * log(inv_logit(ax[x] + kt[t] + gtx[x - t + Y])) + 
        (E[x,t] - d[x,t]) * log(1 - inv_logit(ax[x] + kt[t] + gtx[x - t + Y]));
      
    }
    
  }
  
}
