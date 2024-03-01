data {
  // number of obs
  int<lower = 1> n;
  
  vector[n] E;
  vector[n] d;
  
  int<lower = 1> Y;
  int<lower = 1> A;
  int<lower = 1> P;
  
  int year[n];
  int age[n];
  int product[n];
}

parameters {
  vector[A] a;
  vector[A] b_params;
  vector[Y] k_params;
  vector[P - 1] p_params;
}

transformed parameters {
  
  vector[A] b = b_params / sum(b_params);
  vector[Y] k = k_params - mean(k_params);
  
  vector[P] p;
  p[1] = 0;
  for(i in 1:(P-1)){
    p[i + 1] = p_params[i];
  }
  
  vector[n] logitqx;

  for (i in 1:n) {
    
    logitqx[i] = a[age[i]] + b[age[i]] * k[year[i]] + p[product[i]];
    // logitqx[i] = a[age[i]] + b[age[i]] * k[year[i]];
    
  }  
  
}

model {
  
  for (i in 1:n) {
    
      target += 
        lgamma(E[i] + 1) - lgamma(d[i] + 1) - lgamma(E[i] - d[i] + 1) + 
        d[i] * ( - log(1 + exp(-logitqx[i]))) + 
        (E[i] - d[i]) * ( - log(1 + exp(logitqx[i])));
    
  }
  
}
