
rm(list = ls())

library(Rcpp)
library(RcppArmadillo) 
library(MASS)
library(here)
library(ggplot2)
library(emdbook)
library(optimx)
library(beepr)
library(tidyverse)

sourceCpp("Code/MCMC/code.cpp")
source("Code/MCMC/functions.R")

# SIMULATED DATA ----------

X <- 20 # 30 # ages
Y <- 20 # 36 #  years
P <- 3 # products

ages <- 50 + 1:X
years <- 2000 + 1:Y
products <- as.character(1:P)

E <- array(rpois(X * Y * P, lambda = 5000), dim = c(X, Y, P))

# a
# ax <- seq(-5, -2, length.out = X)
{
  ax <- c(-4.0512213665814, -3.97704124624414, -3.86830606439788, -3.80357367205838,
          
          -3.73200080978762, -3.63738818539699, -3.56804053625705, -3.48365064168728,
          
          -3.39115995197912, -3.31902914229022, -3.22065518527815, -3.15420054660133,
          
          -3.05648688210544, -2.97347316787011, -2.89790411605838, -2.8169553534215,
          
          -2.73895895602344, -2.66346592013173, -2.5943840324644, -2.50680577080775,
          
          -2.3881944793569, -2.32899572754395, -2.24315164400745, -2.16404275881499,
          
          -2.08226515858658, -1.98726563053034, -1.90359661305062, -1.82710159710464,
          
          -1.7396226616441, -1.69860290672302)
  ax <- ax[1:X]
}

# b
# bx <- seq(.8, 0.0, length.out = X)
{
  bx <- c(0.0486060237848308, 0.0447434606869872, 0.0454455937103952,
        
        0.0455343703689019, 0.0450371958284298, 0.045971179464402, 0.0423557607243807,
        
        0.0430873789219571, 0.0407579553985572, 0.0400763401483252, 0.042290002346062,
        
        0.038102592760602, 0.0415188866516725, 0.0371060614373157, 0.0375132576452101,
        
        0.0377232782861309, 0.0352622239546568, 0.0333521889560757, 0.0321508266831311,
        
        0.0321165535952561, 0.0321842233147083, 0.0242184597077586, 0.0237268461855173,
        
        0.0220992202727715, 0.0204801810055311, 0.0164336990044118, 0.0145874922698763,
        
        0.0141527451412956, 0.015237830033942, 0.00812817171090743)
  bx <- bx[1:X]
}

# k
# eps_kt <- rnorm(Y - 1, sd = 1)
# kt <- -11 + c(0, cumsum(eps_kt))
{
  kt <- c(7.30499317745483, 6.73026076213329, 6.37887599533382, 5.97974189114322,

    6.26735910848851, 5.84441416328068, 4.90912390629458, 4.04371301876798,

    3.70435061227499, 3.01335249062068, 2.77082241926598, 1.85918368563587,

    2.03867821280429, 1.44432601169222, 0.97303957101614, 1.08339768889609,

    0.70907265176869, 0.706205461050796, 0.247195664577536, -0.171692466944422,

    -0.343415704243637, -1.27315598510262, -1.83615114239075, -2.32013134814312,

    -2.95894229077233, -2.67850042481684, -3.3179752530899, -3.68705473728354,

    -4.1176429375629, -4.66043301282831, -5.11411187432859, -5.39085145670421,

    -6.17501134761908, -6.79701732830229, -7.21253056233688, -7.95348862003079

  )
  
  kt <- kt[1:Y]
}

# transform the params so that sum b = 1 and sum k = 0
{
  sumbx <- sum(bx)
  bx <- bx / sumbx
  kt <- kt * sumbx
  meankt <- mean(kt)
  kt <- kt - meankt
  ax <- ax + bx * meankt
  
}

a_true <- ax
b_true <- bx
k_true <- kt

# ap
# model NP
{
  ap_true <- matrix(rnorm(X * P, sd = .1), X, P)
  ap_true <- t(apply(ap_true, 1, function(x){
    x - mean(x)
  }))
}
# model A
{
  ap_0 <- rnorm(P, sd = .5)
  ap_0 <- ap_0 - mean(ap_0)
  ap_true <- matrix(ap_0, X, P, byrow = T)
}
# model B
{
  # ap_true <- matrix(0, X, P)
}

# b
bp_true <- matrix(rnorm(X * P, sd = .05), X, P)
bp_true <- t(apply(bp_true, 1, function(x){
x - mean(x)
}))
bp_true <- matrix(0, X, P)

# k
kp_true <- matrix(rnorm(Y * P, sd = 1), Y, P)
kp_true <- t(apply(kp_true, 1, function(x){
x - mean(x)
}))
kp_true <- matrix(0, Y, P)

d <- array(NA, dim = c(X, Y, P))
m <- array(NA, dim = c(X, Y, P))
for (x in 1:X) {
  for (t in 1:Y) {
    for (p in 1:P) {
      m[x,t,p] <- exp(a_true[x] + ap_true[x,p] + 
                        (b_true[x] + bp_true[x,p]) * 
                        (k_true[t] + kp_true[t,p]))
      d[x,t,p] <- rpois(1, E[x,t,p] * m[x,t,p])  
    }
  }
}

# REAL DATA ------

load(here("Data","data_products.rda"))

X <- dim(d)[1]
Y <- dim(d)[2]
P <- dim(d)[3]

ages <- dimnames(d)[[1]]
years <- dimnames(d)[[2]]
products <- dimnames(d)[[3]]

# MCMC -----

nburn <- 500
niter <- 500

# starting values
{
  # true
  {
    a <- a_true
    b <- b_true
    k <- k_true

    # ap <- ap_true
    bp <- bp_true
    kp <- kp_true
    
    gamma_a <- "B"
    gamma_b <- "B"
    gamma_k <- "B"
    
  }
  
  # zero
  {
    # a <- rep(0, X)
    # b <- rep(0, X)
    # k <- rep(0, Y)
    # ap <- matrix(0, X, P)
    # bp <- matrix(0, X, P)
    # kp <- matrix(0, Y, P)
  }
  
  # from MLE
  {
    laplace_fit <- optim(
      par = c(
        rep(0, X),
        rep(1 / X, X),
        rep(0, Y)),
      fn = loglik_r(d, E),
      method = c("BFGS")
    )

    param_star <- laplace_fit$par
    a <- param_star[1:X]
    b <- param_star[X + 1:X]
    k <- param_star[2 *X + 1:Y]

    sumb <- sum(b)
    b <- b / sumb
    k <- k * sumb
    meank <- mean(k)
    k <- k - meank
    a <- a + b * meank

    ap <- matrix(0, X, P)
    bp <- matrix(0, X, P)
    kp <- matrix(0, Y, P)
    
    gamma_a <- "B"
    gamma_b <- "B"
    gamma_k <- "B"
  }
  
}

# prior
{
  a_configs <- c("B","A","NP")
  nConfigs_a <- length(a_configs)
  p_ap <- c(.6, .3, .1)
  p_bp <- .1
  p_kp <- .1
  sd_ap <- .1
  sd_bp <- 1
  sd_kp <- 2
}

# output 
{
  a_output <- matrix(NA, niter, X)
  b_output <- matrix(NA, niter, X)
  k_output <- matrix(NA, niter, Y)
  ap_output <- array(NA, dim = c(niter, X, P))
  bp_output <- array(NA, dim = c(niter, X, P))
  kp_output <- array(NA, dim = c(niter, Y, P))
  gamma_a_output <- rep(NA, niter)
  gamma_b_output <- rep(NA, niter)
  gamma_k_output <- rep(NA, niter)
  loglik_output <- rep(NA, niter)
}

# paramsUpdate
{
  update_k0 <- T
  update_kp <- T
  update_ab <- T
  update_ap <- T
  update_bp <- T
  update_gammaa <- T
  update_gammab <- F
  update_gammak <- F
}

for (iter in 1:(niter + nburn)) {
  
  print(gamma_a)
  
  if(iter %% 10 == 0){
    if(iter < nburn){
      print(paste0("Burn-in Iteration = ",iter))  
    } else {
      print(paste0("Iteration = ",iter - nburn))  
    }  
  }
  
  # update k0
  
  if (update_k0) {
    
    # add product effect to ax and bx
    a2 <- matrix(a, X, P, byrow = F) + ap
    b2 <- matrix(b, X, P, byrow = F) + bp

    # create loglikelihood function of k given ax and bx
    loglik_k_current <- loglik_k0_r(d, E, a2, b2, kp)

    # find maximum over all variables together
    {
      # # find maximum
      # laplace_fit <- optim(
      #   par = k,
      #   fn = loglik_k_current,
      #   method = c("BFGS"),
      #   gr = gr_loglik_k0_r(d, E, a2, b2, kp),
      #   hessian = T
      # )
      # 
      # k_star <- laplace_fit$par
      # H_star <- laplace_fit$hessian
      # 
      # Sigma_star <- 2 * solve(H_star)
    }
    
    # find maximum for each variabe individually (it's the same as independent)
    {
        k_star <- rep(NA, Y)
        Sigma_star <- matrix(0, Y, Y)
        
        for (t in 1:Y) {
          # print(t)
          loglik_kt_current <- loglik_kt_r(t, d, E, a2, b2, k, kp)
          gr_loglik_kt_current <- gr_loglik_kt_r(t, d, E, a2, b2, k, kp)
          
          # find maximum
          
          # microbenchmark({
            laplace_fit <- optim(
              par = k[t],
              lower = k[t] - 1,
              upper = k[t] + 1,
              fn = loglik_kt_current,
              method = c("Brent"),
              # gr = gr_loglik_kt_current,
              hessian = T
            )
          # },{
          #   laplace_fit <- optim(
          #     par = k[t],
          #     # par = c(0, 0),
          #     fn = loglik_kt_current,
          #     method = c("BFGS"),
          #     gr = gr_loglik_kt_current,
          #     hessian = T
          #   )
          # }, times = 2)
          
          
          k_star[t] <- laplace_fit$par
          Sigma_star[t, t] <- solve(laplace_fit$hessian)
          
        }
        
        Sigma_star <- 2 * Sigma_star
      
    }
    
    k_proposed <- as.vector(sampleMVNconstraint_k(k_star, Sigma_star))

    loglik_star <- - loglik_k_current(k_proposed)
    loglik_current <- - loglik_k_current(k)

    # compute proposal for all variables together
    {
      # logproposal_star <- dmvnorm(as.vector(k_proposed), as.vector(k_star), Sigma_star, log = T)
      # logproposal_current <- dmvnorm(as.vector(k), as.vector(k_star), Sigma_star, log = T)
    }
    
    # find proposal for each variable individually and sum
    {
        logproposal_star <- 0
        logproposal_current <- 0
        
        for (t in 1:Y) {
          
          logproposal_star <- logproposal_star + 
            dnorm(k_proposed[t], k_star[t], sqrt(Sigma_star[t,t]), log = T)
       
          logproposal_current <- logproposal_current + 
            dnorm(k[t], k_star[t], sqrt(Sigma_star[t,t]), log = T)
       
        }
      
    }
    
    mh_ratio <- exp(loglik_star - loglik_current + 
                      logproposal_current  - logproposal_star)

    if(runif(1) < mh_ratio){
      k <- as.vector(k_proposed)
    }
    
  }
  
  # update k_p
  
  if (update_kp) {
    
    if(gamma_k == 2){ # if the product effect is present for this variable
      
      # add product effect to ax and bx
      a2 <- matrix(a, X, P, byrow = F) + ap
      b2 <- matrix(b, X, P, byrow = F) + bp
      
      for (t in 1:Y) {
        
        loglik_kpt_current <- loglik_kpt_r(t, a2, b2, k, kp, d, E)
        gr_loglik_kpt_current <- gr_loglik_kpt_r(t, a2, b2, k, kp, d, E)
        
        # find maximum
        laplace_fit <- optim(
          par = kp[t,],
          # par = rep(0, P),
          fn = loglik_kpt_current,
          gr = gr_loglik_kpt_current,
          method = c("BFGS"),
          hessian = T
        )
        
        kp_star <- laplace_fit$par # mle
        H_star <- laplace_fit$hessian # hessian at the mle
        
        # I double the variance of the proposal to increase acceptance rate
        Sigma_star <- 2 * solve(H_star) 
      
        # propose kp subject to the constraint that sum kp = 0
        kp_proposed <- sampleMVNconstraint_k(kp_star, Sigma_star)
        
        loglik_proposal <- - loglik_kpt_current(kp_proposed)
        loglik_current <- - loglik_kpt_current(kp[t,])
        
        logproposal_star <- dmvnorm(as.vector(kp_proposed), as.vector(kp_star), Sigma_star, log = T)
        logproposal_current <- dmvnorm(as.vector(kp[t,]), as.vector(kp_star), Sigma_star, log = T)
        
        logprior_star <- sum(dnorm(as.vector(kp_proposed), 0, sd_kp, log = T))
        logprior_current <- sum(dnorm(as.vector(kp[t,]), 0, sd_kp, log = T))
        
        mh_ratio <- exp(loglik_proposal + logproposal_current + 
                          logprior_star - logprior_current - 
                          loglik_current - logproposal_star)
        
        if(runif(1) < mh_ratio){
          kp[t,] <- kp_proposed
        }
        
      }
      
    }
    
  }
  
  # update a and b jointly
  
  if (update_ab) {
    
    # add product effect to k
    k2 <- matrix(k, Y, P, byrow = F) + kp

    # find maximum over all variables together
    {
      loglik_ab_current <- loglik_ab_r(d, E, ap, bp, k2)
      gr_loglik_ab_current <- gr_loglik_ab_r(d, E, ap, bp, k2)
      # 
      # # find maximum
      # laplace_fit <- optim(
      #   par = c(a, b),
      #   fn = loglik_ab_current,
      #   method = c("BFGS"),
      #   gr = gr_loglik_ab_current,
      #   hessian = T
      # )
      # 
      # ab_star <- laplace_fit$par
      # H_star <- laplace_fit$hessian
      # 
      # Sigma_star <- solve(H_star)
    }
    
    # find maximum for each variabe individually (it's the same as independent)
    {
      ab_star <- rep(NA, 2 * X)
      Sigma_star <- matrix(0, 2 * X, 2 * X)

      for (x in 1:X) {

        loglik_abx_current <- loglik_abx_r(x, d, E, a, b, ap, bp, k2)
        gr_loglik_abx_current <- gr_loglik_abx_r(x, d, E, a, b, ap, bp, k2)

        # find maximum
        laplace_fit <- optim(
          par = c(a[x], b[x]),
          # par = c(0, 0),
          fn = loglik_abx_current,
          method = c("BFGS"),
          gr = gr_loglik_abx_current,
          hessian = T
        )

        ab_star[x + c(0, X)] <- laplace_fit$par
        Sigma_star[x + c(0, X), x + c(0, X)] <- solve(laplace_fit$hessian)

      }
    }
    
    Sigma_star <- 2 * Sigma_star
    
    # loglik_ab_current(c(a,b))
    # gr_loglik_ab_r(d, E, ap, bp, k2)(c(a,b))
    #
    # loglik_ab_current(laplace_fit$par)
    # gr_loglik_ab_r(d, E, ap, bp, k2)(laplace_fit$par)

    ab_proposed <- sampleMVNconstraint_ab(ab_star, Sigma_star)
    a_proposed <- ab_proposed[1:X]
    b_proposed <- ab_proposed[X + 1:X]

    loglik_star <- - loglik_ab_current(ab_proposed)
    loglik_current <- - loglik_ab_current(c(a, b))
    # - loglik_ab_current(c(a_true, b_true))

    # compute proposal for all variables together
    {
      # logproposal_star <- dmvnorm(as.vector(ab_proposed), as.vector(ab_star), Sigma_star, log = T)
      # logproposal_current <- dmvnorm(c(a,b), as.vector(ab_star), Sigma_star, log = T)
    }

    # find maximum for each variabe individually (it's the same as independent)
    {
      logproposal_star <- 0
      logproposal_current <- 0
      
      for (x in 1:X) {
        
        logproposal_star <- logproposal_star + 
          dmvnorm(as.vector(ab_proposed[x + c(0, X)]), 
                  as.vector(ab_star[x + c(0, X)]), 
                  Sigma_star[x + c(0, X), x + c(0, X)], log = T)
        
        logproposal_current <- logproposal_current + 
          dmvnorm(as.vector(c(a,b)[x + c(0, X)]), 
                  as.vector(c(a,b)[x + c(0, X)]), 
                  Sigma_star[x + c(0, X), x + c(0, X)], log = T)
      }
    }
    

    mh_ratio <- exp(loglik_star - loglik_current + 
                      logproposal_current -  logproposal_star)
    
    if(runif(1) < mh_ratio){
      a <- a_proposed
      b <- b_proposed
    }
    
  }
  
  # update a_p
  
  if (update_ap) {
    
    if(gamma_a != "B"){
      
      k2 <- matrix(k, Y, P, byrow = F) + kp
      b2 <- matrix(b, X, P, byrow = F) + bp
      
      if(gamma_a == "A"){ # currently in additive model
        
        loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
        gr_loglik_ap_current <- gr_loglik_ap_r(a, b2, d, E, k2)
        
        # find maximum
        laplace_fit <- optim(
          # par = ap[x,],
          par = rep(0, P),
          fn = loglik_ap_current,
          method = c("BFGS"),
          gr = gr_loglik_ap_current,
          hessian = T
        )
        
        ap_star <- laplace_fit$par
        H_star <- laplace_fit$hessian
        
        Sigma_star <- 2 * solve(H_star)
        
        ap_proposed <- sampleMVNconstraint_k(ap_star, Sigma_star)
        
        loglik_proposal <- - loglik_ap_current(ap_proposed)
        loglik_current <- - loglik_ap_current(ap[x,])
        
        logproposal_star <- dmvnorm(as.vector(ap_proposed), as.vector(ap_star), 
                                    Sigma_star, log = T)
        logproposal_current <- dmvnorm(as.vector(ap[x,]), as.vector(ap_star), 
                                       Sigma_star, log = T)
        
        logprior_star <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
        logprior_current <- sum(dnorm(as.vector(ap[x,]), 0, sd_ap, log = T))
        
        mh_ratio <- exp(loglik_proposal + logproposal_current +
                          logprior_star - logprior_current - 
                          loglik_current - logproposal_star)
        # print(mh_ratio)
        if(runif(1) < mh_ratio){
          ap <- matrix(ap_proposed, X, P, byrow = T)
        }
        
      } else if(gamma_a == "NP"){ # nonparametric model
        
        for (x in 1:X) {
          
          loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k2)
          gr_loglik_apx_current <- gr_loglik_apx_r(x, a, ap, b2, d, E, k2)
          
          # find maximum
          laplace_fit <- optim(
            # par = ap[x,],
            par = rep(0, P),
            fn = loglik_apx_current,
            method = c("BFGS"),
            gr = gr_loglik_apx_current,
            hessian = T
          )
          
          ap_star <- laplace_fit$par
          H_star <- laplace_fit$hessian
          
          Sigma_star <- 2 * solve(H_star)
          
          ap_proposed <- sampleMVNconstraint_k(ap_star, Sigma_star)
          
          loglik_proposal <- - loglik_apx_current(ap_proposed)
          loglik_current <- - loglik_apx_current(ap[x,])
          
          logproposal_star <- dmvnorm(as.vector(ap_proposed), as.vector(ap_star), 
                                      Sigma_star, log = T)
          logproposal_current <- dmvnorm(as.vector(ap[x,]), as.vector(ap_star), 
                                         Sigma_star, log = T)
          
          logprior_star <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
          logprior_current <- sum(dnorm(as.vector(ap[x,]), 0, sd_ap, log = T))
          
          mh_ratio <- exp(loglik_proposal + logproposal_current +
                            logprior_star - logprior_current - 
                            loglik_current - logproposal_star)
          # print(mh_ratio)
          if(runif(1) < mh_ratio){
            ap[x,] <- ap_proposed
          }
          
        }
        
      }
      
    }
    
  }
  
  # update b_p
  
  if (update_bp) {
    
    if(gamma_b == 2){
      
      k2 <- matrix(k, Y, P, byrow = F) + kp
      a2 <- matrix(a, X, P, byrow = F) + ap
      
      for (x in 1:X) {
        
        loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2)
        gr_loglik_bpx_current <- gr_loglik_bpx_r(x, a2, b, bp, d, E, k2)
        
        # find maximum
        laplace_fit <- optim(
          # par = bp[x,],
          par = rep(0, P),
          fn = loglik_bpx_current,
          method = c("BFGS"),
          gr = gr_loglik_bpx_current,
          hessian = T
        )
        
        bp_star <- laplace_fit$par
        H_star <- laplace_fit$hessian
        
        Sigma_star <- 2 * solve(H_star)
        
        bp_proposed <- sampleMVNconstraint_k(bp_star, Sigma_star)
        
        loglik_proposal <- - loglik_bpx_current(bp_proposed)
        loglik_current <- - loglik_bpx_current(bp[x,])
        
        logproposal_star <- dmvnorm(as.vector(bp_proposed), as.vector(bp_star), Sigma_star, log = T)
        logproposal_current <- dmvnorm(as.vector(bp[x,]), as.vector(bp_star), Sigma_star, log = T)
        
        logprior_star <- sum(dnorm(as.vector(bp_proposed), 0, sd_bp, log = T))
        logprior_current <- sum(dnorm(as.vector(bp[x,]), 0, sd_bp, log = T))
        
        mh_ratio <- exp(loglik_proposal + logproposal_current +
                          logprior_star - logprior_current - 
                          loglik_current - logproposal_star)
        # print(mh_ratio)
        if(runif(1) < mh_ratio){
          bp[x,] <- bp_proposed
        }
        
      }
      
    }
    
  }
  
  # propose gamma_a
  
  if (update_gammaa) {
    
    k2 <- matrix(k, Y, P, byrow = F) + kp
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    logPriorCurrent <- log(p_ap[which(a_configs == gamma_a)])
    
    aap <- matrix(a, X, P, byrow = F) + ap
    loglik_current <- loglik_LCp_r(aap, b2, k2, d, E)
    
    if(gamma_a == "B") { # currently in base model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
    } else if(gamma_a == "A"){ # currently in additive model
    
      logProposal_current <- 0
      
      {
      #   loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
      # gr_loglik_ap_current <- gr_loglik_ap_r(a, b2, d, E, k2)
      # 
      # # find maximum
      # laplace_fit <- optim(
      #   par = rep(0, P),
      #   fn = loglik_ap_current,
      #   method = c("BFGS"),
      #   gr = gr_loglik_ap_current,
      #   hessian = T
      # )
      # 
      # ap_star <- laplace_fit$par
      # H_star <- laplace_fit$hessian
      # 
      # Sigma_star <- 2 * solve(H_star)
      }
      
      list_ap_proposal <- findProposalApAdditive(a, b2, d, E, k2)
      ap_star <- list_ap_proposal$ap_star
      Sigma_star <- list_ap_proposal$Sigma_star
      
      ap_proposed_pm1 <- ap[x,1:(P-1)]
      ap_proposed <- c(ap_proposed_pm1, - sum(ap_proposed_pm1))
      
      logProposal_current <- dmvnorm(as.vector(ap_proposed_pm1), 
                                  as.vector(ap_star[1:(P-1)]), 
                                  as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)
      
      logprior_param_current <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
      
      
    } else if(gamma_a == "NP"){ # in nonparametric model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
      for (x in 1:X) {
        
        {
          
          # loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k2)
          # gr_loglik_apx_current <- gr_loglik_apx_r(x, a, ap, b2, d, E, k2)
          # 
          # # find maximum
          # laplace_fit <- optim(
          #   # par = ap[x,],
          #   par = rep(0, P),
          #   fn = loglik_apx_current,
          #   method = c("BFGS"),
          #   gr = gr_loglik_apx_current,
          #   hessian = T
          # )
          # 
          # ap_star <- laplace_fit$par
          # H_star <- laplace_fit$hessian
          # Sigma_star <- 2 * solve(H_star)
        }
      
        list_ap_proposal <- findProposalApXNonparametrics(x, a, ap, b2, d, E, k2)
        ap_star <- list_ap_proposal$ap_star
        Sigma_star <- list_ap_proposal$Sigma_star
        
        logProposal_current <- logProposal_current + 
          dmvnorm(as.vector(ap[x,1:(P-1)]), 
                  as.vector(ap_star[1:(P-1)]), 
                  as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)
        
        logprior_param_current <- logprior_param_current + 
          sum(dnorm(as.vector(ap[x,]), 0, sd_ap, log = T))
        
      }
      
    }

    newConfig <- sample(setdiff(a_configs, gamma_a), 1)
    
    logPriorStar <- log(p_ap[which(a_configs == newConfig)])
    
    if(newConfig == "B"){ # proposing base model
      
      ap_proposed <- matrix(0, X, P)
      
      logProposal_star <- 0
      logprior_param_star <- 0
      
    } else if (newConfig == "A"){ # proposing additive model
      
      list_ap_proposal <- findProposalApAdditive(a, b2, d, E, k2)
      ap_star <- list_ap_proposal$ap_star
      Sigma_star <- list_ap_proposal$Sigma_star
      
      ap_proposed_pm1 <- mvrnorm(
        n = 1,
        ap_star[1:(P-1)],
        as.matrix(Sigma_star[1:(P-1),1:(P-1)])
      )
      
      ap_proposed <- c(ap_proposed_pm1, - sum(ap_proposed_pm1))
      
      logProposal_star <- dmvnorm(as.vector(ap_proposed_pm1), 
                                  as.vector(ap_star[1:(P-1)]), 
                                  as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)
      
      logprior_param_star <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
      
      ap_proposed <- matrix(ap_proposed, X, P, byrow = T)
      
    } else if (newConfig == "NP"){ # proposing nonparametric model
      
      loglik_star <- 0
      logProposal_star <- 0
      logprior_param_star <- 0
      
      ap_proposed <- matrix(NA, X, P)
      
      for (x in 1:X) {
        
        list_ap_proposal <- findProposalApXNonparametrics(x, a, ap, b2, d, E, k2)
        ap_star <- list_ap_proposal$ap_star
        Sigma_star <- list_ap_proposal$Sigma_star
        
        apx_proposed_pm1 <- mvrnorm(
          n = 1,
          ap_star[1:(P-1)],
          as.matrix(Sigma_star[1:(P-1),1:(P-1)])
        )
        
        ap_proposed[x,] <- c(apx_proposed_pm1, - sum(apx_proposed_pm1)) 
        
        loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
        
        logProposal_star <- logProposal_star + 
          dmvnorm(as.vector(apx_proposed_pm1), 
                  as.vector(ap_star[1:(P-1)]), 
                  as.matrix(Sigma_star[1:(P-1),1:(P-1)]), 
                  log = T)
      
      }
      
    }
    
    aap_proposed <- matrix(a, X, P, byrow = F) + ap_proposed
    loglik_proposal <- loglik_LCp_r(aap_proposed, b2, k2, d, E)
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current +  # likelihood 
        logPriorStar - logPriorCurrent + # prior of the model
        logprior_param_star - logprior_param_current + # prior of the parameters
        logProposal_current - logProposal_star # proposal of the parameters
        )
    
    if(runif(1) < mh_ratio){
      gamma_a <- newConfig
      ap <- ap_proposed
    }
    
    #### OLD CODE
    
    {
      # ap_star_all <- matrix(NA, X, P)
      # Sigma_star_all <- array(NA, dim = c(X, P, P))
      # 
      # for (x in 1:X) {
      #   
      #   loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k2)
      #   gr_loglik_apx_current <- gr_loglik_apx_r(x, a, ap, b2, d, E, k2)
      #   
      #   # find maximum
      #   laplace_fit <- optim(
      #     # par = ap[x,],
      #     par = rep(0, P),
      #     fn = loglik_apx_current,
      #     method = c("BFGS"),
      #     gr = gr_loglik_apx_current,
      #     hessian = T
      #   )
      #   
      #   ap_star_all[x,] <- laplace_fit$par
      #   Sigma_star_all[x,,] <- solve(laplace_fit$hessian)
      #   
      # }
      # 
      # if(gamma_a == "NP"){ # currently in nonparametric model
      #   
      #   ap_proposed <- ap
      #   
      #   loglik_proposal <- 0
      #   loglik_current <- 0
      #   logproposal <- 0
      #   
      #   # propose new ap
      #   for (x in 1:X) {
      #     
      #     loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k2)
      #     
      #     loglik_proposal <- loglik_proposal - loglik_apx_current(ap_proposed[x,])
      #     loglik_current <- loglik_current - loglik_apx_current(rep(0, P))
      #     
      #     logproposal <- logproposal + 
      #       dmvnorm(as.vector(ap_proposed[x,]), 
      #               as.vector(ap_star_all[x,]), 
      #               Sigma_star_all[x,,], log = T)
      #     
      #   }
      #   
      #   logprior_proposal <- log(p_ap)
      #   logprior_current <- log(1 - p_ap)
      #   
      #   mh_ratio <- 1 / exp(loglik_proposal - loglik_current + 
      #                         logprior_proposal - logprior_current - 
      #                         logproposal)
      #   
      #   if(runif(1) < mh_ratio){
      #     
      #     gamma_a <- 0
      #     ap <- matrix(0, X, P)
      #     
      #   }
      #   
      # } else { # currently in base model
      #   
      #   ap_proposed <- matrix(NA, X, P)
      #   
      #   loglik_proposal <- 0
      #   loglik_current <- 0
      #   logproposal <- 0
      #   
      #   # propose new ap
      #   for (x in 1:X) {
      #     
      #     apx_proposed_pm1 <- mvrnorm(
      #       n = 1,
      #       ap_star_all[x,1:(P-1)],
      #       as.matrix(Sigma_star_all[x,1:(P-1),1:(P-1)])
      #     )
      #     
      #     ap_proposed[x,] <- c(apx_proposed_pm1, - sum(apx_proposed_pm1)) 
      #     
      #     # ap_proposed[x,] <- sampleMVNconstraint_k(apx_proposed, 
      #     #                                          Sigma_star_all[x,,])
      #     
      #   }
      #   
      #   # compute ap lik
      #   for (x in 1:X) {
      #     
      #     loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k2)
      #     
      #     loglik_proposal <- loglik_proposal - loglik_apx_current(ap_proposed[x,])
      #     loglik_current <- loglik_current - loglik_apx_current(rep(0, P))
      #     
      #     logproposal <- logproposal + 
      #       dmvnorm(as.vector(ap_proposed[x,1:(P-1)]), 
      #               as.vector(ap_star_all[x,1:(P-1)]), 
      #               as.matrix(Sigma_star_all[x,1:(P-1),1:(P-1)]), log = T)
      #     
      #   }
      #   
      #   logprior_proposal <- log(p_ap)
      #   logprior_current <- log(1 - p_ap)
      #   
      #   mh_ratio <- exp(loglik_proposal - loglik_current + 
      #                     logprior_proposal - logprior_current - 
      #                     logproposal)
      #   
      #   if(runif(1) < mh_ratio){
      #     
      #     gamma_a <- 1
      #     ap <- ap_proposed
      #     
      #   }
      #   
      #   
      # } 
    }
    
  }
  
  # propose gamma_b
  
  if (update_gammab) {
    
    k2 <- matrix(k, Y, P, byrow = F) + kp
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    bp_star_all <- matrix(NA, X, P)
    Sigma_star_all <- array(NA, dim = c(X, P, P))
    
    for (x in 1:X) {
      
      loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2)
      gr_loglik_bpx_current <- gr_loglik_bpx_r(x, a2, b, bp, d, E, k2)
      
      # find maximum
      laplace_fit <- optim(
        # par = bp[x,],
        par = rep(0, P),
        fn = loglik_bpx_current,
        method = c("BFGS"),
        gr = gr_loglik_bpx_current,
        hessian = T
      )
   
     # laplace_fit <- optimx(
     #    # par = bp[x,],
     #    par = rep(0, P),
     #    fn = loglik_bpx_current,
     #    method = c("BFGS"),
     #    gr = gr_loglik_bpx_current,
     #    hessian = T
     #  )
      
      bp_star_all[x,] <- laplace_fit$par
      Sigma_star_all[x,,] <- solve(laplace_fit$hessian)
      
      # - loglik_bpx_current(bp_star_all[x,])
      # - loglik_bpx_current(bp_true[x,])
      # - loglik_bpx_current(bp_proposed[x,])
      
    }
    
    if(gamma_b == 1){ # currently in nonparametric model
      
      bp_proposed <- bp
      
      loglik_proposal <- 0
      loglik_current <- 0
      logproposal <- 0
      
      # propose new bp
      for (x in 1:X) {
        
        loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2)
        
        loglik_proposal <- loglik_proposal - loglik_bpx_current(bp_proposed[x,])
        loglik_current <- loglik_current - loglik_bpx_current(rep(0, P))
        
        logproposal <- logproposal + 
          dmvnorm(as.vector(bp_proposed[x,]), 
                  as.vector(bp_star_all[x,]), 
                  Sigma_star_all[x,,], log = T)
        
      }
      
      logprior_proposal <- log(p_bp)
      logprior_current <- log(1 - p_bp)
      
      mh_ratio <- 1 / exp(loglik_proposal - loglik_current + 
                        logprior_proposal - logprior_current - 
                        logproposal)
      
      if(runif(1) < mh_ratio){
        
        gamma_b <- 0
        bp <- matrix(0, X, P)
        
      }
      
    } else { # currently in base model
      
      bp_proposed <- matrix(NA, X, P)
      
      loglik_proposal <- 0
      loglik_current <- 0
      logproposal <- 0
      
      # propose new bp
      for (x in 1:X) {
      
        bpx_proposed_pm1 <- mvrnorm(
          n = 1,
          bp_star_all[x,1:(P-1)],
          as.matrix(Sigma_star_all[x,1:(P-1),1:(P-1)])
        )
        
        bp_proposed[x,] <- c(bpx_proposed_pm1, - sum(bpx_proposed_pm1)) 
        
        # ap_proposed[x,] <- sampleMVNconstraint_k(apx_proposed, 
        #                                          Sigma_star_all[x,,])
        
      }
      
      # compute ap lik
      for (x in 1:X) {
        
        loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2)
        
        loglik_proposal <- loglik_proposal - loglik_bpx_current(bp_proposed[x,])
        loglik_current <- loglik_current - loglik_bpx_current(rep(0, P))
        
        logproposal <- logproposal + 
          dmvnorm(as.vector(bp_proposed[x,1:(P-1)]), 
                  as.vector(bp_star_all[x,1:(P-1)]), 
                  as.matrix(Sigma_star_all[x,1:(P-1),1:(P-1)]), log = T)
        
      }
      
      logprior_proposal <- log(p_bp)
      logprior_current <- log(1 - p_bp)
      
      mh_ratio <- exp(loglik_proposal - loglik_current + 
                        logprior_proposal - logprior_current - 
                        logproposal)
      
      if(runif(1) < mh_ratio){
        
        gamma_b <- 1
        bp <- bp_proposed
        
      }
       
      
    } 
    
  }
  
  # propose gamma_k
  
  if (update_gammak) {
    
    a2 <- matrix(a, X, P, byrow = F) + ap
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    kp_star_all <- matrix(NA, Y, P)
    Sigma_star_all <- array(NA, dim = c(Y, P, P))
    
    for (t in 1:Y) {
      
      loglik_kpt_current <- loglik_kpt_r(t, a2, b2, k, kp, d, E)
      gr_loglik_kpt_current <- gr_loglik_kpt_r(t, a2, b2, k, kp, d, E)
      
      # find maximum
      laplace_fit <- optim(
        # par = kp[t,],
        par = rep(0, P),
        fn = loglik_kpt_current,
        gr = gr_loglik_kpt_current,
        method = c("BFGS"),
        hessian = T
      )
   
      kp_star_all[t,] <- laplace_fit$par
      Sigma_star_all[t,,] <- solve(laplace_fit$hessian)
      
    }
    
    if(gamma_k == 1){ # currently in nonparametric model
      
      kp_proposed <- kp
      
      loglik_proposal <- 0
      loglik_current <- 0
      logproposal <- 0
      
      # propose new kp
      for (t in 1:Y) {
        
        loglik_kpt_current <- loglik_kpt_r(t, a2, b2, k, kp, d, E)
        
        loglik_proposal <- loglik_proposal - loglik_kpt_current(kp_proposed[t,])
        loglik_current <- loglik_current - loglik_kpt_current(rep(0, P))
        
        logproposal <- logproposal + 
          
          dmvnorm(as.vector(kp_proposed[t,1:(P-1)]), 
                  as.vector(kp_star_all[t,1:(P-1)]), 
                  as.matrix(Sigma_star_all[t,1:(P-1),1:(P-1)]), log = T)
        
      }
      
      logprior_proposal <- log(p_kp)
      logprior_current <- log(1 - p_kp)
      
      mh_ratio <- 1 / exp(loglik_proposal - loglik_current + 
                        logprior_proposal - logprior_current - 
                        logproposal)
      
      if(runif(1) < mh_ratio){
        
        gamma_k <- 0
        kp <- matrix(0, Y, P)
        
      }
      
    } else { # currently in base model
      
      kp_proposed <- matrix(NA, Y, P)
      
      loglik_proposal <- 0
      loglik_current <- 0
      logproposal <- 0
      
      # propose new bp
      for (t in 1:Y) {
      
        kpt_proposed_pm1 <- mvrnorm(
          n = 1,
          kp_star_all[t,1:(P-1)],
          as.matrix(Sigma_star_all[t,1:(P-1),1:(P-1)])
        )
        
        kp_proposed[t,] <- c(kpt_proposed_pm1, - sum(kpt_proposed_pm1)) 
        
      }
      
      # compute kp lik
      for (t in 1:Y) {
        
        loglik_kpt_current <- loglik_kpt_r(t, a2, b2, k, kp, d, E)
        
        loglik_proposal <- loglik_proposal - loglik_kpt_current(kp_proposed[t,])
        loglik_current <- loglik_current - loglik_kpt_current(rep(0, P))
        
        logproposal <- logproposal + 
          dmvnorm(as.vector(kp_proposed[t,1:(P-1)]), 
                  as.vector(kp_star_all[t,1:(P-1)]), 
                  as.matrix(Sigma_star_all[t,1:(P-1),1:(P-1)]), log = T)
        
      }
      
      logprior_proposal <- log(p_kp)
      logprior_current <- log(1 - p_kp)
      
      mh_ratio <- exp(loglik_proposal - loglik_current + 
                        logprior_proposal - logprior_current - 
                        logproposal)
      
      if(runif(1) < mh_ratio){
        
        gamma_k <- 1
        kp <- kp_proposed
        
      }
       
      
    } 
    
  }
  
  # output
  
  if(iter > nburn){
    currentIter <- iter - nburn
    a_output[currentIter,] <- a
    b_output[currentIter,] <- b
    k_output[currentIter,] <- k
    ap_output[currentIter,,] <- ap
    bp_output[currentIter,,] <- bp
    kp_output[currentIter,,] <- kp
    
    # loglik
    {
      a2 <- matrix(a, X, P, byrow = F) + ap
      b2 <- matrix(b, X, P, byrow = F) + bp
      k2 <- matrix(k, Y, P, byrow = F) + kp
      
      loglik_output[currentIter] <-  loglik_LCp(a2, b2, k2, d, E)
    }
    
   
  }
  
}

beep()

# OUTPUT ----

# loglik
{
  a2 <- matrix(a, X, P, byrow = F) + ap
  b2 <- matrix(b, X, P, byrow = F) + bp
  k2 <- matrix(k, Y, P, byrow = F) + kp
  
  a2_true <- matrix(a_true, X, P, byrow = F) + ap_true
  b2_true <- matrix(b_true, X, P, byrow = F) + bp_true
  k2_true <- matrix(k_true, Y, P, byrow = F) + kp_true
  
  loglik_est <- loglik_LCp(a2, b2, k2, d, E)
  loglik_true <- loglik_LCp(a2_true, b2_true, k2_true, d, E)
}

qplot(1:niter, loglik_output[1:niter]) + 
  geom_hline(aes(yintercept = loglik_true))

# a
{
  a_CI <- apply(a_output, 2, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  ggplot(data = NULL, aes(x = 1:X,
                          y = a_true,
                          ymin = a_CI[1,],
                          ymax = a_CI[2,])) + geom_errorbar() + 
    geom_point(color = "red")
}

# ap
{
  ap_CI <- apply(ap_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  ap_CI <- t(apply(ap_CI, 1, as.vector))
  
  ggplot(data = NULL, aes(x = 1:(X * P),
                          y = as.vector(ap_true),
                          ymin = ap_CI[1,],
                          ymax = ap_CI[2,])) + geom_errorbar() + 
    geom_point(color = "red")
}

# b
{
  b_CI <- apply(b_output, 2, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  ggplot(data = NULL, aes(x = 1:X,
                          y = b_true,
                          ymin = b_CI[1,],
                          ymax = b_CI[2,])) + geom_errorbar() + 
    geom_point(color = "red")
}

# bp
{
  bp_CI <- apply(bp_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  bp_CI <- t(apply(bp_CI, 1, as.vector))
  
  ggplot(data = NULL, aes(x = 1:(X * P),
                          y = as.vector(bp_true),
                          ymin = bp_CI[1,],
                          ymax = bp_CI[2,])) + geom_errorbar() + 
    geom_point(color = "red")
}

# bp
{
  bbp_output <- bp_output + array(b_output, dim = c(niter, X, P))
  bbp_CI <- apply(bbp_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  bbp_CI <- t(apply(bp_CI, 1, as.vector))
  
  data_plot <- data.frame(Age = rep(1:X, each = P),
                          Product = factor(rep(1:P, X)),
                          True = matrix(b_true, X, P, byrow = T) + bp_true,
                          bbp_CI[1,],
                          bbp_CI[2,])
  
  colnames(data_plot)[3:5] <- c("True","CI1","CI2")
  
  ggplot(data = data_plot, aes(x = Age,
                          color = Product,
                          group = Product,
                          y = True,
                          ymin = CI1,
                          ymax = CI2)) + #geom_errorbar() + 
    geom_line() + 
    geom_point(color = "red")
}

# k
{
  k_CI <- apply(k_output, 2, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  ggplot(data = NULL, aes(x = 1:Y,
                          y = k_true,
                          ymin = k_CI[1,],
                          ymax = k_CI[2,])) + geom_errorbar() + 
    geom_point(color = "red")
}

# kp
{
  kp_CI <- apply(kp_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  kp_CI <- t(apply(kp_CI, 1, as.vector))
  
  ggplot(data = NULL, aes(x = 1:(Y * P),
                          y = as.vector(kp_true),
                          ymin = kp_CI[1,],
                          ymax = kp_CI[2,])) + geom_errorbar() + 
    geom_point(color = "red")
}

# mortality curves
{
  
  data_plot <- expand.grid(1:X, 1:Y, 1:P)
  
  m_CI <- apply(data_plot, 1, function(x){
    
    idx_age <- x[1]
    idx_year <- x[2]
    idx_prod <- x[3]
    
    m_current_output <- 
      a_output[,idx_age] + ap_output[,idx_age, idx_prod] + 
      (b_output[,idx_age] + bp_output[,idx_age, idx_prod]) *
      (k_output[,idx_year] + kp_output[,idx_year, idx_prod])
      
    quantile(m_current_output, probs = c(0.025, 0.975))
    
  })
  
  m_true <- apply(data_plot, 1, function(x){
    
    idx_age <- x[1]
    idx_year <- x[2]
    idx_prod <- x[3]
    
    a_true[idx_age] + ap_true[idx_age, idx_prod] + 
      (b_true[idx_age] + bp_true[idx_age, idx_prod]) *
      (k_true[idx_year] + kp_true[idx_year, idx_prod])
    
    
  })
  
  data_plot$CI_min <- m_CI[1,]
  data_plot$CI_max <- m_CI[2,]
  data_plot$True <- m_true
  
  labeller_years <- as.character(years)
  names(labeller_years) <- 1:Y
  
  data_plot %>% 
    rename(Age = Var1,
           Year = Var2,
           Product = Var3) %>% 
    mutate(Product = factor(Product)) %>% 
    ggplot(aes(x = Age,
               y = True,
               ymin = CI_min,
               ymax = CI_max,
               color = Product,
               group = Product,
               fill = Product)) + geom_line(size = 1) + 
    geom_ribbon(alpha = .3) + theme_bw() + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 14)) + 
    ylab("Log Mortality") + 
    scale_x_continuous(breaks = 1:X, labels = ages) + 
    facet_wrap(~ Year,
               labeller = labeller(Year = labeller_years))
  
  # single year
  {
    data_plot %>% 
      filter(Var2 == 1) %>% 
      rename(Age = Var1,
             Year = Var2,
             Product = Var3) %>% 
      mutate(Product = factor(Product)) %>% 
      ggplot(aes(x = Age,
                 y = True,
                 ymin = CI_min,
                 ymax = CI_max,
                 color = Product,
                 group = Product,
                 fill = Product)) + geom_line(size = 2) + 
      geom_ribbon(alpha = .5) + theme_bw() + 
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14)) + 
      ylab("Log Mortality") + 
      scale_x_continuous(breaks = 1:X, labels = ages)
    }
}
