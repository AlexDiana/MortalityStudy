
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
library(wiqid)
library(coda)

sourceCpp(here("Code/Product","code.cpp"))  
source(here("Code/Product","functions.R"))

# SIMULATED DATA ----------

realData <- F

X <- 20 # 30 # ages
Y <- 20 # 36 #  years
P <- 2 # products

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
# bx <- seq(.04, 0.0, length.out = X)
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
modelA <- "A"

if (modelA == "NP") { # model NP
  ap_true <- matrix(rnorm(X * P, sd = .1), X, P)
  ap_true <- t(apply(ap_true, 1, function(x){
    x - mean(x)
  }))
} else if (modelA == "A") { # model A
  ap_0 <- rnorm(P, sd = .5)
  ap_0 <- ap_0 - mean(ap_0)
  ap_true <- matrix(ap_0, X, P, byrow = T)
} else if (modelA == "M") { # model M
  ap_0 <- rnorm(P, sd = .5)
  ap_0 <- exp(ap_0) / mean(exp(ap_0))
  ap_true <- matrix(ap_0, X, P, byrow = T)
} else if (modelA == "B") { # model B
  ap_true <- matrix(0, X, P)
}

aap_true <- compute_aap(a_true, ap_true, modelA)

# bp
modelB <- "A"

if (modelB == "NP") { # model NP 
  bp_true <- matrix(rnorm(X * P, sd = .05), X, P)
  bp_true <- t(apply(bp_true, 1, function(x){
    x - mean(x)
  }))
} else if (modelB == "A") { # model A
  bp_0 <- rnorm(P, sd = .05)
  bp_0 <- bp_0 - mean(bp_0)
  bp_true <- matrix(bp_0, X, P, byrow = T)
} else if (modelB == "B") { # model B
  bp_true <- matrix(0, X, P)
}

# kp
modelK <- "A"

if (modelK == "NP") { # model NP 
  kp_true <- matrix(rnorm(Y * P, sd = 1), Y, P)
  kp_true <- t(apply(kp_true, 1, function(x){
    x - mean(x)
  }))
} else if (modelK == "A") { # model A
  kp_0 <- rnorm(P, sd = sd(k_true))
  kp_0 <- kp_0 - mean(kp_0)
  kp_true <- matrix(kp_0, Y, P, byrow = T)
} else if (modelK == "B") { # model B
  kp_true <- matrix(0, Y, P)
}

# simulate data
d <- array(NA, dim = c(X, Y, P))
m <- array(NA, dim = c(X, Y, P))
for (x in 1:X) {
  for (t in 1:Y) {
    for (p in 1:P) {
      m[x,t,p] <- exp(aap_true[x,p] + 
                        (b_true[x] + bp_true[x,p]) * 
                        (k_true[t] + kp_true[t,p]))
      d[x,t,p] <- rpois(1, E[x,t,p] * m[x,t,p])  
    }
  }
}

idx <- sample(1:(X*Y*P), replace = T, size = 100)
d[idx] <- NA
E[idx] <- NA

# REAL DATA PRODUCTS ------

realData <- T

load(here("Data","data_products.rda"))

d <- d[10:60,,]
E <- E[10:60,,]

X <- dim(d)[1]
Y <- dim(d)[2]
P <- dim(d)[3]

ages <- dimnames(d)[[1]]
years <- dimnames(d)[[2]]
products <- dimnames(d)[[3]]

# REAL DATA SEX ------

realData <- T

load(here("Data","data_sex_UK.rda"))

X_subset <- which(as.numeric(dimnames(d)[[1]]) > 10 &
                    as.numeric(dimnames(d)[[1]]) < 85)
Y_subset <- which(as.numeric(dimnames(d)[[2]]) > 1961 &
                    as.numeric(dimnames(d)[[2]]) < 2005)

d <- d[X_subset, Y_subset,]
E <- E[X_subset, Y_subset,]

X <- dim(d)[1]
Y <- dim(d)[2]
P <- dim(d)[3]

ages <- dimnames(d)[[1]]
years <- dimnames(d)[[2]]
groups <- dimnames(d)[[3]]

# REAL DATA COUNTRY -------

countries <- c("France","Italy","Japan","Denmark")

P <- length(countries)

age_subset <- c(60,90)
year_subset <- c(1990,)
d <- 

for (p in 1:P) {
  load(here("Data","Countries",paste0(country,"-data.rda")))
  
  X_subset <- which(as.numeric(dimnames(d)[[1]]) > 60 &
                      as.numeric(dimnames(d)[[1]]) < 90)
  Y_subset <- which(as.numeric(dimnames(d)[[2]]) > 1990 &
                      as.numeric(dimnames(d)[[2]]) < 2022)
  
  d <- d[X_subset, Y_subset]
  E <- E[X_subset, Y_subset]
  
  X <- dim(d)[1]
  Y <- dim(d)[2]
  
  ages <- as.numeric(dimnames(d)[[1]])
  years <- as.numeric(dimnames(d)[[2]])  
  
}


# MCMC -----

nburn <- 200
niter <- 200

realData <- T

# starting values
{
  # true
  if(!realData) {
    
    a <- a_true
    b <- b_true
    k <- k_true

    ap <- ap_true
    bp <- bp_true
    kp <- kp_true
    
    gamma_a <- modelA
    gamma_b <- modelB
    gamma_k <- modelK
    
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
    if(F){
      
      a_start <- apply(d / E, 1, function(x){
        mean(log(x + 1), na.rm = T)
      })
      
      laplace_fit <- optim(
        par = c(
          a_start,
          rep(1 / X, X),
          rep(0, Y)),
        fn = loglik_r(d, E),
        method = c("BFGS")
      )
      
      param_star <- laplace_fit$par
      a <- param_star[1:X]
      b <- param_star[X + 1:X]
      k <- param_star[2 *X + 1:Y]
      
      # - loglik_r(d, E)(param_star) 
      # - loglik_r(d, E)(c(a_true, b_true, k_true))
      
      sumb <- sum(b)
      b <- b / sumb
      k <- k * sumb
      meank <- mean(k)
      k <- k - meank
      a <- a + b * meank
      
      ap <- matrix(0, X, P)
      bp <- matrix(0, X, P)
      kp <- matrix(0, Y, P)
      
      gamma_a <- modelA
      gamma_b <- modelB
      gamma_k <- modelK
    }
    
  }
  
}

# prior
{
  # a_configs <- c("B","A")
  # nConfigs_a <- length(a_configs)
  # p_ap <- c(.6, .4)
  a_configs <- c("B","A","NP")
  nConfigs_a <- length(a_configs)
  p_ap <- c(.6, .3, .1)
  
  # b_configs <- c("B","A")
  # nConfigs_b <- length(b_configs)
  # p_bp <- c(.6, .4)
  b_configs <- c("B","A","NP")
  nConfigs_b <- length(b_configs)
  p_bp <- c(.6, .3, .1)
  
  # k_configs <- c("B","A")
  # nConfigs_k <- length(k_configs)
  # p_kp <- c(.6, .4)
  k_configs <- c("B","A","NP")
  nConfigs_k <- length(k_configs)
  p_kp <- c(.6, .3, .1)
  
  sd_ap <- 1
  sd_bp <- 1
  sd_kp <- 1
  
  sd_b <- 3
  sd_k <- 3
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
  update_gammaa <- F
  update_gammab <- F
  update_gammak <- F
}

for (iter in 1:(niter + nburn)) {
  
  # print(a[1:5])
  # print(b[1:5])
  # print(k[1:5])
  # 
  print(paste0("Gamma A = ",gamma_a,
               ", Gamma B = ",gamma_b,
               ", Gamma K = ",gamma_k))
  
  # print(k)
  
  if(iter %% 10 == 0){
    if(iter < nburn){
      print(paste0("Burn-in Iteration = ",iter))  
    } else {
      print(paste0("Iteration = ",iter - nburn))  
    }  
  }
  
  # update k0
  # print("Updatek0")
  if (update_k0) {
    
    # add product effect to ax and bx
    a2 <- matrix(a, X, P, byrow = F) + ap
    b2 <- matrix(b, X, P, byrow = F) + bp

    # find maximum over all variables together
    {
      # loglik_km1_current <- loglik_k0m1_r(d, E, a2, b2, kp)
      # # gr_loglik_km1_current <- gr_loglik_k0m1_r(d, E, a2, b2, kp)
      # 
      # # find maximum
      # laplace_fit <- optim(
      #   # par = k[1:(Y-1)],
      #   par = rep(0, Y-1),
      #   fn = loglik_km1_current,
      #   # gr = gr_loglik_km1_current,
      #   method = c("CG"),
      #   hessian = T
      # )
      # 
      # k_star <- laplace_fit$par
      # H_star <- laplace_fit$hessian
      # 
      # Sigma_star <- solve(H_star)
      # 
      # kpm1_proposed <- mrt2(k_star, Sigma_star, df = 3)
      # 
      # k_proposed <- c(kpm1_proposed, - sum(kpm1_proposed))
      # 
      # logproposal_star <- dmt_cpp(kpm1_proposed, nu = 3, k_star,
      #                             Sigma_star, T)
      # 
      # logproposal_current <- dmt_cpp(k[1:(Y-1)], nu = 3, k_star,
      #                                Sigma_star, T)
    }
    
    # find maximum for each variabe individually (it's the same as independent)
    {
        k_star <- rep(NA, Y)
        Sigma_star <- matrix(0, Y, Y)

        for (t in 1:Y) {
          # print(t)

          list_kstar <- findProposalK0(t, a2, b2, k, kp, d, E)

          loglik_kt_current <- loglik_kt_r(t, d, E, a2, b2, k, kp)

          k_star[t] <- list_kstar$k_star
          Sigma_star[t, t] <- list_kstar$Sigma_star

        }

        k_proposed <- as.vector(sampleMTconstraint_k(k_star, Sigma_star))
        # k_proposed <- as.vector(sampleMVNconstraint_k(k_star, Sigma_star))
        
        logproposal_star <- sum(
          sapply(1:Y, function(t){
            log(dt2(k_proposed[t], k_star[t], sqrt(diag(Sigma_star)[t]), df = 3))
          }))
        # logproposal_star <- sum(dnorm(k_proposed, k_star,
        #                               sqrt(diag(Sigma_star)), log = T))

        logproposal_current <- sum(
          sapply(1:Y, function(t){
            log(dt2(k[t], k_star[t], sqrt(diag(Sigma_star)[t]), df = 3))
          }))
        # logproposal_current <- sum(dnorm(k, k_star,
                                         # sqrt(diag(Sigma_star)), log = T))

    }
    
    # create loglikelihood function of k given ax and bx
    loglik_k_current <- loglik_k0_r(d, E, a2, b2, kp)

    loglik_star <- - loglik_k_current(k_proposed)
    loglik_current <- - loglik_k_current(k)

    logprior_star <- sum(dnorm(k_proposed,
                               0, sd = sd_k, log = T))
    logprior_current <- sum(dnorm(k,
                                  0, sd = sd_k, log = T))

    mh_ratio <- exp(loglik_star - loglik_current +
                      logprior_star - logprior_current +
                      logproposal_current  - logproposal_star)

    if(runif(1) < mh_ratio){
      k <- as.vector(k_proposed)
    }
    
  }
  
  # update k_p
  # print("Updatekp")
  if (update_kp) {
    
    if(gamma_k != "B"){ # if the product effect is present for this variable
      
      # add product effect to ax and bx
      a2 <- matrix(a, X, P, byrow = F) + ap
      b2 <- matrix(b, X, P, byrow = F) + bp
      
      if(gamma_k == "A"){ # additive effect
        
        loglik_kp_current <- loglik_kp_r(a2, b2, d, E, k)
        gr_loglik_kp_current <- gr_loglik_kp_r(a2, b2, d, E, k)
        
        # find maximum
        laplace_fit <- optim(
          # par = ap[x,],
          par = rep(0, P),
          fn = loglik_kp_current,
          method = c("BFGS"),
          # gr = gr_loglik_kp_current,
          hessian = T
        )
        
        kp_star <- laplace_fit$par
        H_star <- laplace_fit$hessian
        
        Sigma_star <- 2 * solve(H_star)
        
        kp_proposed <- sampleMVNconstraint_k(kp_star, Sigma_star)
        
        loglik_proposal <- - loglik_kp_current(kp_proposed)
        loglik_current <- - loglik_kp_current(kp[t,])
        
        logproposal_star <- dmvnorm(as.vector(kp_proposed), as.vector(kp_star), 
                                    Sigma_star, log = T)
        logproposal_current <- dmvnorm(as.vector(kp[t,]), as.vector(kp_star), 
                                       Sigma_star, log = T)
        
        logprior_star <- sum(dnorm(as.vector(kp_proposed), 0, sd_ap, log = T))
        logprior_current <- sum(dnorm(as.vector(kp[t,]), 0, sd_ap, log = T))
        
        mh_ratio <- exp(loglik_proposal + logproposal_current +
                          logprior_star - logprior_current - 
                          loglik_current - logproposal_star)
        # print(mh_ratio)
        if(runif(1) < mh_ratio){
          kp <- matrix(kp_proposed, Y, P, byrow = T)
        }
      
      } else if(gamma_k == "NP"){ 
        
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
          
          # logproposal_star <- dmvnorm(as.vector(kp_proposed), as.vector(kp_star), Sigma_star, log = T)
          # logproposal_current <- dmvnorm(as.vector(kp[t,]), as.vector(kp_star), Sigma_star, log = T)
          
          logproposal_star <- fastdmvnorm(as.vector(kp_proposed),
                                          as.vector(kp_star),
                                          .5 * H_star)
          logproposal_current <- fastdmvnorm(as.vector(kp[t,]),
                                             as.vector(kp_star),
                                             .5 * H_star)
          
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
    
  }
  
  # update a
  # print("Updatea")
  if(update_ab){
    
    # add product effect to k
    k2 <- matrix(k, Y, P, byrow = F) + kp
    b2 <- matrix(b, X, P, byrow = F) + bp
  
    a2_current <- aperm(array(ap, dim = c(X, P, Y)), perm = c(1,3,2))
    b2_current <- aperm(array(b2, dim = c(X, P, Y)), perm = c(1,3,2))
    k2_current <- aperm(array(k2, dim = c(Y, P, X)), perm = c(3,1,2))
    
    m_xtp <- a2_current + b2_current * k2_current + log(E)
    
    for (x in 1:X) {
      
      # cxt <- exp(matrix(ap[x,], Y, P, byrow = T) + 
      #              k2 * matrix(b2[x,], Y, P, byrow = T))
      
      ax_star <- log(
        sum(d[x,,], na.rm = T) / sum(exp(m_xtp[x,,]), na.rm = T)
      )
      
      hessian_ax_star <- derl2der2a(ax_star, m_xtp[x,,])
      sd_star <- sqrt(- 1 / hessian_ax_star)
      
      ax_proposal <- rt2(1, ax_star, sd_star, df = 3)
      
      loglik_proposal <- loglik_LCp_x(x, ax_proposal, ap, b2, k2, d, E)
      
      loglik_current <- loglik_LCp_x(x, a[x], ap, b2, k2, d, E)
      
      logproposal_proposal <- log(dt2(ax_proposal, ax_star, sd_star, df = 3))
      logproposal_current <- log(dt2(a[x], ax_star, sd_star, df = 3))
      
      mh_ratio <- exp(loglik_proposal - loglik_current + 
                        logproposal_current - logproposal_proposal)
      
      if(runif(1) < mh_ratio){
        
        a[x] <- ax_proposal
        
      }
      
    }
    
    
  }
  
  # update b
  # print("Updateb")
  if(update_ab){
    
    # add product effect to k
    k2 <- matrix(k, Y, P, byrow = F) + kp
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    loglik_b_current <- loglik_bm1_r(d, E, a2, bp, k2)
    gr_loglik_b_current <- gr_loglik_bm1_r(d, E, a2, bp, k2)

    # find maximum
    laplace_fit <- optim(
      # par = rep(0, X-1),
      par = b[1:(X-1)],
      fn = loglik_b_current,
      gr = gr_loglik_b_current,
      method = c("BFGS"),
      hessian = T
    )

    loglik_b_current(laplace_fit$par) < loglik_b_current(b[1:(X-1)])
    
    bm1_star <- laplace_fit$par
    Sigma_star <- solve(laplace_fit$hessian)

    bm1_proposed <- mrt2(bm1_star, Sigma_star, df = 3)
    b_proposed <- c(bm1_proposed, 1 - sum(bm1_proposed))

    loglik_b_current <- loglik_b_r(d, E, a2, bp, k2)
    
    loglik_star <- - loglik_b_current(b_proposed)
    loglik_current <- - loglik_b_current(b)
    # loglik_star
    
    logproposal_star <- dmt_cpp(bm1_proposed, nu = 3, bm1_star, Sigma_star, T)
    logproposal_current <- dmt_cpp(b[1:(X-1)], nu = 3, bm1_star, Sigma_star, T)
    
    logprior_star <- sum(dnorm(b_proposed, mean = 0, sd = sd_b, log = T))
    logprior_current <- sum(dnorm(b, mean = 0, sd = sd_b, log = T))
    # loglik_star - loglik_current
    mh_ratio <- exp(
      loglik_star - loglik_current + 
        logprior_star - logprior_current + 
        logproposal_current - logproposal_star
    )
    
    if(runif(1) < mh_ratio){
      
      b <- b_proposed
      
    }
    
  }
  
  # update a and b jointly
  if (update_ab) {
    
    # add product effect to k
    k2 <- matrix(k, Y, P, byrow = F) + kp

    # find maximum over all variables together
    {
      
      # gr_loglik_ab_current <- gr_loglik_ab_r(d, E, ap, bp, k2)
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
    
    # Sigma_star <- 2 * Sigma_star
  
    ab_proposed <- sampleMVNconstraint_ab(ab_star, Sigma_star)
    a_proposed <- ab_proposed[1:X]
    b_proposed <- ab_proposed[X + 1:X]

    loglik_ab_current <- loglik_ab_r(d, E, ap, bp, k2)
    
    loglik_star <- - loglik_ab_current(ab_proposed)
    loglik_current <- - loglik_ab_current(c(a, b))
    
    # compute proposal for all variables together
    {
      # logproposal_star <- dmvnorm(as.vector(ab_proposed), as.vector(ab_star), Sigma_star, log = T)
      # logproposal_current <- dmvnorm(c(a,b), as.vector(ab_star), Sigma_star, log = T)
    }

    # find proposal for each variabe individually (it's the same as independent)
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
    
    logprior_star <- sum(dnorm(ab_proposed[X + 1:X],
                               0, sd = sd_b, log = T))
    logprior_current <- sum(dnorm(b,
                               0, sd = sd_b, log = T))
    
    mh_ratio <- exp(loglik_star - loglik_current + 
                      logprior_star - logprior_current +
                      logproposal_current -  logproposal_star)
    
    if(runif(1) < mh_ratio){
      a <- a_proposed
      b <- b_proposed
    }
    
  }
  
  # update a_p
  # print("Updateap")
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
        
      } else if(gamma_a == "M") { # multiplicative model
        
        # loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
        # gr_loglik_ap_current <- gr_loglik_ap_r(a, b2, d, E, k2)
        # 
        # # find maximum
        # laplace_fit <- optim(
        #   # par = ap[x,],
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
        # 
        # ap_proposed <- sampleMVNconstraint_k(ap_star, Sigma_star)
        # 
        # loglik_proposal <- - loglik_ap_current(ap_proposed)
        # loglik_current <- - loglik_ap_current(ap[x,])
        # 
        # logproposal_star <- dmvnorm(as.vector(ap_proposed), as.vector(ap_star), 
        #                             Sigma_star, log = T)
        # logproposal_current <- dmvnorm(as.vector(ap[x,]), as.vector(ap_star), 
        #                                Sigma_star, log = T)
        # 
        # logprior_star <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
        # logprior_current <- sum(dnorm(as.vector(ap[x,]), 0, sd_ap, log = T))
        # 
        # mh_ratio <- exp(loglik_proposal + logproposal_current +
        #                   logprior_star - logprior_current - 
        #                   loglik_current - logproposal_star)
        # # print(mh_ratio)
        # if(runif(1) < mh_ratio){
        #   ap <- matrix(ap_proposed, X, P, byrow = T)
        # }
        
      } else if(gamma_a == "NP"){ # nonparametric model
        
        for (x in 1:X) {
          
          loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k, sd_ap)
          gr_loglik_apx_current <- gr_loglik_apx_r(x, a, ap, b2, d, E, k2, sd_ap)
          
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
  # print("Updatebp")
  if (update_bp) {
    
    if(gamma_b != "B"){
      
      k2 <- matrix(k, Y, P, byrow = F) + kp
      a2 <- matrix(a, X, P, byrow = F) + ap
      
      if(gamma_b == "A"){ # currently in additive model
        
        # old
        {
          # loglik_bp_current <- loglik_bp_r(a2, b, d, E, k2)
          # gr_loglik_bp_current <- gr_loglik_bp_r(a2, b, d, E, k2)
          # 
          # # find maximum
          # laplace_fit <- optim(
          #   # par = ap[x,],
          #   par = rep(0, P),
          #   fn = loglik_bp_current,
          #   method = c("BFGS"),
          #   gr = gr_loglik_bp_current,
          #   hessian = T
          # )
          # 
          # bp_star <- laplace_fit$par
          # H_star <- laplace_fit$hessian
          # 
          # Sigma_star <- 2 * solve(H_star)
          # # Sigma_star <- 2 * solve(H_star)
          # 
          # bp_proposed <- sampleMVNconstraint_k(bp_star, Sigma_star)
          # 
          # loglik_proposal <- - loglik_bp_current(bp_proposed)
          # loglik_current <- - loglik_bp_current(bp[x,])
          # 
          # logproposal_star <- dmvnorm(as.vector(bp_proposed), 
          #                             as.vector(bp_star), 
          #                             Sigma_star, log = T)
          # 
          # logproposal_current <- dmvnorm(as.vector(bp[x,]), 
          #                                as.vector(bp_star), 
          #                                Sigma_star, log = T)  
        }
        
        # new
        {
          loglik_bpm1_current <- loglik_bp_m1_r(a2, b, d, E, k2)
          # gr_loglik_bp_current <- gr_loglik_bp_r(a2, b, d, E, k2)
          
          # find maximum
          laplace_fit <- optim(
            # par = ap[x,],
            par = rep(0, P - 1),
            fn = loglik_bpm1_current,
            method = c("BFGS"),
            # gr = gr_loglik_bp_current,
            hessian = T
          )
          
          bp_star <- laplace_fit$par
          Sigma_star <- solve(laplace_fit$hessian)
          
          bpm1_proposed <- mrt2(bp_star, Sigma_star, df = 3)
          bp_proposed <- c(bpm1_proposed, - sum(bpm1_proposed))
          
          loglik_bp_current <- loglik_bp_r(a2, b, d, E, k2)
          
          loglik_proposal <- - loglik_bp_current(bp_proposed)
          loglik_current <- - loglik_bp_current(bp[x,])
          
          logproposal_star <- dmt_cpp(as.vector(bpm1_proposed), nu = 3,
                                      as.vector(bp_star), 
                                      Sigma_star, T)
          
          logproposal_current <- dmt_cpp(as.vector(bp[x,1:(P-1)]), nu = 3,
                                         as.vector(bp_star), 
                                         Sigma_star, T)  
        }
        
        logprior_star <- sum(dnorm(as.vector(bp_proposed), 0, sd_bp, log = T))
        logprior_current <- sum(dnorm(as.vector(bp[x,]), 0, sd_bp, log = T))
        
        mh_ratio <- exp(loglik_proposal + logproposal_current +
                          logprior_star - logprior_current - 
                          loglik_current - logproposal_star)
        # print(mh_ratio)
        if(runif(1) < mh_ratio){
          bp <- matrix(bp_proposed, X, P, byrow = T)
        }
        
      } else if(gamma_b == "NP") { # multiplicative model
        
        for (x in 1:X) {
          
          loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
          gr_loglik_bpx_current <- gr_loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
          
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
          
          if(runif(1) < mh_ratio){
            bp[x,] <- bp_proposed
          }
          
        }  
        
      }
      
    }
    
  }
  
  # update a and ap
  if(update_ab & F){
    
    k2 <- matrix(k, Y, P, byrow = F) + kp
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    if(gamma_a == "A"){ # currently in additive model
      
      loglik_aap_current <- loglik_aapm1_r(b2, d, E, k2)
      # gr_loglik_aap_current <- gr_loglik_aap_r(b2, d, E, k2)
      
      aap_current <- c(a, ap[1,1:(P-1)])
      
      # find maximum
      laplace_fit <- optim(
        par = aap_current,
        # par = rep(0, P),
        fn = loglik_aap_current,
        # gr = gr_loglik_aap_current,
        method = c("BFGS"),
        hessian = T
      )
      
      aap_star <- laplace_fit$par
      Sigma_star <- solve(laplace_fit$hessian)
      
      aap_proposed <- mvrnorm(1, aap_star, Sigma_star)
      
      loglik_proposal <- - loglik_aap_current(aap_proposed)
      loglik_current <- - loglik_aap_current(aap_current)
      
      logproposal_star <- dmvnorm(as.vector(aap_proposed), as.vector(aap_star), 
                                  Sigma_star, log = T)
      logproposal_current <- dmvnorm(as.vector(aap_current), as.vector(aap_star), 
                                     Sigma_star, log = T)
      
      ap_proposed <- c(aap_proposed[X + 1:(P-1)], - sum(aap_proposed[X + 1:(P-1)]))
      
      logprior_star <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
      logprior_current <- sum(dnorm(as.vector(ap[1,1:P]), 0, sd_ap, log = T))
      
      mh_ratio <- exp(loglik_proposal - loglik_current + 
                        logproposal_current - logproposal_star +
                        logprior_star - logprior_current)
      
      # print(mh_ratio)
      if(runif(1) < mh_ratio){
        a <- aap_proposed[1:X]
        ap <- matrix(ap_proposed, X, P, byrow = T)
      }
      
    } 
    
  }
  
  # update b and bp
  if(update_ab & F){
    
    k2 <- matrix(k, Y, P, byrow = F) + kp
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    if(gamma_b == "A"){ # currently in additive model
      
      loglik_bbp_current <- loglik_bbp_m1_r(a2, d, E, k2)
      # gr_loglik_aap_current <- gr_loglik_aap_r(b2, d, E, k2)
      
      bbp_current <- c(b, bp[1,1:(P-1)])
      
      # find maximum
      laplace_fit <- optim(
        par = bbp_current,
        # par = rep(0, P),
        fn = loglik_bbp_current,
        # gr = gr_loglik_aap_current,
        method = c("BFGS"),
        hessian = T
      )
      
      bbp_star <- laplace_fit$par
      Sigma_star <- solve(laplace_fit$hessian)
      
      bbp_proposed <- mvrnorm(1, bbp_star, Sigma_star)
      
      loglik_proposal <- - loglik_bbp_current(bbp_proposed)
      loglik_current <- - loglik_bbp_current(bbp_current)
      
      logproposal_star <- dmvnorm(as.vector(bbp_proposed), as.vector(bbp_star), 
                                  Sigma_star, log = T)
      logproposal_current <- dmvnorm(as.vector(bbp_current), as.vector(bbp_star), 
                                     Sigma_star, log = T)
      
      bp_proposed <- c(bbp_proposed[X + 1:(P-1)], - sum(bbp_proposed[X + 1:(P-1)]))
      
      logprior_star <- sum(dnorm(as.vector(bp_proposed), 0, sd_ap, log = T))
      logprior_current <- sum(dnorm(as.vector(bp[1,1:P]), 0, sd_ap, log = T))
      
      mh_ratio <- exp(loglik_proposal - loglik_current + 
                        logproposal_current - logproposal_star +
                        logprior_star - logprior_current)
      
      # print(mh_ratio)
      if(runif(1) < mh_ratio){
        b <- bbp_proposed[1:X]
        bp <- matrix(bp_proposed, X, P, byrow = T)
      }
      
    } 
    
  }
  
  # update k and kp
  if(update_k0 & F){
    
    b2 <- matrix(b, X, P, byrow = F) + bp
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    if(gamma_k == "A"){ # currently in additive model
      
      loglik_kkp_current <- loglik_kppm1_r(a2, b2, d, E)
      # gr_loglik_aap_current <- gr_loglik_aap_r(b2, d, E, k2)
      
      kkp_current <- c(k[1:(Y-1)], kp[1,1:(P-1)])
      
      # find maximum
      laplace_fit <- optim(
        par = kkp_current,
        # par = rep(0, P),
        fn = loglik_kkp_current,
        # gr = gr_loglik_aap_current,
        method = c("BFGS"),
        hessian = T
      )
      
      kkp_star <- laplace_fit$par
      Sigma_star <- solve(laplace_fit$hessian)
      
      kkp_proposed <- mvrnorm(1, kkp_star, Sigma_star)
      
      loglik_proposal <- - loglik_kkp_current(kkp_proposed)
      loglik_current <- - loglik_kkp_current(kkp_current)
      
      logproposal_star <- dmvnorm(as.vector(kkp_proposed), as.vector(kkp_star), 
                                  Sigma_star, log = T)
      logproposal_current <- dmvnorm(as.vector(kkp_current), as.vector(kkp_star), 
                                     Sigma_star, log = T)
      
      kp_proposed <- c(kkp_proposed[Y - 1 + 1:(P-1)], - sum(kkp_proposed[Y - 1 + 1:(P-1)]))
      
      logprior_star <- sum(dnorm(as.vector(kp_proposed), 0, sd_ap, log = T))
      logprior_current <- sum(dnorm(as.vector(kp[1,1:P]), 0, sd_ap, log = T))
      
      mh_ratio <- exp(loglik_proposal - loglik_current + 
                        logproposal_current - logproposal_star +
                        logprior_star - logprior_current)
      
      # print(mh_ratio)
      if(runif(1) < mh_ratio){
        k <- c(kkp_proposed[1:(Y-1)], -sum(kkp_proposed[1:(Y-1)]))
        kp <- matrix(kp_proposed, Y, P, byrow = T)
      }
      
    } 
    
  }
  
  # propose gamma_a
  # print("Updategamma")
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
      
      # old
      {
        # list_ap_proposal <- findProposalApAdditive(a, b2, d, E, k2)
        # ap_star <- list_ap_proposal$ap_star
        # Sigma_star <- list_ap_proposal$Sigma_star
        # 
        # ap_proposed_pm1 <- ap[x,1:(P-1)]
        # ap_proposed <- c(ap_proposed_pm1, - sum(ap_proposed_pm1))
        # 
        # logProposal_current <- dmvnorm(as.vector(ap_proposed_pm1), 
        #                                as.vector(ap_star[1:(P-1)]), 
        #                                as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)
      }
      
      # new
      {
        list_ap_proposal <- findProposalApm1Additive(a, b2, d, E, k2)
        ap_star <- list_ap_proposal$ap_star
        Sigma_star <- list_ap_proposal$Sigma_star
        
        ap_proposed_pm1 <- ap[x,1:(P-1)]
        ap_proposed <- c(ap_proposed_pm1, - sum(ap_proposed_pm1))
        
        logProposal_current <- dmvnorm(as.vector(ap_proposed_pm1), 
                                       as.vector(ap_star), 
                                       as.matrix(Sigma_star), log = T)
      }
      
      
      logprior_param_current <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
      
      
    } else if(gamma_a == "NP"){ # in nonparametric model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
      for (x in 1:X) {
        
        # old
        {
          
          # list_ap_proposal <- findProposalApXNonparametrics(x, a, ap, b2, d, E, k2, sd_ap)
          # ap_star <- list_ap_proposal$ap_star
          # Sigma_star <- list_ap_proposal$Sigma_star
          # 
          # logProposal_current <- logProposal_current + 
          #   dmvnorm(as.vector(ap[x,1:(P-1)]), 
          #           as.vector(ap_star[1:(P-1)]), 
          #           as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4, log = T)
        }
       
         # new
        {
          
          list_ap_proposal <- findProposalApXm1Nonparametrics(x, a, ap, b2, d, E, k2, sd_ap)
          ap_star <- list_ap_proposal$ap_star
          Sigma_star <- list_ap_proposal$Sigma_star
          
          logProposal_current <- logProposal_current + 
            dmvnorm(as.vector(ap[x,1:(P-1)]), 
                    as.vector(ap_star), 
                    as.matrix(Sigma_star), log = T)
        }
      
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
      
      # old
      {
        # list_ap_proposal <- findProposalApAdditive(a, b2, d, E, k2)
        # ap_star <- list_ap_proposal$ap_star
        # Sigma_star <- list_ap_proposal$Sigma_star
        # 
        # ap_proposed_pm1 <- mvrnorm(
        #   n = 1,
        #   ap_star[1:(P-1)],
        #   as.matrix(Sigma_star[1:(P-1),1:(P-1)])
        # )
        # 
        # ap_proposed <- c(ap_proposed_pm1, - sum(ap_proposed_pm1))
        # 
        # logProposal_star <- dmvnorm(as.vector(ap_proposed_pm1), 
        #                             as.vector(ap_star[1:(P-1)]), 
        #                             as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)
        # 
      }
      
      # new
      {
        list_ap_proposal <- findProposalApm1Additive(a, b2, d, E, k2)
        ap_star <- list_ap_proposal$ap_star
        Sigma_star <- list_ap_proposal$Sigma_star
        
        ap_proposed_pm1 <- mvrnorm(
          n = 1,
          ap_star,
          as.matrix(Sigma_star)
        )
        
        ap_proposed <- c(ap_proposed_pm1, - sum(ap_proposed_pm1))
        
        logProposal_star <- dmvnorm(as.vector(ap_proposed_pm1), 
                                    as.vector(ap_star), 
                                    as.matrix(Sigma_star), log = T)
        
      }
      
      logprior_param_star <- sum(dnorm(as.vector(ap_proposed), 0, sd_ap, log = T))
      
      ap_proposed <- matrix(ap_proposed, X, P, byrow = T)
      
    } else if (newConfig == "NP"){ # proposing nonparametric model
      
      # loglik_star <- 0
      logProposal_star <- 0
      logprior_param_star <- 0
      
      ap_proposed <- matrix(NA, X, P)
      
      for (x in 1:X) {
        
        # old
        {
          # list_ap_proposal <- findProposalApXNonparametrics(x, a, ap, b2, d, E, k2, sd_ap)
          # ap_star <- list_ap_proposal$ap_star
          # Sigma_star <- list_ap_proposal$Sigma_star
          # 
          # apx_proposed_pm1 <- mvrnorm(
          #   n = 1,
          #   ap_star[1:(P-1)],
          #   as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4
          # )
          # 
          # ap_proposed[x,] <- c(apx_proposed_pm1, - sum(apx_proposed_pm1)) 
          # 
          # logProposal_star <- logProposal_star + 
          #   dmvnorm(as.vector(apx_proposed_pm1), 
          #           as.vector(ap_star[1:(P-1)]), 
          #           as.matrix(Sigma_star[1:(P-1),1:(P-1)]), 
          #           log = T)  
        }
        
        # new
        {
          list_ap_proposal <- findProposalApXm1Nonparametrics(x, a, ap, b2, d, E, k2, sd_ap)
          ap_star <- list_ap_proposal$ap_star
          Sigma_star <- list_ap_proposal$Sigma_star
          
          #
          
          # aap <- matrix(a, X, P, byrow = F) + ap
          # loglik_current <- loglik_LCp_r(aap, b2, k2, d, E)
          # ap_proposed <- ap
          # ap_proposed[x,] <- c(ap_star, -sum(ap_star))
          # aap_proposed <- matrix(a, X, P, byrow = F) + ap_proposed
          # loglik_proposal <- loglik_LCp_r(aap_proposed, b2, k2, d, E) 
          # loglik_proposal - loglik_current
          #
          
          apx_proposed_pm1 <- 
            mvrnorm(
            n = 1,
            ap_star,
            as.matrix(Sigma_star)
          )
          
          ap_proposed[x,] <- c(apx_proposed_pm1, - sum(apx_proposed_pm1)) 
          
          logProposal_star <- logProposal_star + 
            dmvnorm(as.vector(apx_proposed_pm1), 
                    as.vector(ap_star), 
                    as.matrix(Sigma_star), 
                    log = T)  
        }
       
        
        logprior_param_star <- 
          logprior_param_star + 
          sum(dnorm(as.vector(ap_proposed[x,]), 0, sd_ap, log = T))
        
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
  # print("Updategammab")
  if (update_gammab) {
    
    k2 <- matrix(k, Y, P, byrow = F) + kp
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    logPriorCurrent <- log(p_bp[which(b_configs == gamma_b)])
    
    bbp <- matrix(b, X, P, byrow = F) + bp
    loglik_current <- loglik_LCp_r(a2, bbp, k2, d, E)
    
    if(gamma_b == "B") { # currently in base model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
    } else if(gamma_b == "A"){ # currently in additive model
    
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
      
      # old
      {
        # list_bp_proposal <- findProposalBpAdditive(a2, b, d, E, k2)
        # bp_star <- list_bp_proposal$bp_star
        # Sigma_star <- list_bp_proposal$Sigma_star
        # 
        # bp_proposed_pm1 <- bp[x,1:(P-1)]
        # bp_proposed <- c(bp_proposed_pm1, - sum(bp_proposed_pm1))
        # 
        # logProposal_current <- dmvnorm(as.vector(bp_proposed_pm1), 
        #                                as.vector(bp_star[1:(P-1)]), 
        #                                as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)   
      }
      
      # new
      {
        list_bp_proposal <- findProposalBpm1Additive(a2, b, d, E, k2)
        bp_star <- list_bp_proposal$bp_star
        Sigma_star <- list_bp_proposal$Sigma_star
        
        bp_proposed_pm1 <- bp[x,1:(P-1)]
        bp_proposed <- c(bp_proposed_pm1, - sum(bp_proposed_pm1))
        
        logProposal_current <- dmvnorm(as.vector(bp_proposed_pm1), 
                                       as.vector(bp_star), 
                                       as.matrix(Sigma_star), log = T)
      }
     
      
      logprior_param_current <- sum(dnorm(as.vector(bp_proposed), 0, sd_bp, log = T))
      
      
    } else if(gamma_b == "NP"){ # in nonparametric model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
      for (x in 1:X) {
      
        # old
        {
          # list_bp_proposal <- findProposalBpXNonparametrics(x, a2, b, bp, d, 
          #                                                   E, k2, sd_bp)
          # bp_star <- list_bp_proposal$bp_star
          # Sigma_star <- list_bp_proposal$Sigma_star
          # 
          # logProposal_current <- logProposal_current + 
          #   dmvnorm(as.vector(bp[x,1:(P-1)]), 
          #           as.vector(bp_star[1:(P-1)]), 
          #           as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4, log = T)
        }
        
        # new
        {
          list_bp_proposal <- findProposalBpXM1Nonparametrics(x, a2, b, bp, d, 
                                                            E, k2, sd_bp)
          bp_star <- list_bp_proposal$bp_star
          Sigma_star <- list_bp_proposal$Sigma_star
          
          logProposal_current <- logProposal_current + 
            dmvnorm(as.vector(bp[x,1:(P-1)]), 
                    as.vector(bp_star), 
                    as.matrix(Sigma_star) , log = T)
        }
        
        
        logprior_param_current <- logprior_param_current + 
          sum(dnorm(as.vector(bp[x,]), 0, sd_bp, log = T))
        
      }
      
    }

    newConfig <- sample(setdiff(b_configs, gamma_b), 1)
    
    logPriorStar <- log(p_bp[which(b_configs == newConfig)])
    
    if(newConfig == "B"){ # proposing base model
      
      bp_proposed <- matrix(0, X, P)
      
      logProposal_star <- 0
      logprior_param_star <- 0
      
    } else if (newConfig == "A"){ # proposing additive model
      
      # old
      {
        # list_bp_proposal <- findProposalBpAdditive(a2, b, d, E, k2)
        # bp_star <- list_bp_proposal$bp_star
        # Sigma_star <- list_bp_proposal$Sigma_star
        # 
        # bp_proposed_pm1 <- mvrnorm(
        #   n = 1,
        #   bp_star[1:(P-1)],
        #   as.matrix(Sigma_star[1:(P-1),1:(P-1)])
        # )
        # 
        # bp_proposed <- c(bp_proposed_pm1, - sum(bp_proposed_pm1))
        # 
        # logProposal_star <- dmvnorm(as.vector(bp_proposed_pm1), 
        #                             as.vector(bp_star[1:(P-1)]), 
        #                             as.matrix(Sigma_star[1:(P-1),1:(P-1)]), log = T)
      }
      
      # new
      {
        list_bp_proposal <- findProposalBpm1Additive(a2, b, d, E, k2)
        bp_star <- list_bp_proposal$bp_star
        Sigma_star <- list_bp_proposal$Sigma_star
        
        bp_proposed_pm1 <- mvrnorm(
          n = 1,
          bp_star,
          as.matrix(Sigma_star)
        )
        
        bp_proposed <- c(bp_proposed_pm1, - sum(bp_proposed_pm1))
        
        logProposal_star <- dmvnorm(as.vector(bp_proposed_pm1), 
                                    as.vector(bp_star), 
                                    as.matrix(Sigma_star), log = T)
      }
      
      
      logprior_param_star <- sum(dnorm(as.vector(bp_proposed), 0, sd_bp, log = T))
      
      bp_proposed <- matrix(bp_proposed, X, P, byrow = T)
      
    } else if (newConfig == "NP"){ # proposing nonparametric model
      
      # loglik_star <- 0
      logProposal_star <- 0
      logprior_param_star <- 0
      
      bp_proposed <- matrix(NA, X, P)
      
      for (x in 1:X) {
        
        # old
        {
          # list_bp_proposal <- findProposalBpX(x, a2, b, bp, d, E, k2, sd_bp)
          # bp_star <- list_bp_proposal$bp_star
          # Sigma_star <- list_bp_proposal$Sigma_star
          # 
          # bpx_proposed_pm1 <- mvrnorm(
          #   n = 1,
          #   bp_star[1:(P-1)],
          #   as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4
          # )
          # 
          # bp_proposed[x,] <- c(bpx_proposed_pm1, - sum(bpx_proposed_pm1)) 
          # 
          # # loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
          # 
          # logProposal_star <- logProposal_star + 
          #   dmvnorm(as.vector(bpx_proposed_pm1), 
          #           as.vector(bp_star[1:(P-1)]), 
          #           as.matrix(Sigma_star[1:(P-1),1:(P-1)]), 
          #           log = T)
        }
        
        # new
        {
          list_bp_proposal <- findProposalBpXM1Nonparametrics(x, a2, b, bp, d, E, k2, sd_bp)
          bp_star <- list_bp_proposal$bp_star
          Sigma_star <- list_bp_proposal$Sigma_star
          
          bpx_proposed_pm1 <- mvrnorm(
            n = 1,
            bp_star,
            as.matrix(Sigma_star) 
          )
          
          bp_proposed[x,] <- c(bpx_proposed_pm1, - sum(bpx_proposed_pm1)) 
          
          # loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
          
          logProposal_star <- logProposal_star + 
            dmvnorm(as.vector(bpx_proposed_pm1), 
                    as.vector(bp_star), 
                    as.matrix(Sigma_star), 
                    log = T)
        }
        
        
        logprior_param_star <- 
          logprior_param_star + 
          sum(dnorm(as.vector(bp_proposed[x,]), 0, sd_bp, log = T))
      
      }
      
    }
    
    bbp_proposed <- matrix(b, X, P, byrow = F) + bp_proposed
    loglik_proposal <- loglik_LCp_r(a2, bbp_proposed, k2, d, E)
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current +  # likelihood 
        logPriorStar - logPriorCurrent + # prior of the model
        logprior_param_star - logprior_param_current + # prior of the parameters
        logProposal_current - logProposal_star # proposal of the parameters
        )
    
    if(runif(1) < mh_ratio){
      gamma_b <- newConfig
      bp <- bp_proposed
    }
    
    #### OLD CODE
    
    {
      # k2 <- matrix(k, Y, P, byrow = F) + kp
      # a2 <- matrix(a, X, P, byrow = F) + ap
      # 
      # bp_star_all <- matrix(NA, X, P)
      # Sigma_star_all <- array(NA, dim = c(X, P, P))
      # 
      # for (x in 1:X) {
      #   
      #   loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
      #   gr_loglik_bpx_current <- gr_loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
      #   
      #   # find maximum
      #   laplace_fit <- optim(
      #     # par = bp[x,],
      #     par = rep(0, P),
      #     fn = loglik_bpx_current,
      #     method = c("BFGS"),
      #     gr = gr_loglik_bpx_current,
      #     hessian = T
      #   )
      #   
      #   # laplace_fit <- optimx(
      #   #    # par = bp[x,],
      #   #    par = rep(0, P),
      #   #    fn = loglik_bpx_current,
      #   #    method = c("BFGS"),
      #   #    gr = gr_loglik_bpx_current,
      #   #    hessian = T
      #   #  )
      #   
      #   bp_star_all[x,] <- laplace_fit$par
      #   Sigma_star_all[x,,] <- solve(laplace_fit$hessian)
      #   
      #   # - loglik_bpx_current(bp_star_all[x,])
      #   # - loglik_bpx_current(bp_true[x,])
      #   # - loglik_bpx_current(bp_proposed[x,])
      #   
      # }
      # 
      # if(gamma_b == "NP"){ # currently in nonparametric model
      #   
      #   bp_proposed <- bp
      #   
      #   loglik_proposal <- 0
      #   loglik_current <- 0
      #   logproposal <- 0
      #   
      #   # propose new bp
      #   for (x in 1:X) {
      #     
      #     loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
      #     
      #     loglik_proposal <- loglik_proposal - loglik_bpx_current(bp_proposed[x,])
      #     loglik_current <- loglik_current - loglik_bpx_current(rep(0, P))
      #     
      #     logproposal <- logproposal + 
      #       dmvnorm(as.vector(bp_proposed[x,]), 
      #               as.vector(bp_star_all[x,]), 
      #               Sigma_star_all[x,,], log = T)
      #     
      #   }
      #   
      #   logprior_proposal <- log(p_bp)
      #   logprior_current <- log(1 - p_bp)
      #   
      #   mh_ratio <- 1 / exp(loglik_proposal - loglik_current + 
      #                         logprior_proposal - logprior_current - 
      #                         logproposal)
      #   
      #   if(runif(1) < mh_ratio){
      #     
      #     gamma_b <- "B"
      #     bp <- matrix(0, X, P)
      #     
      #   }
      #   
      # } else { # currently in base model
      #   
      #   bp_proposed <- matrix(NA, X, P)
      #   
      #   loglik_proposal <- 0
      #   loglik_current <- 0
      #   logproposal <- 0
      #   
      #   # propose new bp
      #   for (x in 1:X) {
      #     
      #     bpx_proposed_pm1 <- mvrnorm(
      #       n = 1,
      #       bp_star_all[x,1:(P-1)],
      #       as.matrix(Sigma_star_all[x,1:(P-1),1:(P-1)])
      #     )
      #     
      #     bp_proposed[x,] <- c(bpx_proposed_pm1, - sum(bpx_proposed_pm1)) 
      #     
      #     # ap_proposed[x,] <- sampleMVNconstraint_k(apx_proposed, 
      #     #                                          Sigma_star_all[x,,])
      #     
      #   }
      #   
      #   # compute ap lik
      #   for (x in 1:X) {
      #     
      #     loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
      #     
      #     loglik_proposal <- loglik_proposal - loglik_bpx_current(bp_proposed[x,])
      #     loglik_current <- loglik_current - loglik_bpx_current(rep(0, P))
      #     
      #     logproposal <- logproposal + 
      #       dmvnorm(as.vector(bp_proposed[x,1:(P-1)]), 
      #               as.vector(bp_star_all[x,1:(P-1)]), 
      #               as.matrix(Sigma_star_all[x,1:(P-1),1:(P-1)]), log = T)
      #     
      #   }
      #   
      #   logprior_proposal <- log(p_bp)
      #   logprior_current <- log(1 - p_bp)
      #   
      #   mh_ratio <- exp(loglik_proposal - loglik_current + 
      #                     logprior_proposal - logprior_current - 
      #                     logproposal)
      #   
      #   if(runif(1) < mh_ratio){
      #     
      #     gamma_b <- "NP"
      #     bp <- bp_proposed
      #     
      #   }
      #   
      #   
      # } 
    }
    
  }
  
  # propose gamma_k
  # print("Updategammak")
  if (update_gammak) {
    
    a2 <- matrix(a, X, P, byrow = F) + ap
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    logPriorCurrent <- log(p_kp[which(k_configs == gamma_k)])
    
    kkp <- matrix(k, Y, P, byrow = F) + kp
    loglik_current <- loglik_LCp_r(a2, b2, kkp, d, E)
    
    if(gamma_k == "B") { # currently in base model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
    } else if(gamma_k == "A"){ # currently in additive model
    
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
      
      # list_kp_proposal <- findProposalKpAdditive(a2, b2, d, E, k)
      # kp_star <- list_kp_proposal$kp_star[1:(P-1)]
      # Sigma_star <- list_kp_proposal$Sigma_star
      # Sigma_star <- as.matrix(Sigma_star[1:(P-1),1:(P-1)])
      # kp_proposed_pm1 <- kp[t,1:(P-1)]
      # kp_proposed <- c(kp_proposed_pm1, - sum(kp_proposed_pm1))
      
      list_kp_proposal <- findProposalKpAdditive_m1(a2, b2, d, E, k)
      kp_star <- list_kp_proposal$kp_star
      Sigma_star <- list_kp_proposal$Sigma_star
      
      logProposal_current <- dmvnorm(as.vector(kp[1,1:(P-1)]), 
                                  as.vector(kp_star), 
                                  Sigma_star, log = T)
      
      logprior_param_current <- sum(dnorm(as.vector(kp_proposed), 0, sd_ap, log = T))
      
      
    } else if(gamma_a == "NP"){ # in nonparametric model
      
      logProposal_current <- 0
      logprior_param_current <- 0
      
      for (t in 1:Y) {
        
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
      
        # list_kp_proposal <- findProposalKptNonparametrics(t, a2, b2, d, 
        #                                                   E, k, kp, sd_kp)
        # kp_star <- list_kp_proposal$kp_star
        # Sigma_star <- list_ap_proposal$Sigma_star
        # Sigma_star <- as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4
        # logProposal_current <- logProposal_current + 
        #   dmvnorm(as.vector(kp[t,1:(P-1)]), 
        #           as.vector(kp_star[1:(P-1)]), 
        #           Sigma_star, log = T)
        
        list_kp_proposal <- findProposalKptNonparametrics_m1(t, a2, b2, d, 
                                                             E, k, kp, sd_kp)
        kp_star <- list_kp_proposal$kp_star  
        Sigma_star <- list_kp_proposal$Sigma_star
        
        logProposal_current <- logProposal_current + 
          dmvnorm(as.vector(kp[t,1:(P-1)]), 
                  as.vector(kp_star), 
                  Sigma_star, log = T)
        
        logprior_param_current <- logprior_param_current + 
          sum(dnorm(as.vector(kp[t,]), 0, sd_kp, log = T))
        
      }
      
    }

    newConfig <- sample(setdiff(k_configs, gamma_k), 1)
    
    logPriorStar <- log(p_kp[which(k_configs == newConfig)])
    
    if(newConfig == "B"){ # proposing base model
      
      kp_proposed <- matrix(0, Y, P)
      
      logProposal_star <- 0
      logprior_param_star <- 0
      
    } else if (newConfig == "A"){ # proposing additive model
      
      # old
      {
        # list_kp_proposal <- findProposalKpAdditive(a2, b2, d, E, k)
        # kp_star <- list_kp_proposal$kp_star
        # Sigma_star <- list_kp_proposal$Sigma_star
        # 
        # kp_proposed_pm1 <- mvrnorm(
        #   n = 1,
        #   kp_star[1:(P-1)],
        #   as.matrix(Sigma_star[1:(P-1),1:(P-1)])
        # )
        # 
        # kp_proposed <- c(kp_proposed_pm1, - sum(kp_proposed_pm1))
        # 
        # logProposal_star <- dmvnorm(as.vector(kp_proposed_pm1),
        #                             as.vector(kp_star[1:(P-1)]),
        #                             as.matrix(Sigma_star[1:(P-1),1:(P-1)]),
        #                             log = T)
      }
      
      # new
      {
        list_kp_proposal <- findProposalKpAdditive_m1(a2, b2, d, E, k)
        kp_star <- list_kp_proposal$kp_star
        Sigma_star <- list_kp_proposal$Sigma_star
        
        kp_proposed_pm1 <- mvrnorm(
          n = 1,
          kp_star,
          as.matrix(Sigma_star)
        )
        
        kp_proposed <- c(kp_proposed_pm1, - sum(kp_proposed_pm1))
        
        logProposal_star <- dmvnorm(as.vector(kp_proposed_pm1),
                                    as.vector(kp_star),
                                    as.matrix(Sigma_star),
                                    log = T)
      }
      
      logprior_param_star <- sum(dnorm(as.vector(kp_proposed), 0, 
                                       sd_kp, log = T))
      
      kp_proposed <- matrix(kp_proposed, Y, P, byrow = T)
      
    } else if (newConfig == "NP"){ # proposing nonparametric model
      
      # loglik_star <- 0
      logProposal_star <- 0
      logprior_param_star <- 0
      
      kp_proposed <- matrix(NA, Y, P)
      kp_star_all <- matrix(NA, Y, P)
      
      for (t in 1:Y) {
        
        # list_kp_proposal <- findProposalKptNonparametrics(t, a2, b2, d, 
        #                                                   E, k, kp, sd_kp)
        # kp_star <- list_kp_proposal$kp_star
        # Sigma_star <- list_kp_proposal$Sigma_star
        # kpt_m1 <- (kp_star - mean(kp_star))[1:(P-1)]
        # Sigma_star <- as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4
        list_kp_proposal <- findProposalKptNonparametrics_m1(t, a2, b2, d, 
                                                          E, k, kp, sd_kp)
        kpt_m1 <- list_kp_proposal$kp_star  
        Sigma_star <- list_kp_proposal$Sigma_star
        
        kpt_proposed_pm1 <- mvrnorm(
          n = 1,
          # kp_star[1:(P-1)],
          kpt_m1,
          Sigma_star
          # as.matrix(Sigma_star[1:(P-1),1:(P-1)]) / 4
        )
        
        kp_proposed[t,] <- c(kpt_proposed_pm1, - sum(kpt_proposed_pm1))
        
        logprior_param_star <- logprior_param_star + 
          sum(dnorm(as.vector(kp_proposed[t,]), 0, 
                    sd_kp, log = T))
        
        logProposal_star <- logProposal_star + 
          dmvnorm(as.vector(kpt_proposed_pm1), 
                  as.vector(kpt_m1), 
                  as.matrix(Sigma_star), 
                  log = T)
        # logProposal_star <- logProposal_star + 
        #   dmvnorm(as.vector(kpt_proposed_pm1), 
        #           as.vector(kp_star[1:(P-1)]), 
        #           as.matrix(Sigma_star[1:(P-1),1:(P-1)]), 
        #           log = T)
      
        # kp_star_all[t,] <- kp_star
      }
      
    }
    
    kkp_proposed <- 
      matrix(k, Y, P, byrow = F) + kp_proposed
    loglik_proposal <- loglik_LCp_r(a2, b2, kkp_proposed, d, E)
    
    # kkp_proposed <- matrix(k, Y, P, byrow = F) + kp_star_all
    # loglik_star <- loglik_LCp_r(a2, b2, kkp_proposed, d, E)
    
    # loglik_proposal <- 0
    # for (t in 1:Y) {
    #   loglik_kpt_current <- loglik_kpt_r(t, a2, b2, k, kp, d, E)
    #   loglik_proposal <- loglik_proposal +
    #     (loglik_kpt_current(kp_proposed[t,]) - 
    #        )
    # }
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current +  # likelihood 
        logPriorStar - logPriorCurrent + # prior of the model
        logprior_param_star - logprior_param_current + # prior of the parameters
        logProposal_current - logProposal_star # proposal of the parameters
        )
    
    if(runif(1) < mh_ratio){
      gamma_k <- newConfig
      kp <- kp_proposed
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
    gamma_a_output[currentIter] <- which(a_configs == gamma_a)
    gamma_b_output[currentIter] <- which(b_configs == gamma_b)
    gamma_k_output[currentIter] <- which(k_configs == gamma_k)
    
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

# SAVE OUTPUT -----

setwd(here("Results","Sex"))
# setwd(here("Results","Product"))

save(a_output,
     b_output,
     k_output ,
     ap_output,
     bp_output,
     kp_output,
     gamma_a_output,
     gamma_b_output,
     gamma_k_output,
     loglik_output ,
     file = "results_products.rda"
)

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
                          ymax = kp_CI[2,])) + 
    geom_point(color = "red") +
    geom_errorbar() 
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
  
  realData <- F
  
  if(!realData){
    m_true <- apply(data_plot, 1, function(x){
      
      idx_age <- x[1]
      idx_year <- x[2]
      idx_prod <- x[3]
      
      a_true[idx_age] + ap_true[idx_age, idx_prod] + 
        (b_true[idx_age] + bp_true[idx_age, idx_prod]) *
        (k_true[idx_year] + kp_true[idx_year, idx_prod])
      
      
    })  
    
    data_plot$True <- m_true
  }
  
  data_plot$CI_min <- m_CI[1,]
  data_plot$CI_max <- m_CI[2,]
 
  data_plot <- data_plot %>% 
    rename(Age = Var1,
           Year = Var2,
           Product = Var3) %>% 
    mutate(Product = factor(Product)) 
  
  # levels(data_plot$Product) <- dimnames(d)[[3]]
  levels(data_plot$Product) <- 1:P
  
  data_plot$crudeRate <- as.vector(sapply(1:nrow(data_plot), function(i){
    idx_age <- which(1:X == data_plot$Age[i])
    idx_year <- which(1:Y == data_plot$Year[i])
    idx_product <- which(1:P == as.numeric(data_plot$Product)[i])
    log((d / E)[idx_age, idx_year, idx_product])
  }))
  
  data_plot$zeroDeaths <- is.infinite(data_plot$crudeRate)
    
  data_plot$crudeRate[is.infinite(data_plot$crudeRate)] <- 
    min(data_plot$crudeRate[!is.infinite(data_plot$crudeRate)], na.rm = T) - 1
  
  # by year
  
  labeller_years <- as.character(years)
  names(labeller_years) <- 1:Y
  
  years_subset <- 1:15
  
  data_plot %>% 
    filter(Year %in% years_subset) %>% 
    ggplot(aes(x = Age,
               # y = True,
               y = crudeRate,
               ymin = CI_min,
               ymax = CI_max,
               color = Product,
               group = Product,
               fill = Product)) + 
    # geom_line(size = 1) +
    geom_point(size = 1) + 
    geom_ribbon(alpha = .3) + theme_bw() + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 14)) + 
    ylab("Log Mortality Rate") + 
    scale_x_continuous(breaks = seq(1,X,by = 10), 
                       labels = ages[seq(1,X,by = 10)]) + 
    facet_wrap(~ Year,
               labeller = labeller(Year = labeller_years))
  
  ggsave("Year_MortalityCurves.jpeg")
  
  # by age
  
  labeller_ages <- as.character(ages)
  names(labeller_ages) <- 1:X
  
  data_plot %>% 
    ggplot(aes(x = Year,
               y = True,
               ymin = CI_min,
               ymax = CI_max,
               color = Product,
               group = Product,
               fill = Product)) + 
    # geom_line(size = 1) +
    geom_point(size = 1) +
    geom_ribbon(alpha = .3) + theme_bw() + 
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 14)) + 
    ylab("Log Mortality Rate") + 
    scale_x_continuous(breaks = seq(1,Y,by = 10), 
                       labels = years[seq(1,Y,by = 10)]) + 
    facet_wrap(~ Age,
               labeller = labeller(Age = labeller_ages))
  
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

# DIAGNOSTICS ------

qplot(1:niter, a_output[,2]) + geom_hline(aes(yintercept = a_true[2]))
qplot(1:niter, ap_output[,4,2]) + geom_hline(aes(yintercept = ap_true[4,2]))
qplot(1:niter, b_output[,1]) + geom_hline(aes(yintercept = b_true[1]))
qplot(1:niter, bp_output[,2,1]) + geom_hline(aes(yintercept = bp_true[2,1]))
qplot(1:niter, k_output[,2]) + geom_hline(aes(yintercept = k_true[2]))
qplot(1:niter, kp_output[,1,1]) + geom_hline(aes(yintercept = kp_true[1,1]))
