
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
library(mgcv)
library(coda)
library(wiqid)

sourceCpp(here("Code/ProductSimple","code.cpp"))  
source(here("Code/ProductSimple","functions.R"))

# SIMULATED DATA ----------

realData <- F

X <- 30 # 30 # ages
Y <- 15 # 36 #  years
P <- 5 # products

ages <- 50 + 1:X
years <- 2000 + 1:Y
products <- as.character(1:P)

E <- array(rpois(X * Y * P, lambda = 1000), dim = c(X, Y, P))

# a
modelA <- 3

list_a <- simulateMultModel(X, P, meanVar = -3,
                            stdev = .5, modelA, sum0 = F)
ax_true <- list_a$y
axp_true <- list_a$yip
c1p_a_true <- list_a$c1p
c2xp_a_true <- list_a$c2p

# k
modelK <- 1

list_k <- simulateMultModel(Y, P, meanVar = 0,
                            stdev = .5, modelK, sum0 = T)
kt_true <- list_k$y
ktp_true <- list_k$yip
c1p_k_true <- list_k$c1p
c2tp_k_true <- list_k$c2p

# g
modelG <- 2
list_g <- simulateMultModel(X + Y - 1, P, meanVar = 0,
                            stdev = .5, modelG, sum0 = T, first0 = T)
gc_true <- list_g$y
(gcp_true <- list_g$yip)
c1p_g_true <- list_g$c1p
c2xp_g_true <- list_g$c2p

# p
dp_true <- rnorm(P, sd = .5)
dp_true <- dp_true - mean(dp_true)

# simulate data
d <- array(NA, dim = c(X, Y, P))
m_true <- array(NA, dim = c(X, Y, P))
for (x in 1:X) {
  for (t in 1:Y) {
    for (p in 1:P) {
      m_true[x,t,p] <- exp(
        axp_true[x,p] + ktp_true[t,p] +
          gcp_true[t - x + X,p] +
          dp_true[p]
        )
      d[x,t,p] <- rpois(1, E[x,t,p] * m_true[x,t,p])  
    }
  }
}

# idx <- sample(1:(X*Y*P), replace = T, size = 100)
# d[idx] <- NA
# E[idx] <- NA

# REAL DATA PRODUCTS ------

realData <- T

load(here("Data","data_products.rda"))

d <- d[10:40,,]
E <- E[10:40,,]

d <- d + .0001

X <- dim(d)[1]
Y <- dim(d)[2]
P <- dim(d)[3]

ages <- dimnames(d)[[1]]
years <- dimnames(d)[[2]]
products <- dimnames(d)[[3]]

# SIMULATE GIVEN REAL DATA PRODUCTS ------

realData <- T

load(here("Data","data_products.rda"))

d <- d[10:60,,-3]
E <- E[10:60,,-3]

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

age_subset <- c(50,85)
year_subset <- c(1990,2020)

d_all <- array(NA, dim = c(age_subset[2] - age_subset[1] + 1, 
                           year_subset[2] - year_subset[1] + 1, P))
E_all <- array(NA, dim = c(age_subset[2] - age_subset[1] + 1, 
                           year_subset[2] - year_subset[1] + 1, P))

for (p in 1:P) {
  
  country <- countries[p]
  
  load(here("Data","Countries",paste0(country,"-data.rda")))
  
  X_subset <- which(as.numeric(dimnames(d)[[1]]) >= age_subset[1] &
                      as.numeric(dimnames(d)[[1]]) <= age_subset[2])
  Y_subset <- which(as.numeric(dimnames(d)[[2]]) >= year_subset[1] &
                      as.numeric(dimnames(d)[[2]]) <= year_subset[2])
  
  d <- d[X_subset, Y_subset]
  E <- E[X_subset, Y_subset]
  
  ages <- as.numeric(dimnames(d)[[1]])
  years <- as.numeric(dimnames(d)[[2]])  
  
  d_all[,,p] <- d
  E_all[,,p] <- E
}

d <- d_all
E <- E_all

dimnames(d) <- list(
  ages, years, countries
)

X <- dim(d)[1]
Y <- dim(d)[2]

# MCMC -----

nburn <- 200
niter <- 200

realData <- T

# prior
{
  # a_configs <- c("B","A")
  # nConfigs_a <- length(a_configs)
  # p_ap <- c(.6, .4)
  a_configs <- 0:2
  nConfigs_a <- length(a_configs)
  p_ap <- c(.6, .3, .1)
  
  # b_configs <- c("B","A")
  # nConfigs_b <- length(b_configs)
  # p_bp <- c(.6, .4)
  b_configs <- 0:2
  nConfigs_b <- length(b_configs)
  p_bp <- c(.6, .3, .1)
  
  # k_configs <- c("B","A")
  # nConfigs_k <- length(k_configs)
  # p_kp <- c(.6, .4)
  k_configs <- 0:2
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
  a_output <- array(NA, dim = c(niter, X, P))
  k_output <- array(NA, dim = c(niter, Y, P))
  g_output <- array(NA, dim = c(niter, X + Y - 1, P))
  d_output <- matrix(NA, niter, P)
  
  gamma_a_output <- rep(NA, niter)
  gamma_k_output <- rep(NA, niter)
  gamma_g_output <- rep(NA, niter)
  
  loglik_output <- rep(NA, niter)
  m_output <- array(NA, dim = c(X, Y, P, niter))
}

# paramsUpdate
{
  update_a <- T
  update_k <- T
  update_g <- T
  update_d <- T
  update_gammaa <- T
  update_gammak <- T
  update_gammag <- T
  
  # true params
  {
    param_a_true <- !update_a
    param_k_true <- !update_k
    param_g_true <- !update_g
    param_d_true <- !update_d
    param_gammaa_true <- !update_gammaa
    param_gammak_true <- !update_gammak
    param_gammag_true <- !update_gammag
  }
}

# starting values
{
  # true
  if(!realData) {
    
    param_a <- list_a
    param_k <- list_k
    param_g <- list_g
    dp <- dp_true
    
    gamma_a <- modelA
    gamma_k <- modelK
    gamma_g <- modelG
    
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
      
      # laplace_fit <- optim(
      #   par = c(
      #     a_start,
      #     rep(1 / X, X),
      #     rep(0, Y)),
      #   fn = loglik_r(d, E),
      #   method = c("BFGS")
      # )
      
      # param_star <- laplace_fit$par
      # a <- param_star[1:X]
      # b <- param_star[X + 1:X]
      # k <- param_star[2 *X + 1:Y]
      
      a <- a_start
      k <- rep(0, Y)
      g <- rep(0, X + Y - 1)
      dp <- rep(0, P)
      # dp <- dp_true#rep(0, P)
      
      # - loglik_r(d, E)(param_star) 
      # - loglik_r(d, E)(c(a_true, b_true, k_true))
      
      # sumb <- sum(b)
      # b <- b / sumb
      # k <- k * sumb
      # meank <- mean(k)
      # k <- k - meank
      # a <- a + b * meank
      
      if(gamma_a == 1){
        param_a <- list("y" = a,
                        "yip" = matrix(a, X, P, byrow = F))
      } else if(gamma_a == 2){
        param_a <- list("y" = a,
                        "c1p" = rep(1, P),
                        "yip" = matrix(a, X, P, byrow = F))
      } else if(gamma_a == 3){
        param_a <- list("y" = a,
                        "yip" = matrix(a, X, P, byrow = F))
      }
      
      if(gamma_k == 1){
        param_k <- list("y" = rep(0, Y),
                        "yip" = matrix(k, Y, P, byrow = F))
        
      } else if(gamma_k == 2){
        param_k <- list("y" = rep(0, Y),
                        "c1p" = rep(1, P),
                        "yip" = matrix(k, Y, P, byrow = F))
        
      } else if(gamma_k == 3) {
        param_k <- list("y" = rep(0, Y),
                        "c2p" = matrix(1, Y, P),
                        "yip" = matrix(k, Y, P, byrow = F))
        
      }
      
      # param_k <- list_k
      
      if(gamma_g == 1){
        param_g <- list("y" = rep(0, X + Y - 1),
                        "yip" = matrix(g, X + Y - 1, P, byrow = F))
      } else if(gamma_g == 2){
        param_g <- list("y" = rep(0, X + Y - 1),
                        "c1p" = rep(1, P),
                        "yip" = matrix(g, X + Y - 1, P, byrow = F))
        
      } else if(gamma_g == 3) {
        param_g <- list("y" = rep(0,  X + Y - 1),
                        "c2p" = matrix(1, X + Y - 1, P),
                        "yip" = matrix(g, X + Y - 1, P, byrow = F))
        
      }
      
      # param_g <- list_g
      
    }
    
  }
  
  a <- apply(d / E, 1, function(x){
    mean(log(x + 1), na.rm = T)
  })
  k <- rep(0, Y)
  g <- rep(0, X + Y - 1)
  dp <- rep(0, P)
  
  if(param_gammaa_true){
    
    gamma_a <- modelA
    
  } else {
    
    gamma_a <- 1
    
  }
  
  if(param_gammak_true){
    
    gamma_k <- modelK
    
  } else {
    
    gamma_k <- 1
    
  }
  
  if(param_gammag_true){
    
    gamma_g <- modelG
    
  } else {
    
    gamma_g <- 1
    
  }
  
  if(param_a_true){
    
    param_a <- list_a
    
  } else {
    
    if(gamma_a == 1){
      param_a <- list("y" = a,
                      "yip" = matrix(a, X, P, byrow = F))
    } else if(gamma_a == 2){
      param_a <- list("y" = a,
                      "c1p" = rep(1, P),
                      "yip" = matrix(a, X, P, byrow = F))
    } else if(gamma_a == 3){
      param_a <- list("y" = a,
                      "yip" = matrix(a, X, P, byrow = F))
    }
    
  }
  
  if(param_k_true){
    
    param_k <- list_k
    
  } else {
    
    if(gamma_k == 1){
      param_k <- list("y" = rep(0, Y),
                      "yip" = matrix(k, Y, P, byrow = F))
      
    } else if(gamma_k == 2){
      param_k <- list("y" = rep(0, Y),
                      "c1p" = rep(1, P),
                      "yip" = matrix(k, Y, P, byrow = F))
      
    } else if(gamma_k == 3) {
      param_k <- list("y" = rep(0, Y),
                      "c2p" = matrix(1, Y, P),
                      "yip" = matrix(k, Y, P, byrow = F))
      
    }
    
  }
  
  if(param_g_true){
    
    param_g <- list_g
    
  } else {
    
    if(gamma_g == 1){
      param_g <- list("y" = rep(0, X + Y - 1),
                      "yip" = matrix(g, X + Y - 1, P, byrow = F))
    } else if(gamma_g == 2){
      param_g <- list("y" = rep(0, X + Y - 1),
                      "c1p" = rep(1, P),
                      "yip" = matrix(g, X + Y - 1, P, byrow = F))
      
    } else if(gamma_g == 3) {
      param_g <- list("y" = rep(0,  X + Y - 1),
                      "c2p" = matrix(1, X + Y - 1, P),
                      "yip" = matrix(g, X + Y - 1, P, byrow = F))
      
    }
    
  }
  
  if(param_d_true){
    
    dp <- dp_true 
    
  } else {
    
    dp <- rep(0, P)
    
  }
  
}

for (iter in 1:(niter + nburn)) {
  
  # print(a[1:5])
  # print(b[1:5])
  # print(k[1:5])
  # 
  print(paste0("Gamma A = ",gamma_a,
               ", Gamma K = ",gamma_k,
               ", Gamma G = ",gamma_g))
  
  if(iter %% 10 == 0){
    if(iter < nburn){
      print(paste0("Burn-in Iteration = ",iter))  
    } else {
      print(paste0("Iteration = ",iter - nburn))  
    }  
  }
  
  # update a
  # print("Update a")
  if (update_a) {
    
    m_terms <- computeTerms(1, param_a, param_k, param_g, dp)
    
    if(gamma_a == 1){
      
      ax <- param_a$y
      
      for (x in 1:X) {
        
        param_current <- ax[x]
        
        param_star <- log(sum(d[x,,], na.rm = T) / 
                            sum(exp(m_terms[x,,] + log(E[x,,])), na.rm = T))
        
        hessian_ax_star <- derl2der2a(param_star, m_terms[x,,] + log(E[x,,]))
        sd_star <- sqrt(- 1 / hessian_ax_star)
        
        param_proposed <- rnorm(1, param_star, sd_star)
        
        # create loglikelihood function of k given ax and bx
        loglik_current <- loglik_LC(d[x,,], exp(param_current + m_terms[x,,] + log(E[x,,])))
        loglik_star <- loglik_LC(d[x,,], exp(param_proposed + m_terms[x,,] + log(E[x,,])))
        
        logproposal_current <- dnorm(param_current, param_star, sd_star, log = T)
        logproposal_proposed <- dnorm(param_proposed, param_star, sd_star, log = T)
        
        mh_ratio <- exp(loglik_star - loglik_current +
                          logproposal_current  - logproposal_proposed)
        
        if(runif(1) < mh_ratio){
          ax[x] <- param_proposed
        }  
           
      }
      
      # readjust parameters
      {
        param_a$y <- ax
        param_a$yip <- matrix(ax, X, P, byrow = F)
      }
      
    }
    
    if(gamma_a == 2){
      
      ax <- param_a$y
      c1p <- param_a$c1p
      
      # update ax  
      
      for (x in 1:X) {
        
        param_current <- ax[x]
        
        loglik_ax_current <- loglik_ax_cp(d[x,,], c1p, 
                                        m_terms[x,,] + log(E[x,,]))
        
        modelFit <- optim(-1,
                          lower = -10,
                          upper = 0,
              loglik_ax_current, 
              hessian = T,
              method = "Brent")
        
        param_star <- modelFit$par
        hessian_star <- modelFit$hessian
        
        sd_star <- sqrt( 1 / hessian_star)
        
        param_proposed <- rnorm(1, param_star, sd_star)
        
        loglik_current <- - loglik_ax_current(param_current)
        loglik_star <- - loglik_ax_current(param_star)
        
        logproposal_current <- dnorm(param_current, param_star, sd_star, log = T)
        logproposal_proposed <- dnorm(param_proposed, param_star, sd_star, log = T)
        
        mh_ratio <- exp(loglik_star - loglik_current +
                          logproposal_current  - logproposal_proposed)
        
        if(runif(1) < mh_ratio){
          ax[x] <- param_proposed
        }  
        
      }
      
      # update cp
      
      c1p_star <- rep(0, P)
      
      mu_ax_star <- rep(0, P)
      sd_ax_star <- rep(0, P)
      
      for (p in 1:P) {
        
        param_current <- c1p[p]
        
        loglik_cp_current <- loglik_cp_ax(d[,,p], ax, 
                                        m_terms[,,p] + log(E[,,p]))
        
        modelFit <- optim(0,
                          lower = -3,
                          upper = 3,
              loglik_cp_current, 
              hessian = T,
              method = "Brent")
        
        param_star <- modelFit$par
        hessian_star <- modelFit$hessian
        
        sd_star <- sqrt( 1 / hessian_star)
        
        mu_ax_star[p] <- param_star
        sd_ax_star[p] <- sd_star
        
      }
       
      Sigma_star <- diag(sd_ax_star^2, nrow = P)
      param_proposed <- sampleMTconstraint_c1p(mu_ax_star, Sigma_star)
      param_current <- c1p
      
      loglik_current <- 0
      loglik_star <- 0
      for (p in 1:P) {
        loglik_cp_current <- loglik_cp_ax(d[,,p], ax, 
                                          m_terms[,,p] + log(E[,,p]))
        
        loglik_current <- loglik_current - loglik_cp_current(param_current[p])
        loglik_star <- loglik_star - loglik_cp_current(param_proposed[p])
      }
      
      logproposal_current <- dmt_cpp(param_current, nu = 3, mu_ax_star, Sigma_star, T)
      logproposal_proposed <- dmt_cpp(param_proposed, nu = 3, mu_ax_star, Sigma_star, T)
      
      mh_ratio <- exp(loglik_star - loglik_current +
                        logproposal_current  - logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        c1p <- param_proposed
      }  
      
      # readjust parameters
      {
        param_a$y <- ax
        param_a$c1p <- c1p
        param_a$yip <- matrix(ax, X, P, byrow = F) * matrix(c1p, X, P, byrow = T)
      }
     
      
    }
    
    if(gamma_a == 3){ # update each axp independently
      
      axp <- param_a$yip
      
      for (x in 1:X) {
        for (p in 1:P) {
          
          param_current <- axp[x,p]
          
          param_star <- log(sum(d[x,,p], na.rm = T) / 
                              sum(exp(m_terms[x,,p] + log(E[x,,p])), na.rm = T))
          
          hessian_ax_star <- derl2der2a(param_star, m_terms[x,,p] + log(E[x,,p]))
          sd_star <- sqrt(- 1 / hessian_ax_star)
          
          param_proposed <- rnorm(1, param_star, sd_star)
          
          # create loglikelihood function of k given ax and bx
          loglik_current <- loglik_LC(d[x,,p], exp(param_current + m_terms[x,,p] + log(E[x,,p])))
          loglik_star <- loglik_LC(d[x,,p], exp(param_proposed + m_terms[x,,p] + log(E[x,,p])))
          
          logproposal_current <- dnorm(param_current, param_star, sd_star, log = T)
          logproposal_proposed <- dnorm(param_proposed, param_star, sd_star, log = T)
          
          mh_ratio <- exp(loglik_star - loglik_current +
                            logproposal_current  - logproposal_proposed)
          
          if(runif(1) < mh_ratio){
            axp[x,p] <- param_proposed
          }  
          
        }
      }
      
      # readjust parameters
      {
        param_a$yip <- axp
        
      }
      
    }
    
  }
  
  # update k_p
  # print("Update k")
  if (update_k) {
    
    m_terms <- computeTerms(2, param_a, param_k, param_g, dp)
    
    # update kt
    c2p <- computeCP(param_k, gamma_k, Y)
    
    param_k$y <- updateKT_CP(param_k$y, c2p, d, E, m_terms)
    
    if(gamma_k == 1){ 
      
      param_k$yip <- 
        matrix(param_k$y, Y, P, byrow = F)
      
    } else if(gamma_k == 2){
      
      param_k$c1p <- updateC1P_KT(param_k$c1p, param_k$y, d, E, m_terms)
      
      param_k$yip <- 
        matrix(param_k$y, Y, P, byrow = F) * 
        matrix(param_k$c1p, Y, P, byrow = T)
    
    } else if(gamma_k == 3){
      
      param_k$c2p <- updateC2P_KT(param_k$c2p, param_k$y, d, E, m_terms)
      
      param_k$yip <- 
        matrix(param_k$y, Y, P, byrow = F) * param_k$c2p
    }
    
    
  }
  
  # update g_c
  # print("Update g")
  if (update_g) {
    
    m_terms <- computeTerms(3, param_a, param_k, param_g, dp)
    
    # update kt
    c2p <- computeCP(param_g, gamma_g, X + Y - 1)
    
    param_g$y <- updateGC_CP(param_g$y, c2p, d, E, m_terms)
    
    param_g$yip <- 
      matrix(param_g$y, X + Y - 1, P, F)
    
    if(gamma_g == 2){
      
      param_g$c1p <- updateC1P_GC(param_g$c1p, param_g$y, d, E, m_terms)
      
      param_g$yip <- 
        matrix(param_g$y, X + Y - 1, P, F) * 
        matrix(param_g$c1p, X + Y - 1, P, T)
      
    }
    
    if(gamma_g == 3){
      
      param_g$c2p <- updateC2P_GC(param_g$c2p, param_g$y, d, E, m_terms)
      
      param_g$yip <- 
        matrix(param_g$y, X + Y - 1, P, F) * param_g$c2p
      
    }
    
  }
  
  # update d_p
  # print("Update d")
  if (update_d){
    
    m_terms <- computeTerms(4, param_a, param_k, param_g, dp)
    
    mu_dp_star <- rep(0, P)
    sd_dp_star <- rep(0, P)
    
    for (p in 1:P) {
        
        param_current <- dp[p]
        
        param_star <- log(sum(d[,,p], na.rm = T) / 
                            sum(exp(m_terms[,,p] + log(E[,,p])), na.rm = T))
        
        hessian_star <- derl2der2a(param_star, m_terms[,,p] + log(E[,,p]))
        sd_star <- sqrt(- 1 / hessian_star)
      
        mu_dp_star[p] <-  param_star 
        sd_dp_star[p] <- sd_star
        
    }
    
    Sigma_star <- diag(sd_dp_star^2, nrow = P)
    param_proposed <- sampleMTconstraint_k(mu_dp_star, Sigma_star)
    param_current <- dp
    
    loglik_current <- 0
    loglik_star <- 0
    
    for (p in 1:P) {
      
      loglik_current <- loglik_current + 
        loglik_LC(d[,,p], exp(param_current[p] + m_terms[,,p] + log(E[,,p])))
      loglik_star <- loglik_star + 
        loglik_LC(d[,,p], exp(param_proposed[p] + m_terms[,,p] + log(E[,,p])))
      
    }
    
    logproposal_current <- dmt_cpp(param_current, nu = 3, mu_dp_star, Sigma_star, T)
    logproposal_proposed <- dmt_cpp(param_proposed, nu = 3, mu_dp_star, Sigma_star, T)
    
    mh_ratio <- exp(loglik_star - loglik_current +
                      logproposal_current  - logproposal_proposed)
    
    if(runif(1) < mh_ratio){
      dp <- param_proposed
    }
    
  }
    
  # propose gamma_a
  # print("Update gamma_a")
  if (update_gammaa) {
    
    m_terms <- computeTerms(1, param_a, param_k, param_g, dp)
    
    if(gamma_a == 1){
    
      ax_current <- param_a$y
      
      list_proposal_ax <- findProposalAX(d, E, m_terms)
      mu_ax0_star <- list_proposal_ax$mu_ax_star
      sd_ax0_star <- list_proposal_ax$sd_ax_star
      
      list_proposal_c1p <- findProposalCP_AX(d, E, m_terms)
      mu_ax_star <- list_proposal_c1p$mu_ax_star
      sd_ax_star <- list_proposal_c1p$sd_ax_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_star <- list_proposal_c1p$Sigma_star
      
      ax_proposed <- sapply(1:X, function(x){
        # rt2(1, mu_ax_star[x], sd_ax_star[x], df = 3)
        rnorm(1, mu_ax_star[x], sd_ax_star[x])
      })
      # cpm1_proposed <- mrt2(cpm1_star, Sigma_c1p_star, df = 3)
      cpm1_proposed <- mvrnorm(n = 1, cpm1_star, Sigma_c1p_star)
      c1p_proposed <- c(cpm1_proposed, P - sum(cpm1_proposed))
      
      loglik_current <- loglik_axcp(ax, rep(1, P), d, E, m_terms)
      loglik_proposed <- loglik_axcp(ax_proposed, c1p_proposed, d, E, m_terms)
      
      logprior_proposed <- 0#sum(dgamma(c1p_proposed, 100, 100, log = T))
      
      logproposal_proposed <- 
        sum(sapply(1:X, function(x){
          # log(dt2(ax_proposed[x], mu_ax_star[x], sd_ax_star[x], df = 3))
          dnorm(ax_proposed[x], mu_ax_star[x], sd_ax_star[x], log = T)
        })) + 
        # dmt_cpp(cpm1_proposed, nu = 3, 
                      # cpm1_star, Sigma_c1p_star, returnLog = T)
        dmvnorm(cpm1_proposed, cpm1_star, Sigma_c1p_star, log = T)
      
      logproposal_current <- sum(sapply(1:X, function(x){
        # log(dt2(ax_current[x], mu_ax0_star[x], sd_ax0_star[x], df = 3))
        dnorm(ax_current[x], mu_ax0_star[x], sd_ax0_star[x], log = T)
      }))
      
      logproposal_2 <- log(.5) # from 2 to 1
      logproposal_1 <- log(1) # from 1 to 2
      
      mh_ratio <- exp(loglik_proposed - loglik_current +
                        logprior_proposed +
                        logproposal_2 - logproposal_1 - 
                        logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        
        param_a$c1p <- c1p_proposed
        param_a$yip <- matrix(ax, X, P, byrow = F) * 
          matrix(c1p_proposed, X, P, byrow = T)
        
        gamma_a <- 2
        
      }
      
    } else if (gamma_a == 2){
      
      gamma_a_proposed <- sample(c(1,3), 1)
      
      ax_current <- param_a$y
      c1p_current <- param_a$c1p
      
      list_proposal_c1p <- findProposalCP_AX(d, E, m_terms)
      mu_ax_star <- list_proposal_c1p$mu_ax_star
      sd_ax_star <- list_proposal_c1p$sd_ax_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_star <- list_proposal_c1p$Sigma_star
      
      logproposal_current <- 
        sum(sapply(1:X, function(x){
          # log(dt2(ax_current[x], mu_ax_star[x], sd_ax_star[x], df = 3))
          dnorm(ax_current[x], mu_ax_star[x], sd_ax_star[x], log = T)
        })) + 
        # dmt_cpp(c1p_current[1:(P-1)], nu = 3, 
                      # cpm1_star, Sigma_c1p_star, returnLog = T)
        dmvnorm(c1p_current[1:(P-1)], cpm1_star, Sigma_c1p_star, log = T)
      
      loglik_current <- loglik_axcp(ax_current, c1p_current,d, E, m_terms)
      
      logprior_current <- 0#sum(dgamma(c1p_current, 100, 100, log = T))
      
      if(gamma_a_proposed == 1){
        
        list_proposal_ax <- findProposalAX(d, E, m_terms)
        mu_ax0_star <- list_proposal_ax$mu_ax_star
        sd_ax0_star <- list_proposal_ax$sd_ax_star
        
        ax_proposed <- sapply(1:X, function(x){
          # rt2(1, mu_ax0_star[x], sd_ax0_star[x], df = 3)
          rnorm(1, mu_ax0_star[x], sd_ax0_star[x])
        })
        
        loglik_proposed <- loglik_axcp(ax_proposed, rep(1, P), d, E, m_terms)
        
        logproposal_proposed <- sum(sapply(1:X, function(x){
          # log(dt2(ax_proposed[x], mu_ax0_star[x], sd_ax0_star[x], df = 3))
          dnorm(ax_proposed[x], mu_ax0_star[x], sd_ax0_star[x], log = T)
        }))
        
        logproposal_2 <- log(1) # from 2 to 1
        logproposal_1 <- log(.5) 
        
        mh_ratio <- exp(loglik_proposed - loglik_current -
                          logprior_current +
                          logproposal_2 - logproposal_1 + 
                          logproposal_current - logproposal_proposed)
        
        if(runif(1) < mh_ratio){
          
          param_a$y <- ax_proposed
          param_a$yip <- matrix(ax_proposed, X, P, byrow = F)
          gamma_a <- 1
          param_a$c1p <- NULL
          
        }
        
      } else { 
        
        list_proposal_axp <- findProposalAXP(d, E, m_terms)
        mu_axp_star <- list_proposal_axp$mu_axp_star
        sd_axp_star <- list_proposal_axp$sd_axp_star
        
        axp_proposed <- t(sapply(1:X, function(x){
          sapply(1:P, function(p){
            # rt2(1, mu_axp_star[x,p], sd_axp_star[x,p], df = 3)
            rnorm(1, mu_axp_star[x,p], sd_axp_star[x,p])
          })
        }))
        
        logproposal_proposed <- sum(sapply(1:X, function(x){
          sapply(1:P, function(p){
            # log(dt2(axp_proposed[x,p], mu_axp_star[x,p], sd_axp_star[x,p], df = 3)  )
            dnorm(axp_proposed[x,p], mu_axp_star[x,p], sd_axp_star[x,p], log = T)  
          })
        }))
        
        c1p_proposed <- sapply(1:X, function(x){
          ax_mean <- mean(axp_proposed[x,])
          sapply(1:P, function(p){
            axp_proposed[x,p] / ax_mean
          })
        })
        
        logprior_proposed <- 0#sum(dgamma(c1p_proposed, 100, 100, log = T))
        
        loglik_proposed <- loglik_axp(axp_proposed, d, E, m_terms)
        
        logproposal_2 <- log(1) # going from model 2 to model 1
        logproposal_1 <- log(.5) # going from model 1 to model 2
        
        mh_ratio <- exp(loglik_proposed - loglik_current +
                          logprior_proposed - logprior_current +
                          logproposal_2 - logproposal_1 +
                          logproposal_current - logproposal_proposed)
        
        if(runif(1) < mh_ratio){
          
          param_a$yip <- axp_proposed
          gamma_a <- 3
          param_a$c1p <- NULL
          
        }
        
      }
      
    } else if (gamma_a == 3) { 
      
      axp_current <- param_a$yip
      
      list_proposal_c1p <- findProposalCP_AX(d, E, m_terms)
      mu_ax_star <- list_proposal_c1p$mu_ax_star
      sd_ax_star <- list_proposal_c1p$sd_ax_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_star <- list_proposal_c1p$Sigma_star
      
      list_proposal_axp <- findProposalAXP(d, E, m_terms)
      mu_axp_star <- list_proposal_axp$mu_axp_star
      sd_axp_star <- list_proposal_axp$sd_axp_star
      
      ax_proposed <- sapply(1:X, function(x){
        # rt2(1, mu_ax_star[x], sd_ax_star[x], df = 3)
        rnorm(1, mu_ax_star[x], sd_ax_star[x])
      })
      # cpm1_proposed <- mrt2(cpm1_star, Sigma_c1p_star, df = 3)
      cpm1_proposed <- mvrnorm(n = 1, cpm1_star, Sigma_c1p_star)
      c1p_proposed <- c(cpm1_proposed, P - sum(cpm1_proposed))
      
      c1p_current <- sapply(1:X, function(x){
        ax_mean <- mean(axp_current[x,])
        sapply(1:P, function(p){
          axp_current[x,p] / ax_mean
        })
      })
      
      logprior_current <- 0 #sum(dgamma(c1p_current, 100, 100, log = T))
      logprior_proposed <- 0#sum(dgamma(c1p_proposed, 100, 100, log = T))
      
      loglik_proposed <- loglik_axcp(ax_proposed, c1p_proposed, d, E, m_terms)
      loglik_current <- loglik_axp(axp_current, d, E, m_terms)
      
      logproposal_proposed <- 
        sum(sapply(1:X, function(x){
          # log(dt2(ax_proposed[x], mu_ax_star[x], sd_ax_star[x], df = 3))
          dnorm(ax_proposed[x], mu_ax_star[x], sd_ax_star[x], log = T)
        })) + 
        dmvnorm(cpm1_proposed, cpm1_star, Sigma_c1p_star, log = T)
        # dmt_cpp(cpm1_proposed, nu = 3, 
        #               cpm1_star, Sigma_c1p_star, returnLog = T)
      
      logproposal_current <- sum(sapply(1:X, function(x){
        sapply(1:P, function(p){
          # log(dt2(axp_current[x,p], mu_axp_star[x,p], sd_axp_star[x,p], df = 3)  )
          dnorm(axp_current[x,p], mu_axp_star[x,p], sd_axp_star[x,p], log = T)
        })
      }))
      
      logproposal_2 <- log(1) # going from model 2 to model 1
      logproposal_1 <- log(.5) # going from model 1 to model 2
      
      mh_ratio <- exp(loglik_proposed - loglik_current +
                        logproposal_2 - logproposal_1 +
                        logproposal_current - logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        
        param_a$y <- ax_proposed
        param_a$c1p <- c1p_proposed
        param_a$yip <- 
          matrix(ax_proposed, X, P, byrow = F) * matrix(c1p_proposed, X, P, byrow = T)
        gamma_a <- 2
        
      }
      
    }
  
  }
    
  # propose gamma_k
  # print("Update gamma_k")
  if (update_gammak) {
    
    m_terms <- computeTerms(2, param_a, param_k, param_g, dp)
    
    if(gamma_k == 1){
    
      kt_current <- param_k$y
      
      list_proposal_kt <- findProposalKT(d, E, m_terms)
      mu_ktm1_0_star <- list_proposal_kt$ktm1_star
      Sigma_ktm1_0_star <- list_proposal_kt$Sigma_ktm1_star
      
      list_proposal_c1p <- findProposalC1P_KT(d, E, m_terms)
      mu_ktm1_star <- list_proposal_c1p$mu_ktm1_star
      Sigma_ktm1_star <- list_proposal_c1p$Sigma_ktm1_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_star <- list_proposal_c1p$Sigma_star
      
      # ktm1_proposed <- mrt2(mu_ktm1_star, Sigma_ktm1_star, df = 3)
      ktm1_proposed <- mvrnorm(n = 1, mu_ktm1_star, Sigma_ktm1_star)
      kt_proposed <- c(ktm1_proposed, - sum(ktm1_proposed))
      # cpm1_proposed <- mrt2(cpm1_star, Sigma_c1p_star, df = 3)
      cpm1_proposed <- mvrnorm(n = 1, cpm1_star, Sigma_c1p_star)
      c1p_proposed <- c(cpm1_proposed, P - sum(cpm1_proposed))
      
      loglik_current <- loglik_ktcp(kt_current, rep(1, P), d, E, m_terms)
      loglik_proposed <- loglik_ktcp(kt_proposed, c1p_proposed, d, E, m_terms)
      
      logproposal_proposed <- 
        # dmt_cpp(ktm1_proposed, nu = 3, mu_ktm1_star, Sigma_ktm1_star, T) + 
        # dmt_cpp(cpm1_proposed, nu = 3, cpm1_star, Sigma_c1p_star, returnLog = T)
        dmvnorm(ktm1_proposed, mu_ktm1_star, Sigma_ktm1_star, T) + 
        dmvnorm(cpm1_proposed, cpm1_star, Sigma_c1p_star, T)
      
      logproposal_current <- 
        # dmt_cpp(kt_current[1:(Y-1)], nu = 3, mu_kt0_star, Sigma_kt0_star, T) 
        dmvnorm(kt_current[1:(Y-1)], mu_ktm1_0_star, Sigma_ktm1_0_star, T)
    
      logproposal_2 <- log(.5) # from 2 to 1
      logproposal_1 <- log(1) # from 1 to 2
      
      mh_ratio <- exp(loglik_proposed - loglik_current +
                        logproposal_2 - logproposal_1 - 
                        logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        
        param_k$y <- kt_proposed
        param_k$c1p <- c1p_proposed
        param_k$yip <- matrix(kt_proposed, Y, P, byrow = F) * 
          matrix(c1p_proposed, Y, P, byrow = T)
        
        gamma_k <- 2
        
      }
      
    } else if (gamma_k == 2){
      
      gamma_k_proposed <- sample(c(1,3), 1)
      
      kt_current <- param_k$y
      c1p_current <- param_k$c1p
      
      list_proposal_c1p <- findProposalC1P_KT(d, E, m_terms)
      mu_ktm1_star <- list_proposal_c1p$mu_ktm1_star
      Sigma_ktm1_star <- list_proposal_c1p$Sigma_ktm1_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_star <- list_proposal_c1p$Sigma_star
      
      logproposal_current <- 
        dmvnorm(kt_current[1:(Y-1)], mu_ktm1_star, Sigma_ktm1_star, T) +
        # sum(sapply(1:Y, function(t){
        #   # log(dt2(kt_current[t], mu_kt_star[t], sd_kt_star[t], df = 3))
        #   dnorm(kt_current[t], mu_kt_star[t], sd_kt_star[t], log = T)
        # })) + 
        dmvnorm(c1p_current[1:(P-1)], cpm1_star, Sigma_c1p_star, T)
        # dmt_cpp(c1p_current[1:(P-1)], nu = 3, 
        #         cpm1_star, Sigma_c1p_star, returnLog = T)
      
      loglik_current <- loglik_ktcp(kt_current, c1p_current, d, E, m_terms)
      
      if(gamma_k_proposed == 1){
        
        list_proposal_kt <- findProposalKT(d, E, m_terms)
        mu_ktm1_star <- list_proposal_kt$ktm1_star
        Sigma_ktm1_star <- list_proposal_kt$Sigma_ktm1_star
        
        kt_m1_proposed <- 
          mvrnorm(n = 1, mu_ktm1_star, Sigma_ktm1_star)
        kt_proposed <- c(kt_m1_proposed, - sum(kt_m1_proposed))
        
        loglik_proposed <- loglik_ktcp(kt_proposed, rep(1, P), d, E, m_terms)
       
        logproposal_proposed <- 
          dmvnorm(kt_m1_proposed, mu_ktm1_star, Sigma_ktm1_star, log = T)
        #   sum(sapply(1:Y, function(x){
        #   # log(dt2(kt_proposed[t], mu_kt0_star[t], sd_kt0_star[t], df = 3))
        #   dnorm(kt_proposed[t], mu_kt0_star[t], sd_kt0_star[t], log = T)
        # }))
        
        logproposal_2 <- log(1) # from 2 to 1
        logproposal_1 <- log(.5) 
        
        mh_ratio <- exp(loglik_proposed - loglik_current +
                          logproposal_2 - logproposal_1 + 
                          logproposal_current - logproposal_proposed)
        
        if(runif(1) < mh_ratio){
          
          ktp_proposed <- matrix(kt_proposed, Y, P, byrow = F)
          
          param_k$y <- kt_proposed
          param_k$c1p <- c1p_proposed
          param_k$yip <- ktp_proposed
          gamma_k <- 1
          param_k$c2p <- NULL
          
        }
        
      } else { 
        
        list_proposal_c2pktp <- findProposalC2P_KT(d, E, m_terms)
        mu_ktm1_star <- list_proposal_c2pktp$mu_ktm1_star
        Sigma_ktm1_star <- list_proposal_c2pktp$Sigma_ktm1_star
        mu_c2pm1_star <- list_proposal_c2pktp$mu_c2pm1_star
        Sigma_c2pm1_star <- list_proposal_c2pktp$Sigma_c2pm1_star
        
        kt_m1_proposed <- 
          mvrnorm(n = 1, mu_ktm1_star, Sigma_ktm1_star)
        kt_proposed <- c(kt_m1_proposed, - sum(kt_m1_proposed))
        
        c2pm1_proposed <- t(sapply(1:Y, function(t){
          # sapply(1:P, function(p){
          #   # rt2(1, mu_ktp_star[t,p], sd_ktp_star[t,p], df = 3)
          #   rnorm(1, mu_ktp_star[t,p], sd_ktp_star[t,p])
          # })
          mvrnorm(n = 1, mu_c2pm1_star[t,], Sigma_c2pm1_star[t,,])
        }))
        c2p_proposed <- cbind(c2pm1_proposed, 
                              P - apply(c2pm1_proposed, 1, sum))
        
        logproposal_proposed <- 
          dmvnorm(kt_m1_proposed, mu_ktm1_star, Sigma_ktm1_star, log = T) + 
          sum(
            sapply(1:Y, function(t){
              dmvnorm(c2pm1_proposed[t,], mu_c2pm1_star[t,], Sigma_c2pm1_star[t,,], log = T)
            })
          )
        #   sum(sapply(1:Y, function(t){
        #   sapply(1:P, function(p){
        #     # log(dt2(ktp_proposed[t,p], mu_ktp_star[t,p], sd_ktp_star[t,p], df = 3)  )
        #     dnorm(ktp_proposed[t,p], mu_ktp_star[t,p], sd_ktp_star[t,p], log = T)
        #   })
        # }))
        
        ktp_proposed <- t(sapply(1:Y, function(t){
          kt_proposed[t] * c2p_proposed[t,]
        }))
        
        loglik_proposed <- loglik_ktp(ktp_proposed, d, E, m_terms)
        
        logproposal_2 <- log(1) # going from model 2 to model 1
        logproposal_1 <- log(.5) # going from model 1 to model 2
        
        mh_ratio <- exp(loglik_proposed - loglik_current +
                          logproposal_2 - logproposal_1 +
                          logproposal_current - logproposal_proposed)
        
        if(runif(1) < mh_ratio){
          
          param_k$y <- kt_proposed
          param_k$c2p <- c2p_proposed
          param_k$yip <- ktp_proposed
          gamma_k <- 3
          param_k$c1p <- NULL
          
        }
        
      }
      
    } else if (gamma_k == 3) { 
      
      ktp_current <- param_k$yip
      kt_current <- param_k$y
      c2p_current <- param_k$c2p
      
      list_proposal_c1p <- findProposalC1P_KT(d, E, m_terms)
      mu_ktm1_star <- list_proposal_c1p$mu_ktm1_star
      Sigma_ktm1_star <- list_proposal_c1p$Sigma_ktm1_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_cpm1_star <- list_proposal_c1p$Sigma_star
      
      list_proposal_ktp <- findProposalC2P_KT(d, E, m_terms)
      mu_ktm1_0_star <- list_proposal_ktp$mu_ktm1_star
      Sigma_ktm1_0_star <- list_proposal_ktp$Sigma_ktm1_star
      mu_c2pm1_star <- list_proposal_ktp$mu_c2pm1_star
      Sigma_c2pm1_star <- list_proposal_ktp$Sigma_c2pm1_star
      
      ktm1_proposed <- mvrnorm(n = 1, mu_ktm1_star, Sigma_ktm1_star)
      kt_proposed <- c(ktm1_proposed, - sum(ktm1_proposed))
      cpm1_proposed <- mvrnorm(n = 1, cpm1_star, Sigma_cpm1_star)
      c1p_proposed <- c(cpm1_proposed, P - sum(cpm1_proposed))
      
      loglik_current <- loglik_ktp(ktp_current, d, E, m_terms)
      loglik_proposed <- loglik_ktcp(kt_proposed, c1p_proposed, d, E, m_terms)
      # loglik_proposed <- loglik_ktp(param_k$yip, d, E, m_terms)
      
      logproposal_current <- 
        dmvnorm(kt_current[1:(Y-1)], mu_ktm1_0_star, Sigma_ktm1_0_star, log = T) + 
        sum(
          sapply(1:Y, function(t){
            dmvnorm(c2p_current[t,1:(P-1)], mu_c2pm1_star[t,], Sigma_c2pm1_star[t,,], log = T)
          })
        )
      
      logproposal_proposed <- 
        dmvnorm(ktm1_proposed, mu_ktm1_star, Sigma_ktm1_star, T) + 
        dmvnorm(cpm1_proposed, cpm1_star, Sigma_cpm1_star, T)
      
      # kt_proposed <- sapply(1:Y, function(t){
      #   # rt2(1, mu_kt_star[t], sd_kt_star[t], df = 3)
      #   rnorm(1, mu_kt_star[t], sd_kt_star[t], df = 3)
      # })
      # cpm1_proposed <- mrt2(cpm1_star, Sigma_c1p_star, df = 3)
      # c1p_proposed <- c(cpm1_proposed, P - sum(cpm1_proposed))
      # 
      # loglik_proposed <- loglik_ktcp(kt_proposed, c1p_proposed, d, E, m_terms)
      # loglik_current <- loglik_ktp(ktp_current, d, E, m_terms)
      # 
      # logproposal_proposed <- 
      #   sum(sapply(1:Y, function(t){
      #     log(dt2(kt_proposed[t], mu_kt_star[t], sd_kt_star[t], df = 3))
      #   })) + dmt_cpp(cpm1_proposed, nu = 3, 
      #                 cpm1_star, Sigma_c1p_star, returnLog = T)
      
      logproposal_2 <- log(1) # going from model 2 to model 1
      logproposal_1 <- log(.5) # going from model 1 to model 2
      
      mh_ratio <- exp(loglik_proposed - loglik_current +
                        logproposal_2 - logproposal_1 +
                        logproposal_current - logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        
        param_k$y <- kt_proposed
        param_k$c1p <- c1p_proposed
        param_k$yip <- 
          matrix(kt_proposed, Y, P, byrow = F) * matrix(c1p_proposed, Y, P, byrow = T)
        gamma_k <- 2
        
      }
      
    }
    
  }
  
  # propose gamma_g
  # print("Update gamma_g")
  if (update_gammag) {
    
    m_terms <- computeTerms(3, param_a, param_k, param_g, dp)
    
    if(gamma_g == 1){
      
      gc_current <- param_g$y
      
      # list_proposal_gtx <- findProposalGTX(d, E, m_terms)
      # mu_gtm2_0_star <- list_proposal_gtx$gtm2_star
      # Sigma_gtm2_0_star <- list_proposal_gtx$Sigma_gtm2_star
      
      list_proposal_c1p <- findProposalC1P_GC(d, E, m_terms)
      mu_gtm2_star <- list_proposal_c1p$mu_gtm2_star
      Sigma_gtm2_star <- list_proposal_c1p$Sigma_gtm2_star
      cpm1_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_star <- list_proposal_c1p$Sigma_star
      
      mu_gtm2_0_star <- mu_gtm2_star
      Sigma_gtm2_0_star <- Sigma_gtm2_star
      
      gcm2_proposed <- mvrnorm(n = 1, mu_gtm2_star, Sigma_gtm2_star)
      gc_proposed <- c(0, gcm2_proposed, - sum(gcm2_proposed))
      cpm1_proposed <- mvrnorm(n = 1, cpm1_star, Sigma_c1p_star)
      c1p_proposed <- c(cpm1_proposed, P - sum(cpm1_proposed))
      c2p_proposed <- matrix(c1p_proposed, X + Y - 1, P, byrow = T)
      
      loglik_current <- - loglik_gccp(gc_current, 
                                    matrix(1, X + Y - 1, P), d, E, m_terms)
      loglik_proposed <- - loglik_gccp(gc_proposed, 
                                     c2p_proposed, d, E, m_terms)
      
      logproposal_current <- 
        dmvnorm(gc_current[-c(1,X+Y-1)], mu_gtm2_0_star, Sigma_gtm2_0_star, T)
      
      logproposal_proposed <- 
        dmvnorm(gcm2_proposed, mu_gtm2_star, Sigma_gtm2_star, T) + 
        dmvnorm(cpm1_proposed, cpm1_star, Sigma_c1p_star, T)
      
      logproposal_2 <- log(1) # from 2 to 1
      logproposal_1 <- log(1) # from 1 to 2
      
      mh_ratio <- exp(loglik_proposed - loglik_current +
                        logproposal_2 - logproposal_1 - 
                        logproposal_current - logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        
        param_g$y <- gc_proposed
        param_g$c1p <- c1p_proposed
        param_g$yip <- 
          matrix(gc_proposed, X + Y - 1, P, byrow = F) * 
          c2p_proposed
        
        gamma_g <- 2
        
      }
      
    } else if (gamma_g == 2){
      
      gc_current <- param_g$y
      c1p_current <- param_g$c1p
      
      list_proposal_c1p <- findProposalC1P_GC(d, E, m_terms)
      mu_gtm2_0_star <- list_proposal_c1p$mu_gtm2_star
      Sigma_gtm2_0_star <- list_proposal_c1p$Sigma_gtm2_star
      cpm1_0_star <- list_proposal_c1p$cpm1_star
      Sigma_c1p_0_star <- list_proposal_c1p$Sigma_star
      
      logproposal_current <- 
        dmvnorm(gc_current[-c(1,X+Y-1)], mu_gtm2_0_star, Sigma_gtm2_0_star, T) + 
        dmvnorm(c1p_current[1:(P-1)], cpm1_0_star, Sigma_c1p_0_star, T)
        
      c2p_current <- matrix(c1p_current, X + Y - 1, P, byrow = T)
      
      loglik_current <- - loglik_gccp(gc_current, c2p_current, d, E, m_terms)
      
      # list_proposal_gtx <- findProposalGTX(d, E, m_terms)
      mu_gtm2_star <- mu_gtm2_0_star
      Sigma_gtm2_star <- Sigma_gtm2_0_star
      
      gc_m2_proposed <- 
        mvrnorm(n = 1, mu_gtm2_star, Sigma_gtm2_star)
      gc_proposed <- c(0, gc_m2_proposed, - sum(gc_m2_proposed))
      
      loglik_proposed <- - loglik_gccp(gc_proposed, 
                                     matrix(1, X + Y - 1, P), 
                                     d, E, m_terms)
      
      logproposal_proposed <- 
        dmvnorm(gc_m2_proposed, mu_gtm2_star, Sigma_gtm2_star, log = T)
      
      logproposal_2 <- log(1) # from 2 to 1
      logproposal_1 <- log(1)
      
      # loglik_gtxm2_current <- loglik_gcm2(d, E, m_terms)
      # 
      # loglik_gtxm2_current(gc_true[-c(1,X+Y-1)])
      # loglik_gtxm2_current(mu_gtm2_star)
      
      mh_ratio <- exp(loglik_proposed - loglik_current +
                        logproposal_2 - logproposal_1 + 
                        logproposal_current - logproposal_proposed)
      
      if(runif(1) < mh_ratio){
        
        param_g$y <- gc_proposed
        param_g$c1p <- NULL
        param_g$yip <- matrix(gc_proposed, X + Y - 1, P, F)
        
        gamma_g <- 1
        
      }
      
    } 
    
  }
  
  # output
  
  if(iter > nburn){
    currentIter <- iter - nburn
    
    a_output[currentIter,,] <- param_a$yip
    k_output[currentIter,,] <- param_k$yip
    g_output[currentIter,,] <- param_g$yip
    d_output[currentIter,] <- dp
    # ap_output[currentIter,,] <- ap
    # bp_output[currentIter,,] <- bp
    # kp_output[currentIter,,] <- kp
    # gamma_a_output[currentIter] <- which(a_configs == gamma_a)
    # gamma_k_output[currentIter] <- which(b_configs == gamma_b)
    # gamma_g_output[currentIter] <- which(k_configs == gamma_k)
    
    # loglik
    {
      # a2 <- matrix(a, X, P, byrow = F) + ap
      # b2 <- matrix(b, X, P, byrow = F) + bp
      # k2 <- matrix(k, Y, P, byrow = F) + kp
      
      # loglik_output[currentIter] <-  loglik_LCp(a2, b2, k2, d, E)
    }
    
    m_current <- computeTerms(0, param_a, param_k, param_g, dp)
    
    m_output[,,,iter - nburn] <- m_current
  }
  
}

beep()

qplot(1:niter, m_output[5,2,3,]) #+ geom_hline(aes(yintercept = m_true[1,1,1]))
qplot(1:niter, a_output[,1,1])
qplot(1:niter, g_output[,2,3])

# TRUE DATA COMPARISON -----

data_plot <- expand.grid(
  x = 1:X,
  t = 1:Y,
  p = 1:P
)

data_plot_CI <- t(apply(data_plot, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  p <- dat[3]
  quantile(m_output[x,t,p,], probs = c(0.025, 0.975))
}))

data_plot$CI1 <- data_plot_CI[,1]
data_plot$CI2 <- data_plot_CI[,2]

data_plot$CrudeRates <- apply(data_plot, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  p <- dat[3]
  log(d[x,t,p] / E[x,t,p] + .0000000001)
})

idxBreaks <- floor(seq(1, X, by = 5))

data_plot %>%
  filter(
    x %in% c(1:10),
         t %in% c(1, 5, 10, 15)
         # p %in% c(1,2)
         ) %>%
ggplot(aes(x = x,
           y = CrudeRates,
           color = factor(t),
           group = factor(t),
           ymin = CI1,
           ymax = CI2)) + 
  geom_errorbar() +
  scale_x_continuous(breaks = idxBreaks,
                     labels = dimnames(d)[[1]][idxBreaks]) + 
  # coord_cartesian(ylim = c(0, - 10)) +
  geom_point() + 
  geom_line() + 
  facet_grid(vars(p), scales = "free")

data_plot %>%
  filter(
    x %in% c(1,5, 10)#,
         # t %in% c(1, 5, 10, 15)
         # p %in% c(1,2)
         ) %>%
ggplot(aes(x = t,
           y = CrudeRates,
           color = factor(x),
           group = factor(x),
           ymin = CI1,
           ymax = CI2)) + 
  geom_errorbar() +
  scale_x_continuous(breaks = idxBreaks,
                     labels = dimnames(d)[[1]][idxBreaks]) + 
  # coord_cartesian(ylim = c(0, - 10)) +
  geom_point() + 
  geom_line() + 
  facet_grid(vars(p), scales = "free")


# SIMULATED DATA COMPARISON -----

data_plot <- expand.grid(
  x = 1:X,
  t = 1:Y,
  p = 1:P
)

data_plot_CI <- t(apply(data_plot, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  p <- dat[3]
  quantile(m_output[x,t,p,], probs = c(0.025, 0.975))
}))

data_plot$CI1 <- data_plot_CI[,1]
data_plot$CI2 <- data_plot_CI[,2]

data_plot$True <- apply(data_plot, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  p <- dat[3]
  log(m_true[x,t,p])
})

data_plot %>% 
  filter(
    # x %in% c(1, 3, 6),
         t %in% c(1, 4)
         # p == 1
         ) %>%
ggplot(aes(x = x,
           y = True,
           color = factor(t),
           group = factor(t),
           ymin = CI1,
           ymax = CI2)) + 
  geom_errorbar() +
  geom_point() + 
  facet_grid(vars(p))

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
