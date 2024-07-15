
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

sourceCpp("Code/8Models/code.cpp")
source("Code/8Models/functions.R")

# SIMULATED DATA ----------

realData <- F

X <- 20 # 30 # ages
Y <- 10 # 36 #  years

ages <- 50 + 1:X
years <- 2000 + 1:Y

E <- array(rpois(X * Y, lambda = 5000), dim = c(X, Y))

c2x <- (ages - mean(ages)) / (- ages[1] + mean(ages))

# model choices
{
  choice1 <- 1 # 1 unconstrained, 2 line
  choice2 <- 2 # 1 not additive , 2 additive
  choice3 <- 3 # 1 no effect, 2 varying line, 3 lee carter
  choice4 <- 1 # 1 no cohort effect, 2 cohort effect
}

list_d <- simulateData(X, Y, E, c2x,
                       choice1, choice2, choice3, choice4)
d <- list_d$d
term1 <- list_d$term1
term2 <- list_d$term2
term3 <- list_d$term3
term4 <- list_d$term4

# adding missing values
{
  # idx <- sample(1:(X*Y), replace = T, size = 100)
  # d[idx] <- NA
  # E[idx] <- NA
  
}

# REAL DATA SEX ------

realData <- T

load(here("Data","data_sex_UK.rda"))

d <- apply(d, c(1,2), sum)
E <- apply(E, c(1,2), sum)

X_subset <- which(as.numeric(dimnames(d)[[1]]) > 59 &
                    as.numeric(dimnames(d)[[1]]) < 90)
Y_subset <- which(as.numeric(dimnames(d)[[2]]) > 1963 &
                    as.numeric(dimnames(d)[[2]]) < 2004)

d <- d[X_subset, Y_subset]
E <- E[X_subset, Y_subset]

X <- dim(d)[1]
Y <- dim(d)[2]

ages <- as.numeric(dimnames(d)[[1]])
years <- as.numeric(dimnames(d)[[2]])

# MCMC -----

nburn <- 3000
niter <- 3000

trueStartingVal <- F

# starting values
{
  
  c2x <- (ages - mean(ages)) / (- ages[1] + mean(ages))
  
  # true
  if(trueStartingVal) {
    
    term1 <- term1
    param1 <- list_d$param1
    term2 <- term2
    param2 <- list_d$param2
    term3 <- term3
    param3 <- list_d$param3
    term4 <- term4
    param4 <- list_d$param4
    
    delta1 <- choice1
    delta2 <- choice2
    delta3 <- choice3
    delta4 <- choice4
    
    if(delta1 == 1){
      
      ax <- param1
      
    } else if(delta1 == 2){
      
      ab <- param1
      
    }
    
    if(delta2 == 2){
      k1t <- param2
    } 
    
    if(delta3 == 2){
      
      k2t <- param3$k2t
      
    } else if(delta3 == 3){
      
      k2t <- param3$k2t
      bx <- param3$bx
      
    }
    
    if(delta4 == 2){
      
      gtx <- param4
      
    }
    
  } else { # from MLE
    
    delta1 <- 1
    delta2 <- 2
    delta3 <- 3
    delta4 <- 1
    
    if(delta1 == 1){
      
      ax <- apply(d / E, 1, function(x){
        mean(log(x + .001), na.rm = T)
      })
      
      term1 <- ax
      
    } else {
      
      ab <- c(
        log(sum(d) / sum(E)), 0
      )
      
      term1 <- ab[1] + ab[2] * c2x
      
    }
    
    term1 <- matrix(term1, X, Y, byrow = F)
    
    # k1t <- rep(0, Y)
    # k2t <- rep(0, Y)
    # bx <- c2x
    
    term2 <- matrix(0, X, Y)
    
    if(delta2 == 2){
      k1t <- rep(0, Y)
    }
    
    if(delta3 == 3){
      if(delta2 == 1){
        bx <- seq(-1, 0, length.out = X)
      } else if (delta2 == 2){
        bx <- c2x  
      }
      
      k2t <- rnorm(Y)
    }
    
    if(delta4 == 2){
      gtx <- rep(0, X + Y - 1)
    }
    
    term3 <- matrix(0, X, Y)
    term4 <- matrix(0, X, Y)
    
  }
  
}

# prior
{
  sd_b <- 2
  sd_k <- 2
  sd_g <- 2
  
  sd_bx <- 4
  
  p_k1t <- .1
  p_gtx <- .1
}

# output 
{
  term1_output <- array(NA, dim = c(niter, X, Y))
  term2_output <- array(NA, dim = c(niter, X, Y))
  term3_output <- array(NA, dim = c(niter, X, Y))
  term4_output <- array(NA, dim = c(niter, X, Y))
  
  delta1_output <- rep(NA, niter)
  delta2_output <- rep(NA, niter)
  delta3_output <- rep(NA, niter)
  delta4_output <- rep(NA, niter)
  
  loglik_output <- rep(NA, niter)
  
  ax_output <- matrix(NA, niter, X)
  ab_output <- matrix(NA,niter, 2)
  k1t_output <- matrix(NA,niter, Y)
  k2t_output <- matrix(NA,niter, Y)
  bx_output <- matrix(NA,niter, X)
  
}

# paramsUpdate
{
  update_term1 <- T
  update_term2 <- T
  update_term3 <- T
  update_term4 <- F
  
  update_delta1 <- F
  update_delta2 <- F
  update_delta3 <- F
  update_delta4 <- F
}

for (iter in 1:(niter + nburn)) {
  
  if(delta3 == 3) print(bx[1:5])
  # if(delta4 == 2) print(gtx[1:5])
  
  if(iter %% 10 == 0){
    if(iter < nburn){
      print(paste0("Burn-in Iteration = ",iter))  
    } else {
      print(paste0("Iteration = ",iter - nburn))  
    }  
    
    print(paste0(
      "delta1 = ",delta1,
      " / delta2 = ",delta2,
      " / delta3 = ",delta3,
      " / delta4 = ",delta4
    ))
    
  }
  
  # update term1 ------
  
  if(update_term1){
    
    if(delta1 == 1){
      param1 <- ax
    } else {
      param1 <- ab
    }
    
    list_term1 <- updateTerm1(param, delta1, term1, term2, term3, term4)
    
    if(delta1 == 1){
      ax <- list_term1$param
    } else {
      ab <- list_term1$param
    }
    term1 <- list_term1$term1
    
  }
  
  # update term2 ----
  
  if(update_term2){
    
    if(delta2 == 2){ 
      
      term134 <- computeTerms(idx = 2, term1, term2, term3, term4)
      
      k1t <- updateK1T(k1t, d, E, term134)
      term2 <- matrix(k1t, X, Y, byrow = T)
      
    } 
    
  }
  
  # update term3 ----
  
  if(update_term3){
    
    term124 <- computeTerms(idx = 3, term1, term2, term3, term4)
    
    if(delta3 == 2){ # varying line 
      
      list_k2t <- update_k2t(k2t, d, E, term124, c2x)
      k2t <- list_k2t$k2t
      term3 <- list_k2t$term3
      
    } else if(delta3 == 3){ # lee carter
      
      # list_k2t <- update_k2t(k2t, d, E, term124, bx)
      # k2t <- list_k2t$k2t
      # term3 <- list_k2t$term3
      # 
      # list_bx <- update_bx(bx, d, E, term124, k2t)
      # bx <- list_bx$bx
      # term3 <- list_bx$term3
      
      list_bxk2t <- update_bxk2t_twosteps(bx, k2t, d, E, term124, delta2, sd_bx)
      # list_bxk2t <- update_bxk2t(bx, k2t, d, E, term124)
      bx <- list_bxk2t$bx
      k2t <- list_bxk2t$k2t
      term3 <- list_bxk2t$term3
      
    }
    
  }
  
  # update term4 ----
  
  if(update_term4){
    
    term123 <- computeTerms(idx = 4, term1, term2, term3, term4)
    
    if(delta4 == 2){ # cohort effect present
      
      if(delta2 == 1){ 
        
        list_gtx <- updateGTX(gtx, d, E, term123)
        
      } else if(delta2 == 2){
        # if the additive effect on kt is present we 
        # need an additional constraint
        list_gtx <- updateGTX_kt(gtx, d, E, term123)
        
      }
      
      gtx <- list_gtx$gtx
      term4 <- list_gtx$term4
      
    } 
    
  }
  
  # update delta1 ----
  
  if(update_delta1){
    
    term234 <- computeTerms(idx = 1, term1, term2, term3, term4)
    cxt <- term234
    
    # find proposal for line model
    list_proposal_abx <- findProposalABX(c2x, d, E, cxt)
    ab_star <- list_proposal_abx$ab_star
    Sigma_ab_star <- list_proposal_abx$Sigma_star
    
    # find proposal for unconstrained model
    list_ax_proposal <- buildProposalAX(d, E, cxt)
    ax_star_all <- list_ax_proposal$ax_star_all
    sd_ax_star <- list_ax_proposal$sd_ax_star
    
    if(delta1 == 1){ # unconstrained 
      
      # propose new ab
      ab_proposal <- mvrnorm(1, ab_star, Sigma_ab_star)
      
      loglik_proposal <- - loglik_term1_ab_fun(ab_proposal, c2x, d, E, cxt)
      loglik_current <- - loglik_term1_ax_fun(ax, d, E, cxt)
      
      logproposal_current <- sum(dnorm(ax, 
                                       ax_star_all, 
                                       sd_ax_star,
                                       log = T))
      logproposal_proposal <- 
        dmvnorm(ab_proposal, ab_star, Sigma_ab_star, log = T)
      
      # logprior_param_current <- 
      #   sum(dnorm(ax, sd = sd_b, log = T))
      # 
      # logprior_param_star <- 
      #   sum(dnorm(ax, sd = sd_b, log = T))
      
      mh_ratio <- exp(
        loglik_proposal - loglik_current +
          logproposal_current - logproposal_proposal
      )
      
      if(runif(1) < mh_ratio){
        
        delta1 <- 2
        ab <- ab_proposal
        term1 <- matrix(ab[1] + ab[2] * c2x, X, Y, byrow = F)
        ax <- NULL
        
      }
      
    } else if(delta1 == 2){ # line
      
      ax_proposed <- rnorm(X, ax_star_all, 
                           sd_ax_star)
      
      loglik_proposal <- - loglik_term1_ax_fun(ax_proposed, d, E, cxt)
      loglik_current <- - loglik_term1_ab_fun(ab, c2x, d, E, cxt)
      
      logproposal_proposal <- sum(dnorm(ax_proposed, 
                                        ax_star_all, 
                                        sd_ax_star,
                                        log = T))
      
      logproposal_current <- dmvnorm(ab, 
                                     ab_star, 
                                     Sigma_ab_star,
                                     log = T)
      
      mh_ratio <- exp(
        loglik_proposal - loglik_current +
          logproposal_current - logproposal_proposal
      )
      
      if(runif(1) < mh_ratio){
        
        delta1 <- 1
        ax <- ax_proposed
        term1 <- matrix(ax, X, Y, byrow = F)
        ab <- NULL
        
      }
      
    }
    
  }
  
  # update delta2 ----
  
  if(update_delta2){
    
    term134 <- computeTerms(idx = 2, term1, term2, term3, term4)
    cxt <- term134
    
    if(delta3 == 1){ 
      # no lee carter effect present, in this case
      # we can simply propose to add or remove k1t
      
      list_proposal <- buildProposalK1t_m1(d, E, cxt)
      km1_star <- list_proposal$km1_star
      Sigma_km1_star <- list_proposal$Sigma_star
      
      if(delta2 == 1){ # no effect
        
        # k1t_proposed_pm1 <- rnorm(Y - 1, k_star[1:(Y-1)], 
        #                           sd = sqrt(Sigma_star[1:(Y-1),1:(Y-1)]))
        # k1t_proposed <- c(k1t_proposed_pm1, - sum(k1t_proposed_pm1))
        # 
        # k1t_proposed_2 <- sampleMVNconstraint_k(k_star, Sigma_star)
        
        k1t_proposed <- mvrnorm(n = 1,
                                km1_star,
                                Sigma_km1_star)
        k1t_proposed <- c(k1t_proposed, - sum(k1t_proposed))
        
      } else {
        
        k1t_proposed <- k1t
        
      }
      
      logprior_param_star <- sum(dnorm(k1t_proposed, sd = sd_k, log = T))
      
      k1t_current <- rep(0, Y)
      
      loglik_proposal <- loglik_term1_k(k1t_proposed, d, E, cxt)
      loglik_current <- loglik_term1_k(k1t_current, d, E, cxt)
      
      logprior_star <- log(p_k1t)
      logprior_current <- log(1 - p_k1t)
      
      logProposal_proposal <- 
        dmvnorm(
          k1t_proposed[1:(Y-1)], 
          km1_star[1:(Y-1)], 
          Sigma_km1_star, 
          log = T)
      
      mh_ratio <- exp(
        loglik_proposal - loglik_current + 
          logprior_param_star +
          logprior_star - logprior_current - 
          logProposal_proposal
      )
      
      if(delta2 == 2){ 
        
        mh_ratio <- 1 / mh_ratio
        
      }
      
      if(runif(1) < mh_ratio){ 
        
        if(delta2 == 1){ # additive effect not present and adding it
          
          k1t <- k1t_proposed
          delta2 <- 2
          
          term2 <- matrix(k1t, X, Y, byrow = F)
          
        } else { # additive effect present and removing it
          
          k1t <- NULL
          delta2 <- 1
          
          term2 <- matrix(0, X, Y)
          
        }
        
        
      }
      
      
    } else if(delta3 == 3){
      # lee carter effect present, therefore removing k1t 
      # now needs to take into account that bx
      # doesn't have the constraint sum(bx) = 0
      
      # compute base model with k1t removed
      
      if(delta2 == 1) { # additive effect not present
        
        list_proposal <- buildProposalK1t_m1(d, E, cxt)
        km1_star <- list_proposal$km1_star
        Sigma_km1_star <- list_proposal$Sigma_star
        
        k1t_proposed <- mvrnorm(n = 1,
                                km1_star,
                                Sigma_km1_star)
        k1t_proposed <- c(k1t_proposed, - sum(k1t_proposed))
        
        k1t_current <- rep(0, Y)
        
        loglik_proposal <- loglik_term1_k(k1t_proposed, d, E, cxt)
        loglik_current <- loglik_term1_k(k1t_current, d, E, cxt)
        
        logProposal_proposal <- 
          dmvnorm(
            k1t_proposed[1:(Y-1)], 
            km1_star[1:(Y-1)], 
            Sigma_km1_star, 
            log = T)
        
        
      } else if (delta2 == 2){
        
        # approximate with base model
        param1 <- list("ax"= ax)
        param2 <- list("k1t" = k1t)
        param3 <- list("k2t" = k2t, "bx" = bx)
        param4 <- NULL
        list_params <- mapModels(delta1, delta2, delta3, delta4,
                                 delta1_new = 1, delta2_new = 1, 
                                 delta3_new = 3, delta4_new = 1,
                                 param1, param2, param3, param4)
        k2t_base1 <- list_params$k2t
        bx_base1 <- list_params$bx
        
        # reconvert the approximated model to a model with k1t, k2t, bx
        param1 <- list("ax"= ax)
        param3 <- list("k2t" = k2t_base1, "bx" = bx_base1)
        param4 <- NULL
        list_params <- mapModels(delta1, delta2 = 1, delta3, delta4,
                                 delta1_new = 1, delta2_new = 2, 
                                 delta3_new = 3, delta4_new = 1,
                                 param1, param2, param3, param4)
        k1t_base <- list_params$k1t
        k2t_base <- list_params$k2t
        bx_base <- list_params$bx
        
        term2_base <- matrix(k1t_base, X, Y, byrow = T)
        term3_base <- 
          matrix(k2t_base, X, Y, byrow = T) * 
          matrix(bx_base, X, Y, byrow = F)
        cxt_base <- computeTerms(idx = 0, term1, term2_base, 
                                 term3_base, term4)
        
        list_proposal <- buildProposalK1t_m1(d, E, cxt_base)
        km1_star <- list_proposal$km1_star
        Sigma_km1_star <- list_proposal$Sigma_star
        
        k1t_proposed <- k1t - k1t_base
        
        loglik_proposal <- loglik_term1_k(k1t_proposed, d, E, cxt_base)
        loglik_current <- loglik_term1_k(rep(0, Y), d, E, cxt_base)
        
        logProposal_proposal <- 
          dmvnorm(
            k1t_proposed[1:(Y-1)], 
            km1_star[1:(Y-1)], 
            Sigma_km1_star, 
            log = T)
        
      }
      
      logprior_param_star <- sum(dnorm(k1t_proposed, sd = sd_k, log = T))
      
      logprior_star <- log(p_k1t)
      logprior_current <- log(1 - p_k1t)
      
      mh_ratio <- exp(
        loglik_proposal - loglik_current + 
          logprior_param_star +
          logprior_star - logprior_current - 
          logProposal_proposal
      )
      
      if(delta2 == 2){
        
        mh_ratio <- 1 / mh_ratio
        
      }
      
      if(runif(1) < mh_ratio){
        
        if(delta2 == 1){
          # adding the additive effect to lee carter
          # therefore reparametrising to add the constraint sum(b_x) = 0
          
          delta2 <- 2
          
          list_params <- reparamK1tK2txBx(k1t_proposed, k2t, bx)
          k1t <- list_params$k1t_tilde
          k2t <- list_params$k2t_tilde
          bx <- list_params$bx_tilde
          
          term2 <- matrix(k1t, X, Y, byrow = T)
          
          term3 <- matrix(k2t, X, Y, byrow = T) * 
            matrix(bx, X, Y, byrow = T)
          
        } else if (delta2 == 2){
          # removing the effect
          
          delta2 <- 1
          
          k2t <- k2t_base1
          bx <- bx_base1
          
          k1t <- NULL
          
          term2 <- matrix(0, X, Y, byrow = T)
          
          term3 <- matrix(k2t, X, Y, byrow = T) * 
            matrix(bx, X, Y, byrow = T)
          
        }
        
      }
      
    }
    
  }
  
  # update delta3 ----
  
  if(update_delta3){
    
    term124 <- computeTerms(idx = 3, term1, term2, term3, term4)
    cxt <- term124
    
    if(delta3 == 1) { # no effect
      
      k2t_current <- rep(0, Y)
      bx_current <- rep(0, X)
      
      delta3_proposed <- 2
      
      logproposal_current <- 0
      
      logmove_current <- log(.5)
      
      logprior_param_current <- 0
      
    } else if(delta3 == 2) { # varying line
      
      k2t_current <- k2t
      bx_current <- c2x
      
      delta3_proposed <- sample(c(1,3), 1)
      
      list_proposal_k2t <- findProposalk2t_m1(c2x, d, E, cxt)
      k2t_m1_star <- list_proposal_k2t$k2t_star
      Sigma_k2t_m1_star <- list_proposal_k2t$Sigma_star
      
      logprior_param_current <- 
        sum(dnorm(k2t[1:(Y-1)], sd = sd_k, log = T))
      
      logproposal_current <- 
        dmvnorm(k2t[1:(Y-1)], 
                k2t_m1_star, 
                Sigma_k2t_m1_star, 
                log = T)
      
      logmove_current <- log(1)
      
    } else if(delta3 == 3) { # lee carter
      
      k2t_current <- k2t
      bx_current <- bx
      
      delta3_proposed <- 2
      
      # find proposal for lee carter
      {
        # list_proposal_bxk2t <- findProposalBXk2t(d, E, cxt, c2x)
        # bxl2t_star <- list_proposal_bxk2t$bxl2t_star
        # Sigma_bxl2t_star <- list_proposal_bxk2t$Sigma_star
        # 
        # k2t_bx_m_current <- c(bx[-c(1,X)], k2t[-Y])
        # 
        # logproposal_current <- 
        #   dmvnorm(k2t_bx_m_current,
        #           bxl2t_star,
        #           Sigma_bxl2t_star, log = T)
        # 
      }
      
      # find proposal for lee carter in two steps
      {
        if(delta2 == 1){
          bx_start <- seq(-1, 0, length.out = X)
        } else if(delta2 == 2){
          bx_start <- c2x  
        }
        
        list_proposal_kt <- findProposalk2t_given_bx(k2t, bx_start, d, E, cxt)
        kt_star <- list_proposal_kt$bxl2t_star
        Sigma_kt_star <- list_proposal_kt$Sigma_star
        
        if(delta2 == 1){
          list_proposal_bx <- findProposalBX_given_k2t_m1(bx_start, k2t, d, E, cxt, c2x)  
        } else if(delta2 == 2){
          list_proposal_bx <- findProposalBX_given_k2t_m2(bx_start, k2t, d, E, cxt, c2x)  
        }
        
        bx_star <- list_proposal_bx$bxl2t_star
        Sigma_bx_star <- list_proposal_bx$Sigma_star
        
        if(delta2 == 1){
          
          logproposal_current <- 
            dmvnorm(bx[-1], 
                    bx_star,
                    Sigma_bx_star,
                    log = T) +
            dmvnorm(k2t[-Y], 
                    kt_star,
                    Sigma_kt_star,
                    log = T) 
          
        } else if (delta2 == 2){
          
          logproposal_current <- 
            dmvnorm(bx[-c(1,X)], 
                    bx_star,
                    Sigma_bx_star,
                    log = T) +
            dmvnorm(k2t[-Y], 
                    kt_star,
                    Sigma_kt_star,
                    log = T) 
          
        }
        
      }
      
      logprior_param_current <- 
        sum(dnorm(k2t[1:(Y-1)], sd = sd_k, log = T)) +
        sum(dnorm(bx, sd = sd_b, log = T))
      
      logmove_current <- log(.5)
      
    }
    
    loglik_current <- - loglik_term3(k2t_current, bx_current, d, E, cxt)
    
    if(delta3_proposed == 1){
      
      k2t_proposed <- rep(0, Y)
      bx_proposed <- rep(0, X)
      
      logproposal_star <- 0
      
      logmove_proposal <- log(.5)
      
    } else if(delta3_proposed == 2){
      
      list_proposal_k2t <- findProposalk2t_m1(c2x, d, E, cxt)
      k2t_m1_star <- list_proposal_k2t$k2t_star
      Sigma_k2t_m1_star <- list_proposal_k2t$Sigma_star
      
      k2t_proposed_m1 <- mvrnorm(1,
                                 k2t_m1_star,
                                 Sigma_k2t_m1_star)
      k2t_proposed <- c(k2t_proposed_m1, - sum(k2t_proposed_m1))
      
      bx_proposed <- c2x
      
      logproposal_star <- 
        dmvnorm(k2t_proposed_m1, 
                k2t_m1_star, 
                Sigma_k2t_m1_star, 
                log = T)
      
      logprior_param_star <- 
        sum(dnorm(k2t_proposed_m1, sd = sd_k, log = T)) 
      
      logmove_proposal <- log(1)
      
    } else if(delta3_proposed == 3){
      
      # find proposal for lee carter
      {
        # list_proposal_bxk2t <- findProposalBXk2t(d, E, cxt, c2x)
        # bxl2t_star <- list_proposal_bxk2t$bxl2t_star
        # Sigma_bxl2t_star <- list_proposal_bxk2t$Sigma_star
        # 
        # k2tbx_proposed_m1 <- mvrnorm(n = 1,
        #                              bxl2t_star,
        #                              Sigma_bxl2t_star)
        # 
        # bx_m1X <- k2tbx_proposed_m1[1:(X-2)]
        # k2t_mY <- k2tbx_proposed_m1[X-2 + 1:(Y-1)]
        # 
        # bx_proposed <- c(-1, bx_m1X, 1 - sum(bx_m1X) + 1)
        # k2t_proposed <- c(k2t_mY, - sum(k2t_mY))
        # 
        # logproposal_star <- 
        #   dmvnorm(k2tbx_proposed_m1, 
        #           bxl2t_star,
        #           Sigma_bxl2t_star,
        #           log = T)  
      }
      
      # find proposal for lee carter in two steps
      {
        if(delta2 == 1){
          bx_start <- seq(-1, 0, length.out = X)
        } else if(delta2 == 2){
          bx_start <- c2x  
        }
        
        list_proposal_kt <- findProposalk2t_given_bx(k2t, bx_start, d, E, cxt)
        kt_star <- list_proposal_kt$bxl2t_star
        Sigma_kt_star <- list_proposal_kt$Sigma_star
        
        kt_proposed_m1 <- mvrnorm(n = 1,
                                  kt_star,
                                  Sigma_kt_star)
        
        k2t_proposed <- c(kt_proposed_m1, - sum(kt_proposed_m1))
        
        if (delta2 == 1){
          
          list_proposal_bx <- findProposalBX_given_k2t_m1(bx_start, k2t_proposed, 
                                                          d, E, cxt, bx_start)
          
        } else  if (delta2 == 2){
          
          list_proposal_bx <- findProposalBX_given_k2t_m2(bx_start, k2t_proposed, 
                                                          d, E, cxt, bx_start)
          
        }
        
        
        bx_star <- list_proposal_bx$bxl2t_star
        Sigma_bx_star <- list_proposal_bx$Sigma_star
        
        bx_proposed_m <- mvrnorm(n = 1,
                                 bx_star,
                                 Sigma_bx_star)
        
        if (delta2 == 1){
          
          bx_proposed <- c(-1, bx_proposed_m)
          
        } else  if (delta2 == 2){
          
          bx_proposed <- c(-1, bx_proposed_m, - sum(bx_proposed_m) + 1)
          
        }
        
        logproposal_star <- 
          dmvnorm(bx_proposed_m, 
                  bx_star,
                  Sigma_bx_star,
                  log = T) +
          dmvnorm(kt_proposed_m1, 
                  kt_star,
                  Sigma_kt_star,
                  log = T)  
      }
      
      logprior_param_star <- 
        sum(dnorm(k2t_proposed_m1, sd = sd_k, log = T)) +
        sum(dnorm(bx_proposed, sd = sd_b, log = T)) 
      
      logmove_proposal <- log(.5)
      
    }
    
    loglik_proposal <- - loglik_term3(k2t_proposed, 
                                      bx_proposed, 
                                      d, E, cxt)
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current + 
        logmove_current - logmove_proposal +
        logprior_param_star - logprior_param_current +
        logproposal_current - logproposal_star
    )
    
    if(runif(1) < mh_ratio){
      
      delta3 <- delta3_proposed
      
      if(delta3_proposed == 1){
        
        bx <- NULL
        k2t <- NULL
        
      } else if(delta3_proposed == 2){
        
        bx <- NULL
        k2t <- k2t_proposed
        
      } else {
        
        bx <- bx_proposed
        k2t <- k2t_proposed
        
      }
      
    }
    
    if(delta3 == 1){
      
      term3 <- matrix(0, X, Y)  
      
    } else if (delta3 == 1){
      
      term3 <- 
        matrix(c2x, X, Y, byrow = F) *   
        matrix(k2t, X, Y, byrow = T)
      
    } else if (delta3 == 1){
      
      term3 <- 
        matrix(bx, X, Y, byrow = F) *   
        matrix(k2t, X, Y, byrow = T)
      
    } 
    
  }
  
  # update delta4 ----
  
  if(update_delta4){
    
    term123 <- computeTerms(idx = 4, term1, term2, term3, term4)
    cxt <- term123
    
    if(delta2 == 1){
      
      # list_proposal <- buildProposalGtx(d, E, cxt)
      list_proposal <- buildProposalGtx_m1(d, E, cxt)
      gtx_star <- list_proposal$gtx_star
      Sigma_star <- list_proposal$Sigma_star
      
      idxElemProposed <- setdiff(1:(X+Y-1), X+Y-1)
      
    } else if(delta2 == 2){
      
      # list_proposal <- buildProposalGtx(d, E, cxt)
      list_proposal <- buildProposalGtx_m2(d, E, cxt)
      gtx_star <- list_proposal$gtx_star
      Sigma_star <- list_proposal$Sigma_star
      
      idxElemProposed <- setdiff(1:(X+Y-1), c(1, X+Y-1))
      
    }
    
    numElemProposed <- length(idxElemProposed)
    
    if(delta4 == 1){
      
      gtx_proposed_pm1 <- mvrnorm(n = 1, gtx_star, Sigma_star)
      
      if(delta2 == 1){
        
        gtx_proposed <- c(gtx_proposed_pm1, - sum(gtx_proposed_pm1))
        
      } else if(delta2 == 2){
        
        gtx_proposed <- c(0, gtx_proposed_pm1, 0)
        term4_proposed <- createTerm4(gtx_proposed, X, Y)
        gtx_proposed[length(gtx_proposed)] <- - sum(term4_proposed)
        
      }
      
    } else {
      
      gtx_proposed <- gtx
      
    }
    
    gtx_current <- rep(0, X + Y - 1)
    
    loglik_proposal <- loglik_term_gtx(gtx_proposed, d, E, cxt)
    loglik_current <- loglik_term_gtx(gtx_current, d, E, cxt)
    
    logprior_star <- log(p_gtx)
    logprior_current <- log(1 - p_gtx)
    
    logprior_param_star <- sum(dnorm(gtx_proposed, sd = sd_g, log = T))
    
    logProposal_proposal <- 
      dmvnorm(gtx_proposed[idxElemProposed], gtx_star, Sigma_star, log = T)
    
    # sum(
    #   dnorm(gtx_proposed[1:numElemProposed], gtx_star[1:numElemProposed], 
    #         sd = diag(sqrt(Sigma_star[1:numElemProposed,
    #                                   1:numElemProposed])), log = T)
    # )
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current + 
        logprior_param_star +
        logprior_star - logprior_current - 
        logProposal_proposal
    )
    
    if(delta4 == 2){ # 
      
      mh_ratio <- 1 / mh_ratio
      
    }
    
    if(delta4 == 1){ # additive effect not present
      
      if(runif(1) < mh_ratio){
        
        gtx <- gtx_proposed
        delta4 <- 2
        
      }
      
    } else { # additive effect not present
      
      if(runif(1) < mh_ratio){
        
        gtx <- NULL
        delta4 <- 1
        
      }
      
    }
    
    if(delta4 == 1){
      
      term4 <- matrix(0, X, Y)
      
    } else {
      
      term4 <- createTerm4(gtx, X, Y)
      
    }
    
    
    
  }
  
  # output -----
  
  if(iter > nburn){
    
    currentIter <- iter - nburn
    
    term1_output[currentIter,,] <- term1
    term2_output[currentIter,,] <- term2
    term3_output[currentIter,,] <- term3
    term4_output[currentIter,,] <- term4
    
    delta1_output[currentIter] <- delta1
    delta2_output[currentIter] <- delta2
    delta3_output[currentIter] <- delta3
    delta4_output[currentIter] <- delta4
    
    # loglik
    {
      # a2 <- matrix(a, X, P, byrow = F) + ap
      # b2 <- matrix(b, X, P, byrow = F) + bp
      # k2 <- matrix(k, Y, P, byrow = F) + kp
      # 
      loglik_output[currentIter] <-  loglik_m(d, E, term1, term2,
                                              term3, term4)
    }
    
    # additional params
    {
      if(delta1 == 1){
        ax_output[currentIter,] <- ax  
      } else if(delta1 == 1){
        ab_output[currentIter,] <- ab
      }
      
      if (delta2 == 2){
        k1t_output[currentIter,] <- k1t
      }
      
      if (delta3 == 2){
        k2t_output[currentIter,] <- k2t
      } else if (delta3 == 3){
        k2t_output[currentIter,] <- k2t
        bx_output[currentIter,] <- bx
      }
      
    }
    
  }
  
}

# beep()

# SAVE OUTPUT -----

# setwd(here("Results","Sex"))
# # setwd(here("Results","Product"))
# 
# save(a_output,
#      b_output,
#      k_output ,
#      ap_output,
#      bp_output,
#      kp_output,
#      delta_a_output,
#      delta_b_output,
#      delta_k_output,
#      loglik_output ,
#      file = "results_products.rda"
# )

# OUTPUT 
# DIAGNOSTICS ------

ess_term1 <- diagnosticsCheck(term1_output)
ess_term2 <- diagnosticsCheck(term2_output)
ess_term3 <- diagnosticsCheck(term3_output)
ess_term4 <- diagnosticsCheck(term4_output)

x <- 1
t <- 1
qplot(1:niter, term1_output[,x,t]) + 
  geom_hline(aes(yintercept = list_d$term1[x]))
qplot(1:niter, term2_output[,x,t]) + 
  geom_hline(aes(yintercept = list_d$term2[t]))
qplot(1:niter, term3_output[,x,t]) + 
  geom_hline(aes(yintercept = list_d$term3[x,t]))
qplot(1:niter, term4_output[,x,t]) + 
  geom_hline(aes(yintercept = list_d$term4[x,t]))

qplot(1:niter, term1_output[,x] + 
        term2_output[,t] + term4_output[,x,t]) + 
  geom_hline(aes(yintercept = 
                   list_d$term1[x] + 
                   list_d$term2[t] + 
                   list_d$term4[x,t]))

# FITTED CURVES VS TRUE CURVES ---------

# true curves
{
  curvesTrue <- plotCurves(list_d$term1, list_d$term2,
                           list_d$term3, list_d$term4)
}

# fitted curves CI
{
  
  mxt_CI <- matrix(NA, nrow(curvesTrue), 2)
  
  for (i in seq_len(nrow(curvesTrue))) {
    
    x_current <- curvesTrue$Age[i]
    t_current <- curvesTrue$Year[i]
    
    mxt_current <- 
      term1_output[,x_current,t_current] + 
      term2_output[,x_current,t_current] + 
      term3_output[,x_current,t_current] +
      term4_output[,x_current,t_current] 
    
    mxt_CI[i,] <- 
      quantile(mxt_current, probs = c(0.025, 0.975))
    
  }
  
  curvesTrue$CI_min <- mxt_CI[,1]
  curvesTrue$CI_max <- mxt_CI[,2]
  
}

curvesTrue %>% 
  filter(Year %in% 1:5) %>% 
  ggplot(aes(x = Age,
             y = m,
             ymin = CI_min,
             ymax = CI_max,
             color = Year,
             group = Year)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = .25) +
  theme_bw()

# CURRENT CURVES VS TRUE CURVES ---------

# true curves
{
  curvesTrue <- plotCurves(list_d$term1, list_d$term2,
                           list_d$term3, list_d$term4)
}

# fitted curves CI
{
  
  mxt_current <- matrix(NA, nrow(curvesTrue), 1)
  
  for (i in seq_len(nrow(curvesTrue))) {
    
    x_current <- curvesTrue$Age[i]
    t_current <- curvesTrue$Year[i]
    
    mxt_current[i] <- 
      term1[x_current,t_current] + 
      term2[x_current,t_current] + 
      term3[x_current,t_current] +
      term4[x_current,t_current] 
    
  }
  
  curvesTrue$current <- mxt_current
  
}

curvesTrue %>% 
  filter(Year %in% 1:5) %>% 
  ggplot(aes(x = Age,
             color = Year,
             group = Year)) + 
  geom_line(aes(y = current)) + 
  geom_line(aes(y = m), size = 1) +
  # geom_ribbon(alpha = .25) +
  theme_bw()

# TRANSFORM PARAMS ----

ax_LC_output <- ax_output
sum_bx_output <- apply(bx_output, 1, sum)
bx_LC_output <- t(sapply(1:niter, function(i) {
  bx_output[i,] / sum_bx_output[i]
}))
k2t_LC_output <- t(sapply(1:niter, function(i) {
  k2t_output[i,] * sum_bx_output[i]
}))

ax_mean <- apply(ax_LC_output, 2, mean)
bx_mean <- apply(bx_LC_output, 2, mean)
kt_mean <- apply(k2t_LC_output, 2, mean)

kt_true <- list_d$param3$k2t * sum(list_d$param3$bx)
bx_true <- list_d$param3$bx / sum(list_d$param3$bx)

cbind(ax_mean, list_d$param1)
cbind(bx_mean, bx_true)
cbind(kt_mean, kt_true)

save(ax_mean, bx_mean, kt_mean, file = "output_UKdata.rda")

abx_results <- data.frame(
  ax = ax_mean,
  bx = bx_mean
)
row.names(abx_results) <- ages

kt_results <- data.frame(
  kt = kt_mean
)
row.names(kt_results) <- years

write.csv(abx_results, file = "abx_results.csv")
write.csv(kt_results, file = "kt_results.csv")
