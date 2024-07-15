# PACKAGES --------

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

sourceCpp(here("Code","8Models - New/code.cpp"))
source(here("Code","8Models - New/functions.R"))

realData <- T

# SIMULATED DATA ----------

if(!realData){
  
  X <- 30 # 30 # ages
  Y <- 20 # 36 #  years
  
  ages <- 50 + 1:X
  years <- 2000 + 1:Y
  
  E <- array(rpois(X * Y, lambda = 5000), dim = c(X, Y))
  
  cx <- (ages - mean(ages)) / (- ages[1] + mean(ages))
  
  # model choices
  {
    choice1 <- 1 # 1 unconstrained, 2 line
    choice2 <- 1 # 1 0, 2 1, 3 k2t
    choice3 <- 1 # 1 cx, 2 bx
    choice4 <- 1 # 1 no cohort effect, 2 cohort effect
  }
  
  list_d <- simulateData(X, Y, E, cx,
                         choice1, choice2, choice3, choice4)
  d <- list_d$d
  
  # adding missing values
  {
    # idx <- sample(1:(X*Y), replace = T, size = 100)
    # d[idx] <- NA
    # E[idx] <- NA
    
  }
    
}

# REAL DATA SEX ------

if(F){
  
  load(here("Data","data_sex_UK.rda"))
  
  d <- apply(d, c(1,2), sum)
  E <- apply(E, c(1,2), sum)
  
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

# REAL DATA COUNTRY ------

country <- "France"

if(realData){
  
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

nburn <- 500
niter <- 500

trueStartingVal <- F

# starting values
{
  
  cx <- (ages - mean(ages)) / (- ages[1] + mean(ages))
  
  # true
  if(trueStartingVal) {
    
    delta <- c(choice1, choice2, choice3, choice4)
    
    param <- list("param1" = list_d$param1,
                  "param2" = list_d$param2,
                  "param3" = list_d$param3,
                  "param4" = list_d$param4)
    
  } else { # from MLE

    delta1 <- 2
    delta2 <- 1
    delta3 <- 1
    delta4 <- 1
    
    delta <- c(delta1, delta2, delta3, delta4)
    # delta <- c(choice1, choice2, choice3, choice4)
    
    ax <- apply(d / E, 1, function(x){
      mean(log(x + .001), na.rm = T)
    })
    
    param <- list()
    
    if(delta[1] == 1){
      
      param$param1 <- list("ax" = ax)
      
    } else {
      
      ab <- lm(ax ~ cx)
      ab <- as.numeric(ab$coefficients)
      
      param$param1 <- list("ab" = ab)
      
    }
    
    # term2
    
    param$param2 <- list("k1t" = rep(0,Y))
    
    # term 3
    
    param$param3 <- list()
    
    if(delta[2] == 2){
      
      if(delta[3] == 1){
        
        param$param3$b_bar <- 0
        
      } else if(delta[3] == 2){
        
        param$param3$bx <- cx
        
      }
      
      
    } else if (delta[2] == 3){
      
      param$param3$k2t <- rep(0, Y)
      
      if(delta[3] == 2){
        
        param$param3$bx <- cx
        
      }
      
    }
    
    # term 4
    
    if(delta[4] == 2){
      gtx <- rep(0, X + Y - 1)
      param$param4 <- list("gtx" = gtx)
    } else {
      param$param4 <- NULL
    }
     
  }
  
}

# prior
{
  sd_a <- 2
  sd_b <- 2
  sd_k <- 2
  sd_g <- 2
  
  sd_bx <- 4
  
  p_k1t <- .1
  p_gtx <- .1
  
  priorParams <- list("sd_a" = sd_a,
                      "sd_k" = sd_k,
                      "sd_b" = sd_b,
                      "sd_g" = sd_g)
}

# output 
{
  term1_output <- array(NA, dim = c(niter, X, Y))
  term2_output <- array(NA, dim = c(niter, X, Y))
  term3_output <- array(NA, dim = c(niter, X, Y))
  term4_output <- array(NA, dim = c(niter, X, Y))
  
  m_output <- array(NA, dim = c(niter, X, Y))
  
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
  gtx_output <- matrix(NA,niter, X+Y-1)
  
}

# paramsUpdate
{
  update_term1 <- T
  update_term2 <- T
  update_term3_1 <- T
  update_term3_2 <- T
  update_term4 <- T
  
  update_delta1 <- T
  update_delta2 <- T
  update_delta3 <- T
  update_delta4 <- T
}

for (iter in 1:(niter + nburn)) {

  # if(delta3 == 3) print(bx[1:5])
  # if(delta4 == 2) print(gtx[1:5])
  
  if(iter %% 2 == 0){
    if(iter < nburn){
      print(paste0("Burn-in Iteration = ",iter))  
    } else {
      print(paste0("Iteration = ",iter - nburn))  
    }  
    
    print(paste0(
      "delta1 = ",delta[1],
      " / delta2 = ",delta[2],
      " / delta3 = ",delta[3],
      " / delta4 = ",delta[4]
    ))
    
  }
  
  # update term1 ------
  
  if(update_term1){
    
    param$param1 <- updateTerm1(param, delta, d, E)
    
  }
  
  # update term2 ----
  
  if(update_term2){
    
    param$param2$k1t <- updateK1(param, delta, d, E)
    
  }
  
  # update term3_1 ----
  
  if(update_term3_1){
    
      param$param3$k2t <- updateK2(param, delta, d, E)
      
  }
  
  # update term3_2 ----
  
  if(update_term3_2){
    
    param$param3 <- updateBXB(param, delta, d, E)
    
  }
  
  # update term4 ----
  
  if(update_term4){
    
    param$param4 <- updateTerm4(param, delta, d, E)
    
  }
  
  # update delta1  -----
  
  if(update_delta1){
    
    list_params <- proposeNewState(idx_delta = 1, delta, 
                                   param, d, E)
    param <- list_params$param
    delta <- list_params$delta
    
  }
  
  # update delta2  -----
  
  if(update_delta2){
    
    list_params <- proposeNewState(idx_delta = 2, delta, 
                                   param, d, E)
    param <- list_params$param
    delta <- list_params$delta
    
  }
  
  # update delta3  -----
  
  if(update_delta3){
    
    if(delta[2] != 1){
      
      list_params <- proposeNewState(idx_delta = 3, delta, 
                                     param, d, E)
      param <- list_params$param
      delta <- list_params$delta
      
    }
    
  }
  
  # update delta4  -----

  if(update_delta4){
    
    list_params <- proposeNewState(idx_delta = 4, delta, 
                                   param, d, E)
    param <- list_params$param
    delta <- list_params$delta
    
  }
  
  # output -----
  
  if(iter > nburn){
    
    currentIter <- iter - nburn
    
    m_output[currentIter,,] <- computemxt(idx_delta = 0, delta, 
                      param$param1, param$param2,
                      param$param3, param$param4)
    
    # term1_output[currentIter,,] <- term1
    # term2_output[currentIter,,] <- term2
    # term3_output[currentIter,,] <- term3
    # term4_output[currentIter,,] <- term4
    
    delta1_output[currentIter] <- delta[1]
    delta2_output[currentIter] <- delta[2]
    delta3_output[currentIter] <- delta[3]
    delta4_output[currentIter] <- delta[4]
    
    # loglik
    {
      # a2 <- matrix(a, X, P, byrow = F) + ap
      # b2 <- matrix(b, X, P, byrow = F) + bp
      # k2 <- matrix(k, Y, P, byrow = F) + kp
      # 
      mxt <- computemxt(idx_delta = 0, delta,
                        param$param1, param$param2,
                        param$param3, param$param4)
      loglik_output[currentIter] <-  - loglik_term1(d, E, mxt,
                                                    matrix(0, X, Y))
    }
    
    # additional params
    {
      if(delta[1] == 1){
        ax_output[currentIter,] <- param$param1$ax
      } else if(delta[1] == 2){
        ab_output[currentIter,] <- param$param1$ab
      }
      # 
      # if (delta2 == 2){
      k1t_output[currentIter,] <- param$param2$k1t
      # }
      # 
      # if (delta3 == 2){
      #   k2t_output[currentIter,] <- k2t
      # } else if (delta3 == 3){
      #   k2t_output[currentIter,] <- k2t
      #   bx_output[currentIter,] <- bx
      # }
      
      if(delta[4] == 2){
        gtx_output[currentIter,] <- param$param4$gtx
      }
      
    }
   
  }
  
}

beep()

# FITTED CURVES VS TRUE CURVES ---------

# true curves
{
  mortalityCurves <- list_d$m %>% 
    as.data.frame() %>% 
    `colnames<-` (c(1:Y)) %>%
    rownames_to_column() %>%
    pivot_longer(-1, names_to = "Column", values_to = "Value") %>% 
    rename("Age" = "rowname", "Year" = "Column", "True" = "Value") %>% 
    mutate(Age = as.numeric(Age),
           Year = as.factor(as.numeric(Year)))
}

# fitted curves CI
{

  mxt_CI <- matrix(NA, nrow(mortalityCurves), 2)
  
  for (i in seq_len(nrow(mortalityCurves))) {
    
    x_current <- as.numeric(mortalityCurves$Age[i])
    t_current <- as.numeric(mortalityCurves$Year[i])
    
    mxt_CI[i,] <- 
      quantile(m_output[,x_current, t_current], probs = c(0.025, 0.975))
    
  }
  
  mortalityCurves$CI_min <- mxt_CI[,1]
  mortalityCurves$CI_max <- mxt_CI[,2]
  
}

mortalityCurves %>% 
  filter(Year %in% 1:10) %>% 
ggplot(aes(x = Age,
           y = True,
           ymin = CI_min,
           ymax = CI_max,
           color = Year,
           group = Year)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = .25) +
  theme_bw()


# FITTED CURVES VS CRUDE RATES ---------

# true curves
{
  mortalityCurves <- expand.grid(X = 1:X, Y = 1:Y)
  
  mortalityCurves$crudeRates <- sapply(1:nrow(mortalityCurves), function(i){
    
    idx_age <- which(1:X == mortalityCurves$X[i])
    idx_year <- which(1:Y == mortalityCurves$Y[i])
    log((d / E)[idx_age, idx_year])
    
  })
  
  # mortalityCurves <- list_d$m %>% 
  #   as.data.frame() %>% 
  #   `colnames<-` (c(1:Y)) %>%
  #   rownames_to_column() %>%
  #   pivot_longer(-1, names_to = "Column", values_to = "Value") %>% 
  #   rename("Age" = "rowname", "Year" = "Column", "True" = "Value") %>% 
  #   mutate(Age = as.numeric(Age),
  #          Year = as.factor(as.numeric(Year)))
}

# fitted curves CI
{

  mxt_CI <- matrix(NA, nrow(mortalityCurves), 2)
  
  for (i in seq_len(nrow(mortalityCurves))) {
    
    x_current <- as.numeric(mortalityCurves$X[i])
    t_current <- as.numeric(mortalityCurves$Y[i])
    
    mxt_CI[i,] <- 
      quantile(m_output[,x_current, t_current], probs = c(0.025, 0.975))
    
  }
  
  mortalityCurves$CI_min <- mxt_CI[,1]
  mortalityCurves$CI_max <- mxt_CI[,2]
  
}

mortalityCurves %>% 
  rename("Age" = "X", "Year" = "Y") %>% 
  mutate(Age = ages[Age],
         Year = years[Year]) %>% 
  mutate(Year = factor(Year)) %>% 
  filter(Year %in% (1980 + 5*(0:6))) %>% 
ggplot(aes(x = Age,
           y = crudeRates,
           ymin = CI_min,
           ymax = CI_max,
           color = Year,
           group = Year)) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  geom_ribbon(alpha = .25) +
  ylab("Log-Mortality") +
  theme_bw()

# SAVE OUTPUT -----

# setwd(here("Results","Sex"))
# setwd(here("Results","Product"))
setwd(here("Results","Country HMD"))
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

ess_mxt <- diagnosticsCheck(m_output)
which(ess_mxt == min(ess_mxt), arr.ind = T) 

qplot(1:niter, m_output[,4,30])

#

mxt <- computemxt(idx_delta = 0,
                  c(choice1, choice2, choice3, choice4),
                  list_d$param1, list_d$param2,
                  list_d$param3, list_d$param4)
loglik_true <-  - loglik_term1(d, E, mxt,
                               matrix(0, X, Y))

qplot(1:niter, loglik_output) + 
  geom_hline(aes(yintercept = loglik_true))
  
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

qplot(1:niter, m_output[,X,1]) + 
  geom_hline(aes(yintercept = list_d$m[X,1]))

