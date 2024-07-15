
simulateData <- function(X, Y, E,
                         cx,
                         choice1,
                         choice2,
                         choice3,
                         choice4){
  
  
  # Choice 1 (Base effect)
  
  {
    
    if(choice1 == 2){
      param1 <- c(-3,1)
      term1 <- param1[1] + param1[2] * cx
      param1 <- list("ab" = param1)
    } else {
      ax <- seq(-3, -1, length.out = X) + rnorm(X, sd = .03)
      term1 <- ax
      param1 <- list("ax" = ax)
    }
    
  }
  
  k1t <- seq(-1.5, 1.5, length.out = Y)
  if(!(choice1 == 2 & choice2 == 2 & choice3 == 2)){
    k1t <- k1t - mean(k1t)  
  }
  
  param2 <- list("k1t" = k1t)
  
  # Choice 2 
  
  {
     
    param3 <- list()
    
    if(choice2 == 1){
      term2 <- rep(0, Y)
    } else if(choice2 == 2){
      k2t <- k1t
      term2 <- k2t
    } else if(choice2 == 3){
      k2t <- seq(0, 1.5, length.out = Y) + rnorm(Y, sd = .1)
      if(!(choice1 == 2 & choice3 == 2)){
        k2t <- k2t - mean(k2t)  
      }
      term2 <- k2t
      param3$k2t <- k2t
    }
    
  }
  
  # Choice 3
  
  {
    
    if (choice3 == 1){
      
      if(choice2 != 2){
        param3$b_bar <- 1
      } else {
        b_bar <- .5
        param3$b_bar <- b_bar
      }
      
      term3 <- cx * param3$b_bar
      
    } else if (choice3 == 2){
      
      bx <- seq(.3, 0.0, length.out = X) + rnorm(X, sd = .2)
      bx <- bx - mean(bx)  
      if(choice2 == 3){
        bx <- bx / (-bx[1])  
      }
      
      term3 <- bx
      
      param3$bx <- bx
      
    }
    
 
  }
  
  # Choice 4 (Cohort effect)
  
  {
    
    # if(choice2 == 1){
      
      gtx <- rnorm(X + Y - 1, sd = .2)
      
      term4 <- createTerm4_cpp(gtx, X, Y)  
    
      gtx <- gtx - mean(term4)
      
    # } else if(choice2 == 2){
    #   
    #   gtx_m1 <- rnorm(X + Y - 2)
    #   
    #   ik <- create_ik(X, Y)
    #   
    #   gtx_m1 <- gtx_m1 - sum(ik[-1] * gtx_m1) / (X*Y - 1)
    #   
    #   gtx <- c(0, gtx_m1)
    #   
    #   term4 <- createTerm4_cpp(gtx, X, Y)
    #   # (sumParams <- sum(term4))
    #   
    #   # gtx[length(gtx)] <- - sumParams
    #   
    # }
    
    
    if(choice4 == 1){
      term4 <- matrix(0, X, Y)
      param4 <- NULL
    } else {
      term4 <- createTerm4_cpp(gtx, X, Y)
      param4 <- list("gtx" = gtx)
    }
    
  }
  
  # simulate data
  d <- array(NA, dim = c(X, Y))
  m <- array(NA, dim = c(X, Y))
  for (x in 1:X) {
    for (t in 1:Y) {
      m[x,t] <- 
        term1[x] + 
        # k1t[t] * (1 + term2[t] * term3[x]) + 
        k1t[t] + term2[t] * term3[x] +
        term4[t - x + X + 1]
      d[x,t] <- rpois(1, E[x,t] * exp(m[x,t]))  
      
    }
  }
  
  list("d" = d,
       "param1" = param1,
       "param2" = param2,
       "param3" = param3,
       "param4" = param4,
       "m" = m)
  
}

computemxt <- function(idx_delta, # which term to remove
                       delta, param1, param2, param3, param4, ...){
  
  # term1 
  if(delta[1] == 1){
    
    term1 <- matrix(param1$ax, X, Y, byrow = F)
    
  } else {
    
    term1 <- matrix(param1$ab[1] + param1$ab[2] * cx, X, Y, byrow = F)
    
  }
  
  # term2
  term2 <- matrix(param2$k1t, X, Y, byrow = T)  
  
  # term3
  if(delta[2] == 1){ # no additive term
    term3 <- matrix(0, X, Y)
  } else if(delta[2] == 2){
    if(delta[3] == 1){
      term3 <- matrix(param2$k1t, X, Y, byrow = T) *
        param3$b_bar * matrix(cx, X, Y, byrow = F)
    } else if(delta[3] == 2){
      term3 <- matrix(param2$k1t, X, Y, byrow = T) *
        matrix(param3$bx, X, Y, byrow = F)
    }
  } else if(delta[2] == 3){
    if(delta[3] == 1){
      term3 <- matrix(param3$k2t, X, Y, byrow = T) *
        matrix(cx, X, Y, byrow = F)
    } else if(delta[3] == 2){
      term3 <- matrix(param3$k2t, X, Y, byrow = T) *
        matrix(param3$bx, X, Y, byrow = F)
    }

  }
  
  if(delta[4] == 1){
    term4 <- matrix(0, X, Y, byrow = F)  
  } else if(delta[4] == 2){
    term4 <- createTerm4_cpp(param4$gtx, X, Y)
  }
  
  sumTerms <- term1 + term2 + term3 + term4
  
  if(1 %in% idx_delta) sumTerms <- sumTerms - term1
  
  if(2 %in% idx_delta) sumTerms <- sumTerms - term2
  
  if(3 %in% idx_delta) sumTerms <- sumTerms - term3

  if(4 %in% idx_delta) sumTerms <- sumTerms - term4
  
  sumTerms
}

loglik_term1 <- function(d, E, abx, cxt){
  
  mxt <- abx + cxt + log(E)
  
  - sum(d * mxt - exp(mxt))
  
}

sampleMVNconstraint_sum0 <- function(mu, Sigma){
  
  d <- length(mu)
  
  A <- matrix(1, 1, d)
  b <- 0
  
  P <- matrix(0, d, d)
  diag(P)[1:d - 1] <- 1
  P <- P[-d,]
  
  C <- rbind(P, A)
  
  nrow1 <- 1:(d - nrow(A))
  nrow2 <- d - nrow(A) + 1:nrow(A)
  
  nu <- C %*% mu
  Omega <- C %*% Sigma %*% t(C)
  Omega11 <- Omega[nrow1, nrow1]
  Omega22 <- Omega[nrow2, nrow2]
  Omega12 <- Omega[nrow1, nrow2]
  Omega21 <- Omega[nrow2, nrow1]
  
  mu1 <- nu[nrow1] + Omega12 %*% solve(Omega22) %*% (b - nu[nrow2])
  Sigma1 <- Omega11 - Omega12 %*% solve(Omega22) %*% Omega21
  
  pxcd <- mvrnorm(1, mu1, Sigma1)
  
  x <- solve(C) %*% c(pxcd, b)
  
  x
}

createTerm4_cpp <- function(gtx, X, Y){
  
  term4 <- matrix(0, X, Y)
  for (g in 1:(X+Y-1)) {
    sdiag(term4, k = g - X) <- gtx[g]
  }
  
  term4
}

create_ik <- function(X,Y){
  c(
    seq(1, (min(X,Y) - 1)),
    rep(min(X,Y), max(X,Y) - min(X,Y) + 1),
    X + Y - seq(max(X,Y) + 1, X + Y - 1)
  )
}

# term1

# second derivative of the loglikelihood with respect to a

derl2der2a <- function(a, cx){
  
  - exp(a) * sum(exp(cx))
  
}

buildProposalAx <- function(d, E, cxt){
  
  ax_star <- rep(0, X)
  sd_star <- rep(0, X)
  
  for (x in 1:X) {
    
    ax_star[x] <- log(
      sum(d[x,]) / sum(exp(cxt[x,] + log(E[x,])))
    )
    
    hessian_ax_star <- derl2der2a(ax_star[x], cxt[x,] + log(E[x,]))
    sd_star[x] <- sqrt(- 1 / hessian_ax_star)
    
  }
  
  list("ax_star" = ax_star,
       "sd_star" = sd_star)
}

updateAX <- function(ax, d, E, cxt) {
  
  for (x in 1:X) {
    
    ax_star <- log(
      sum(d[x,]) / sum(exp(cxt[x,] + log(E[x,])))
    )
    
    hessian_ax_star <- derl2der2a(ax_star, cxt[x,] + log(E[x,]))
    sd_star <- sqrt(- 1 / hessian_ax_star)
    
    ax_proposal <- rnorm(1, ax_star, sd_star)
    
    loglik_proposal <- loglik_term1_a(ax_proposal, x, d, E, cxt)
    loglik_current <- loglik_term1_a(ax[x], x, d, E, cxt)
    
    logproposal_proposal <- dnorm(ax_proposal, ax_star, sd_star, log = T)
    logproposal_current <- dnorm(ax[x], ax_star, sd_star, log = T)
    
    mh_ratio <- exp(loglik_proposal - loglik_current + 
                      logproposal_proposal - logproposal_current)
    
    if(runif(1) < mh_ratio){
      
      ax[x] <- ax_proposal
      
    }
    
  }
  
  ax
}

updateABX <- function(ab, cx, d, E, cxt){
  
  list_proposal_abx <- findProposalABX(cx, d, E, cxt)
  ab_star <- list_proposal_abx$ab_star
  Sigma_star <- list_proposal_abx$Sigma_star
  
  ab_proposal <- mrt2(ab_star, Sigma_star, df = 3)
  
  loglik_ab_current <- loglik_term1_ab(cx, d, E, cxt)
  
  loglik_proposal <- - loglik_ab_current(ab_proposal)
  loglik_current <- - loglik_ab_current(ab)
  
  logproposal_proposal <- dmt_cpp(ab_proposal, nu = 3, ab_star, Sigma_star, T)
  logproposal_current <- dmt_cpp(ab,  nu = 3, ab_star, Sigma_star, T)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    ab <- ab_proposal
    
  }
  
  ab
  
}

updateTerm1 <- function(param, delta, d, E){
  
  term234 <- computemxt(idx_delta = 1, delta, 
                        param$param1, param$param2,
                        param$param3, param$param4)
  
  if(delta[1] == 1){ # unconstrained model
    
    ax <- updateAX(param$param1$ax, d, E, term234)  
    term1 <- matrix(ax, X, Y, byrow = F)
    
    param1 <- list("ax" = ax)
    
    return(
      param1
    )
    
  } else if(delta[1] == 2){ # linear model
    
    ab <- updateABX(param$param1$ab, cx, d, E, term234)
    
    param1 <- list("ab" = ab)
    
    return(
      param1
    )
    
  }
  
}

# term2

gr_loglik_k1t_add_m1 <- function(kt, d, E, cxt){
  
  function(param){
    
    kt <- c(param, -sum(param))
    
    kt_mat <- matrix(kt, X, Y, byrow = T)
    
    # priorTerm <- sum(
    #   dnorm(bx, 0, sd = sd_bx, log = T)
    # )
    
    loglik_term1(d, E, kt_mat, cxt) #- priorTerm
    
  }
  
}

loglik_k1t_add_m1 <- function(kt, d, E, cxt){
  
  function(param){
    
    kt <- c(param, -sum(param))
    
    kt_mat <- matrix(kt, X, Y, byrow = T)
    
    # priorTerm <- sum(
    #   dnorm(bx, 0, sd = sd_bx, log = T)
    # )
    
    loglik_term1(d, E, kt_mat, cxt) #- priorTerm
    
  }
  
}

buildProposalK1t_m1 <- function(kt, d, E, cxt){
  
  loglik_kt_current <- loglik_k1t_add_m1(kt, d, E, cxt)
  # gr_loglik_bx_current <- gr_loglik_term3_bx_k2t_m1(kt, d, E, cxt, sd_bx)
  
  # find maximum
  laplace_fit <- optim(
    par = kt[-Y],
    fn = loglik_kt_current,
    # gr = gr_loglik_bx_k2t_current,
    method = c("BFGS"),
    hessian = T
  )
  
  ktm1_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  list("ktm1_star" = ktm1_star,
       "Sigma_star" = Sigma_star)
  
}

buildProposalK1t <- function(d, E, cxt){
  
  k_star <- rep(NA, Y)
  Sigma_star <- matrix(0, Y, Y)
  
  for (t in 1:Y) {
    
    k_star[t] <- log(
      sum(d[,t]) / sum(exp(cxt[,t] + log(E[,t])))
    )
    
    hessian_kt_star <- derl2der2a(k_star[t], cxt[,t] + log(E[,t]))
    Sigma_star[t,t] <- - 1 / hessian_kt_star
    
  }
  
  Sigma_star <- 2 * Sigma_star
  
  list("k_star" = k_star,
       "Sigma_star" = Sigma_star)
  
}

loglik_term1_k <- function(k, d, E, cxt){
  
  mxt <- matrix(k, X, Y, byrow = T) + cxt + log(E)
  
  sum(d * mxt - exp(mxt))
  
}

updateK1T_add <- function(kt, d, E, cxt) {
  
  # list_proposal <- buildProposalK1t(d, E, cxt)
  list_proposal <- buildProposalK1t_m1(kt, d, E, cxt)
  km1_star <- list_proposal$ktm1_star
  Sigma_star <- list_proposal$Sigma_star
  
  k_proposed_m1 <- mrt2(km1_star, Sigma_star, df = 3)
  k_proposed <- c(k_proposed_m1, - sum(k_proposed_m1))
  
  # k_proposed <- as.vector(sampleMVNconstraint_sum0(k_star, Sigma_star))
  
  loglik_star <- loglik_term1_k(k_proposed, d, E, cxt)
  loglik_current <- loglik_term1_k(kt, d, E, cxt)
  
  logproposal_star <- dmt_cpp(k_proposed_m1, nu = 3, km1_star, 
                                  Sigma_star, T)
  
  logproposal_current <-  dmt_cpp(kt[-Y], nu = 3, km1_star, 
                                  Sigma_star, T)
  
  # logproposal_star <- sum(dnorm(k_proposed, k_star, 
  #                               sqrt(diag(Sigma_star)), log = T))
  # 
  # logproposal_current <- sum(dnorm(kt, k_star, 
  #                                  sqrt(diag(Sigma_star)), log = T))
  
  mh_ratio <- exp(loglik_star - loglik_current + 
                    logproposal_current  - logproposal_star)
  
  if(runif(1) < mh_ratio){
    kt <- as.vector(k_proposed)
  }
  
  kt
  
}

loglik_kt_bx_m1 <- function(bx, d, E, cxt){
  
  function(param){
    
    k2tm1 <- param[1:(Y-1)]
    
    k2t <- c(k2tm1, - sum(k2tm1))
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

gr_loglik_kt_bx_m1 <- function(bx, d, E, cxt){
  
  function(param){
    
    k2tm1 <- param[1:(Y-1)]
    k2t <- c(k2tm1, - sum(k2tm1))
    
    k2t_mat <- matrix(k2t, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term1_kt <- d * bx_mat - exp(lambda_xt) * bx_mat
    
    term2_kt <- (- d[,Y] * bx_mat[,Y] ) - exp(lambda_xt[,Y]) * (- bx_mat[,Y])
    
    gr_kt <- apply(term1_kt[,1:(Y-1)], 2, sum) + sum(term2_kt)
    
    - gr_kt
    
  }
  
}

loglik_kt_bx <- function(bx, d, E, cxt){
  
  function(param){
    
    k2t <- param
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

gr_loglik_kt_bx <- function(bx, d, E, cxt){
  
  function(param){
    
    k2t <- param
    
    k2t_mat <- matrix(k2t, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term_kt <- d * bx_mat - exp(lambda_xt) * bx_mat
    
    gr_kt <- apply(term_kt, 2, sum)
    
    - gr_kt
    
  }
  
}

findProposalkt_given_bx <- function(k2t, bx, d, E, cxt, sum0){
  
  if(sum0){
    loglik_kt_current <- loglik_kt_bx_m1(bx, d, E, cxt)
    gr_loglik_kt_current <- gr_loglik_kt_bx_m1(bx, d, E, cxt)
    startPar <-  k2t[-Y]
  } else {
    loglik_kt_current <- loglik_kt_bx(bx, d, E, cxt)
    gr_loglik_kt_current <- gr_loglik_kt_bx(bx, d, E, cxt)
    startPar <-  k2t
  }
  
  # find maximum
  laplace_fit <- optim(
    par = startPar,
    fn = loglik_kt_current,
    gr = gr_loglik_kt_current,
    method = c("CG"),
    hessian = T
  )
  
  kt_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  list("kt_star" = kt_star,
       "Sigma_star" = Sigma_star)
  
}

loglik_term3 <- function(k2t, bx, d, E, cxt){
  
  k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
  
  loglik_term1(d, E, k2tbx, cxt)
  
}

updateMultTermKt <- function(kt, bx, d, E, cxt, sum0){
  
  list_proposal <- findProposalkt_given_bx(kt, bx, d, E, cxt, sum0)
  kt_star <- list_proposal$kt_star
  Sigma_kt_star <- list_proposal$Sigma_star
  
  if(sum0){
    
    kt_proposed_m1 <- mrt2(kt_star,
                           Sigma_kt_star,
                           df = 3)
    
    kt_proposed <- c(kt_proposed_m1, - sum(kt_proposed_m1))
    
    logproposal_proposal <- dmt_cpp(kt_proposed_m1, nu = 3, kt_star, 
                                    Sigma_kt_star, returnLog = T)
    logproposal_current <- dmt_cpp(kt[-Y], nu = 3, kt_star, 
                                   Sigma_kt_star, returnLog = T)
    
  } else {
    
    kt_proposed <- mrt2(kt_star,
                        Sigma_kt_star,
                        df = 3)
    
    logproposal_proposal <- dmt_cpp(kt_proposed, nu = 3, kt_star, 
                                    Sigma_kt_star, returnLog = T)
    logproposal_current <- dmt_cpp(kt, nu = 3, kt_star, 
                                   Sigma_kt_star, returnLog = T)
    
  }
  
  loglik_proposal <- - loglik_term3(kt_proposed, bx, d, E, cxt)
  loglik_current <- - loglik_term3(kt, bx, d, E, cxt)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    kt <- kt_proposed
    
  }
  
  kt
}

updateK1 <- function(param, delta, d, E){
  
  # compute remaining additive terms of mxt
  
  if(delta[2] != 3){
    
    cxt <- computemxt(idx_delta = c(2,3), delta, 
                          param$param1, param$param2,
                          param$param3, param$param4)
    
  } else {
    
    cxt <- computemxt(idx_delta = 2, delta, 
                         param$param1, param$param2,
                         param$param3, param$param4)
    
  }
  
  # update
  
  # additive update
  if(delta[2] != 2){
    
    param$param2$k1t <- updateK1T_add(param$param2$k1t, d, E, cxt)
    
  } else { # multiplicative update
    
    sum0 <- !(delta[1] == 2 & 
                delta[3] == 2)
    if(delta[3] == 1){
      bx_current <- (1 + param$param3$b_bar * cx)
    } else if(delta[3] == 2){
      bx_current <- (1 + param$param3$bx)
    }
    
    param$param2$k1t <- updateMultTermKt(param$param2$k1t, bx_current, 
                                         d, E, cxt, sum0)
    
  }
  
  param$param2$k1t
}

# term3.1

updateK2 <- function(param, delta, d, E){
  
  if(delta[2] == 3){
  # compute remaining additive terms of mxt
  
  cxt <- computemxt(idx_delta = 3, delta, 
                    param$param1, param$param2,
                    param$param3, param$param4)
  
  # compute multiplicative term
  if(delta[3] == 1){
    bx_current <- cx
  } else if(delta[3] == 2){
    bx_current <- param$param3$bx
  }
  
  sum0 <- !(delta[1] == 2 & 
              delta[3] == 2)
  
  param$param3$k2t <- updateMultTermKt(param$param3$k2t, 
                          bx_current, d, E, cxt, sum0)
  
  }
  
  param$param3$k2t
}

# term3.2

loglik_bbar_fun <- function(b_bar, ktcx, cxt, d, E){
  
  ktcxbbar <- ktcx * b_bar
  
  loglik_term1(d, E, ktcxbbar, cxt)
  
}

loglik_bbar <- function(ktcx, cxt, d, E){
  
  function(param){
    
    loglik_bbar_fun(param, ktcx, cxt, d, E)
    
  }
  
}
  
findProposalBBar <- function(bbar, ktcx, cxt, d, E){
  
  loglik_bbar_current <- loglik_bbar(ktcx, cxt, d, E)
  # gr_loglik_bbar_current <- gr_loglik_kt_bx_m1(bx, d, E, cxt)
  
  # find maximum
  laplace_fit <- optim(
    par = bbar,
    fn = loglik_bbar_current,
    # gr = gr_loglik_kt_current,
    method = c("BFGS"),
    hessian = T
  )
  
  bbar_star <- laplace_fit$par
  sd_star <- sqrt( 1 / laplace_fit$hessian)
  
  list("b_bar" = bbar_star,
       "sd_star" = sd_star)
  
}

updateBbar <- function(b_bar, ktcx, cxt, d, E){
  
  list_proposal <- findProposalBBar(b_bar, ktcx, cxt, d, E)
  bbar_star <- list_proposal$b_bar
  sd_star <- list_proposal$sd_star
  
  bbar_proposal <- as.numeric(rt2(1, bbar_star, sd_star, df = 3))
  
  loglik_bbar_current <- loglik_bbar(ktcx, cxt, d, E)
  
  loglik_proposal <- - loglik_bbar_current(bbar_proposal)
  loglik_current <- - loglik_bbar_current(b_bar)
  
  logproposal_proposal <- log(dt2(bbar_proposal, bbar_star, 
                                sd_star, df = 3))
  logproposal_current <- log(dt2(b_bar, bbar_star, 
                               sd_star, df = 3))
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    b_bar <- bbar_proposal
    
  }
  
  b_bar
}

loglik_term3_bx_m2 <- function(kt, d, E, cxt, sd_bx){
  
  function(param){
    
    bxm1X <- param[1:(X-2)]
    
    bx <- c(-1, bxm1X, - sum(bxm1X) + 1)
    
    ktbx <- matrix(kt, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    priorTerm <- sum(
      dnorm(bx, 0, sd = sd_bx, log = T)
    )
    
    loglik_term1(d, E, ktbx, cxt) - priorTerm
    
  }
  
}

gr_loglik_term3_bx_m2 <- function(kt, d, E, cxt, sd_bx){
  
  function(param){
    
    bxm1X <- param[1:(X-2)]
    bx <- c(-1, bxm1X, - sum(bxm1X) + 1)
    
    k2t_mat <- matrix(kt, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term1_bx <- d * k2t_mat - exp(lambda_xt) * k2t_mat
    
    term2_bx <- d[X,] * (-k2t_mat[X,]) - exp(lambda_xt[X,]) * (- k2t_mat[X,])
    
    gr_bx <- apply(term1_bx[2:(X-1),], 1, sum) + sum(term2_bx)
    
    - gr_bx + 1 / sd_bx^2
    
  }
  
}

findProposalBX_m2 <- function(bx, kt, d, E, cxt, sd_bx){
  
  loglik_bx_current <- loglik_term3_bx_m2(kt, d, E, cxt, sd_bx)
  # gr_loglik_bx_current <- gr_loglik_term3_bx_m2(kt, d, E, cxt, sd_bx)
  
  # find maximum
  laplace_fit <- optim(
    par = bx[-c(1,X)],
    # par = rep(0, X - 2),
    fn = loglik_bx_current,
    # gr = gr_loglik_bx_k2t_current,
    method = c("CG"),
    hessian = T
  )
  
  bx_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  list("bx_star" = bx_star,
       "Sigma_star" = Sigma_star)
  
}

loglik_term3_bx_m1 <- function(kt, d, E, cxt, sd_bx){
  
  function(param){
    
    bxmX <- param
    
    bx <- c(bxmX, -sum(bxmX))
    
    ktbx <- matrix(kt, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    priorTerm <- sum(
      dnorm(bx, 0, sd = sd_bx, log = T)
    )
    
    loglik_term1(d, E, ktbx, cxt) - priorTerm
    
  }
  
}

gr_loglik_term3_bx_m1 <- function(k2t, d, E, cxt, sd_bx){
  
  function(param){
    
    bxmX <- param
    
    bx <- c(bxmX, - sum(bxmX))
    
    k2t_mat <- matrix(k2t, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term1_bx <- d * k2t_mat - exp(lambda_xt) * k2t_mat
    
    # gr_bx <- apply(term1_bx[1:X,], 1, sum)
    
    term2_bx <- d[X,] * (-k2t_mat[X,]) - exp(lambda_xt[X,]) * (- k2t_mat[X,])
    
    gr_bx <- apply(term1_bx[1:(X-1),], 1, sum) + sum(term2_bx)
    
    - gr_bx + 1 / sd_bx^2
    
  }
  
}

findProposalBX_m1 <- function(bx, kt, d, E, cxt, sd_bx){
  
  loglik_bx_current <- loglik_term3_bx_m1(kt, d, E, cxt, sd_bx)
  gr_loglik_bx_current <- gr_loglik_term3_bx_m1(kt, d, E, cxt, sd_bx)
  
  # find maximum
  laplace_fit <- optim(
    par = bx[-X],
    fn = loglik_bx_current,
    gr = gr_loglik_bx_current,
    method = c("BFGS"),
    hessian = T
  )
  
  bx_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  list("bx_star" = bx_star,
       "Sigma_star" = Sigma_star)
  
}

updateBX <- function(bx, kt, d, E, cxt, delta, sd_bx){
  
  if(delta[2] == 2){ # no constraint that b_1 = -1
    
    list_proposal <- findProposalBX_m1(bx, kt, d, E, cxt, sd_bx)
    bx_star <- list_proposal$bx_star
    Sigma_bx_star <- list_proposal$Sigma_star
    
    bx_proposed_m1 <- mrt2(bx_star, Sigma_bx_star, df = 3)
    bx_proposed <- c(bx_proposed_m1, - sum(bx_proposed_m1))
    
    loglik_bx_current <- loglik_term3_bx_m1(kt, d, E, cxt, sd_bx)
    
    loglik_proposal <- - loglik_bx_current(bx_proposed[-X])
    loglik_current <- - loglik_bx_current(bx[-X])
    
    logproposal_proposal <- dmt_cpp(bx_proposed_m1, nu = 3, 
                                    bx_star, Sigma_bx_star, T)
    logproposal_current <- dmt_cpp(bx[-X], nu = 3,
                                   bx_star, Sigma_bx_star, T)
    
    logprior_proposal <- sum(dnorm(bx_proposed, sd = sd_bx, log = T))
    logprior_current <- sum(dnorm(bx, sd = sd_bx, log = T))
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current + 
        logprior_proposal - logprior_current + 
        logprior_current - logprior_proposal
    )
    
    if(runif(1) < mh_ratio){
      bx <- bx_proposed
    }
    
    
  } else if(delta[2] == 3){
    
    list_proposal <- findProposalBX_m2(bx, kt, d, E, cxt, sd_bx)
    bx_star <- list_proposal$bx_star
    Sigma_bx_star <- list_proposal$Sigma_star
    
    bx_proposed_m2 <- mrt2(bx_star, Sigma_bx_star, df = 3)
    
    bx_proposed <- c(-1, bx_proposed_m2, - sum(bx_proposed_m2) + 1)
    
    loglik_bx_current <- loglik_term3_bx_m2(kt, d, E, cxt, sd_bx)
    
    loglik_proposal <- - loglik_bx_current(bx_proposed[-c(1,X)])
    loglik_current <- - loglik_bx_current(bx[-c(1,X)])
    
    logproposal_proposal <- dmt_cpp(bx_proposed_m2, nu = 3,
                                    bx_star, Sigma_bx_star, T)
    logproposal_current <- dmt_cpp(bx[-c(1,X)], nu = 3,
                                   bx_star, Sigma_bx_star, T)
    
    logprior_proposal <- sum(dnorm(bx_proposed, sd = sd_bx, log = T))
    logprior_current <- sum(dnorm(bx, sd = sd_bx, log = T))
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current + 
        logprior_proposal - logprior_current + 
        logproposal_current - logproposal_proposal
    )
    
    if(runif(1) < mh_ratio){
      
      bx <- bx_proposed
      
    }
    
  }
  
  bx
}

updateBXB <- function(param, delta, d, E, ...){
  
  # compute remaining additive terms of mxt
  
  cxt <- computemxt(idx_delta = 3, delta, 
                    param$param1, param$param2,
                    param$param3, param$param4)
  
  if(delta[3] == 1 & delta[2] == 2){ # cx
    
    ktcx <- 
      matrix(param$param2$k1t, X, Y, byrow = T) *
      matrix(cx, X, Y, byrow = F)
    
    param$param3$b_bar <- updateBbar(param$param3$b_bar, ktcx, cxt, d, E)
    
  } else if(delta[3] == 2){ # bx
    
    if(delta[2] == 2){
      kt_current <- param$param2$k1t  
    } else if(delta[2] == 3){
      kt_current <- param$param3$k2t
    }
    
    param$param3$bx <- updateBX(param$param3$bx, kt_current,
                                d, E, cxt, delta, sd_bx)
    
    
  } 
  
  return(param$param3)
}

# term 4

loglik_term4_gtx_m2 <- function(d, E, cxt){
  
  function(param){
    
    param_matrix <- createTerm4_cpp(c(0, param, 0), X, Y)
    
    sumParams <- sum(param_matrix)
    
    gtx <- c(0, param, - sumParams)
    
    gtx_mat <- createTerm4_cpp(gtx, X, Y)
    
    loglik_term1(d, E, gtx_mat, cxt)
    
  }
  
}

gr_loglik_term4_gtx_m2 <- function(d, E, cxt){
  
  ik <- create_ik(X,Y)
  
  function(param){
    
    param_matrix <- createTerm4_cpp(c(0, param, 0), X, Y)
    
    sumParams <- sum(param_matrix)
    
    gtx <- c(0, param, - sumParams)
    
    gtx_mat <- createTerm4_cpp(gtx, X, Y)
    
    term1 <- sapply(1:(X+Y-2), function(g){
      sum(sdiag(d, k = g - X))
    })
    
    term2 <- sapply(1:(X+Y-2), function(g){
      - sum(sdiag(exp(cxt + gtx_mat + log(E)), k = g - X))
    })
    
    term3 <- d[1,Y] * (- ik[-c(X+Y-1)])
    
    term4 <- - exp(cxt[1,Y] + gtx_mat[1,Y] + log(E[1,Y])) * (-ik[-c(X+Y-1)])
    
    - (term1 + term2 + term3 + term4)[-1]
    
  }
  
}

buildProposalGtx_m2 <- function(gtx, d, E, cxt){
  
  loglik_gtx_current <- loglik_term4_gtx_m2(d, E, cxt)
  gr_loglik_gtx_current <- gr_loglik_term4_gtx_m2(d, E, cxt)
  
  # find maximum
  laplace_fit <- optim(
    par = gtx[-c(1,X+Y-1)],
    fn = loglik_gtx_current,
    gr = gr_loglik_gtx_current,
    method = c("BFGS"),
    hessian = T
  )
  
  gtx_star_m2 <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  list("gtx_star" = gtx_star_m2,
       "Sigma_star" = Sigma_star)
  
}

updateGTX_m2 <- function(gtx, d, E, cxt) {
  
  list_proposal <- buildProposalGtx_m2(gtx, d, E, cxt)
  gtxm2_star <- list_proposal$gtx_star
  Sigma_star <- list_proposal$Sigma_star
  
  gtx_m2_proposed <- mrt2(gtxm2_star, Sigma_star, df = 3)
  
  gtx_proposed <- c(0, gtx_m2_proposed, 0)
  term4_proposed <- createTerm4_cpp(gtx_proposed, X, Y)
  gtx_proposed[length(gtx_proposed)] <- - sum(term4_proposed)
  
  loglik_star <- loglik_term_gtx(gtx_proposed, d, E, cxt)
  loglik_current <- loglik_term_gtx(gtx, d, E, cxt)
  
  logproposal_star <- 
    dmt_cpp(gtx_proposed[-c(1,length(gtx_proposed))], 
            nu = 3,
            gtxm2_star, 
            Sigma_star, T)
  
  logproposal_current <- 
    dmt_cpp(gtx[-c(1,length(gtx_proposed))], 
            nu = 3,
            gtxm2_star, 
            Sigma_star, T)
  
  mh_ratio <- exp(loglik_star - loglik_current + 
                    logproposal_current  - logproposal_star)
  
  if(runif(1) < mh_ratio){
    gtx <- as.vector(gtx_proposed)
  }
  
  gtx
  
}

loglik_term4_gtx_m1 <- function(d, E, cxt){
  
  ik <- create_ik(X,Y)
  
  function(param){
    
    gtx <- c(param, - sum(ik[-length(ik)] * param))
    
    gtx_mat <- createTerm4_cpp(gtx, X, Y)
    
    loglik_term1(d, E, gtx_mat, cxt)
    
  }
  
  
}

gr_loglik_term4_gtx_m1 <- function(d, E, cxt){
  
  ik <- create_ik(X,Y)
  
  function(param){
  
    gr_loglik_term4_gtx_m1_fun_cpp(param,
                                   d, E, cxt, 
                                   create_ik(X,Y)[-c(X+Y-1)])
  
  #   
  #   param_matrix <- createTerm4_cpp(c(param, 0), X, Y)
  #   
  #   sumParams <- sum(param_matrix)
  #   
  #   gtx <- c(param, - sumParams)
  #   
  #   gtx_mat <- createTerm4_cpp(gtx, X, Y)
  #   
  #   term1 <- sapply(1:(X+Y-2), function(g){
  #     sum(sdiag(d, k = g - X))
  #   })
  #   
  #   term2 <- sapply(1:(X+Y-2), function(g){
  #     - sum(sdiag(exp(cxt + gtx_mat + log(E)), k = g - X))
  #   })
  #   
  #   term3 <- d[1,Y] * (- ik[-c(X+Y-1)])
  #   
  #   term4 <- - exp(cxt[1,Y] + gtx_mat[1,Y] + log(E[1,Y])) * (-ik[-c(X+Y-1)])
  # 
  #   - (term1 + term2 + term3 + term4)
  #   
  }
  
}

hessian_loglik_term4_gtx_m1 <- function(d, E, cxt){
  
  ik <- create_ik(X,Y)
  
  function(param){
    
    hess_loglik_term4_gtx_m1_fun_cpp(param,
                                     d, E, cxt, 
                                     ik[-c(X+Y-1)])
    
    # param_matrix <- createTerm4_cpp(c(param, 0), X, Y)
    # 
    # sumParams <- sum(param_matrix)
    # 
    # gtx <- c(param, - sumParams)
    # 
    # gtx_mat <- createTerm4_cpp(gtx, X, Y)
    # 
    # hess <- matrix(0, X + Y - 2, X + Y - 2)
    # 
    # for (g1 in 1:(X+Y-2)) {
    #   
    #   # term2
    #   
    #   hess[g1,g1] <- hess[g1,g1] - 
    #     sum(sdiag(exp(cxt + gtx_mat + log(E)), k = g1 - X))
    #     
    #   for (g2 in 1:(X+Y-2)) {
    #     
    #     # term3
    #     
    #     hess[g1,g2] <- hess[g1,g2] -
    #       ik[g1] * ik[g2] * exp(cxt[1,Y] + gtx_mat[1,Y] + log(E[1,Y]))
    #     
    #   }  
    # }
    
    hess
  }
  
}

buildProposalGtx_m1 <- function(gtx, d, E, cxt){
  
  loglik_gtx_current <- loglik_term4_gtx_m1(d, E, cxt)
  gr_loglik_gtx_current <- gr_loglik_term4_gtx_m1(d, E, cxt)
  # hess_loglik_gtx_current <- hessian_loglik_term4_gtx_m1(d, E, cxt)
  
  # find maximum
  laplace_fit <- optim(
    par = gtx[-c(X+Y-1)],
    fn = loglik_gtx_current,
    gr = gr_loglik_gtx_current,
    method = c("CG"),
    hessian = T
  )
  
  # f <- function(param) {
  #   out <- loglik_gtx_current(param)
  #   attr(out, 'gradient') <- gr_loglik_gtx_current(param)
  #   attr(out, 'hessian') <-  hess_loglik_gtx_current(param)
  #   return(out)
  # }
  # 
  # out_nlm <- nlm(f, gtx[-c(X+Y-1)], hessian = TRUE, check.analyticals = TRUE)
  
  gtx_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  # hess1 <- hess_loglik_gtx_current(gtx_star)
  # hess2 <- laplace_fit$hessian
  
  list("gtx_star" = gtx_star,
       "Sigma_star" = Sigma_star)
  
}

loglik_term_gtx <- function(gtx, d, E, cxt){
  
  sum(
    sapply(1:(X+Y-1), function(i){
      
      d_current <- sdiag(d, k = i - X)
      E_current <- sdiag(E, k = i - X)
      cxt_current <- sdiag(cxt, k = i - X)
      
      mxt_current <- gtx[i] + cxt_current + log(E_current)
      
      sum(d_current * mxt_current - exp(mxt_current))
      
    })
  )
  
}

buildProposalGtx <- function(d, E, cxt){
  
  gtx_star <- rep(NA, X + Y - 1)
  Sigma_star <- matrix(0, X + Y - 1, X + Y - 1)
  
  for (i in 1:(X+Y-1)) {
    
    d_current <- sdiag(d, k = i - X)
    E_current <- sdiag(E, k = i - X)
    cxt_current <- sdiag(cxt, k = i - X)
    
    gtx_star[i] <- log(
      sum(d_current) / sum(exp(cxt_current + log(E_current)))
    )
    
    hessian_gxt_star <- derl2der2a(gtx_star[i], 
                                   cxt_current + log(E_current))
    
    Sigma_star[i,i] <- - 1 / hessian_gxt_star
    
  }
  
  list("gtx_star" = gtx_star,
       "Sigma_star" = Sigma_star)
  
}

updateGTX <- function(gtx, d, E, cxt) {
  
  list_proposal <- buildProposalGtx_m1(gtx, d, E, cxt)
  gtx_star <- list_proposal$gtx_star
  Sigma_star <- list_proposal$Sigma_star
  
  # gtx_proposed <- as.vector(sampleMVNconstraint_sum0(gtx_star, Sigma_star))
  
  gtx_proposed_m1 <- mrt2(gtx_star, Sigma_star, df = 3)
  ik <- create_ik(X,Y)
  term4_sum <- sum(ik[-length(ik)] * gtx_proposed_m1)
                     
  gtx_proposed <- c(gtx_proposed_m1, -term4_sum)

  logproposal_star <- dmt_cpp(gtx_proposed_m1, nu = 3,
                              gtx_star, Sigma_star, T)  
  logproposal_current <- dmt_cpp(gtx[-c(X+Y-1)], nu = 3,
                                 gtx_star, Sigma_star, T)  
  
  # find proposal for each variable individually and sum
  # logproposal_star <- sum(dnorm(gtx_proposed, gtx_star, 
  #                               sqrt(diag(Sigma_star)), log = T))
  # 
  # logproposal_current <- sum(dnorm(gtx, gtx_star, 
  #                                  sqrt(diag(Sigma_star)), log = T))
  
  loglik_star <- loglik_term_gtx(gtx_proposed, d, E, cxt)
  loglik_current <- loglik_term_gtx(gtx, d, E, cxt)
  
  mh_ratio <- exp(loglik_star - loglik_current + 
                    logproposal_current  - logproposal_star)
  
  if(runif(1) < mh_ratio){
    gtx <- as.vector(gtx_proposed)
  }
  
  gtx
  
}

updateTerm4 <- function(param, delta, d, E){
  
  # compute remaining additive terms of mxt
  
  cxt <- computemxt(idx_delta = 4, delta, 
                    param$param1, param$param2,
                    param$param3, param$param4)
  
  if(delta[4] == 2){ # cohort effect present
    
    if(delta[2] == 2){ 
      
      gtx <- updateGTX(param$param4$gtx, d, E, cxt)
      
    } else {
      # in this case we need the additional constraint gamma_1 = 0
      
      gtx <- updateGTX_m2(param$param4$gtx, d, E, cxt)
      
    }
    
    param$param4$gtx <- gtx
    
  } 
  
  param$param4
}

# moving

proposeNewDelta <- function(idx_delta, delta){
  
  if(idx_delta == 1){
    delta1_star <- ifelse(delta[1] == 1, 2, 1)
    delta_star <- c(delta1_star, delta[2], delta[3], delta[4])
  } else if(idx_delta == 2){
    if(delta[2] == 1){
      delta2_star <- 2
      delta3_star <- 1#sample(c(1,2), 1)
    } else if(delta[2] == 2){
      if(delta[3] == 2){ # bx
        delta2_star <- 3#sample(c(1,3), 1)
        delta3_star <- delta[3]
      } else {
        delta2_star <- sample(c(1,3), 1)
        delta3_star <- delta[3]
      }
    } else if(delta[2] == 3){
      delta2_star <- 2
      delta3_star <- delta[3]
    }
    delta_star <- c(delta[1], delta2_star, delta3_star, delta[4])
    
  } else if(idx_delta == 3){
    delta3_star <- ifelse(delta[3] == 1, 2, 1)
    delta_star <- c(delta[1], delta[2], delta3_star, delta[4])
  } else if(idx_delta == 4){
    delta4_star <- ifelse(delta[4] == 1, 2, 1)
    delta_star <- c(delta[1], delta[2], delta[3], delta4_star)
  } 
  
  delta_star
}

isModelBeingExtended <- function(delta, delta_star){
  
  delta[1] <- 3 - delta[1]
  delta_star[1] <- 3 - delta_star[1]
  
  isModelExtended <- sum(delta) < sum(delta_star)
  
  return(isModelExtended)
  
}

mapModels <- function(delta, delta_star,
                      param1, param2, param3, param4,
                      ...){
  
  param_new <- list("param1" = param1,
                    "param2" = param2,
                    "param3" = param3,
                    "param4" = param4)
  
  # changing the first term
  if(delta[1] != delta_star[1]){
    
    if(delta[1] == 1){ # going from unconstrained model to linear
      
      if(delta[2] == 2 & delta[3] == 2){ 
        # in the lee carter, here you should in theory
        # disregard the constraint sum k1t = 0 and therefore match
        # a + b c_x + k_bar b_x to a_x
        
        # TO CHECK AGAIN
        
        ab_tilde <- lm(param1$ax ~ cx)
        ab_tilde <- as.numeric(ab_tilde$coefficients)
        
        param_new$param1$ab <- ab_tilde
        param_new$param1$ax <- NULL
        
      } else if(delta[2] == 3 & delta[3] == 2){ 
        # here you should in theory
        # disregard the constraint sum kt2 = 0 and therefore match
        # a + b c_x + k_bar b_x to a_x
        
        # TO CHECK AGAIN
        
        ab_tilde <- lm(param1$ax ~ cx)
        ab_tilde <- as.numeric(ab_tilde$coefficients)
        
        param_new$param1$ab <- ab_tilde
        param_new$param1$ax <- NULL
        
      } else {
        
        ab_tilde <- lm(param1$ax ~ cx)
        ab_tilde <- as.numeric(ab_tilde$coefficients)
        
        param_new$param1$ab <- ab_tilde
        param_new$param1$ax <- NULL
        
      }
      
    } else { # going from linear to unconstrainted
      
      if(delta[2] == 2 & delta[3] == 2){ 
        # in the lee carter, so now enforcing the constraint sum kt = 0
        
        ax_tilde <- param1$ab[1] + param1$ab[2] * cx
        
        k1t_current <- param2$k1t
        kbar <- mean(k1t_current)
        k1t_tilde <- k1t_current - kbar
        
        param_new$param2$k1t <- k1t_tilde

        param_new$param1$ax <- ax_tilde + kbar
        param_new$param1$ab <- NULL
        
      } else if(delta[2] == 3 & delta[3] == 2){ 
       
        # we now need to enforce sum k2t = 0
        
        ax_tilde <- param1$ab[1] + param1$ab[2] * cx
        
        k2t_current <- param3$k2t
        kbar <- mean(k2t_current)
        k2t_tilde <- k2t_current - kbar
        
        param_new$param3$k2t <- k2t_tilde
        
        param_new$param1$ax <- ax_tilde + kbar
        param_new$param1$ab <- NULL
        
      } else {
        
        ax_tilde <- param1$ab[1] + param1$ab[2] * cx
        
        param_new$param1$ax <- ax_tilde
        param_new$param1$ab <- NULL
        
      }
      
    }
    
  }
  
  # changing the second term 
  if(delta[2] != delta_star[2]){
    
    # extending to k1t (1 + ..)
    if(delta[2] == 1 & delta_star[2] == 2){ 
      
      if(delta[3] == 1){ # cx
        param_new$param3$b_bar <- 0
      } else { # bx
        param_new$param3$bx <- rep(0, X)
      }
      
    } 
    
    # from k1t (1 + ..) to k1t
    if(delta[2] == 2 & delta_star[2] == 1){ 
        
      if(delta[3] == 1){ # cx
        baseTerm <- 
          matrix(param2$k1t, X, Y, byrow = T) * 
          (1 + matrix(param3$b_bar * cx, X, Y, byrow = F))
      } else { # bx
        baseTerm <- 
          matrix(param2$k1t, X, Y, byrow = T) * 
          (1 + matrix(param3$bx, X, Y, byrow = F))
      }
      
      param_new$param2$k1t <- apply(baseTerm, 2, mean)
      
      if(delta[4] == 2){
        # we now need to enforce the constraint gamma_1 = 0
        
        # TO DO
        
      }
      
    } 
    
    # adding k2t 
    if(delta[2] == 2 & delta_star[2] == 3){
      
      if(delta[3] == 1) { # cx
        
        param_new$param3$k2t <- param2$k1t * param3$b_bar   
        
      } else if(delta[3] == 2) { # bx
        
        if(delta[1] == 2){
          
          # sum k1t is not 0 so we need to enforce that now
          
          param_new$param2$k1t <- param$param2$k1t
          meank1t <- mean(param$param2$k1t)
          param_new$param2$k1t <- param_new$param2$k1t - meank1t
          param_new$param1$ab[1] <- param_new$param1$ab[1] + meank1t
          
        } 
        
        # enforcing b_1 = -1
        
        param_new$param3$k2t <- param2$k1t
        param_new$param3$bx <- param$param3$bx
        bm1 <- param_new$param3$bx[1]
        param_new$param3$bx <- param_new$param3$bx / (-bm1)
        param_new$param3$k2t <- param_new$param3$k2t * (-bm1)
        
      }
      
      if(delta[4] == 2){
        # we now need to enforce the constraint gamma_1 = 0
        
        # TO DO 
        
        
      }
      
    } 
    
    if(delta[2] == 3 & delta_star[2] == 2){
      # removing k2t
      
      c_bar <- mean(param$param3$k2t  / param$param2$k1t)
      param_new$param2$k1t <-  (param$param2$k1t + param$param3$k2t / c_bar) / 2
      
      if(delta[3] == 2){
        # bx
        
        param_new$param3$bx <- param$param3$bx * c_bar
        
      } else {
        # cx
        
        param_new$param3$b_bar <- c_bar
        param_new$param3$k2t <- NULL
         
      }
      
      if(delta[4] == 2){
        # well in theory now you don't have anymore the constraint gamma1 = 0 
      }   
    }
    
  }
  
  # changing the third term
  if(delta[3] != delta_star[3]){  
    
    if(delta_star[3] == 2){ # going from cx to bx 
      
      if(delta[2] == 2){
        bx_tilde <- param3$b_bar * cx  
        param_new$param3$b_bar <- NULL # would be better to set it to NULL when needed
      } else {
        bx_tilde <- cx
      }
      
      param_new$param3$bx <- bx_tilde
      
    } else if (delta_star[3] == 1){ # going from bx to cx
      
      ab_tilde <- lm(param3$bx ~ cx)
      ab_tilde <- as.numeric(ab_tilde$coefficients)
      
      if(delta[2] == 2){ 
        # approximate k1t bx with k1t (1 + b cx)
        
        param_new$param2$k1t <- param$param2$k1t * (1 + ab_tilde[1])
        param_new$param3$b_bar <- ab_tilde[2] / (1 + ab_tilde[1])
        
        if(delta[1] == 2){ 
          # we need to enforce sum k1t = 0
          
          k1t_bar <- mean(param$param2$k1t)
          
          param_new$param2$k1t <- param$param2$k1t - k1t_bar
          param_new$param1$ab[1] <- param_new$param1$ab[1] + k1t_bar
          param_new$param1$ab[2] <- 
            param_new$param1$ab[2] + k1t_bar * param_new$param3$b_bar
          
        }
        
      } else if (delta[2] == 3){
        # approximate k2t bx with k2t cx
        
        param_new$param2$k1t <- param_new$param2$k1t + 
          param_new$param3$k2t * ab_tilde[1]
        
        param_new$param3$k2t <- param_new$param3$k2t * ab_tilde[2]
        
        if(delta[1] == 2){
          # in this case we also need to enforce sum k2t = 0
          
          k2t_bar <- mean(param_new$param3$k2t)
          
          param_new$param3$k2t <- param_new$param3$k2t - k2t_bar
          param_new$param1$ab[2] <- param_new$param1$ab[2] + k2t_bar
          
          
        }
        
      }
      
    } 
    
    # param_new$param3 <- param3
    
  }
  
  # changing the fouth term
  if(delta[4] != delta_star[4]){
    
    if(delta_star[4] == 2){ # extending the model
      
      param_new$param4 <- list("gtx" = rep(0, X + Y - 1)   )
      
    } else if(delta_star[4] == 1){ # reducing the model
      
      param_new$param4 <- NULL
      
    }
    
  }
  
  return(param_new)
  
}

proposeNewParams <- function(param, delta, idx_delta, 
                             d, E){
  
  if(idx_delta == 1){
    
    # compute current value of mxt
    cxt <- computemxt(idx_delta = 1,
                      delta, 
                      param$param1, 
                      param$param2, 
                      param$param3, 
                      param$param4)
    
    # proposing ax from a + b c_x
    
    list_proposal <- buildProposalAx(d, E, cxt)
    ax_star <- list_proposal$ax_star
    sd_star <- list_proposal$sd_star
    
    ax_proposal <- rnorm(X, ax_star, sd_star)
      
    logproposal <- sum(dnorm(ax_proposal, ax_star, 
                             sd_star, log = T))
    
    param_proposed_param1 <- list("ax" = ax_proposal)
    
    param_proposed <- list("param1" = param_proposed_param1,
                           "param2" = param$param2,
                           "param3" = param$param3,
                           "param4" = param$param4)
    
    list_return <- list("param" = param_proposed,
                        "logproposal" = logproposal)
    
  }
  
  if(idx_delta == 2){
    
    if(delta[2] == 2 & delta[3] == 1){
      # proposing b_bar and k1t
      
      cxt <- computemxt(idx_delta = 3,
                        delta, 
                        param$param1, 
                        param$param2, 
                        param$param3, 
                        param$param4)
      
      ktcx <- 
        matrix(param$param2$k1t, X, Y, byrow = T) *
        matrix(cx, X, Y, byrow = F)
      
      if(any(ktcx != 0)){
        
        list_proposal <- findProposalBBar(0, ktcx, cxt, d, E)
        bbar_star <- list_proposal$b_bar
        sd_star <- list_proposal$sd_star
        
      } else {
        
        bbar_star <- 0
        sd_star <- .0001
        
      }
      
      bbar_proposal <- rnorm(1, bbar_star, sd_star)
      
      logproposal <- dnorm(bbar_proposal, bbar_star, sd_star)
      bx_current <- (1 + bbar_proposal * cx)
      
      param$param3$b_bar <- bbar_proposal
      
      cxt <- computemxt(idx_delta = c(2,3),
                        delta, 
                        param$param1, 
                        param$param2, 
                        param$param3, 
                        param$param4)
      
      sum0 <- !(delta[1] == 2 & 
                  delta[3] == 2)
      
      list_proposal <- findProposalkt_given_bx(param$param2$k1t, 
                                               bx_current, 
                                               d, E, cxt, sum0)
      kt_star <- list_proposal$kt_star
      Sigma_kt_star <- list_proposal$Sigma_star
      
      if(sum0){
        
        kt_proposed_m1 <- mvrnorm(1, kt_star,
                               Sigma_kt_star)
        
        k1t_proposed <- c(kt_proposed_m1, - sum(kt_proposed_m1))
        
        logproposal <- logproposal + dmvnorm(kt_proposed_m1, kt_star, 
                                        Sigma_kt_star, T)
        
      } else {
        
        k1t_proposed <- mvrnorm(1, kt_star,
                            Sigma_kt_star)
        
        logproposal <- logproposal + dmvnorm(k1t_proposed, kt_star, 
                               Sigma_kt_star, T)
        
      }
      
      param_proposed <- list("param1" = param$param1,
                             "param2" = list("k1t" = k1t_proposed),
                             "param3" = param$param3,
                             "param4" = param$param4)

      list_return <- list("param" = param_proposed,
                          "logproposal" = logproposal)
      
    } 
    
    if(delta[2] == 3){
      
      # propose k2t and k1t
      
      cxt <- computemxt(idx_delta = 3,
                        delta, 
                        param$param1, 
                        param$param2, 
                        param$param3, 
                        param$param4)
      
      if(delta[3] == 1){
        bx_current <- cx
      } else if(delta[3] == 2){
        bx_current <- param$param3$bx
      }
      
      sum0 <- !(delta[1] == 2 & 
                  delta[3] == 2)
      
      list_proposal <- findProposalkt_given_bx(param$param3$k2t, 
                                               bx_current, d, E, cxt, sum0)
      kt_star <- list_proposal$kt_star
      Sigma_kt_star <- list_proposal$Sigma_star
      
      if(sum0){
        
        k2t_proposed_m1 <- mrt2(kt_star,
                                Sigma_kt_star,
                                df = 3)
        
        k2t_proposed <- c(k2t_proposed_m1, - sum(k2t_proposed_m1))
        
        logproposal <- dmt_cpp(k2t_proposed_m1, nu = 3, kt_star, 
                               Sigma_kt_star, returnLog = T)
        
      } else {
        
        k2t_proposed <- mrt2(kt_star,
                             Sigma_kt_star,
                             df = 3)
        
        logproposal <- dmt_cpp(k2t_proposed, nu = 3, kt_star, 
                               Sigma_kt_star, returnLog = T)
        
      }
      
      param$param3$k2t <- k2t_proposed
        
      cxt <- computemxt(idx_delta = 2,
                        delta, 
                        param$param1, 
                        param$param2, 
                        param$param3, 
                        param$param4)
      
      list_proposal <- buildProposalK1t_m1(param$param2$k1t, d, E, cxt)
      km1_star <- list_proposal$ktm1_star
      Sigma_star <- list_proposal$Sigma_star
      
      k1t_proposed_m1 <- mrt2(km1_star, Sigma_star, df = 3)
      k1t_proposed <- c(k1t_proposed_m1, - sum(k1t_proposed_m1))
      
      logproposal <- logproposal + dmt_cpp(k1t_proposed_m1, nu = 3, km1_star, 
                                  Sigma_star, T)
      
      param_proposed <- list("param1" = param$param1,
                             "param2" = list("k1t" = k1t_proposed),
                             "param3" = param$param3,
                             "param4" = param$param4)
      
      list_return <- list("param" = param_proposed,
                          "logproposal" = logproposal)
      
      
    }
    
  } 
  
  if(idx_delta == 3){
    
    cxt <- computemxt(idx_delta = 3, delta, 
                      param$param1, param$param2,
                      param$param3, param$param4)
    
    if(delta[2] == 2){
      kt_current <- param$param2$k1t  
    } else if(delta[2] == 3){
      kt_current <- param$param3$k2t
    }
    
    if(delta[2] == 2){ # no constraint that b_1 = -1
      
      startVal <- param$param3$bx
      list_proposal <- findProposalBX_m1(startVal, kt_current, d, E, cxt, sd_bx)
      bx_star <- list_proposal$bx_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      bx_proposed_m <- mrt2(bx_star, Sigma_bx_star, df = 3)
      
      bx_proposed <- c(bx_proposed_m, - sum(bx_proposed_m))
      
    } else if(delta[2] == 3){
      
      startVal <- cx
      list_proposal <- findProposalBX_m2(startVal, kt_current, d, E, cxt, sd_bx)
      bx_star <- list_proposal$bx_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      bx_proposed_m <- mrt2(bx_star, Sigma_bx_star, df = 3)
      
      bx_proposed <- c(-1, bx_proposed_m, - sum(bx_proposed_m) + 1)
      
    }

    logproposal <- dmt_cpp(bx_proposed_m, nu = 3,
                                    bx_star, Sigma_bx_star, T)
    
    param$param3$bx <- bx_proposed
    
    param_proposed <- list("param1" = param$param1,
                           "param2" = param$param2,
                           "param3" = param$param3,
                           "param4" = param$param4)
    
    list_return <- list("param" = param_proposed,
                        "logproposal" = logproposal)
        
  } 
  
  if(idx_delta == 4){
    
    cxt <- computemxt(idx_delta = 4, delta, 
                      param$param1, param$param2,
                      param$param3, param$param4)
    
    if(delta[2] == 2){ 
      # no additive effect on k1t and therefore
      # no constraint gamma1 = 0
      
      list_proposal <- buildProposalGtx_m1(rep(0, X+ Y - 1),
                                           d, E, cxt)
      gtx_star <- list_proposal$gtx_star
      Sigma_star <- list_proposal$Sigma_star
      
      gtx_proposed_pm1 <- mvrnorm(n = 1, gtx_star, Sigma_star)
      
      ik <- create_ik(X, Y)
      
      gtx_proposed <- c(gtx_proposed_pm1, 
                        - sum(ik[-length(ik)] * gtx_proposed_pm1))
      
    } else { 
      # additive effect on k1t and therefore
      # constraint gamma1 = 0
      
      list_proposal <- buildProposalGtx_m2(rep(0, X+ Y - 1),
                                           d, E, cxt)
      gtx_star <- list_proposal$gtx_star
      Sigma_star <- list_proposal$Sigma_star
      
      gtx_proposed_pm1 <- mvrnorm(n = 1, gtx_star, Sigma_star)
      
      ik <- create_ik(X, Y)
      
      gtx_proposed <- c(0, gtx_proposed_pm1, 
                        - sum(ik[-c(1,length(ik))] * gtx_proposed_pm1))
      
    }
    
    logproposal <- dmvnorm(gtx_proposed_pm1, gtx_star, 
                           Sigma_star, log = T)
    
    ik <- create_ik(X, Y)
    
    param4 <- list("gtx" = gtx_proposed)
    
    param_proposed <- list("param1" = param$param1,
                           "param2" = param$param2,
                           "param3" = param$param3,
                           "param4" = param4)
    
    list_return <- list("param" = param_proposed,
                        "logproposal" = logproposal)
    
  }
  
  return(list_return)
  
}

computeProposalProb <- function(param_prop, param_base,
                                d, E, idx_delta, delta){

  if(idx_delta == 1){
    
    cxt <- computemxt(idx_delta = 1, delta, 
                      param_base$param1, 
                      param_base$param2, 
                      param_base$param3, 
                      param_base$param4)
    
    # proposing ax from a + b c_x
    
    list_proposal <- buildProposalAx(d, E, cxt)
    ax_star <- list_proposal$ax_star
    sd_star <- list_proposal$sd_star
    
    logproposal <- sum(dnorm(param_prop$param1$ax, ax_star, 
                             sd_star, log = T))
    
  } 
  
  if(idx_delta == 2){
    
    if(delta[2] == 2 & delta[3] == 1){
      # proposing b_bar and k1t
      
      cxt <- computemxt(idx_delta = 3,
                        delta, 
                        param_base$param1, 
                        param_base$param2, 
                        param_base$param3, 
                        param_base$param4)
      
      ktcx <- 
        matrix(param_base$param2$k1t, X, Y, byrow = T) *
        matrix(cx, X, Y, byrow = F)
      
      list_proposal <- findProposalBBar(0, ktcx, cxt, d, E)
      bbar_star <- list_proposal$b_bar
      sd_star <- list_proposal$sd_star
      
      logproposal <- dnorm(param_prop$param3$b_bar, bbar_star, sd_star)
      bx_current <- (1 + param_prop$param3$b_bar * cx)
      
      param_base$param3$b_bar <- param_prop$param3$b_bar
      
      cxt <- computemxt(idx_delta = c(2,3),
                        delta, 
                        param_base$param1, 
                        param_base$param2, 
                        param_base$param3, 
                        param_base$param4)
      
      sum0 <- !(delta[1] == 2 & 
                  delta[3] == 2)
      
      list_proposal <- findProposalkt_given_bx(rep(0, Y), 
                                               bx_current, 
                                               d, E, cxt, sum0)
      kt_star <- list_proposal$kt_star
      Sigma_kt_star <- list_proposal$Sigma_star
      
      if(sum0){
        
        logproposal <- logproposal + dmvnorm(param_prop$param2$k1t[-Y], 
                                             kt_star, 
                                             Sigma_kt_star, T)
        
      } else {
        
        logproposal <- logproposal + dmvnorm(param_prop$param2$k1t, 
                                             kt_star, 
                                             Sigma_kt_star, T)
        
      }
      
      
    }
    
    if(delta[2] == 3){
      
      # propose k2t and k1t
      
      cxt <- computemxt(idx_delta = 3,
                        delta, 
                        param_base$param1, 
                        param_base$param2, 
                        param_base$param3, 
                        param_base$param4)
      
      if(delta[3] == 1){
        bx_current <- cx
      } else if(delta[3] == 2){
        bx_current <- param_base$param3$bx
      }
      
      sum0 <- !(delta[1] == 2 & 
                  delta[3] == 2)
      
      list_proposal <- findProposalkt_given_bx(param_base$param3$k2t, 
                                               bx_current, d, E, cxt, sum0)
      kt_star <- list_proposal$kt_star
      Sigma_kt_star <- list_proposal$Sigma_star
      
      if(sum0){
        
        logproposal <- dmt_cpp(param_prop$param3$k2t[-Y], nu = 3, kt_star, 
                               Sigma_kt_star, returnLog = T)
        
      } else {
        
        logproposal <- dmt_cpp(param_prop$param3$k2t, nu = 3, kt_star, 
                               Sigma_kt_star, returnLog = T)
        
      }
      
      cxt <- computemxt(idx_delta = 2,
                        delta, 
                        param_prop$param1, 
                        param_base$param2, 
                        param_prop$param3, 
                        param_prop$param4)
      
      list_proposal <- buildProposalK1t_m1(param_base$param2$k1t, d, E, cxt)
      km1_star <- list_proposal$ktm1_star
      Sigma_star <- list_proposal$Sigma_star
      
      logproposal <- logproposal + dmt_cpp(param_prop$param2$k1t[-Y], nu = 3, 
                                           km1_star, Sigma_star, T)
      
      
    }
    
  } 
  
  if(idx_delta == 3){
    
    cxt <- computemxt(idx_delta = 3, delta, 
                      param$param1, param$param2,
                      param$param3, param$param4)
    
    if(delta[2] == 2){
      kt_current <- param$param2$k1t  
    } else if(delta[2] == 3){
      kt_current <- param$param3$k2t
    }
    
    if(delta[2] == 2){ # no constraint that b_1 = -1
      
      # startVal <- param$param3$b_bar * cx
      list_proposal <- findProposalBX_m1(rep(0,X), kt_current, d, E, cxt, sd_bx)
      bx_star <- list_proposal$bx_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      logproposal <- dmt_cpp(param_prop$param3$bx[-X], nu = 3,
                             bx_star, Sigma_bx_star, T)
      
    } else if(delta[2] == 3){
      
      list_proposal <- findProposalBX_m2(rep(0,X), kt_current, d, E, cxt, sd_bx)
      bx_star <- list_proposal$bx_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      logproposal <- dmt_cpp(param_prop$param3$bx[-c(1,X)], nu = 3,
                             bx_star, Sigma_bx_star, T)
      
    }
    
  } 
  
  if(idx_delta == 4){
    
    # compute current value of mxt
    cxt <- computemxt(idx_delta = 4, 
                      delta, 
                      param_base$param1, 
                      param_base$param2, 
                      param_base$param3, 
                      param_base$param4)
    
    proposedParams <- 
      param_prop$param4$gtx - 
      param_base$param4$gtx
    
    if(delta[2] == 2){ 
      # no additive effect on k1t and therefore
      # no constraint gamma1 = 0
      
      list_proposal <- buildProposalGtx_m1(param_base$param4$gtx,
                                           d, E, cxt)
      gtx_star <- list_proposal$gtx_star
      Sigma_star <- list_proposal$Sigma_star
      
      logproposal <- dmvnorm(proposedParams[-c(X+Y-1)], 
                             gtx_star, 
                             Sigma_star, log = T)
      
      # gtx_proposed_pm1 <- mvrnorm(n = 1, gtx_star, Sigma_star)
      # 
      # ik <- create_ik(X, Y)
      # 
      # gtx_proposed <- c(gtx_proposed_pm1, 
      #                   - sum(ik[-length(ik)] * gtx_proposed_pm1))
      
    } else { 
      # additive effect on k1t and therefore
      # constraint gamma1 = 0
      
      list_proposal <- buildProposalGtx_m2(param_base$param4$gtx, d, E, cxt)
      gtx_star <- list_proposal$gtx_star
      Sigma_star <- list_proposal$Sigma_star
      
      logproposal <- dmvnorm(proposedParams[-c(1,X+Y-1)], 
                             gtx_star, 
                             Sigma_star, log = T)
      
    }
    
    
  }
  
  return(logproposal)
}

computePrior <- function(param, delta, priorParams){
  
  logprior <- 0
  
  if(delta[1] == 1){
    
    logprior <- logprior + sum(dnorm(
      param$param1$ax, 0, sd = priorParams$sd_a, log = T
    ))
    
  }
  
  logprior <- logprior + sum(dnorm(
    param$param2$k1t, 0, sd = priorParams$sd_k, log = T
  ))
  
  if(delta[2] == 2 & delta[3] == 1){
    
    logprior <- logprior + dnorm(
      param$param3$b_bar, 0, sd = priorParams$sd_b, log = T
    )
    
  }
  
  if(delta[2] == 3){
    
    logprior <- logprior + sum(dnorm(
      param$param3$k2t, 0, sd = priorParams$sd_k, log = T
    ))
    
  }
  
  if(delta[3] == 2){
    
    logprior <- logprior + sum(dnorm(
      param$param3$bx, 0, sd = priorParams$sd_b, log = T
    ))
    
  }
  
  if(delta[4] == 2){
    
    logprior <- logprior + sum(dnorm(
      param$param4$gtx, 0, sd = priorParams$sd_g, log = T
    ))
    
  }
  
  logprior
}

proposeNewState <- function(idx_delta, delta,
                            param, d, E){
  
  # propose new index
  delta_star <- proposeNewDelta(idx_delta, delta)
  
  extendingModel <- isModelBeingExtended(delta, delta_star)
  
  # map current params to the new state
  param_tilde <- mapModels(delta, delta_star, 
                           param$param1, param$param2, 
                           param$param3, param$param4)
  
  if(extendingModel){ # proposing a more complex model
    
    # propose new params
    list_proposal <- proposeNewParams(
      param_tilde, 
      delta_star,# propose under the extended model
      idx_delta, d, E)
    param_star <- list_proposal$param
    logproposal <- list_proposal$logproposal
    
    param_base <- param
    
    m_base <- computemxt(idx_delta = 0 ,
                         delta, 
                         param$param1, param$param2,
                         param$param3, param$param4)
    
    m_star <- computemxt(idx_delta = 0 ,
                         delta_star, 
                         param_star$param1, param_star$param2,
                         param_star$param3, param_star$param4)
    
    logprior_base <- computePrior(param, delta, priorParams)
    logprior_star <- computePrior(param_star, delta_star, priorParams)
    
  } else { # proposing a simpler model
    
    param_tilde2 <- mapModels(delta_star, delta, 
                              param_tilde$param1, param_tilde$param2,
                              param_tilde$param3, param_tilde$param4)
    
    logproposal <- computeProposalProb(
      param, param_tilde2,
      d, E, idx_delta, 
      delta # the extended model is the current model
    )
    
    param_star <- param # the larger model is the current model
    param_base <- param_tilde # the smaller model is the new model
    
    m_star <- computemxt(idx_delta = 0 ,
                         delta, 
                         param_star$param1, param_star$param2,
                         param_star$param3, param_star$param4)
    
    m_base <- computemxt(idx_delta = 0,
                         delta_star, 
                         param_base$param1, param_base$param2,
                         param_base$param3, param_base$param4)
    
    logprior_star <- computePrior(param_star, delta, priorParams)
    logprior_base <- computePrior(param_base, delta_star, priorParams)
    
  }
  
  loglik_star <- loglik(d, E, m_star)
  loglik_base <- loglik(d, E, m_base)
  
  mh_ratio <- exp(
    loglik_star - loglik_base +
      logprior_star - logprior_base -
      logproposal
  )
  
  if(!extendingModel){
    
    mh_ratio <- 1 / mh_ratio
    
  }
  
  if(runif(1) < mh_ratio){
    
    delta <- delta_star
    
    if(extendingModel){
      
      param <- param_star
      
    } else {
      
      param <- param_tilde
      
    }
    
  } 
  
  list("param" = param,
       "delta" = delta)
}

############ OLD -------

plotCurves <- function(term1, term2, term3, term4){
  
  data <- expand.grid(x = 1:X, t = 1:Y)

  mxtp <- apply(data, 1, function(dat){
    x <- dat[1]
    t <- dat[2]
    term1[x,t] + term2[x,t] + term3[x,t] + term4[x,t]
  })
  
  df <- data.frame(Age = data$x,
                   Year = factor(data$t),
                   m = mxtp)
  
}


computeTerms <- function(idx, term1, term2, term3, term4){
  
  sumterms <- term1 + term2 + term3 + term4
  
  if(idx == 1){
    sumterms <- sumterms - term1 
  } else if(idx == 2){
    sumterms <- sumterms - term2 
  } else if(idx == 3){
    sumterms <- sumterms - term3 
  } else if(idx == 4){
    sumterms <- sumterms - term4
  } 
  
  sumterms
}

convertK2tBxtoK1tK2txBx <- function(k2t, bx){
  
  meanb <- mean(bx)
  k1t_tilde <- k2t * meanb
  
  bx_hat <- bx - mean(bx)
  
  bx_tilde <- bx_hat / (-bx_hat[1])
  
  k2t_tilde <- k2t * (-bx_hat[1])
  
  list("k1t_tilde" = k1t_tilde,
       "k2t_tilde" = k2t_tilde,
       "bx_tilde" = bx_tilde)  
}

convertK1tK2txBxtoK2tBx <- function(k1t, k2t, bx){
  
  # kbar <-  - sum(bx) / (X - sum(bx))
  
  # bx_tilde_new <- kbar + bx * (1 + kbar)
  
  # mean(bx_old) + bx * (1 + mean(bx_old))
  
  b_mean <- mean( - k1t / (k1t - k2t) , na.rm = T)
  
  k1t_hat <- k1t / b_mean
  k2t_hat <- k2t / (1 + b_mean)
  
  k2t_tilde <- (k1t_hat + k2t_hat) / 2
  
  bx_tilde <- b_mean + bx * (1 + b_mean)
  
  list("k2t_tilde" = k2t_tilde,
       "bx_tilde" = bx_tilde)  
}

reparamK1tK2txBx <- function(k1t, k2t, bx){
  
  b_bar <- mean(bx)
  k1t_tilde <- k1t + k2t * b_bar
  
  k2t_tilde <- k2t * (- (bx[1]- b_bar))
  bx_tilde <- (bx - b_bar) / (- (bx[1]- b_bar))
  
  list("bx_tilde" = bx_tilde,
       "k1t_tilde" = k1t_tilde,
       "k2t_tilde" = k2t_tilde)
}

diagnosticsCheck <- function(params_output){
  
  apply(params_output, c(2,3), function(x){
    
    if(any(!is.na(x))){
      x <- x[!is.na(x)]
      x_current <- mcmc(x)
      return(as.numeric(effectiveSize(x_current)))
      
    } else {
      return(NA)
    }
  })
  
}

loglik <- function(d, E, m){
  
  mxt <- m + log(E)
  
  sum(d * mxt - exp(mxt))
  
}

loglik_m <- function(d, E, term1, term2, term3, term4){
  
  mxt <- term1 + term2 + term3 + term4 + log(E)
  
  sum(d * mxt - exp(mxt))
  
}

# UPDATE TERM 1 ------------

loglik_term1_a <- function(a, x, d, E, cxt){
  
  mxt <- a + cxt[x,] + log(E[x,])
  
  sum(d[x,] * mxt - exp(mxt))
  
}

buildProposalAX <- function(d, E, cxt) {
  
  ax_star_all <- rep(0, X)
  sd_ax_star <- rep(0, X)
  
  for (x in 1:X) {
    
    ax_star <- log(
      sum(d[x,]) / sum(exp(cxt[x,] + log(E[x,])))
    )
    
    hessian_ax_star <- derl2der2a(ax_star, cxt[x,] + log(E[x,]))
    sd_star <- sqrt(- 1 / hessian_ax_star)
    
    ax_star_all[x] <- ax_star
    sd_ax_star[x] <- sd_star
    
  }
  
  list("ax_star_all" = ax_star_all,
       "sd_ax_star" = sd_ax_star)
  
}



loglik_term1_ab_fun <- function(ab, c2x, d, E, cxt){

    a <- ab[1]
    b <- ab[2]
    
    abx <- matrix(a + b * c2x, X, Y, byrow = F)
    
    loglik_term1(d, E, abx, cxt)
    
    # mxt <- abx + cxt + log(E)
    # 
    # - sum(d * mxt - exp(mxt))

}

loglik_term1_ax_fun <- function(ax, d, E, cxt){
  
  abx <- matrix(ax, X, Y, byrow = F)
  
  loglik_term1(d, E, abx, cxt)
  
  # mxt <- abx + cxt + log(E)
  # 
  # sum(d * mxt - exp(mxt))
  
}

loglik_term1_ab <- function(c2x, d, E, cxt){
  
  function(param){
    
    a <- param[1]
    b <- param[2]
    
    abx <- matrix(a + b * c2x, X, Y, byrow = F)
    
    loglik_term1(d, E, abx, cxt)
    
  }
  
}

gr_loglik_term1_ab <- function(c2x, d, E, cxt){
  
  function(param){
    
    a <- param[1]
    b <- param[2]
    
    abx <- matrix(a + b * c2x, X, Y, byrow = F)
    
    mxt <- abx + cxt + log(E)
    
    gr_a <- sum(d - exp(mxt))
    
    term2 <- (d - exp(mxt)) * matrix(c2x, X, Y, byrow = F)
    
    gr_b <- sum(term2)
    
    - c(gr_a, gr_b)
    
  }
  
}

findProposalABX <- function(c2x, d, E, cxt){
  
  loglik_ab_current <- loglik_term1_ab(c2x, d, E, cxt)
  gr_loglik_ab_current <- gr_loglik_term1_ab(c2x, d, E, cxt)
  
  # find maximum
  laplace_fit <- optim(
    # par = ab,
    par = c(0, 0),
    fn = loglik_ab_current,
    method = c("BFGS"),
    gr = gr_loglik_ab_current,
    hessian = T
  )
  
  ab_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  list("ab_star" = ab_star,
       "Sigma_star" = Sigma_star)
  
}


# UPDATE TERM 2 -----
  

updateTerm2 <- function(param, delta2, term1, term2, term3, term4){
  
  if(delta2 == 2){
    
    term134 <- computeTerms(idx = 2, term1, term2, term3, term4)
    
    k1t <- updateK1T(param$k1t, d, E, term134)
    term2 <- matrix(k1t, X, Y, byrow = T)
    
    param <- list("k1t" = k1t)
    
  } else {
   
    param <- NULL
     
  }
  
  return(
    list("param" = param,
         "term" = term2)
  )
  
}

# UPDATE TERM 3 -------

loglik_term3_k2t_m1 <- function(bx, d, E, cxt){
  
  function(param){
    
    k2t <- c(param, - sum(param))
    
    k2tx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tx, cxt)
    
  }
  
}

findProposalk2t_m1 <- function(bx, d, E, cxt){
  
  loglik_k2t_current <- loglik_term3_k2t_m1(bx, d, E, cxt)
  # gr_loglik_ab_current <- gr_loglik_abx_r(x, d, E, a, b, ap, bp, k2)
  
  # find maximum
  laplace_fit <- optim(
    par = rep(0, Y - 1),
    # par = k2t,
    fn = loglik_k2t_current,
    method = c("BFGS"),
    # gr = gr_loglik_ab_current,
    hessian = T
  )
  
  k2t_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  Sigma_star <- Sigma_star
  
  list("k2t_star" = k2t_star,
       "Sigma_star" = Sigma_star)
  
}

loglik_term3_k2t <- function(bx, d, E, cxt){
  
  function(param){
    
    k2t <- param
    
    k2tx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tx, cxt)
    
  }
  
}

findProposalk2t <- function(bx, d, E, cxt){
  
  loglik_k2t_current <- loglik_term3_k2t(bx, d, E, cxt)
  # gr_loglik_ab_current <- gr_loglik_abx_r(x, d, E, a, b, ap, bp, k2)
  
  # find maximum
  laplace_fit <- optim(
    par = rep(0, Y),
    # par = k2t,
    fn = loglik_k2t_current,
    method = c("BFGS"),
    # gr = gr_loglik_ab_current,
    hessian = T
  )
  
  k2t_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  Sigma_star <- Sigma_star
  
  list("k2t_star" = k2t_star,
       "Sigma_star" = Sigma_star)
  
}

update_k2t <- function(k2t, d, E, cxt, bx){
  
  list_proposal_k2t <- findProposalk2t(bx, d, E, cxt)
  k2t_star <- list_proposal_k2t$k2t_star
  Sigma_star <- list_proposal_k2t$Sigma_star
  
  k2t_proposal <- sampleMVNconstraint_k(k2t_star, Sigma_star)
  
  loglik_k2t_current <- loglik_term3_k2t(bx, d, E, cxt)
  
  loglik_proposal <- - loglik_k2t_current(k2t_proposal)
  loglik_current <- - loglik_k2t_current(k2t)
  
  logproposal_proposal <- dmvnorm(as.vector(k2t_proposal), k2t_star, Sigma_star, log = T)
  logproposal_current <- dmvnorm(as.vector(k2t), k2t_star, Sigma_star, log = T)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  mh_ratio
  
  if(runif(1) < mh_ratio){
    
    k2t <- k2t_proposal
    
  }
  
  term3 <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
  
  list("k2t" = k2t,
       "term3" = term3)
}

sampleMVNconstraint_b <- function(mu, Sigma){
  
  d <- length(mu)
  
  A <- matrix(1, 1, d)
  b <- 1
  
  P <- matrix(0, d, d)
  diag(P)[1:d - 1] <- 1
  P <- P[-d,]
  
  C <- rbind(P, A)
  
  nrow1 <- 1:(d - nrow(A))
  nrow2 <- d - nrow(A) + 1:nrow(A)
  
  nu <- C %*% mu
  Omega <- C %*% Sigma %*% t(C)
  Omega11 <- Omega[nrow1, nrow1]
  Omega22 <- Omega[nrow2, nrow2]
  Omega12 <- Omega[nrow1, nrow2]
  Omega21 <- Omega[nrow2, nrow1]
  
  mu1 <- nu[nrow1] + Omega12 %*% solve(Omega22) %*% (b - nu[nrow2])
  Sigma1 <- Omega11 - Omega12 %*% solve(Omega22) %*% Omega21
  
  pxcd <- mvrnorm(1, mu1, Sigma1)
  
  x <- solve(C) %*% c(pxcd, b)
  
  x
}

loglik_term3_bx <- function(k2t, d, E, cxt){
  
  function(param){
    
    bx <- param
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

findProposalbx <- function(k2t, d, E, cxt){
  
  loglik_bx_current <- loglik_term3_bx(k2t, d, E, cxt)
  # gr_loglik_ab_current <- gr_loglik_abx_r(x, d, E, a, b, ap, bp, k2)
  
  # find maximum
  laplace_fit <- optim(
    par = rep(0, X),
    # par = k2t,
    fn = loglik_bx_current,
    method = c("BFGS"),
    # gr = gr_loglik_ab_current,
    hessian = T
  )
  
  bx_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  # Sigma_star <- 2 * Sigma_star
  
  list("bx_star" = bx_star,
       "Sigma_star" = Sigma_star)
  
}

update_bx <- function(bx, d, E, cxt, k2t){
  
  list_proposal_bx <- findProposalbx(k2t, d, E, cxt)
  bx_star <- list_proposal_bx$bx_star
  Sigma_star <- list_proposal_bx$Sigma_star
  
  bx_proposal <- sampleMVNconstraint_b(bx_star, Sigma_star)
  
  loglik_bx_current <- loglik_term3_bx(k2t, d, E, cxt)
  
  loglik_proposal <- - loglik_bx_current(bx_proposal)
  loglik_current <- - loglik_bx_current(bx)
  
  logproposal_proposal <- dmvnorm(as.vector(bx_proposal), bx_star, Sigma_star, log = T)
  logproposal_current <- dmvnorm(as.vector(bx), bx_star, Sigma_star, log = T)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    bx <- bx_proposal
    
  }
  
  term3 <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
  
  list("bx" = bx,
       "term3" = term3)
}

loglik_term3_bxk2t <- function(d, E, cxt){
  
  function(param){
    
    bx <- param[1:X]
    k2t <- param[X + 1:Y]
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

loglik_term3_bxk2t_m1 <- function(d, E, cxt){
  
  function(param){
    
    bxm1X <- param[1:(X-2)]
    k2tm1 <- param[X - 2 + 1:(Y-1)]

    bx <- c(-1, bxm1X, - sum(bxm1X) + 1)
    k2t <- c(k2tm1, - sum(k2tm1))
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

gr_loglik_term3_bxk2t_m1 <- function(d, E, cxt){
  
  function(param){
    
    bxm1X <- param[1:(X-2)]
    k2tm1 <- param[X-2 + 1:(Y-1)]

    bx <- c(-1, bxm1X, - sum(bxm1X) + 1)
    k2t <- c(k2tm1, - sum(k2tm1))
    
    k2t_mat <- matrix(k2t, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term1_bx <- d * k2t_mat - exp(lambda_xt) * k2t_mat
    
    term2_bx <- d[X,] * (-k2t_mat[X,]) - exp(lambda_xt[X,]) * (- k2t_mat[X,])
    
    gr_bx <- apply(term1_bx[2:(X-1),], 1, sum) + sum(term2_bx)
    
    term1_kt <- d * bx_mat - exp(lambda_xt) * bx_mat
    
    term2_kt <- (- d[,Y] * bx_mat[,Y] ) - exp(lambda_xt[,Y]) * (- bx_mat[,Y])
    
    gr_kt <- apply(term1_kt[,1:(Y-1)], 2, sum) + sum(term2_kt)
    
    - c(gr_bx, gr_kt)
  
  }
  
}

loglik_term3_bxk2t <- function(d, E, cxt){
  
  function(param){
    
    bx <- param[1:X]
    k2t <- param[X + 1:Y]
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

findProposalBXk2t <- function(d, E, cxt, c2x){
  
  loglik_bxk2t_current <- loglik_term3_bxk2t_m1(d, E, cxt)
  gr_loglik_bxk2t_current <- gr_loglik_term3_bxk2t_m1(d, E, cxt)
  
  startVal_b <- c2x
  
  # find maximum
  laplace_fit <- optim(
    par = c(startVal_b[-c(1,X)], rep(0, Y - 1)),
    # par = k2t,
    fn = loglik_bxk2t_current,
    method = c("BFGS"),
    gr = gr_loglik_bxk2t_current,
    hessian = T
  )
  
  bxl2t_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian + diag(1, nrow = length(bxl2t_star)))
  
  # Sigma_star <- Sigma_star + diag(.001, nrow = length(bxl2t_star))
  # Sigma_star <- 2 * Sigma_star
  
  list("bxl2t_star" = bxl2t_star,
       "Sigma_star" = Sigma_star)
  
}

update_bxk2t <- function(bx, k2t, d, E, cxt){
  
  list_proposal <- findProposalBXk2t(d, E, cxt, c2x)
  bxl2t_star <- list_proposal$bxl2t_star
  Sigma_bxl2t_star <- list_proposal$Sigma_star
  
  k2tbx_proposed_m1 <- mvrnorm(n = 1,
                               bxl2t_star,
                               Sigma_bxl2t_star)
  
  bx_m1X <- k2tbx_proposed_m1[1:(X-2)]
  k2t_mY <- k2tbx_proposed_m1[X-2 + 1:(Y-1)]
  
  bx_proposed <- c(-1, bx_m1X, - sum(bx_m1X) + 1)
  k2t_proposed <- c(k2t_mY, - sum(k2t_mY))

  loglik_proposal <- - loglik_term3(k2t_proposed, bx_proposed, d, E, cxt)
  loglik_current <- - loglik_term3(k2t, bx, d, E, cxt)
  
  k2tbx_current_m1 <- c(bx[-c(1,X)], k2t[-Y])
  
  logproposal_proposal <- dmvnorm(k2tbx_proposed_m1, bxl2t_star, Sigma_bxl2t_star, log = T)
  logproposal_current <- dmvnorm(k2tbx_current_m1, bxl2t_star, Sigma_bxl2t_star, log = T)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    bx <- bx_proposed
    k2t <- k2t_proposed
    
  }
  
  term3 <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
  
  list("bx" = bx,
       "k2t" = k2t,
       "term3" = term3)
    
}

loglik_term3_bxk2t_m1_delta2u1 <- function(d, E, cxt){
  
  function(param){
    
    bxm1X <- param[1:(X-1)]
    k2tm1 <- param[X - 1 + 1:(Y-1)]
    
    bx <- c(-1, bxm1X)
    k2t <- c(k2tm1, - sum(k2tm1))
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    loglik_term1(d, E, k2tbx, cxt)
    
  }
  
}

gr_loglik_term3_bxk2t_m1_delta2u1 <- function(d, E, cxt){
  
  function(param){
    
    bxm1X <- param[1:(X-1)]
    k2tm1 <- param[X-1 + 1:(Y-1)]
    
    bx <- c(-1, bxm1X)
    k2t <- c(k2tm1, - sum(k2tm1))
    
    k2t_mat <- matrix(k2t, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term1_bx <- d * k2t_mat - exp(lambda_xt) * k2t_mat
    
    gr_bx <- apply(term1_bx[2:X,], 1, sum)
    
    term1_kt <- d * bx_mat - exp(lambda_xt) * bx_mat
    
    term2_kt <- (- d[,Y] * bx_mat[,Y] ) - exp(lambda_xt[,Y]) * (- bx_mat[,Y])
    
    gr_kt <- apply(term1_kt[,1:(Y-1)], 2, sum) + sum(term2_kt)
    
    - c(gr_bx, gr_kt)
    
  }
  
}

findProposalBXk2t_delta2u1 <- function(d, E, cxt, c2x){
  
  loglik_bxk2t_current <- loglik_term3_bxk2t_m1_delta2u1(d, E, cxt)
  gr_loglik_bxk2t_current <- gr_loglik_term3_bxk2t_m1_delta2u1(d, E, cxt)
  
  startVal_b <- c2x
  
  # find maximum
  laplace_fit <- optim(
    par = c(startVal_b[-1], rep(0, Y - 1)),
    # par = k2t,
    fn = loglik_bxk2t_current,
    method = c("BFGS"),
    gr = gr_loglik_bxk2t_current,
    hessian = T
  )
  
  bxl2t_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian + diag(1, nrow = length(bxl2t_star)))
  
  # Sigma_star <- Sigma_star + diag(.001, nrow = length(bxl2t_star))
  # Sigma_star <- 2 * Sigma_star
  
  list("bxl2t_star" = bxl2t_star,
       "Sigma_star" = Sigma_star)
  
}

update_bxk2t_delta2u1 <- function(bx, k2t, d, E, cxt){
  
  list_proposal <- findProposalBXk2t(d, E, cxt, c2x)
  bxl2t_star <- list_proposal$bxl2t_star
  Sigma_bxl2t_star <- list_proposal$Sigma_star
  
  k2tbx_proposed_m1 <- mvrnorm(n = 1,
                               bxl2t_star,
                               Sigma_bxl2t_star)
  
  bx_m1X <- k2tbx_proposed_m1[1:(X-2)]
  k2t_mY <- k2tbx_proposed_m1[X-2 + 1:(Y-1)]
  
  bx_proposed <- c(-1, bx_m1X, - sum(bx_m1X) + 1)
  k2t_proposed <- c(k2t_mY, - sum(k2t_mY))

  loglik_proposal <- - loglik_term3(k2t_proposed, bx_proposed, d, E, cxt)
  loglik_current <- - loglik_term3(k2t, bx, d, E, cxt)
  
  k2tbx_current_m1 <- c(bx[-c(1,X)], k2t[-Y])
  
  logproposal_proposal <- dmvnorm(k2tbx_proposed_m1, bxl2t_star, Sigma_bxl2t_star, log = T)
  logproposal_current <- dmvnorm(k2tbx_current_m1, bxl2t_star, Sigma_bxl2t_star, log = T)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    bx <- bx_proposed
    k2t <- k2t_proposed
    
  }
  
  term3 <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
  
  list("bx" = bx,
       "k2t" = k2t,
       "term3" = term3)
    
}

# doing two steps

# for bx

loglik_term3_bx_k2t_m2 <- function(k2t, d, E, cxt, sd_bx){
  
  function(param){
    
    bxm1X <- param[1:(X-2)]
    
    bx <- c(-1, bxm1X, - sum(bxm1X) + 1)
    
    k2tbx <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
    
    priorTerm <- sum(
      dnorm(bxm1X, 0, sd = sd_bx, log = T)
      )
    
    loglik_term1(d, E, k2tbx, cxt) - priorTerm
    
  }
  
}

gr_loglik_term3_bx_k2t_m2 <- function(k2t, d, E, cxt, sd_bx){
  
  function(param){
    
    bxm1X <- param[1:(X-2)]
    bx <- c(-1, bxm1X, - sum(bxm1X) + 1)
    
    k2t_mat <- matrix(k2t, X, Y, byrow = T)
    bx_mat <- matrix(bx, X, Y, byrow = F)
    
    k2tbx <- k2t_mat * bx_mat
    
    cxtlogE <- cxt + log(E)
    
    lambda_xt <- cxtlogE + k2t_mat * bx_mat
    
    term1_bx <- d * k2t_mat - exp(lambda_xt) * k2t_mat
    
    term2_bx <- d[X,] * (-k2t_mat[X,]) - exp(lambda_xt[X,]) * (- k2t_mat[X,])
    
    gr_bx <- apply(term1_bx[2:(X-1),], 1, sum) + sum(term2_bx)
    
    - gr_bx + 1 / sd_bx^2
    
  }
  
}

findProposalBX_given_k2t_m2 <- function(bx, k2t, d, E, cxt, c2x, sd_bx){
  
  loglik_bx_k2t_current <- loglik_term3_bx_k2t_m2(k2t, d, E, cxt, sd_bx)
  gr_loglik_bx_k2t_current <- gr_loglik_term3_bx_k2t_m2(k2t, d, E, cxt, sd_bx)
  
  startVal_b <- c2x
  
  # find maximum
  laplace_fit <- optim(
    par = bx[-c(1,X)],
    # par = rep(0, X - 2),
    fn = loglik_bx_k2t_current,
    method = c("CG"),
    # gr = gr_loglik_bx_k2t_current,
    hessian = T
  )
  
  bxl2t_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  # Sigma_star <- Sigma_star + diag(.001, nrow = length(bxl2t_star))
  # Sigma_star <- 2 * Sigma_star
  
  list("bxl2t_star" = bxl2t_star,
       "Sigma_star" = Sigma_star)
  
}



findProposalBX_given_k2t_m1 <- function(bx, k2t, d, E, cxt, c2x, sd_bx){
  
  loglik_bx_k2t_current <- loglik_term3_bx_k2t_m1(k2t, d, E, cxt, sd_bx)
  gr_loglik_bx_k2t_current <- gr_loglik_term3_bx_k2t_m1(k2t, d, E, cxt, sd_bx)
  
  startVal_b <- c2x
  
  # find maximum
  laplace_fit <- optim(
    par = bx[-1],
    # par = rep(0, X - 2),
    fn = loglik_bx_k2t_current,
    method = c("CG"),
    gr = gr_loglik_bx_k2t_current,
    hessian = T
  )
  
  bxl2t_star <- laplace_fit$par
  Sigma_star <- solve(laplace_fit$hessian)
  
  # Sigma_star <- Sigma_star + diag(.001, nrow = length(bxl2t_star))
  # Sigma_star <- 2 * Sigma_star
  
  list("bxl2t_star" = bxl2t_star,
       "Sigma_star" = Sigma_star)
  
}

# for kt

update_bx_given_k2t_joint <- function(bx, k2t, d, E, cxt, delta2, sd_bx){
  
  if(!(all(k2t == 0))){ # no point in updating bx if k2t is 0
    
    if(delta2 == 1){
      
      list_proposal <- findProposalBX_given_k2t_m1(bx, k2t, d, E, cxt, c2x, sd_bx)
      bx_star <- list_proposal$bxl2t_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      Sigma_bx_star <- Sigma_bx_star
      
      bx_proposed_m1 <- mvrnorm(n = 1,
                                bx_star,
                                Sigma_bx_star)
      
      bx_proposed <- c(-1, bx_proposed_m1)
      
      loglik_proposal <- - loglik_term3(k2t, bx_proposed, d, E, cxt)
      loglik_current <- - loglik_term3(k2t, bx, d, E, cxt)
      
      logproposal_proposal <- dmvnorm(bx_proposed_m1, bx_star, Sigma_bx_star, log = T)
      logproposal_current <- dmvnorm(bx[-1], bx_star, Sigma_bx_star, log = T)
      
    } else if(delta2 == 2){
      
      list_proposal <- findProposalBX_given_k2t_m2(bx, k2t, d, E, cxt, c2x, sd_bx)
      bx_star <- list_proposal$bxl2t_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      Sigma_bx_star <- Sigma_bx_star
      
      bx_proposed_m2 <- mvrnorm(n = 1,
                                bx_star,
                                Sigma_bx_star)
      
      bx_proposed <- c(-1, bx_proposed_m2, - sum(bx_proposed_m2) + 1)
      
      loglik_proposal <- - loglik_term3(k2t, bx_proposed, d, E, cxt)
      loglik_current <- - loglik_term3(k2t, bx, d, E, cxt)
      
      logproposal_proposal <- dmvnorm(bx_proposed_m2, bx_star, Sigma_bx_star, log = T)
      logproposal_current <- dmvnorm(bx[-c(1,X)], bx_star, Sigma_bx_star, log = T)
      
    }
    
    logprior_proposal <- sum(dnorm(bx_proposed, sd = sd_bx, log = T))
    logprior_current <- sum(dnorm(bx, sd = sd_bx, log = T))
    
    mh_ratio <- exp(
      loglik_proposal - loglik_current + 
        logprior_proposal - logprior_current + 
        logproposal_current - logproposal_proposal
    )
    
    if(runif(1) < mh_ratio){
      
      bx <- bx_proposed
      
    }
    
  }
  
  bx
}

loglik_term3_x <- function(x, d, E, bx_x, k2t, cxt){
  
  mxt_x <- k2t * bx_x + cxt[x,] + log(E[x,])
  
  sum(d[x,] * mxt_x - exp(mxt_x))
  
}

update_bx_given_k2t <- function(bx, k2t, d, E, cxt, delta2, sd_bx){
  
  if(!(all(k2t == 0))){ # no point in updating bx if k2t is 0
    
    if(delta2 == 1){ # no constraint that sum bx = 0
      
      list_proposal <- findProposalBX_given_k2t_m1(bx, k2t, d, E, cxt, c2x, sd_bx)
      bx_star <- list_proposal$bxl2t_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      Sigma_bx_star <- Sigma_bx_star
      
      for (x in 2:X) {
        
        bx_proposed_x <- rt2(1, bx_star[x-1], 
                             sqrt(Sigma_bx_star[x-1,x-1]), df = 3)
        
        loglik_proposal <- loglik_term3_x(x, d, E, bx_proposed_x, k2t, cxt)
        loglik_current <- loglik_term3_x(x, d, E, bx[x], k2t, cxt)
        
        logproposal_proposal <- log(dt2(bx_proposed_x, bx_star, 
                                    sqrt(Sigma_bx_star[x-1,x-1]), df = 3))
        logproposal_current <- log(dt2(bx[x], bx_star, 
                                   sqrt(Sigma_bx_star[x-1,x-1]), df = 3))
        
        logprior_proposal <- dnorm(bx_proposed_x, sd = sd_bx, log = T)
        logprior_current <- dnorm(bx[x], sd = sd_bx, log = T)
        
        mh_ratio <- exp(
          loglik_proposal - loglik_current + 
            logprior_proposal - logprior_current + 
            logprior_current - logprior_proposal
        )
        
        if(runif(1) < mh_ratio){
          bx[x] <- bx_proposed_x
        }
        
      }
      
    } else if(delta2 == 2){
      
      list_proposal <- findProposalBX_given_k2t_m2(bx, k2t, d, E, cxt, c2x, sd_bx)
      bx_star <- list_proposal$bxl2t_star
      Sigma_bx_star <- list_proposal$Sigma_star
      
      Sigma_bx_star <- Sigma_bx_star
      
      bx_proposed_m2 <- mvrnorm(n = 1,
                                bx_star,
                                Sigma_bx_star)
      
      bx_proposed <- c(-1, bx_proposed_m2, - sum(bx_proposed_m2) + 1)
      
      loglik_proposal <- - loglik_term3(k2t, bx_proposed, d, E, cxt)
      loglik_current <- - loglik_term3(k2t, bx, d, E, cxt)
      
      logproposal_proposal <- dmvnorm(bx_proposed_m2, bx_star, Sigma_bx_star, log = T)
      logproposal_current <- dmvnorm(bx[-c(1,X)], bx_star, Sigma_bx_star, log = T)
      
      logprior_proposal <- sum(dnorm(bx_proposed, sd = sd_bx, log = T))
      logprior_current <- sum(dnorm(bx, sd = sd_bx, log = T))
      
      mh_ratio <- exp(
        loglik_proposal - loglik_current + 
          logprior_proposal - logprior_current + 
          logproposal_current - logproposal_proposal
      )
      
      if(runif(1) < mh_ratio){
        
        bx <- bx_proposed
        
      }
      
    }
    
  }
  
  bx
}

update_bxk2t_twosteps <- function(bx, k2t, d, E, cxt, delta1, delta2, sd_bx){
  
  # update kt
  
  list_proposal <- findProposalk2t_given_bx(k2t, bx, d, E, cxt)
  kt_star <- list_proposal$bxl2t_star
  Sigma_kt_star <- list_proposal$Sigma_star
  
  Sigma_kt_star <- Sigma_kt_star
  
  if(delta1 == 2){
   
    kt_proposed <- mrt2(kt_star,
                           Sigma_kt_star,
                           df = 3)
    
    logproposal_proposal <- dmt_cpp(kt_proposed, nu = 3, kt_star, 
                                    Sigma_kt_star, returnLog = T)
    logproposal_current <- dmt_cpp(k2t, nu = 3, kt_star, 
                                   Sigma_kt_star, returnLog = T)
    
  } else if (delta1 == 1){
    
    kt_proposed_m1 <- mrt2(kt_star,
                           Sigma_kt_star,
                           df = 3)
    
    kt_proposed <- c(kt_proposed_m1, - sum(kt_proposed_m1))
    
    logproposal_proposal <- dmt_cpp(kt_proposed_m1, nu = 3, kt_star, 
                                    Sigma_kt_star, returnLog = T)
    logproposal_current <- dmt_cpp(k2t[-Y], nu = 3, kt_star, 
                                   Sigma_kt_star, returnLog = T)
    
  }
  
  loglik_proposal <- - loglik_term3(kt_proposed, bx, d, E, cxt)
  loglik_current <- - loglik_term3(k2t, bx, d, E, cxt)
  
  # logproposal_proposal <- dmvnorm(kt_proposed_m1, kt_star, Sigma_kt_star, log = T)
  # logproposal_current <- dmvnorm(k2t[-Y], kt_star, Sigma_kt_star, log = T)
  
  mh_ratio <- exp(
    loglik_proposal - loglik_current + 
      logproposal_current - logproposal_proposal
  )
  
  if(runif(1) < mh_ratio){
    
    k2t <- kt_proposed
    
  }
  
  # update bx
  
  bx <- update_bx_given_k2t(bx, k2t, d, E, cxt, delta2, sd_bx)
  
  term3 <- matrix(k2t, X, Y, byrow = T) * matrix(bx, X, Y, byrow = F)
  
  list("bx" = bx,
       "k2t" = k2t,
       "term3" = term3)
    
}

computeTerm3 <- function(param3, delta3, ...){
  
  if(delta3 == 1){
    term3 <- matrix(0, X, Y, byrow = F)  
  } else if(delta3 == 2){
    term3 <- matrix(param3$k2t, X, Y, byrow = T) * c2x   
  } else if(delta3 == 3){
    term3 <- 
      matrix(param3$k2t, X, Y, byrow = T) * 
      matrix(param3$bx, X, Y, byrow = F)   
  }
  
  term3
}

# UPDATE TERM 4 -------


gr_loglik_term4_gtx <- function(d, E, cxt){
  
  function(param){
    
    param_matrix <- createTerm4_cpp(c(0, param, 0), X, Y)
    
    sumParams <- sum(param_matrix)
    
    gtx <- c(0, param, - sumParams)
    
    gtx_mat <- createTerm4_cpp(gtx, X, Y)
    
    term1 <- sapply(1:(X+Y-2), function(g){
      sum(sdiag(d, k = g - X))
    })
    
    term2 <- sapply(1:(X+Y-2), function(g){
      - sum(sdiag(exp(cxt + gtx_mat + log(E)), k = g - X))
    })
    
    term3 <- d[1,Y] * (- ik[-c(X+Y-1)])
    
    term4 <- - exp(cxt[1,Y] + gtx_mat[1,Y] + log(E[1,Y])) * (-ik[-c(X+Y-1)])
    
    - (term1 + term2 + term3 + term4)[-1]
    
  }
  
}






# OLD ----------

fastdmvnorm <- function(x, mu, Lambda){
  - .5 * t(x - mu) %*% Lambda %*% (x- mu)
}

compute_aap <- function(a, ap, model){
  
  if(model == "B"){
    
    aap <- matrix(a, X, P, byrow = F)
    
  } else if(model == "A"){
    
    aap <- matrix(a, X, P, byrow = F) + ap
    
  } else if(model == "M"){
    
    aap <- matrix(a, X, P, byrow = F) * ap
    
  } else if(model == "NP"){
    
    aap <- matrix(a, X, P, byrow = F) + ap
    
  }
  
  
  aap
}

loglik_r_cpp <- function(d, E){

  function(param){

    a <- param[1:X]
    b <- param[X + 1:X]
    k <- param[2 *X + 1:Y]
    - loglik_LC(a, b, k, d, E)

  }

}

loglik_r <- function(d, E){

  function(param){

    a <- param[1:X]
    b <- param[X + 1:X]
    k <- param[2 *X + 1:Y]
    
    mxt <- matrix(a, X, Y, byrow = F) + 
      matrix(b, X, Y, byrow = F) * 
      matrix(k, X, Y, byrow = T) 
    
    mxtp <- array(mxt, dim = c(X, Y, P)) + log(E)
    
    - sum(d * mxtp - exp(mxtp), na.rm = T)
    
    # - loglik_LC(a, b, k, d, E)

  }

}

# R loglikelihood 

loglik_LCp_r <- function(a, b, k, d, E){
  
  a <- aperm(array(a, dim = c(X, P, Y)), perm = c(1,3,2))
  b <- aperm(array(b, dim = c(X, P, Y)), perm = c(1,3,2))
  k <- aperm(array(k, dim = c(Y, P, X)), perm = c(3,1,2))

  mxtp <- a + b * k + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T)
}

# R loglikelihood for a specific x

loglik_LCp_x_r <- function(x, a, b, k, d, E){
  
  mxtp <- matrix(a[x,], Y, P, byrow = T) + 
    matrix(b[x,], Y, P, byrow = T) * k + log(E[x,,])
  
  term1 <- d[x,,] * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T) 
}

# R loglikelihood for a specific t

loglik_LCp_t_r <- function(t, a, b, k, d, E){
  
  mxtp <- a + b * matrix(k[t,], X, P, byrow = T) + log(E[,t,])
  
  term1 <- d[,t,] * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T)
}

# gradient of loglikelihood of ab for a specific x

gr_loglik_ab_x_r <- function(x, a, b, ab, bp, k, d, E){
  
  mxtp <- matrix(a[x] + ap[x,], Y, P, byrow = T) + 
    (matrix(b[x] + bp[x,], Y, P, byrow = T)) * k + log(E[x,,])
  
  grad1 <- sum( d[x,,]  - exp(mxtp), na.rm = T)
  grad2 <- sum( k * d[x,,] - exp(mxtp) * k, na.rm = T)
  
  c(grad1, grad2)
}

# UPDATE OF AB

{
  loglik_abx_r <- function(x, d, E, a, b, ap, bp, k){
    
    function(param){
      
      a[x] <- param[1]
      b[x] <- param[2]
      
      a2 <- matrix(a, X, P, byrow = F) + ap
      b2 <- matrix(b, X, P, byrow = F) + bp
      
      # - loglik_LCp_x(x - 1, a2, b2, k, d, E)
      - loglik_LCp_x_r(x, a2, b2, k, d, E)
        
      
      
    }
    
  }
  
  gr_loglik_abx_r <- function(x, d, E, a, b, ap, bp, k){
    
    function(param){
      
      a[x] <- param[1]
      b[x] <- param[2]
      
      # - gr_loglik_ab_x(x - 1, a, b, ap, bp, k, d, E)  
      - gr_loglik_ab_x_r(x, a, b, ap, bp, k, d, E)  
      
    }
    
  }
  
  findabMax <- function(d, E, k2, a, b){
    
    ab_star <- rep(NA, 2 * X)
    Sigma_star <- matrix(0, 2 * X, 2 * X)
    
    for (x in 1:X) {
      
      loglik_ab_current <- loglik_abx_r(x, d, E, a, b, ap, bp, k2)
      gr_loglik_ab_current <- gr_loglik_abx_r(x, d, E, a, b, ap, bp, k2)
      
      # find maximum
      laplace_fit <- optim(
        # par = c(a[x], b[x]),
        par = c(0, 0),
        fn = loglik_ab_current,
        method = c("BFGS"),
        gr = gr_loglik_ab_current,
        hessian = T
      )
      
      # laplace_fit <- optimx(
      #   par = c(a[x], b[x]),
      #   # par = c(0, 0),
      #   fn = loglik_ab_current,
      #   method = c("BFGS"),
      #   gr = gr_loglik_ab_current
      #   # hessian = T
      # )
      
      ab_star[x + c(0, X)] <- laplace_fit$par
      Sigma_star[x + c(0, X), x + c(0, X)] <- solve(laplace_fit$hessian)
      
    }
    
  }
  
  # loglikelihood of ab 
  
  loglik_ab_r <- function(d, E, ap, bp, k){
    
    function(param){
      
      a <- param[1:X]
      b <- param[X + 1:X]
      
      a2 <- matrix(a, X, P, byrow = F) + ap
      b2 <- matrix(b, X, P, byrow = F) + bp
      
      - loglik_LCp(a2, b2, k, d, E)  
      
    }
    
  }
  
  # loglik_ab_r <- function(d, E, ap, bp, k){
  #   
  #   function(param){
  #     
  #     a <- param[1:X]
  #     b <- param[X + 1:X]
  #     
  #     # a2 <- matrix(a, X, P, byrow = F) + ap
  #     # b2 <- matrix(b, X, P, byrow = F) + bp
  #     
  #     - loglik_LC(a, b, k, d, E)  
  #     
  #   }
  #   
  # }
  
  gr_loglik_ab_r <- function(d, E, ap, bp, k){
    
    function(param){
      
      a <- param[1:X]
      b <- param[X + 1:X]
      
      - gr_loglik_ab(a, b, ap, bp, k, d, E)  
      
    }
    
  }
  
      
}

# UPDATE OF K

{
  
  # loglikelihood of k given the rest for a specific t
  
  loglik_kt_r <- function(t, d, E, a, b, k, kp){
    
    function(param){
      
      k[t] <- param
      # k2 <- matrix(k, Y, P, byrow = F) + kp
      
      # - loglik_LCp_t(t - 1, a, b, k2, d, E)
      
      mxp <- a + b * matrix(k[t] + kp[t,], X, P, byrow = T) + log(E[,t,])

      term1 <- sum(d[,t,] * mxp, na.rm = T)
      term2 <- sum(exp(mxp), na.rm = T)

      - (term1 - term2)
      
    }
    
  }
  
  # gradient of loglikelihood of k given the rest for a specific t
  
  gr_loglik_kt_r <- function(t, d, E, ap, bp, k, kp){
    
    function(param){
      
      k[t] <- param
      
      # - gr_loglik_k_t(t - 1, ap, bp, k, kp, d, E)
      
      mxp <- ap + bp * matrix(k[t] + kp[t,], X, P, byrow = T) + log(E[,t,])
      
      term1 <- sum(d[,t,] * bp, na.rm = T)
      term2 <- sum(exp(mxp) * bp, na.rm = T)
      
      - (term1 - term2)
      
    }
    
  }
  
  # loglikelihood of k conditional on kp
  
  loglik_k0_r <- function(d, E, a, b, kp){
    
    function(param){
      
      k2 <- matrix(param, Y, P, byrow = F) + kp
      
      - loglik_LCp(a, b, k2, d, E)
      
    }
    
  }
  
  gr_loglik_k0_r <- function(d, E, a, b, kp){
    
    function(param){
      
      k <- param
      - gr_loglik_k(a, b, k, kp, d, E)  
      
    }
    
  }  
  
  # find proposal values for k0
  
  findProposalK0 <- function(t, a2, b2, k, kp, d, E){
    
    loglik_kt_current <- loglik_kt_r(t, d, E, a2, b2, k, kp)
    gr_loglik_kt_current <- gr_loglik_kt_r(t, d, E, a2, b2, k, kp)
    
    # find maximum
    
    laplace_fit <- optim(
      par = k[t],
      fn = loglik_kt_current,
      method = c("BFGS"),
      gr = gr_loglik_kt_current,
      hessian = T
    )
    
    # Sigma_star <- 2 * solve(laplace_fit$hessian)  
    Sigma_star <- solve(laplace_fit$hessian)  
    
    list("k_star" = laplace_fit$par,
         "Sigma_star" = Sigma_star)
  }
  
}

# UPDATE OF AP

{
  # loglikelihood of ap under the nonparametric model for a specific x
  
  loglik_apx_r <- function(x, a, ap, bp, d, E, k, sd_ap){
    
    function(param){
      
      ap[x,] <- param
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - (loglik_LCp_x_r(x, a + ap, bp, k, d, E) + 
           sum(dnorm(param, mean = 0, sd = sd_ap, log = T)))
      
    }
    
  }
  
  gr_loglik_ap_r_fun <- function(x, a, ap, bp, k, d, E){
    
    m_xtp <- matrix(a[x] + ap[x,], Y, P, byrow = T) + 
      matrix(bp[x,], Y, P, byrow = T) * k + log(E[x,,])
    
    term2 <- exp(m_xtp)
    
    gr_loglik <- d[x,,] - term2
    
    apply(gr_loglik, 2, function(x){
      sum(x, na.rm = T)
    })
    
  }
  
  gr_loglik_apx_r <- function(x, a, ap, bp, d, E, k, sd_ap){
    
    function(param){
      
      ap[x,] <- param
      # - gr_loglik_ap_x(x - 1, a, ap, bp, k, d, E)
      - (gr_loglik_ap_r_fun(x, a, ap, bp, k, d, E) - param / sd_ap^2)
      
    }
    
  }  
  
  # loglikelihood under the additive model for a
  
  loglik_ap_r <- function(a, bp, d, E, k){
    
    function(param){
      
      aap <- matrix(a, X, P, byrow = F) + matrix(param, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(aap, bp, k, d, E)
      
    }
    
  }
  
  # function: gradient of loglikelihood of ap under the additive model for a
  
  gr_loglik_ap_add_r_fun <- function(a, ap, b, k, d, E){
    
    aap <- matrix(a, X, P, byrow = F) + ap
    
    a <- aperm(array(aap, dim = c(X, P, Y)), perm = c(1,3,2))
    b <- aperm(array(b, dim = c(X, P, Y)), perm = c(1,3,2))
    k <- aperm(array(k, dim = c(Y, P, X)), perm = c(3,1,2))
    
    m_xtp <- a + b * k + log(E)
    
    term2 <- exp(m_xtp)
    
    gr_loglik <- d - term2
    
    apply(gr_loglik, 3, function(x){
      sum(x, na.rm = T)
    })
    
  }
  
  # gradient of loglikelihood of ap under the additive model for a
  
  gr_loglik_ap_r <- function(a, bp, d, E, k){
    
    function(param){
      
      ap <- matrix(param, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - gr_loglik_ap_add_r_fun(a, ap, bp, k, d, E)
      
    }
    
  }
  
  # find proposal values in the additive model
  
  findProposalApAdditive <- function(a, b2, d, E, k2){
    
    loglik_ap_current <- loglik_ap_r(a, b2, d, E, k2)
    gr_loglik_ap_current <- gr_loglik_ap_r(a, b2, d, E, k2)
    
    # find maximum
    laplace_fit <- optim(
      par = rep(0, P),
      fn = loglik_ap_current,
      method = c("BFGS"),
      gr = gr_loglik_ap_current,
      hessian = T
    )
    
    ap_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- 2 * solve(H_star)  
    
    list("ap_star" = ap_star,
         "Sigma_star" = Sigma_star)
  }
  
  # find proposal values in the nonparametric model for a specific x
  
  findProposalApXNonparametrics <- function(x, a, ap, b2, d, E, k2, sd_ap){
    
    loglik_apx_current <- loglik_apx_r(x, a, ap, b2, d, E, k2, sd_ap)
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
    
    list("ap_star" = ap_star,
         "Sigma_star" = Sigma_star)
  }
  
  # loglikelihood under the multiplicative model for a
  
  loglik_ap_m_r <- function(a, bp, d, E, k){
    
    function(param){
      
      ap <- matrix(param, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(a + ap, bp, k, d, E)
      
    }
    
  }
}

# UPDATE OF BP

{
  
  loglik_bpx_r <- function(x, ap, b, bp, d, E, k, sd_bp){
    
    function(param){
      
      bp[x,] <- param
      bbp <- matrix(b, X, P, byrow = F) + bp
      # - loglik_LCp_x(x - 1, ap, bbp, k, d, E)
      - (loglik_LCp_x_r(x, ap, bbp, k, d, E) + 
           sum(dnorm(param, mean = 0, sd = sd_bp, log = T)))
      
    }
    
  }
  
  gr_loglik_bp_r_fun <- function(x, ap, b, bp, k, d, E){
    
    m_xtp <- matrix(ap[x,], Y, P, byrow = T) + 
      matrix(b[x] + bp[x,], Y, P, byrow = T) * k + log(E[x,,])

    term2 <- exp(m_xtp)
    
    gr_loglik <- d[x,,] * k - term2 * k
    
    apply(gr_loglik, 2, function(x){
      sum(x, na.rm = T)
    })
    
  }
  
  gr_loglik_bpx_r <- function(x, ap, b, bp, d, E, k, sd_bp){
    
    function(param){
      
      bp[x,] <- param
      # - gr_loglik_bp_x(x - 1, ap, b, bp, k, d, E)
      - (gr_loglik_bp_r_fun(x, ap, b, bp, k, d, E) - param / sd_bp^2)
      
    }
    
  }  
  
  # loglikelihood under the additive model for b
  
  loglik_bp_r <- function(ap, b, d, E, kp){
    
    function(param){
      
      bbp <- matrix(b, X, P, byrow = F) + matrix(param, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(ap, bbp, kp, d, E)
      
    }
    
  }
  
  # function: gradient of loglikelihood of ap under the additive model for a
  
  gr_loglik_bp_add_r_fun <- function(a, b, bp, k, d, E){
    
    bbp <- matrix(b, X, P, byrow = F) + bp
    
    a <- aperm(array(a, dim = c(X, P, Y)), perm = c(1,3,2))
    b <- aperm(array(bbp, dim = c(X, P, Y)), perm = c(1,3,2))
    k <- aperm(array(k, dim = c(Y, P, X)), perm = c(3,1,2))
    
    m_xtp <- a + b * k + log(E)
    
    term2 <- exp(m_xtp)
    
    gr_loglik <- (d - term2) * k
    
    apply(gr_loglik, 3, function(x){
      sum(x, na.rm = T)
    })
    
  }
  
  # gradient of loglikelihood of ap under the additive model for a
  
  gr_loglik_bp_r <- function(ap, b, d, E, k){
    
    function(param){
      
      bp <- matrix(param, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - gr_loglik_bp_add_r_fun(ap, b, bp, k, d, E)
      
    }
    
  }
  
  # find proposal values in the additive model
  
  findProposalBpAdditive <- function(a2, b, d, E, k2){
    
    loglik_bp_current <- loglik_bp_r(a2, b, d, E, k2)
    gr_loglik_bp_current <- gr_loglik_bp_r(a2, b, d, E, k2)
    
    # find maximum
    laplace_fit <- optim(
      par = rep(0, P),
      fn = loglik_bp_current,
      method = c("BFGS"),
      gr = gr_loglik_bp_current,
      hessian = T
    )
    
    bp_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- 2 * solve(H_star)  
    
    list("bp_star" = bp_star,
         "Sigma_star" = Sigma_star)
  }
  
  # find proposal values in the nonparametric model for a specific x
  
  findProposalBpXNonparametrics <- function(x, a2, b, bp, d, E, k2, sd_bp){
    
    loglik_bpx_current <- loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
    gr_loglik_bpx_current <- gr_loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
    
    # find maximum
    laplace_fit <- optim(
      # par = ap[x,],
      par = rep(0, P),
      fn = loglik_bpx_current,
      method = c("BFGS"),
      gr = gr_loglik_bpx_current,
      hessian = T
    )
    
    bp_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- 2 * solve(H_star)
    
    list("bp_star" = bp_star,
         "Sigma_star" = Sigma_star)
  }
  
}

# UPDATE OF KP

{
  
  # loglikelihood of kp for a specific t
  
  loglik_kpt_r <- function(t, ap, bp, k, kp, d, E){
    
    function(param){
      
      kp[t,] <- param
      # - loglik_LCp_t(t - 1, ap, bp, k + kp, d, E)
      - loglik_LCp_t_r(t, ap, bp, k + kp, d, E)
      
    }
    
  }
  
  # gradient of loglikelihood of kp for a specific t
  
  gr_loglik_kpt_r <- function(t, ap, bp, k, kp, d, E){
    
    function(param){
      
      kp[t,] <- param
      - gr_loglik_kpt_r_fun(t, ap, bp, k, kp, d, E)
      # - gr_loglik_kp_t(t - 1, ap, bp, k, kp, d, E)
      
    }
    
  }  
  
  # function: gradient of loglikelihood of kp for a specific t
  
  gr_loglik_kpt_r_fun <- function(t, ap, bp, k, kp, d, E){
    
    mxtp <- ap + bp * matrix(k[t] + kp[t,], X, P, byrow = T) + log(E[,t,])
    
    term1 <- d[,t,] * bp
    term2 <- exp(mxtp) * bp
    
    apply(term1 - term2, 2, function(x){
      sum(x, na.rm = T)
    })
    
  }  
  
  # loglikelihood under the additive model for k
  
  loglik_kp_r <- function(ap, bp, d, E, k){
    
    function(param){
      
      kkp <- 
        matrix(k, Y, P, byrow = F) + 
        matrix(param, Y, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(ap, bp, kkp, d, E)
      
    }
    
  }
  
  # function: gradient of loglikelihood of kp under the additive model for k
  
  gr_loglik_kp_add_r_fun <- function(a, b, k, kp, d, E){
    
    kkp <- matrix(k, Y, P, byrow = F) + kp
    
    a <- aperm(array(p, dim = c(X, P, Y)), perm = c(1,3,2))
    b <- aperm(array(b, dim = c(X, P, Y)), perm = c(1,3,2))
    k <- aperm(array(kp, dim = c(Y, P, X)), perm = c(3,1,2))
    
    m_xtp <- a + b * k + log(E)
    
    term2 <- exp(m_xtp)
    
    gr_loglik <- d - term2
    
    apply(gr_loglik, 3, function(x){
      sum(x, na.rm = T)
    })
    
  }
  
  # gradient of loglikelihood of ap under the additive model for a
  
  gr_loglik_kp_r <- function(a, b, d, E, k){
    
    function(param){
      
      kp <- matrix(param, Y, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - gr_loglik_kp_add_r_fun(a, b, k, kp, d, E)
      
    }
    
  }
   
  # find proposal values in the additive model
  
  findProposalKpAdditive <- function(a2, b2, d, E, k){
    
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
    
    list("kp_star" = kp_star,
         "Sigma_star" = Sigma_star)
  }
  
  # find proposal values in the nonparametric model for a specific x
  
  findProposalKptNonparametrics <- function(t, a2, b2, d, E, 
                                            k, kp, sd_kt){
    
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
    
    kp_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- 2 * solve(H_star)
    
    list("kp_star" = kp_star,
         "Sigma_star" = Sigma_star)
  }
  
}

#

sampleMVNconstraint <- function(mu, Sigma){
  
  d <- length(mu)
  
  A <- matrix(0, 2, d)
  A[1, 1:(d / 2)] <- 1
  A[2, d / 2 + 1:(d / 2)] <- 1
  b <- c(0, 1)
  
  P <- matrix(0, d, d)
  diag(P)[1:(d/2-1)] <- 1
  diag(P)[d/2 + 1:(d/2-1)] <- 1
  P <- P[-c(d/2, d),]
  
  C <- rbind(P, A)
  
  nrow1 <- 1:(d - nrow(A))
  nrow2 <- d - nrow(A) + 1:nrow(A)
  
  nu <- C %*% mu
  Omega <- C %*% Sigma %*% t(C)
  Omega11 <- Omega[nrow1, nrow1]
  Omega22 <- Omega[nrow2, nrow2]
  Omega12 <- Omega[nrow1, nrow2]
  Omega21 <- Omega[nrow2, nrow1]
  
  mu1 <- nu[nrow1] + Omega12 %*% solve(Omega22) %*% (b - nu[nrow2])
  Sigma1 <- Omega11 - Omega12 %*% solve(Omega22) %*% Omega21
  
  pxcd <- mvrnorm(1, mu1, Sigma1)
  
  x <- solve(C) %*% c(pxcd, b)
  
  x
}


sampleMVNconstraint_ab <- function(mu, Sigma){
  
  d <- length(mu)
  
  A <- matrix(0, 1, d)
  A[1, d / 2 + 1:(d / 2)] <- 1
  b <- 1
  
  P <- matrix(0, d, d)
  diag(P)[1:d - 1] <- 1
  P <- P[-d,]
  
  C <- rbind(P, A)
  
  nrow1 <- 1:(d - nrow(A))
  nrow2 <- d - nrow(A) + 1:nrow(A)
  
  nu <- C %*% mu
  Omega <- C %*% Sigma %*% t(C)
  Omega11 <- Omega[nrow1, nrow1]
  Omega22 <- Omega[nrow2, nrow2]
  Omega12 <- Omega[nrow1, nrow2]
  Omega21 <- Omega[nrow2, nrow1]
  
  mu1 <- nu[nrow1] + Omega12 %*% solve(Omega22) %*% (b - nu[nrow2])
  Sigma1 <- Omega11 - Omega12 %*% solve(Omega22) %*% Omega21
  
  pxcd <- mvrnorm(1, mu1, Sigma1)
  
  x <- solve(C) %*% c(pxcd, b)
  
  x
}





plotCurves_old <- function(delta1, delta2, delta3, delta4,
                       param1, param2, param3, param4,
                       c2x){
  
  data <- expand.grid(x = 1:X, t = 1:Y)
  
  {
    
    if(delta1 == 2){
      mu0 <- param1[1]
      b0 <- param1[2]
      term1 <- mu0 + b0 * c2x
    } else {
      ax <- param1
      term1 <- ax
    }
  }
  
  # Choice 2 (Additive year effect)
  
  {
    
    if(delta2 == 2){
      k1t <- param2
      term2 <- k1t
    } else {
      param2 <- rep(0, Y)
    }
    
  }
  
  # Choice 3 (Varying year-age effect)
  
  {
    
    if(delta3 == 2){
      k2t <- param3$k2t
      term3 <- matrix(c2x, X, Y, byrow = F) * matrix(k2t, X, Y, byrow = T)
    } else if(delta3 == 3){
      bx <- param3$bx
      k2t <- param3$k2t
      term3 <- matrix(bx, X, Y, byrow = F) * matrix(k2t, X, Y, byrow = T)
    } else if(delta3 == 1){
      term3 <- matrix(0, X, Y)  
    }
    
  }
  
  # Choice 4 (Cohort effect)
  
  {
    
    
    if(delta4 == 2){
      gtx <- param4
    } else {
      gtx <- rep(0, X + Y - 1)
    }
    
    term4 <- createTerm4_cpp(gtx, X, Y)
    
  }
  
  mxtp <- apply(data, 1, function(dat){
    x <- dat[1]
    t <- dat[2]
    term1[x] + term2[t] + term3[x,t] + term4[x,t]
  })
  
  df <- data.frame(Age = data$x,
                   Year = factor(data$t),
                   m = mxtp)
  
}