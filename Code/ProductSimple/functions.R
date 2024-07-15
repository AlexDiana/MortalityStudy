
computeTerms <- function(idx_term = 0, 
                         param_a, 
                         param_k, 
                         param_g, 
                         dp){
  
  a_term <- aperm(array(param_a$yip, dim = c(X, P, Y)), perm = c(1,3,2))
  k_term <- aperm(array(param_k$yip, dim = c(Y, P, X)), perm = c(3,1,2))
  g_term <- createTermG(param_g$yip, X, Y)
  p_term <- aperm(array(dp, dim = c(P, X, Y)), perm = c(2,3,1))
  
  allTerms <- a_term + k_term + g_term + p_term
  
  if(idx_term == 1){
    return(allTerms - a_term)
  } else if(idx_term == 2){
    return(allTerms - k_term)
  } else if(idx_term == 3){
    return(allTerms - g_term)
  } else if(idx_term == 4){
    return(allTerms - p_term)
  } 
  
  return(allTerms)
}

simulateMultModelOld <- function(N, P, meanVar, stdev,
                              modelType, sum0){
  
  y <- rnorm(N, meanVar, stdev)
  
  if(sum0) y <- y - mean(y)
  
  if (modelType >= 1) { 
    yip <- matrix(y, N, P, byrow = F)
    c1p <- NULL
    c2p <- NULL
  } 
  
  if (modelType >= 2) { 
    cpm1 <- rgamma(P - 1, 100, 100)
    c1p <- c(cpm1, P - sum(cpm1))
    yip <- t(sapply(1:N, function(i){
      y[i] * c1p
    }))
    
    c2p <- NULL
  } 
  
  if (modelType == 3) { 
    
    c2p <- sapply(1:(P-1), function(p){
      rgamma(N - 1, 10, 10)
    })
    
    c2P <- sapply(1:(P-1), function(p){
      (c1p[p] * sum(y) - sum(y[1:(N-1)] * c2p[1:(N-1),p])) / y[N]
    })
    
    c2p <- rbind(c2p, c2P)
    # c2p <- rbind(c2p, X - apply(c2p, 2, sum))
    
    c2P <- sapply(1:N, function(i){
      (P - sum(c1p[1:(P-1)] * c2p[i,1:(P-1)])) / c1p[P]
    })
    c2p <- cbind(c2p, c2P)
    # c2xp[X,P] <- X - sum(c2xp[1:(X-1),P])
    
    # mean(c1p * c2xp[x,])
    
    y_ip <- 
      matrix(y, N, P, byrow = F) *
      matrix(c1p, N, P, byrow = T) * c2p
    
    # first check
    apply(y_ip, 1, mean) - y
    
    # second check
    sapply(1:P, function(p){
      sum(y_ip[,p])
    }) - sapply(1:P, function(p){
      sum(y * c1p[p])
    })
    
    # inverse calcs
    # y_2 <- apply(y_ip, 1, mean)
    # 
    # c1p_2 <- sapply(1:P, function(p){
    #   sum(yip[,p]) / sum(y_2)
    # })
    
    
      
    # apply(axp_true, 1, mean)
    # ax
    
  } 
  
  list("y" = y,
       "yip" = yip,
       "c1p" = c1p,
       "c2p" = c2p)
  
}

simulateMultModel <- function(N, P, meanVar, stdev,
                              modelType, sum0 = F, first0 = F){
  
  y <- rnorm(N, meanVar, stdev)
  
  if(first0) y[1] <- 0
  
  if(sum0) {
    if(!first0){
      y <- y - mean(y)  
    } else {
      y[2:N] <- y[2:N] - mean(y[2:N])
    }
    
  }
  
  if (modelType >= 1) { 
    yip <- matrix(y, N, P, byrow = F)
    c1p <- NULL
    c2p <- NULL
  } 
  
  if (modelType == 2) { 
    cpm1 <- rgamma(P - 1, 100, 100)
    c1p <- c(cpm1, P - sum(cpm1))
    yip <- t(sapply(1:N, function(i){
      y[i] * c1p
    }))
    
    c2p <- NULL
  } 
  
  if (modelType == 3) { 
    
    c2p <- t(sapply(1:N, function(p){
      rgamma(P - 1, 100, 100)
    }))
    
    c2P <- P - apply(c2p, 1, sum)
    
    c2p <- cbind(c2p, c2P)
    
    yip <- 
      matrix(y, N, P, byrow = F)* c2p
    
    c1p <- NULL
    
    # first check
    # apply(y_ip, 1, mean) - y
    # 
    # # second check
    # sapply(1:P, function(p){
    #   sum(y_ip[,p])
    # }) - sapply(1:P, function(p){
    #   sum(y * c1p[p])
    # })
    # 
    # inverse calcs
    # y_2 <- apply(y_ip, 1, mean)
    # 
    # c1p_2 <- sapply(1:P, function(p){
    #   sum(yip[,p]) / sum(y_2)
    # })
    
    
      
    # apply(axp_true, 1, mean)
    # ax
    
  } 
  
  list("y" = y,
       "yip" = yip,
       "c1p" = c1p,
       "c2p" = c2p)
  
}

createTermG <- function(gcp, X, Y){
  
  term4 <- array(0, dim = c(X, Y, P))
  for (g in 1:(X+Y-1)) {
    for (p in 1:P) {
      sdiag(term4[,,p], k = g - X) <- gcp[g,p] 
    }
  }
  
  term4
}

createTerm4 <- function(gc, X, Y){
  
  term4 <- matrix(0, X, Y)
  for (g in 1:(X+Y-1)) {
    sdiag(term4, k = g - X) <- gc[g] 
  }
  
  term4
}

updateParameters <- function(param_current){
  
  
  
}

computeCP <- function(param, gamma_param, N){
  
  if(gamma_param == 1){
    
    matrix(1, N, P)
    
  } else if(gamma_param == 2){
    
    matrix(param$c1p, N, P, T)
    
  } else {
    
    param$c2p
    
  }
  
}

create_ik <- function(X,Y){
  c(
    seq(1, (min(X,Y) - 1)),
    rep(min(X,Y), max(X,Y) - min(X,Y) + 1),
    X + Y - seq(max(X,Y) + 1, X + Y - 1)
  )
}

#

loglik_LC <- function(d, lambda){
  
  # sum(dpois(as.vector(d), as.vector(lambda), log = T), na.rm  = T)
  
  loglambda <- log(lambda)
  
  sum(
    as.vector(d) * as.vector(loglambda) - exp(as.vector(loglambda)),
  na.rm = T)
  
}

# A PARAMETER UPDATE ---

{
  
  loglik_ax_cp <- function(d, cp, cxtp){
    
    function(ax){
      terms2 <- ax * matrix(cp, Y, P, byrow = T)
      # - sum(dpois(d, lambda = exp(terms2 + cxtp), log = T), na.rm = T)  
      
      - loglik_LC(d, exp(terms2 + cxtp))
      
    }
    
  }
  
  loglik_cp_ax <- function(d, ax, cxtp){
    
    function(cp){
      terms2 <- cp * matrix(ax, X, Y, byrow = F)
      
      loglik <- loglik_LC(d, exp(terms2 + cxtp))
      logprior <- sum(dgamma(cp, 100, 100, log = T))
                          
      # - sum(dpois(d, lambda = exp(terms2 + cxtp), log = T), na.rm = T)  
      - (loglik + logprior)
    }
    
  }
  
}

# A MODEL UPDATE ----

{
  findProposalAX <- function(d, E, m_terms){
    
    mu_ax_star <- rep(NA, X)
    sd_ax_star <- rep(NA, X)
    
    for (x in 1:X) {
      
      mu_ax_star[x] <- log(sum(d[x,,], na.rm = T) / 
                             sum(exp(m_terms[x,,] + log(E[x,,])), na.rm = T))
      
      hessian_ax_star <- derl2der2a(mu_ax_star[x], m_terms[x,,] + log(E[x,,]))
      sd_ax_star[x] <- sqrt(- 1 / hessian_ax_star)
      
      
    }
    
    
    list("mu_ax_star" = mu_ax_star,
         "sd_ax_star" = sd_ax_star)
    
  }
  
  loglik_cpm1_ax <- function(d, ax, E, cxtp){
    
    function(cpm1){
      
      cp <- c(cpm1,P - sum(cpm1))
      
      ax_mat <- array(ax, dim = c(X, Y, P))
      cp_mat <- aperm(array(cp, dim = c(P, X, Y)), perm = c(2,3,1))

      mxtp <- ax_mat * cp_mat + cxtp + log(E)

      term1 <- d * mxtp
      term2 <- exp(mxtp)

      loglik <- sum(term1 - term2, na.rm = T)
      
      # loglik <- loglik_cpm1_ax_cpp(cp, d, E, ax, cxtp)
      
      logprior <- sum(dgamma(cp, 100, 100, log = T))
      
      - (logprior + loglik)
      
    }
    
  }
  
  findProposalCP_AX <- function(d, E, cxt){
    
    # find ax first
    
    list_proposalax <- findProposalAX(d, E, cxt)
    mu_ax_star <- list_proposalax$mu_ax_star
    sd_ax_star <- list_proposalax$sd_ax_star
    
    # next c1p
    loglik_cp_current <- loglik_cpm1_ax(d, mu_ax_star, E, cxt)
    
    fitModel <- optim(rep(1, P - 1),
                      loglik_cp_current,
                      hessian = T,
                      method = "BFGS")
    cpm1_star <- fitModel$par
    Sigma_star <- solve(fitModel$hessian)
    
    list(
      "mu_ax_star" = mu_ax_star,
      "sd_ax_star" = sd_ax_star,
      "cpm1_star" = cpm1_star,
      "Sigma_star" = Sigma_star)
    
  }
  
  findProposalAXP <- function(d, E, cxt){
    
    mu_axp_star <- matrix(NA, X, P)
    sd_axp_star <- matrix(NA, X, P)
    
    for (x in 1:X) {
      for (p in 1:P) {
        
        mu_axp_star[x,p] <- log(sum(d[x,,p], na.rm = T) / 
                                  sum(exp(cxt[x,,p] + log(E[x,,p])), na.rm = T))
        
        hessian_ax_star <- derl2der2a(mu_axp_star[x,p], cxt[x,,p] + log(E[x,,p]))
        sd_axp_star[x,p] <- sqrt(- 1 / hessian_ax_star)
        
      }
    }
    
    list("mu_axp_star" = mu_axp_star,
         "sd_axp_star" = sd_axp_star)
  }
}

# K PARAMETER UPDATE ----

{
  loglik_kt_cxp <- function(d, E, cxp, m_terms){
    
    function(kt){
      
      ktcp <- matrix(kt, Y, P, byrow = F) * cxp
      
      ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
      
    }  
    
  }
  
  loglik_ktm1_cp_cxp <- function(d, E, cxp, m_terms){
    
    function(ktm1){
      
      kt <- c(ktm1, - sum(ktm1))

      ktcp <- matrix(kt, Y, P, byrow = F) * cxp

      ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))

      mxtp <- ktcp_mat + m_terms + log(E)

      term1 <- d * mxtp
      term2 <- exp(mxtp)

      - sum(term1 - term2, na.rm = T)
      
      # - loglik_ktm1_cp_cxp_cpp(kt, d, E, cxp, m_terms)
      
      
    }  
    
  }
  
  gr_loglik_ktm1_cp_cxp <- function(d, E, cxp, m_terms){
    
    function(ktm1){
      
      kt <- c(ktm1, - sum(ktm1))
      
      ktcp <- matrix(kt, Y, P, byrow = F) * cxp
      
      ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      cxp_array <- aperm(array(cxp, dim = c(Y, P, X)), perm = c(3,1,2))
      
      term1_grad <- apply(d * cxp_array, 2, function(x){
        sum(x, na.rm = T)
      })
      
      term2_grad <- apply(exp(mxtp) * cxp_array, 2, function(x){
        sum(x, na.rm = T)
      })
      
      term1_1_grad <- term1_grad[1:(Y-1)]
      term1_2_grad <- term1_grad[Y]
      
      term2_1_grad <- term2_grad[1:(Y-1)]
      term2_2_grad <- term2_grad[Y]
      
      - (term1_1_grad - term2_1_grad - term1_2_grad + term2_2_grad)
      
      
    }  
    
  }
  
  updateKT_CP <- function(kt_current, c2xp, d, E, m_terms){
    
    ktm1_current <- kt_current[1:(Y-1)]
    
    # loglik_kt_current <- loglik_kt_cxp(d, E, c2p, m_terms)
    loglik_kt_current <- loglik_ktm1_cp_cxp(d, E, c2xp, m_terms)
    gr_loglik_kt_current <- gr_loglik_ktm1_cp_cxp(d, E, c2xp, m_terms)
    
    fitModel <- optim(
      # rep(0, Y),
      rep(0, Y - 1),
      fn = loglik_kt_current,
      gr = gr_loglik_kt_current,
      hessian = T,
      method = "BFGS")
    # kt_star <- fitModel$par
    ktm1_star <- fitModel$par
    Sigma_star <- solve(fitModel$hessian)
    
    # kt_proposed <- sampleMTconstraint_k(kt_star, Sigma_star)
    ktm1_proposed <- mrt2(ktm1_star, Sigma_star, df = 3)
    kt_proposed <- c(ktm1_proposed, - sum(ktm1_proposed))
    
    loglik_current <- - loglik_kt_current(ktm1_current)
    loglik_proposed <- - loglik_kt_current(ktm1_proposed)
    
    logproposal_current <- dmt_cpp(ktm1_current, nu = 3, ktm1_star, Sigma_star, T)
    logproposal_proposal <- dmt_cpp(ktm1_proposed, nu = 3, ktm1_star, Sigma_star, T)
    
    mh_ratio <- exp(loglik_proposed - loglik_current + 
                      logproposal_current - logproposal_proposal)
    
    if(runif(1) < mh_ratio){
      
      return(kt_proposed)
      
    }
    
    return(kt_current)
    
  }
  
  loglik_c1p_kt <- function(d, E, kt, m_terms){
    
    function(c1p){
      
      ktcp <- matrix(kt, Y, P, byrow = F) * matrix(c1p, Y, P, byrow = T)
      
      ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
      
    }  
    
  }
  
  updateC1P_KT <- function(c1p_current, kt, d, E, m_terms){
    
    loglik_c1p_current <- loglik_c1p_kt(d, E, kt, m_terms)
    
    fitModel <- optim(rep(1, P),
                      loglik_c1p_current,
                      hessian = T,
                      method = "BFGS")
    c1p_star <- fitModel$par
    Sigma_star <- solve(fitModel$hessian)
    
    c1p_proposed <- sampleMTconstraint_c1p(c1p_star, Sigma_star)
    
    loglik_current <- - loglik_c1p_current(c1p_current)
    loglik_proposed <- - loglik_c1p_current(c1p_proposed)
    
    logproposal_current <- dmt_cpp(c1p_current, nu = 3, c1p_star, Sigma_star, T)
    logproposal_proposal <- dmt_cpp(c1p_proposed, nu = 3, c1p_star, Sigma_star, T)
    
    mh_ratio <- exp(loglik_proposed - loglik_current + 
                      logproposal_current - logproposal_proposal)
    
    if(runif(1) < mh_ratio){
      
      return(c1p_proposed)
      
    }
    
    return(c1p_current)
    
  }
  
  loglik_c2p_kt <- function(d, E, kt, m_terms){
    
    function(c2p_t){
      
      ktcp_mat <- matrix(kt, X, P) * matrix(c2p_t, X, P, byrow = T)
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - (sum(term1 - term2, na.rm = T) + sum(dnorm(c2p_t, 0, sd = 1, log = T)))
      
      
    }  
    
  }
  
  updateC2P_KT <- function(c2p_current, kt, d, E, m_terms){
    
    c2p_new <- matrix(NA, Y, P)
    
    for (t in 1:Y) {
      
      loglik_c2p_current <- loglik_c2p_kt(d[,t,], E[,t,], kt[t], m_terms[,t,])
      
      fitModel <- optim(rep(1, P),
                        loglik_c2p_current,
                        hessian = T,
                        method = "BFGS")
      c2p_star <- fitModel$par
      Sigma_star <- solve(fitModel$hessian)
      
      c2p_proposed <- sampleMTconstraint_c1p(c2p_star, Sigma_star)
      
      loglik_current <- - loglik_c2p_current(c2p_current[t,])
      loglik_proposed <- - loglik_c2p_current(c2p_proposed)
      
      logproposal_current <- dmt_cpp(c2p_current[t,], nu = 3, c2p_star, Sigma_star, T)
      logproposal_proposal <- dmt_cpp(as.vector(c2p_proposed), nu = 3, c2p_star, Sigma_star, T)
      
      mh_ratio <- exp(loglik_proposed - loglik_current + 
                        logproposal_current - logproposal_proposal)
      
      if(runif(1) < mh_ratio){
        
        c2p_new[t,] <- c2p_proposed
        
      } else {
        
        c2p_new[t,] <- c2p_current[t,]
        
        
      }
      
    }
    
    c2p_new
  }
}

# K MODEL UPDATE ------

{
  loglik_ktm1_cxp <- function(d, E, m_terms){
    
    function(ktm1){
      
      kt <- c(ktm1, - sum(ktm1))
      
      ktcp <- matrix(kt, Y, P, byrow = F)
      
      ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
      
    }  
    
  }
  
  gr_loglik_ktm1_cxp <- function(d, E, m_terms){
    
    function(ktm1){
      
      kt <- c(ktm1, - sum(ktm1))
      
      ktcp <- matrix(kt, Y, P, byrow = F)
      
      ktcp_mat <- aperm(array(ktcp, dim = c(Y, P, X)), perm = c(3, 1, 2))
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1_grad <- apply(d, 2, function(x){
        sum(x, na.rm = T)
      })
      
      term2_grad <- apply(exp(mxtp), 2, function(x){
        sum(x, na.rm = T)
      })
      
      term1_1_grad <- term1_grad[1:(Y-1)]
      term1_2_grad <- term1_grad[Y]
      
      term2_1_grad <- term2_grad[1:(Y-1)]
      term2_2_grad <- term2_grad[Y]
      
      - (term1_1_grad - term2_1_grad - term1_2_grad + term2_2_grad)
      
      
    }  
    
  }
  
  findProposalKT <- function(d, E, m_terms){
    
    loglik_kt_current <- loglik_ktm1_cxp(d, E, m_terms)
    gr_loglik_kt_current <- gr_loglik_ktm1_cxp(d, E, m_terms)
    
    fitModel <- optim(rep(0, Y - 1),
                      fn = loglik_kt_current,
                      gr = gr_loglik_kt_current,
                      hessian = T,
                      method = "BFGS")
    ktm1_star <- fitModel$par
    Sigma_star <- solve(fitModel$hessian)
    
    list("ktm1_star" = ktm1_star,
         "Sigma_ktm1_star" = Sigma_star)
    
  }
  
  loglik_cpm1_kt <- function(d, kt, cxtp){
    
    function(cpm1){
      
      cp <- c(cpm1,P - sum(cpm1))
      
      kt_mat <- aperm(array(kt, dim = c(Y, X, P)), perm = c(2,1,3))
      cp_mat <- aperm(array(cp, dim = c(P, X, Y)), perm = c(2,3,1))
      
      mxtp <- kt_mat * cp_mat + cxtp + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
    }
    
  }
  
  findProposalC1P_KT <- function(d, E, cxt){
    
    # find kt first
    
    list_proposalkt <- findProposalKT(d, E, cxt)
    mu_ktm1_star <- list_proposalkt$ktm1_star
    Sigma_ktm1_star <- list_proposalkt$Sigma_ktm1_star
    
    mu_kt_star <- c(mu_ktm1_star, - sum(mu_ktm1_star))
    
    # next c1p
    loglik_cp_current <- loglik_cpm1_kt(d, mu_kt_star, cxt)
    
    fitModel <- optim(rep(1, P - 1),
                      loglik_cp_current,
                      hessian = T,
                      method = "BFGS")
    cpm1_star <- fitModel$par
    Sigma_star <- solve(fitModel$hessian)
    
    list(
      "mu_ktm1_star" = mu_ktm1_star,
      "Sigma_ktm1_star" = Sigma_ktm1_star,
      "cpm1_star" = cpm1_star,
      "Sigma_star" = Sigma_star)
    
  }
  
  loglik_c2pm1_kt <- function(d, E, kt, m_terms){
    
    function(c2pm1_t){
      
      c2p_t <- c(c2pm1_t, P - sum(c2pm1_t))
      
      ktcp_mat <- matrix(kt, X, P) * matrix(c2p_t, X, P, byrow = T)
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
      
    }  
    
  }
  
  gr_loglik_c2pm1_kt <- function(d, E, kt, m_terms){
    
    function(c2pm1_t){
      
      c2p_t <- c(c2pm1_t, P - sum(c2pm1_t))
      
      kt_mat <- matrix(kt, X, P)
      c2p_mat <- matrix(c2p_t, X, P, byrow = T)
      
      ktcp_mat <- kt_mat * c2p_mat
      
      mxtp <- ktcp_mat + m_terms + log(E)
      
      term1_1 <- (d * kt_mat)[,1:(P-1)]
      term2_1 <- exp(mxtp[,1:(P-1)]) * kt_mat[,1:(P-1)]
      
      term1_2 <- (d * kt_mat)[,P]
      term2_2 <- exp(mxtp[,P]) * kt_mat[,P]
      
      grad_term1_1 <- apply(term1_1, 2, function(x){
        sum(x, na.rm = T)
      })
      
      grad_term2_1 <- apply(term2_1, 2, function(x){
        sum(x, na.rm = T)
      })
      
      grad_term1_2 <- sum(term1_2, na.rm = T)
      
      grad_term2_2 <- sum(term2_2, na.rm = T)
      
      - ((grad_term1_1 - grad_term2_1) -
        (grad_term1_2 - grad_term2_2))
      
      
    }  
    
  }
  
  findProposalC2P_KT <- function(d, E, cxt){
    
    # find kt first
    
    list_proposalkt <- findProposalKT(d, E, cxt)
    mu_ktm1_star <- list_proposalkt$ktm1_star
    Sigma_ktm1_star <- list_proposalkt$Sigma_ktm1_star
    
    kt <- c(mu_ktm1_star, - sum(mu_ktm1_star))
    
    mu_c2pm1_star <- matrix(NA, Y, P - 1)
    Sigma_c2pm1_star <- array(NA, dim = c(Y, P - 1, P - 1))
    
    for (t in 1:Y) {
      
      loglik_c2p_current <- loglik_c2pm1_kt(d[,t,], E[,t,], kt[t], m_terms[,t,])
      gr_loglik_c2p_current <- gr_loglik_c2pm1_kt(d[,t,], E[,t,], kt[t], m_terms[,t,])
      
      fitModel <- optim(rep(1, P - 1),
                        loglik_c2p_current,
                        gr_loglik_c2p_current,
                        hessian = T,
                        method = "BFGS")
      mu_c2pm1_star[t,] <- fitModel$par
      Sigma_c2pm1_star[t,,] <- solve(fitModel$hessian)
      
    }
    
    list( "mu_ktm1_star" = mu_ktm1_star,
          "Sigma_ktm1_star" = Sigma_ktm1_star,
          "mu_c2pm1_star" = mu_c2pm1_star,
          "Sigma_c2pm1_star" = Sigma_c2pm1_star)
  }
  
}

# G PARAMETER UPDATE -----

{
  loglik_gc_cxp <- function(d, cxp, m_terms){
    
    function(gc){
      
      numElems <- nrow(d)
      
      gccp <- matrix(gc, numElems, P, byrow = F) * matrix(cxp, numElems, P, byrow = T)
      
      mxtp <- gccp + m_terms
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
      
    }  
    
  }
  
  updateGC_CP <- function(gc_current, c2p, d, E, m_terms){
    
    gc_star <- rep(NA, X + Y - 2)
    sd_star <- rep(NA, X + Y - 2)
    
    for (cc in 2:(X+Y-1)) {
      
      d_current <- sapply(1:P, function(p){
        sdiag(d[,,p], k = cc - X)
      })
      
      m_terms_current <- sapply(1:P, function(p){
        sdiag(m_terms[,,p] + log(E[,,p]), k = cc - X)
      })
      
      if(is.null(nrow(d_current))){
        d_current <- matrix(d_current, 1, P)
        m_terms_current <- matrix(m_terms_current, 1, P)
      }
      
      loglik_gc_current <- loglik_gc_cxp(d_current, c2p[cc,], m_terms_current)
      
      fitModel <- optim(0,
                        loglik_gc_current,
                        hessian = T,
                        method = "BFGS")
      gc_star[cc - 1] <- fitModel$par
      sd_star[cc - 1] <- sqrt(1 / (fitModel$hessian))
      
      loglik_gc_current(gc_star[cc-1]) <   loglik_gc_current(gc_current[cc])
      
      
    }
    
    gc_proposed <- sampleMTconstraint_k(gc_star, diag(sd_star^2, X + Y - 2))
    gc_proposed <- c(0, gc_proposed)
    # gc_proposed <- c(0, gc_star)
    
    # gc_proposed <- rep(0, X + Y - 2)
    # 
    # for (cc in 2:(X+Y-1)) {
    # 
    #   gc_proposed[cc] <- rt2(1, gc_star[cc - 1], sd_star[cc - 1], df = 3)
    # 
    # }
    
    loglik_current <- 0
    loglik_proposal <- 0
    logproposal_current <- 0
    logproposal_proposal <- 0
    for (cc in 2:(X+Y-1)) {
      
      d_current <- sapply(1:P, function(p){
        sdiag(d[,,p], k = cc - X)
      })
      
      m_terms_current <- sapply(1:P, function(p){
        sdiag(m_terms[,,p] + log(E[,,p]), k = cc - X)
      })
      
      if(is.null(nrow(d_current))){
        d_current <- matrix(d_current, 1, P)
        m_terms_current <- matrix(m_terms_current, 1, P)
      }
      
      loglik_gc_current <- loglik_gc_cxp(d_current, c2p[cc,], m_terms_current)
      
      loglik_current <- loglik_current - 
        loglik_gc_current(gc_current[cc])
      loglik_proposal <- loglik_proposal - 
        loglik_gc_current(gc_proposed[cc])
      # print(loglik_gc_current(gc_proposed[cc]) - loglik_gc_current(gc_current[cc]))
      
      logproposal_current <- logproposal_current + 
        log(dt2(gc_current[cc], gc_star[cc - 1], sd_star[cc - 1], df = 3))
      logproposal_proposal <- logproposal_proposal + 
        log(dt2(gc_proposed[cc], gc_star[cc - 1], sd_star[cc - 1], df = 3))
    }
    
    (mh_ratio <- exp(loglik_proposal - loglik_current + 
                      logproposal_current - logproposal_proposal))
    
    if(runif(1) < mh_ratio){
      
      # gc_proposed <- c(0, gc_proposed)
      return(gc_proposed)
      
    }
    
    return(gc_current)
    
  }
  
  loglik_c1p_gc <- function(d, E, gc, m_terms){
    
    function(c1p){
      
      gccp <- array(createTerm4(gc, X, Y), dim = c(X, Y, P)) * 
        aperm(array(c1p, dim = c(P, X, Y)), perm = c(2,3,1))
      
      mxtp <- gccp + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - sum(term1 - term2, na.rm = T)
      
      
    }  
    
  }
  
  updateC1P_GC <- function(c1p_current, gc, d, E, m_terms){
    
    loglik_c1p_current <- loglik_c1p_gc(d, E, gc, m_terms)
    
    fitModel <- optim(rep(1, P),
                      loglik_c1p_current,
                      hessian = T,
                      method = "BFGS")
    c1p_star <- fitModel$par
    Sigma_star <- solve(fitModel$hessian)
    
    c1p_proposed <- sampleMTconstraint_c1p(c1p_star, Sigma_star)
    
    loglik_current <- - loglik_c1p_current(c1p_current)
    loglik_proposed <- - loglik_c1p_current(c1p_star)
    
    logproposal_current <- dmt_cpp(c1p_current, nu = 3, c1p_star, Sigma_star, T)
    logproposal_proposal <- dmt_cpp(c1p_proposed, nu = 3, c1p_star, Sigma_star, T)
    
    mh_ratio <- exp(loglik_proposed - loglik_current + 
                      logproposal_current - logproposal_proposal)
    
    if(runif(1) < mh_ratio){
      
      return(c1p_proposed)
      
    }
    
    return(c1p_current)
    
  }
  
  loglik_c2p_gc <- function(d, gc, m_terms){
    
    function(c2p_c){
      
      numElems <- nrow(d)
      
      gccp_mat <- 
        matrix(gc, numElems, P, byrow = F) * matrix(c2p_c, numElems, P, byrow = T)
      
      mxtp <- gccp_mat + m_terms + log(E)
      
      term1 <- d * mxtp
      term2 <- exp(mxtp)
      
      - (sum(term1 - term2, na.rm = T) + sum(dnorm(c2p_c, 0, sd = 1, log = T)))
      
    }  
    
  }
  
  updateC2P_GC <- function(c2p_current, gc, d, E, m_terms){
    
    c2p_new <- matrix(NA, X + Y - 1, P)
    
    for (cc in 2:(X+Y-1)) {
      
      d_current <- sapply(1:P, function(p){
        sdiag(d[,,p], k = cc - X)
      })
      
      m_terms_current <- sapply(1:P, function(p){
        sdiag(m_terms[,,p] + log(E[,,p]), k = cc - X)
      })
      
      if(is.null(nrow(d_current))){
        d_current <- matrix(d_current, 1, P)
        m_terms_current <- matrix(m_terms_current, 1, P)
      }
      
      loglik_c2p_current <- loglik_c2p_gc(d_current, gc[cc], m_terms_current)
      
      fitModel <- optim(rep(1, P),
                        loglik_c2p_current,
                        hessian = T,
                        method = "BFGS")
      c2p_star <- fitModel$par
      Sigma_star <- solve(fitModel$hessian)
      
      c2p_proposed <- sampleMTconstraint_c1p(c2p_star, Sigma_star)
      
      loglik_current <- - loglik_c2p_current(c2p_current[cc,])
      loglik_proposed <- - loglik_c2p_current(c2p_proposed)
      
      logproposal_current <- dmt_cpp(c2p_current[cc,], nu = 3, c2p_star, Sigma_star, T)
      logproposal_proposal <- dmt_cpp(as.vector(c2p_proposed), nu = 3, c2p_star, Sigma_star, T)
      
      mh_ratio <- exp(loglik_proposed - loglik_current + 
                        logproposal_current - logproposal_proposal)
      
      if(runif(1) < mh_ratio){
        
        c2p_new[cc,] <- c2p_proposed
        
      } else {
        
        c2p_new[cc,] <- c2p_current[cc,]
        
        
      }
      
    }
    
    c2p_new
  }
}

# G MODEL UPDATE ----

loglik_gccp <- function(gc, cxp, d,  E, m_terms){
  
  cxp_term <- createTermG(cxp, X, Y)
  
  gc_term <- createTerm4_cpp(gc, X, Y)
  
  gccp <- 
    array(gc_term, dim = c(X, Y, P)) * cxp_term
  
  mxtp <- gccp + m_terms + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  - sum(term1 - term2, na.rm = T)
  
}

loglik_cpm1_gc <- function(d, gc_term, cxtp, E){
  
  function(cpm1){
    
    cp <- c(cpm1,P - sum(cpm1))
    
    gc_mat <- array(gc_term, dim = c(X, Y, P))
    cp_mat <- aperm(array(cp, dim = c(P, X, Y)), perm = c(2,3,1))
    
    mxtp <- gc_mat * cp_mat + cxtp + log(E)
    
    term1 <- d * mxtp
    term2 <- exp(mxtp)
    
    - sum(term1 - term2, na.rm = T)
    
  }
  
}

findProposalC1P_GC <- function(d, E, cxt){
  
  # find gc first
  
  list_proposalgc <- findProposalGTX(d, E, cxt)
  mu_gtm2_star <- list_proposalgc$gtm2_star
  Sigma_gtm2_star <- list_proposalgc$Sigma_gtm2_star
  
  mu_gc_star <- c(0, mu_gtm2_star, - sum(mu_gtm2_star))
  mu_gc_term <- createTerm4_cpp(mu_gc_star, X, Y)
    
  loglik_cp_current <- loglik_cpm1_gc(d, mu_gc_term, cxt, E)
  
  fitModel <- optim(rep(1, P - 1),
                    loglik_cp_current,
                    hessian = T,
                    method = "BFGS")
  cpm1_star <- fitModel$par
  Sigma_star <- solve(fitModel$hessian)
  
  list(
    "mu_gtm2_star" = mu_gtm2_star,
    "Sigma_gtm2_star" = Sigma_gtm2_star,
    "cpm1_star" = cpm1_star,
    "Sigma_star" = Sigma_star)
  
}

loglik_c2pm1_gc <- function(d, E, kt, m_terms){
  
  function(c2pm1_t){
    
    c2p_t <- c(c2pm1_t, P - sum(c2pm1_t))
    
    ktcp_mat <- matrix(kt, X, P) * matrix(c2p_t, X, P, byrow = T)
    
    mxtp <- ktcp_mat + m_terms + log(E)
    
    term1 <- d * mxtp
    term2 <- exp(mxtp)
    
    - sum(term1 - term2, na.rm = T)
    
    
  }  
  
}

findProposalC2P_GC <- function(d, E, cxt){
  
  # find kt first
  
  list_proposalkt <- findProposalKT(d, E, cxt)
  mu_ktm1_star <- list_proposalkt$ktm1_star
  Sigma_ktm1_star <- list_proposalkt$Sigma_ktm1_star
  
  kt <- c(mu_ktm1_star, - sum(mu_ktm1_star))
  
  mu_c2pm1_star <- matrix(NA, Y, P - 1)
  Sigma_c2pm1_star <- array(NA, dim = c(Y, P - 1, P - 1))
  
  for (t in 1:Y) {
    
    loglik_c2p_current <- loglik_c2pm1_kt(d[,t,], E[,t,], kt[t], m_terms[,t,])
    
    fitModel <- optim(rep(1, P - 1),
                      loglik_c2p_current,
                      hessian = T,
                      method = "BFGS")
    mu_c2pm1_star[t,] <- fitModel$par
    Sigma_c2pm1_star[t,,] <- solve(fitModel$hessian)
    
  }
  
  list( "mu_ktm1_star" = mu_ktm1_star,
        "Sigma_ktm1_star" = Sigma_ktm1_star,
        "mu_c2pm1_star" = mu_c2pm1_star,
        "Sigma_c2pm1_star" = Sigma_c2pm1_star)
}

loglik_gcm2 <- function(d, E, m_terms){
  
  function(param){
    
    gc <- c(0, param, - sum(param))
    
    # loglik_gcm2_cpp(gc, d, E, m_terms)
    gc_term <- createTerm4_cpp(gc, X, Y)

    gccp <- array(gc_term, dim = c(X, Y, P))

    mxtp <- gccp + m_terms + log(E)

    term1 <- d * mxtp
    term2 <- exp(mxtp)

    - sum(term1 - term2, na.rm = T)
    
    
  }  
  
}

gr_loglik_gcm2 <- function(d, E, m_terms){
  
  # ik <- create_ik(X,Y)
  
  function(param){
    
    gtx <- c(0, param, - sum(param))
    
    # gr_loglik_gcm2_cpp(gtx, d, E, m_terms)
    
    param_matrix <- createTerm4_cpp(c(0, param, 0), X, Y)

    sumParams <- sum(param)

    gtx_mat <- createTerm4_cpp(gtx, X, Y)

    gtx_array <- array(gtx_mat, dim = c(X, Y, P))

    # term1 <- sapply(1:P, function(p){
    #   sapply(2:(X+Y-2), function(g){
    #     sum(sdiag(d[,,p], k = g - X), na.rm = T)
    #   })
    # })
    # 
    # term1 <- apply(term1, 1, sum)
    # 
    # term2 <-  sapply(1:P, function(p){
    #   sapply(2:(X+Y-2), function(g){
    #     - sum(sdiag(exp(m_terms[,,p] + gtx_array[,,p] + log(E[,,p])), k = g - X), na.rm = T)
    #   })
    # })
    # 
    # term2 <- apply(term2, 1, sum)

    term12 <- createTerm1(gtx, d, E, m_terms)
    term12 <- apply(term12, 1, sum)
    
    term3 <- sum(d[1,Y,] * (-1), na.rm = T)

    term4 <- - sum(exp(m_terms[1,Y,] + gtx_array[1,Y,] + log(E[1,Y,])) * (-1), na.rm = T)

    - (term12 + term3 + term4)
    # - (term1 + term2 + term3 + term4)
    
  }
  
}

findProposalGTX <- function(d, E, m_terms){
  
  loglik_gtxm2_current <- loglik_gcm2(d, E, m_terms)
  gr_loglik_gtxm2_current <- gr_loglik_gcm2(d, E, m_terms)
  
  fitModel <- optim(rep(0, X + Y - 3),
                    loglik_gtxm2_current,
                    gr_loglik_gtxm2_current,
                    hessian = T,
                    method = "BFGS")
  
  gtm2_star <- fitModel$par
  Sigma_star <- solve(fitModel$hessian)
  
  list("gtm2_star" = gtm2_star,
       "Sigma_gtm2_star" = Sigma_star)
  
}

###

fastdmvnorm <- function(x, mu, Lambda){
  - .5 * t(x - mu) %*% Lambda %*% (x- mu)
}

derl2der2a <- function(a, cx){
  
  - exp(a) * sum(exp(cx), na.rm = T)
  
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

loglik_axcp <- function(ax, cp, d, E, cxt){
  
  axp <- matrix(ax, X, P, F) * matrix(cp, X, P, T) 
  
  axp_mat <- aperm(array(axp, dim = c(X, P, Y)), perm = c(1,3,2))
  
  mxtp <- axp_mat + cxt + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T)
  
}

loglik_axcp2 <- function(ax, cp, d, E, cxt){
  
  axp <- matrix(ax, X, P, F) * matrix(cp, X, P, T) 
  
  axp_mat <- aperm(array(axp, dim = c(X, P, Y)), perm = c(1,3,2))
  
  mxtp <- axp_mat + cxt + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  term1 - term2
  
}

loglik_ktcp <- function(kt, cp, d, E, cxt){
  
  ktp <- matrix(kt, Y, P, F) * matrix(cp, Y, P, T) 
  
  ktp_mat <- aperm(array(ktp, dim = c(Y, P, X)), perm = c(3,1,2))
  
  mxtp <- ktp_mat + cxt + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T)
  
}

# NONPARAMETRIC UPDATE ------



# UPDATE ---------

loglik_axp <- function(axp, d, E, cxtp){
  
  axp_mat <- aperm(array(axp, dim = c(X, P, Y)), perm = c(1,3,2))
  
  mxtp <- axp_mat + cxtp + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T)
  
}



loglik_ktp <- function(ktp, d, E, cxtp){
  
  ktp_mat <- aperm(array(ktp, dim = c(Y, P, X)), perm = c(3,1,2))
  
  mxtp <- ktp_mat + cxtp + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T)
  
}





loglik_c2xp_axc1 <- function(d, E, ax, c1p, cxtp){
  
  function(c2p){
    
    c2p <- matrix(c2p, X - 1, P - 1, byrow = F)
    
    c2P <- sapply(1:(P-1), function(p){
      (c1p[p] * sum(ax) - sum(ax[1:(X-1)] * c2p[1:(X-1),p])) / ax[X]
    })
    
    c2p <- rbind(c2p, c2P)
    c2P <- sapply(1:X, function(i){
      (P - sum(c1p[1:(P-1)] * c2p[i,1:(P-1)])) / c1p[P]
    })
    c2p <- cbind(c2p, c2P)

    axp <- 
      matrix(ax, X, P, byrow = F) *
      matrix(c1p, X, P, byrow = T) * c2p
    
    axp_mat <- aperm(array(axp, dim = c(X, P, Y)), perm = c(1, 3, 2))
    
    mxtp <- axp_mat + cxtp + log(E)
    
    term1 <- d * mxtp
    term2 <- exp(mxtp)
    
    - sum(term1 - term2, na.rm = T)
    
  }
  
}

loglik_kt_cp <- function(d, cp, cxtp){
  
  function(kt){
    terms2 <- kt * matrix(cp, X, P, byrow = T)
    - sum(dpois(d, lambda = exp(terms2 + cxtp), log = T), na.rm = T)  
  }
  
}

loglik_cp_kt <- function(d, kt, cxtp){
  
  function(cp){
    terms2 <- cp * matrix(kt, X, Y, byrow = T)
    - sum(dpois(d, lambda = exp(terms2 + cxtp), log = T), na.rm = T)  
  }
  
}

loglik_cp_gc <- function(d, gc, cxtp){
  
  function(cp){
    terms2 <- cp * matrix(kt, X, P, byrow = F)
    - sum(dpois(d, lambda = exp(terms2 + cxtp), log = T), na.rm = T)  
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
  # print(a[x,])
  term1 <- d[x,,] * mxtp
  term2 <- exp(mxtp)
  
  sum(term1 - term2, na.rm = T) 
}

loglik_LCp_x <- function(x, a, ap, b, k, d, E){
  
  mxtp <- a + matrix(ap[x,], Y, P, byrow = T) + 
    matrix(b[x,], Y, P, byrow = T) * k + log(E[x,,])
  # print((a + matrix(ap[x,], Y, P, byrow = T))[x,])
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
  
  loglik_apm1x_r <- function(x, a, ap, bp, d, E, k, sd_ap){
    
    function(param){
      
      ap[x,] <- c(param, -sum(param))
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - (loglik_LCp_x_r(x, a + ap, bp, k, d, E) + 
           sum(dnorm(ap[x,], mean = 0, sd = sd_ap, log = T)))
      
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
  
  loglik_apm1_r <- function(a, bp, d, E, k){
    
    function(param){
      
      ap <- c(param, -sum(param))
      
      aap <- matrix(a, X, P, byrow = F) + matrix(ap, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(aap, bp, k, d, E)
      
    }
    
  }
  
  loglik_aapm1_r <- function(bp, d, E, k){
    
    function(param){
      
      a <- param[1:X]
      ap <- c(param[X + 1:(P-1)], -sum(param[X + 1:(P-1)]))
      
      aap <- matrix(a, X, P, byrow = F) + matrix(ap, X, P, byrow = T)
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
  
  findProposalApm1Additive <- function(a, b2, d, E, k2){
    
    loglik_ap_current <- loglik_apm1_r(a, b2, d, E, k2)
    # gr_loglik_ap_current <- gr_loglik_ap_r(a, b2, d, E, k2)
    
    # find maximum
    laplace_fit <- optim(
      par = rep(0, P - 1),
      fn = loglik_ap_current,
      method = c("BFGS"),
      # gr = gr_loglik_ap_current,
      hessian = T
    )
    
    ap_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- solve(H_star)  
    
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
  
  findProposalApXm1Nonparametrics <- function(x, a, ap, b2, d, E, k2, sd_ap){
    
    loglik_apx_current <- loglik_apm1x_r(x, a, ap, b2, d, E, k2, sd_ap)
    # gr_loglik_apx_current <- gr_loglik_apx_r(x, a, ap, b2, d, E, k2, sd_ap)
    
    # find maximum
    laplace_fit <- optim(
      # par = ap[x,],
      par = rep(0, P - 1),
      fn = loglik_apx_current,
      method = c("BFGS"),
      # gr = gr_loglik_apx_current,
      hessian = T
    )
    
    ap_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- solve(H_star)
    
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
  
  loglik_bpxm1_r <- function(x, ap, b, bp, d, E, k, sd_bp){
    
    function(param){
      
      bp[x,] <- c(param, - sum(param))
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
  
  loglik_bp_m1_r <- function(ap, b, d, E, kp){
    
    function(param){
      
      bp <- c(param, - sum(param))
      
      bbp <- matrix(b, X, P, byrow = F) + matrix(bp, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(ap, bbp, kp, d, E)
      
    }
    
  }
  
  loglik_bbp_m1_r <- function(ap, d, E, kp){
    
    function(param){
      
      b <- c(param[1:(X-1)], 1 - sum(param[1:(X-1)]))
      bp <- c(param[X + 1:(P-1)], - sum(param[X + 1:(P-1)]))
      
      bbp <- matrix(b, X, P, byrow = F) + matrix(bp, X, P, byrow = T)
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
  
  findProposalBpm1Additive <- function(a2, b, d, E, k2){
    
    loglik_bp_current <- loglik_bp_m1_r(a2, b, d, E, k2)
    # gr_loglik_bp_current <- gr_loglik_bp_r(a2, b, d, E, k2)
    
    # find maximum
    laplace_fit <- optim(
      par = rep(0, P - 1),
      fn = loglik_bp_current,
      method = c("BFGS"),
      # gr = gr_loglik_bp_current,
      hessian = T
    )
    
    bp_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- solve(H_star)  
    
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
  
  findProposalBpXM1Nonparametrics <- function(x, a2, b, bp, d, E, k2, sd_bp){
    
    loglik_bpx_current <- loglik_bpxm1_r(x, a2, b, bp, d, E, k2, sd_bp)
    # gr_loglik_bpx_current <- gr_loglik_bpx_r(x, a2, b, bp, d, E, k2, sd_bp)
    
    # find maximum
    laplace_fit <- optim(
      # par = ap[x,],
      par = rep(0, P - 1),
      fn = loglik_bpx_current,
      method = c("BFGS"),
      # gr = gr_loglik_bpx_current,
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
 
   # loglikelihood of kp for a specific t
  
  loglik_kpt_m1_r <- function(t, ap, bp, k, kp, d, E){
    
    function(param){
      
      kp[t,1:(P-1)] <- param
      kp[t,P] <- - sum(param)
      # - loglik_LCp_t(t - 1, ap, bp, k + kp, d, E)
      - loglik_LCp_t_r(t, ap, bp, k + kp, d, E)
      
    }
    
  }
  
  # gradient of loglikelihood of kp for a specific t
  
  gr_loglik_kpt_m1_r <- function(t, ap, bp, k, kp, d, E){
    
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
  
  loglik_kppm1_r <- function(ap, bp, d, E){
    
    function(param){
      
      k <- c(param[1:(Y-1)], - sum(param[1:(Y-1)]))
      kp <- c(param[Y - 1 + 1:(P-1)], - sum(param[Y - 1 + 1:(P-1)]))
      
      kkp <- 
        matrix(k, Y, P, byrow = F) + 
        matrix(kp, Y, P, byrow = F)
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
  
  # loglikelihood under the additive model for k
  
  loglik_kp_m1_r <- function(ap, bp, d, E, k){
    
    function(param){
      
      kp <- c(param,-sum(param))
      
      kkp <- 
        matrix(k, Y, P, byrow = F) + 
        matrix(kp, Y, P, byrow = T)
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
   
  findProposalKpAdditive_m1 <- function(a2, b2, d, E, k){
    
    loglik_kp_current <- loglik_kp_m1_r(a2, b2, d, E, k)
    # gr_loglik_kp_current <- gr_loglik_kp_r(a2, b2, d, E, k)
    
    # find maximum
    laplace_fit <- optim(
      # par = ap[x,],
      par = rep(0, P - 1),
      fn = loglik_kp_current,
      method = c("BFGS"),
      # gr = gr_loglik_kp_current,
      hessian = T
    )
    
    kp_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- solve(H_star)  
    
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
  
  findProposalKptNonparametrics_m1 <- function(t, a2, b2, d, E, 
                                            k, kp, sd_kt){
    
    loglik_kpt_current <- loglik_kpt_m1_r(t, a2, b2, k, kp, d, E)
    # gr_loglik_kpt_current <- gr_loglik_kpt_m1_r(t, a2, b2, k, kp, d, E)
    
    # find maximum
    laplace_fit <- optim(
      par = kp[t,1:(P-1)],
      # par = rep(0, P),
      fn = loglik_kpt_current,
      # gr = gr_loglik_kpt_current,
      method = c("BFGS"),
      hessian = T
    )
    
    kp_star <- laplace_fit$par
    H_star <- laplace_fit$hessian
    
    Sigma_star <- solve(H_star)
    
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


sampleMVNconstraint_k <- function(mu, Sigma){
  
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

sampleMTconstraint_k <- function(mu, Sigma){
  
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
  
  pxcd <- mrt2(mu1, Sigma1, df = 3)
  
  x <- solve(C) %*% c(pxcd, b)
  
  x
}

sampleMTconstraint_c1p <- function(mu, Sigma){
  
  d <- length(mu)
  
  A <- matrix(1, 1, d)
  b <- d
  
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
  
  pxcd <- mrt2(mu1, Sigma1, df = 3)
  
  x <- solve(C) %*% c(pxcd, b)
  
  x
}



# OLD ----


createC2p <- function(ax, c1p, c2p){
  
  c2P <- sapply(1:(P-1), function(p){
    (c1p[p] * sum(ax) - sum(ax[1:(X-1)] * c2p[1:(X-1),p])) / ax[X]
  })
  
  c2p <- rbind(c2p, c2P)
  c2P <- sapply(1:X, function(i){
    (P - sum(c1p[1:(P-1)] * c2p[i,1:(P-1)])) / c1p[P]
  })
  c2p <- cbind(c2p, c2P)
  
  c2p
  
}

loglik_axp2 <- function(axp, d, E, cxtp){
  
  axp_mat <- aperm(array(axp, dim = c(X, P, Y)), perm = c(1,3,2))
  
  mxtp <- axp_mat + cxtp + log(E)
  
  term1 <- d * mxtp
  term2 <- exp(mxtp)
  
  term1 - term2
  
}

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

loglik_a_r <- function(d, E, ap, b2, k2){
  
  function(param){
    
    a <- param
    
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    - loglik_LCp(a2, b2, k2, d, E)  
    
  }
  
}

gr_loglik_a_r <- function(d, E, ap, b2, k2){
  
  function(param){
    
    a <- param
    
    a2 <- matrix(a, X, P, byrow = F) + ap
    
    - loglik_LCp(a2, b2, k2, d, E)  
    
  }
  
}

loglik_b_r <- function(d, E, a2, bp, k2){
  
  function(param){
    
    b <- param
    
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    - loglik_LCp(a2, b2, k2, d, E)  
    
    # priorTerm <- sum(
    #   dnorm(bx, 0, sd = sd_bx, log = T)
    # )
    
  }
  
}

loglik_bm1_r <- function(d, E, a2, bp, k2){
  
  function(param){
    
    bmx <- param
    
    b <- c(bmx, 1 - sum(bmx))
    
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    priorTerm <- sum(
      dnorm(param, 0, sd = sd_b, log = T)
    )
    # priorTerm <- 0
    
    - (loglik_LCp_r(a2, b2, k2, d, E) + priorTerm) 
    
    
  }
  
}

gr_loglik_bm1_r <- function(d, E, a2, bp, k2){
  
  function(param){
    
    bmx <- param
    
    b <- c(bmx, 1 - sum(bmx))
    
    b2 <- matrix(b, X, P, byrow = F) + bp
    
    a2_current <- aperm(array(a2, dim = c(X, P, Y)), perm = c(1,3,2))
    b2_current <- aperm(array(b2, dim = c(X, P, Y)), perm = c(1,3,2))
    k2_current <- aperm(array(k2, dim = c(Y, P, X)), perm = c(3,1,2))
    
    m_xtp <- a2_current + b2_current * k2_current + log(E)
    
    term2 <- exp(m_xtp)
    
    gr_loglik1 <- (d - term2) * k2_current
    
    gr_loglik1 <- apply(gr_loglik1[1:(X-1),,], 1, function(x){
      sum(x, na.rm = T)
    })
    
    gr_loglik2 <- - sum(((d - term2) * k2_current)[X,,], na.rm = T)
    
    # - (gr_loglik1 + gr_loglik2)
    - (gr_loglik1 + gr_loglik2 - param / sd_b^2)
    
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

loglik_k0m1_r <- function(d, E, a, b, kp){
  
  function(param){
    
    k <- c(param, -sum(param))
    
    k2 <- matrix(k, Y, P, byrow = F) + kp
    
    - loglik_LCp(a, b, k2, d, E)
    
  }
  
}

gr_loglik_k0m1_r <- function(d, E, a, b, kp){
  
  function(param){
    
    k <- c(param, -sum(param))
    
    gr_loglik_km1(a, b, k, kp, d, E)
    
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


findProposalKT_old <- function(d, E, m_terms){
  
  mu_kt_star <- rep(NA, Y)
  sd_kt_star <- rep(NA, Y)
  
  for (t in 1:Y) {
    
    mu_kt_star[t] <- log(sum(d[,t,], na.rm = T) / 
                           sum(exp(m_terms[,t,] + log(E[,t,])), na.rm = T))
    
    hessian_kt_star <- derl2der2a(mu_kt_star[t], m_terms[,t,] + log(E[,t,]))
    sd_kt_star[t] <- sqrt(- 1 / hessian_kt_star)
    
    
  }
  
  
  list("mu_kt_star" = mu_kt_star,
       "sd_kt_star" = sd_kt_star)
  
}



findProposalKTP <- function(d, E, cxt){
  
  mu_ktp_star <- matrix(NA, Y, P)
  sd_ktp_star <- matrix(NA, Y, P)
  
  for (t in 1:Y) {
    for (p in 1:P) {
      
      mu_ktp_star[t,p] <- log(sum(d[,t,p], na.rm = T) / 
                                sum(exp(cxt[,t,p] + log(E[,t,p])), na.rm = T))
      
      hessian_kt_star <- derl2der2a(mu_ktp_star[t,p], cxt[,t,p] + log(E[,t,p]))
      sd_ktp_star[t,p] <- sqrt(- 1 / hessian_kt_star)
      
    }
  }
  
  list("mu_ktp_star" = mu_ktp_star,
       "sd_ktp_star" = sd_ktp_star)
}
