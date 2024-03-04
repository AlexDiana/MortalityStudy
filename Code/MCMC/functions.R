
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

loglik_LCp_r <- function(a, b, k, d, E){
  
  a <- aperm(array(a, dim = c(X, P, Y)), perm = c(1,3,2))
  b <- aperm(array(b, dim = c(X, P, Y)), perm = c(1,3,2))
  k <- aperm(array(k, dim = c(Y, P, X)), perm = c(3,1,2))

  mxtp <- a + b * k
  
  term1 <- d * (mxtp + log(E))
  term2 <- exp(mxtp + log(E))
  
  sum(term1 - term2)
}

loglik_LCp_x_r <- function(x, a, b, k, d, E){
  
  mxtp <- matrix(a[x,], Y, P, byrow = T) + 
    matrix(b[x,], Y, P, byrow = T) * k
  
  term1 <- d[x,,] * (mxtp + log(E[x,,]))
  term2 <- exp(mxtp + log(E[x,,]))
  
  sum(term1 - term2)
}

loglik_LCp_t_r <- function(t, a, b, k, d, E){
  
  mxtp <- a + b * matrix(k[t,], X, P, byrow = T)
  
  term1 <- d[,t,] * (mxtp + log(E[,t,]))
  term2 <- exp(mxtp + log(E[,t,]))
  
  sum(term1 - term2)
}

gr_loglik_ab_x_r <- function(x, a, b, ab, bp, k, d, E){
  
  mxtp <- matrix(a[x] + ap[x,], Y, P, byrow = T) + 
    (matrix(b[x] + bp[x,], Y, P, byrow = T)) * k
  
  grad1 <- sum( d[x,,]  - exp(mxtp + log(E[x,,])))
  grad2 <- sum( k * d[x,,] - exp(mxtp + log(E[x,,])) * k)
  
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
  
  loglik_kt_r <- function(t, d, E, a, b, k, kp){
    
    function(param){
      
      k[t] <- param
      # k2 <- matrix(k, Y, P, byrow = F) + kp
      
      # - loglik_LCp_t(t - 1, a, b, k2, d, E)
      
      mxp <- a + b * matrix(k[t] + kp[t,], X, P, byrow = T)

      term1 <- sum(d[,t,] * (mxp + log(E[,t,])))
      term2 <- sum(exp(mxp + log(E[,t,])))

      - (term1 - term2)
      
    }
    
  }
  
  gr_loglik_kt_r <- function(t, d, E, ap, bp, k, kp){
    
    function(param){
      
      k[t] <- param
      
      - gr_loglik_k_t(t - 1, ap, bp, k, kp, d, E)  
      
    }
    
  }
  
  
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
  
}

# UPDATE OF AP

{
  
  loglik_apx_r <- function(x, a, ap, bp, d, E, k){
    
    function(param){
      
      ap[x,] <- param
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_x_r(x, a + ap, bp, k, d, E)
      
    }
    
  }
  
  gr_loglik_ap_r_fun <- function(x, a, ap, bp, k, d, E){
    
    m_xtp <- matrix(a[x] + ap[x,], Y, P, byrow = T) + 
      matrix(bp[x,], Y, P, byrow = T) * k + log(E[x,,])
    
    term2 <- exp(m_xtp)
    
    gr_loglik <- d[x,,] - term2
    
    apply(gr_loglik, 2, sum)
    
  }
  
  gr_loglik_apx_r <- function(x, a, ap, bp, d, E, k){
    
    function(param){
      
      ap[x,] <- param
      # - gr_loglik_ap_x(x - 1, a, ap, bp, k, d, E)
      - gr_loglik_ap_r_fun(x, a, ap, bp, k, d, E)
      
    }
    
  }  
  
  # loglikelihood under the additive model for a
  
  loglik_ap_r <- function(a, bp, d, E, k){
    
    function(param){
      
      ap <- matrix(param, X, P, byrow = T)
      # - loglik_LCp_x(x - 1, a + ap, bp, k, d, E)
      - loglik_LCp_r(a + ap, bp, k, d, E)
      
    }
    
  }
  
  # gradient of loglikelihood under the additive model for a
  
  gr_loglik_ap_add_r_fun <- function(a, ap, b, k, d, E){
    
    aap <- matrix(a, X, P, byrow = F) + ap
    
    a <- aperm(array(aap, dim = c(X, P, Y)), perm = c(1,3,2))
    b <- aperm(array(b, dim = c(X, P, Y)), perm = c(1,3,2))
    k <- aperm(array(k, dim = c(Y, P, X)), perm = c(3,1,2))
    
    m_xtp <- a + b * k + log(E)
    
    term2 <- exp(m_xtp)
    
    gr_loglik <- d - term2
    
    apply(gr_loglik, 3, sum)
    
  }
  
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
  
  findProposalApXNonparametrics <- function(x, a, ap, b2, d, E, k2){
    
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
    
    list("ap_star" = ap_star,
         "Sigma_star" = Sigma_star)
  }
  
}

# UPDATE OF BP

{
  
  loglik_bpx_r <- function(x, ap, b, bp, d, E, k){
    
    function(param){
      
      bp[x,] <- param
      bbp <- matrix(b, X, P, byrow = F) + bp
      # - loglik_LCp_x(x - 1, ap, bbp, k, d, E)
      - loglik_LCp_x_r(x, ap, bbp, k, d, E)
      
    }
    
  }
  
  gr_loglik_bp_r_fun <- function(x, ap, b, bp, k, d, E){
    
    m_xtp <- matrix(ap[x,], Y, P, byrow = T) + 
      matrix(b[x] + bp[x,], Y, P, byrow = T) * k + log(E[x,,])

    term2 <- exp(m_xtp)
    
    gr_loglik <- d[x,,] * k - term2 * k
    
    apply(gr_loglik, 2, sum)
    
  }
  
  gr_loglik_bpx_r <- function(x, ap, b, bp, d, E, k){
    
    function(param){
      
      bp[x,] <- param
      # - gr_loglik_bp_x(x - 1, ap, b, bp, k, d, E)
      - gr_loglik_bp_r_fun(x, ap, b, bp, k, d, E)
      
    }
    
  }  
  
    
}

# UPDATE OF KP

{
  
  loglik_kpt_r <- function(t, ap, bp, k, kp, d, E){
    
    function(param){
      
      kp[t,] <- param
      # - loglik_LCp_t(t - 1, ap, bp, k + kp, d, E)
      - loglik_LCp_t_r(t, ap, bp, k + kp, d, E)
      
    }
    
  }
  
  gr_loglik_kpt_r_fun <- function(t, ap, bp, k, kp, d, E){
    
    mxtp <- ap + bp * matrix(k[t] + kp[t,], X, P, byrow = T) + log(E[,t,])
    
    term1 <- d[,t,] * bp
    term2 <- exp(mxtp) * bp
    
    apply(term1 - term2, 2, sum)
    
  }  
  
  gr_loglik_kpt_r <- function(t, ap, bp, k, kp, d, E){
    
    function(param){
      
      kp[t,] <- param
      - gr_loglik_kpt_r_fun(t, ap, bp, k, kp, d, E)
      # - gr_loglik_kp_t(t - 1, ap, bp, k, kp, d, E)
      
    }
    
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

