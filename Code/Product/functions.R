
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

