logistic <- function(x){
  1 / (1 + exp(-x))
}

createTerm4 <- function(gtx, X, Y){
  
  term4 <- matrix(0, X, Y)
  for (g in 1:(X+Y-1)) {
    sdiag(term4, k = g - X) <- gtx[g]
  }
  
  term4
}

simulateData <- function(qx, Ex){
  
  dx <- sapply(1:X, function(x){
    sapply(1:Yd, function(t){
      rbinom(1, 1, prob = qx[x,t])
    })
  })
  
  dx
}

# old
{
  
  LC_loglik <- function(d, E){
    
    function(param){
      
      ax <- param[1:X]
      bm1 <- param[X + 1:(X-1)]
      bx <- c(bm1, 1 - sum(bm1))
      km1 <- param[X + X - 1 + 1:(Y-1)]
      kt <- c(km1, - sum(km1))
      
      q <- exp( matrix(ax, X, Y, byrow = F) + 
                  matrix(bx, X, Y, byrow = F) * 
                  matrix(kt, X, Y, byrow = T))
      
      - sum(d * log(q) + (E - d) * log(1 - q))
      
    }
    
  }
  
  
  fitModel_LC <- function(d, E){
    
    ax <- sapply(1:X, function(x){
      mean(log(d[x,] / E[x,]))
    })
    
    modelFit <- optim(
      par = c(ax, rep(0, X + Y - 2)),
      fn = LC_loglik(d, E),
      hessian = T,
      method = "BFGS"
    )
    
    modelFit
    
  }
  
  
  fitRW <- function(kt){
    
    laggedKt <- kt[-1] - kt[-length(kt)]
    alpha <- mean(laggedKt)
    sigma <- sd(laggedKt - alpha)
    
    c("alpha" = alpha,
      "sigma" = sigma)
  }
  
  predictKt <- function(ktY, rwParam, tstar){
    
    alpha <- rwParam[1]
    sigma <- rwParam[2]
    
    kt_star <- rep(0, tstar)
    kt_star[1] <- ktY + alpha
    
    for (t in 1:t_star) {
      
      kt_star[t] <- kt[t]
      
    }
    
  }
  
  modelPred_LC <- function(modelFit, ages){
    
    ax <- modelFit$par[1:X]
    bxm1 <- modelFit$par[X + 1:(X-1)]
    bx <- c(bxm1, 1 - sum(bxm1))
    ktm1 <- modelFit$par[X + X - 1 + 1:(Y-1)]
    kt <- c(ktm1, - sum(ktm1))
    
    # predict kt
    rwParam <- fitRW(kt)
    
    
    qplot(1:Y, kt)
    
    # kt_future <- 
      
  }
}

fitModel_LCstan <- function(d_current, E_current){
  
  X <- nrow(d_current)
  Y <- ncol(d_current)
  
  LC_dat <- list(
    X = X,
    Y = Y,
    d = d_current,
    E = E_current
  )
  
  startAx <- sapply(1:X, function(x){
    mean(log(d_current[x,] / E_current[x,]))
  })
  
  init_fun <- function(...) list(
    ax = startAx,
    bxm1 = rep(0 / (X-1), X - 1), 
    ktm1 = rep(0, Y - 1)
  )
  
  modelLC <- stan_model(file = here("Code","CaseStudy","LC.stan"), verbose = T)
  modelLCrw <- stan_model(file = here("Code","CaseStudy","LC_rw.stan"), verbose = T)
  
  # mcmc_fit <- vb(modelLC, data = LC_dat,
  #                algorithm = "meanfield",
  #                pars = c("ax","bx","kt","alpha","sigmak"),
  #                elbo_samples = 500,
  #                init = init_fun,
  #                tol_rel_obj = 0.00001,
  #                output_samples = 10000)
  mcmc_fit <- sampling(
    # modelLC,
    modelLCrw,
    data = LC_dat,
    init = init_fun,
    pars = c("ax","bx","kt","alpha","sigmak"),
    # pars = c("ax","bx","kt"),
    chains = 1,
    iter = 5000,
    verbose = T)
  
  matrix_of_draws <- as.matrix(mcmc_fit)
  
  matrix_of_draws
  
}

modelPred_LC_stan <- function(modelFit, Ycurrent, tstar){
  
  # modelFit <- matrix_of_draws
  
  ax_output <- modelFit[,1:X]
  bx_output <- modelFit[,X + 1:X]
  kt_output <- modelFit[,2 * X + 1:Ycurrent]
  alpha_output <- modelFit[,2 * X + Ycurrent + 1]
  sigmak_output <- modelFit[,2 * X + Ycurrent + 2]
  
  #
  {
    # cbind(t(apply(ax_output, 2, function(x){
    #        quantile(x, probs = c(0.025, 0.975))
    #    })), ax_true) %>% 
    #   ggplot(aes(x = 1:X, 
    #              ymin = `2.5%`,
    #              ymax = `97.5%`,
    #              y = ax_true)) + 
    #   geom_errorbar() + geom_point()
  }
  
  niter <- nrow(ax_output)
  
  kt_star_sims <- matrix(NA, niter, tstar)
  qx_star <- array(NA, dim = c(niter, X, tstar))
  
  for (iter in 1:niter) {
    
    kt_star_sims[iter,1] <- kt_output[iter,Ycurrent] + 
      alpha_output[iter] + rnorm(1, sd = sigmak_output[iter])
    
    for (x in 1:X) {
      qx_star[iter,x,1] <- 
        logistic(ax_output[iter,x] + bx_output[iter,x] * kt_star_sims[iter,1])
    }
    
    for (t in 2:tstar) {
      kt_star_sims[iter,t] <- kt_star_sims[iter,t - 1] + 
        alpha_output[iter] + rnorm(1, sd = sigmak_output[iter])
      
      for (x in 1:X) {
        qx_star[iter,x,t] <- 
          logistic(ax_output[iter,x] + bx_output[iter,x] * kt_star_sims[iter,t])
      }
      
    }
    
  }
  
  qx_star
  
}

fitModel_M5stan <- function(dx_current, Ex_current){
  
  X <- nrow(dx_current)
  Y <- ncol(dx_current)
  
  LC_dat <- list(
    X = X,
    Y = Y,
    d = dx_current,
    E = Ex_current
  )
  
  startAx <- sapply(1:X, function(x){
    mean(log(d[x,] / E[x,]))
  })
  
  init_fun <- function(...) list(
    ax = startAx,
    ktm1 = rep(0, Y - 1),
    gtxm1XY = rep(0, X + Y - 3)
  )
  
  modelM5rw <- stan_model(file = here("Code","CaseStudy","M5_rw.stan"), verbose = T)
  
  mcmc_fit <- sampling(modelM5rw, data = LC_dat,
                       init = init_fun,
                       pars = c("ax","kt","gtx","alpha","sigmak"),
                       # pars = c("ax","kt","alpha","sigmak"),
                       # pars = c("ax","bx","kt"),
                       chains = 1,
                       iter = 5000,
                       verbose = T)
  
  matrix_of_draws <- as.matrix(mcmc_fit)
  
  matrix_of_draws
    
}

modelPred_M5_stan <- function(modelFit, Ycurrent, tstar){
  
  ax_output <- modelFit[,1:X]
  kt_output <- modelFit[,X + 1:Ycurrent]
  gtx_output <- modelFit[,X + Ycurrent + 1:(X+Ycurrent-1)]
  alpha_output <- modelFit[,X + Ycurrent + X + Ycurrent - 1 + 1]
  sigmak_output <- modelFit[,X + Ycurrent + X + Ycurrent - 1 + 2]
  
  # cbind(t(apply(kt_output, 2, function(x){
  #   quantile(x, probs = c(0.025, 0.975))
  # })), True = kt_true[1:ncol(kt_output)]) %>%
  #   ggplot(aes(x = 1:ncol(kt_output),
  #              ymin = `2.5%`,
  #              ymax = `97.5%`,
  #              y = True)) +
  #   geom_errorbar() + geom_point()
  
  niter <- nrow(ax_output)
  
  kt_star_sims <- matrix(NA, niter, tstar)
  qx_star <- array(NA, dim = c(niter, X, tstar))
  
  for (iter in 1:niter) {
    
    kt_star_sims[iter,1] <- kt_output[iter,Yd] + 
      alpha_output[iter] + rnorm(1, sd = sigmak_output[iter])
    
    for (x in 1:X) {
      qx_star[iter,x,1] <- 
        logistic(ax_output[iter,x] + 
                   kt_star_sims[iter,1] + gtx_output[iter,1 - t + Ycurrent])
    }
    
    for (t in 2:tstar) {
      kt_star_sims[iter,t] <- kt_star_sims[iter,t - 1] + 
        alpha_output[iter] + rnorm(1, sd = sigmak_output[iter])
      
      for (x in 1:X) {
        qx_star[iter,x,t] <- 
          logistic(ax_output[iter,x] + 
                     kt_star_sims[iter,t] + gtx_output[iter,x - t + Ycurrent])
      }
      
    }
    
  }
  
  qx_star
  
}
