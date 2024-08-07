library(here)
library("rstan") 
library(tidyverse)
library(mgcv)
library(coda)
library(wiqid)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source(here("Code","CaseStudy","function.R"))

# SIMULATE DATA -----

Y <- 15
Yd <- 10
Ys <- Y - Yd

X <- 10

# model LC

{
  
  # a
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
    
    param_rw <- fitRW(kt)
    
    for (t in 2:Y) {
      kt[t] <- kt[t - 1] + param_rw[1] + param_rw[2] * rnorm(1)
    }
    
    # qplot(1:Y, kt)
  }
  
  meankt <- mean(kt)
  kt <- kt - meankt
  ax <- ax + bx * meankt
  sumbx <- sum(bx)
  bx <- bx / sumbx
  kt <- kt * sumbx
  
  ax_true <- ax
  bx_true <- bx
  kt_true <- kt
  
  q_true <- logistic( matrix(ax, X, Y, byrow = F) + 
                        matrix(bx, X, Y, byrow = F) * 
                        matrix(kt, X, Y, byrow = T))  
}

# model 5

{
  # a
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
  
  # k
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
  
  kt <- kt  / 5
  
  kt <- kt[1:Y]
  
  param_rw <- fitRW(kt)
  
  for (t in 2:Y) {
    kt[t] <- kt[t - 1] + param_rw[1] + param_rw[2] * rnorm(1)
  }
  
  kt <- kt - mean(kt)
  
  # g
  
  gtx <- rnorm(X + Y - 1)
  gtx[1] <- 0
  gtx[2:(X+Y-1)] <- gtx[2:(X+Y-1)] - mean(gtx[2:(X+Y-1)])
  
  gtx <- rep(0, X + Y - 1)
  
  ax_true <- ax
  kt_true <- kt
  gtx_true <- gtx
  
  gtx_term <- createTerm4(gtx_true, X, Y)
  
  q_true <- logistic( matrix(ax, X, Y, byrow = F) +  
                        matrix(kt, X, Y, byrow = T) + 
                        gtx_term) 
  
}

E <- matrix(rpois(X * Y, lambda = 5000), X, Y)

d <- t(sapply(1:X, function(x){
  sapply(1:Y, function(t){
    rbinom(1, E[x,t], prob = q_true[x,t])
  })
}))

# LOOP ------

npred <- 1000
qx_pred1 <- matrix(NA, Ys, npred)
qx_pred2 <- matrix(NA, Ys, npred)

for (t in 1:Ys) {
  
  # fit model on previous years
  {
     d_current <- d[,1:(Yd + t - 1)]
     E_current <- E[,1:(Yd + t - 1)]
     
     modelFit_LC <- fitModel_LCstan(d_current, E_current)
     modelFit_M5 <- fitModel_M5stan(d_current, E_current)
     
     # apply(modelFit_LC[,1:X], 2, function(x){
     #   quantile(x, probs = c(0.025, 0.975))
     # })
     # ax_true
     # 
     
     # cbind(t(apply(kt_output, 2, function(x){
     #   quantile(x, probs = c(0.025, 0.975))
     # })), kt_true[1:ncol(kt_output)]) %>%
     #   ggplot(aes(x = 1:ncol(kt_output),
     #              ymin = `2.5%`,
     #              ymax = `97.5%`,
     #              y = kt_true)) +
     #   geom_errorbar() + geom_point()
     
  }
  
  # predict qx on future years
  {
    qx_star_LC <- modelPred_LC_stan(modelFit_LC, ncol(d_current), Ys - t + 1)
    qx_star_M5 <- modelPred_M5_stan(modelFit_M5, ncol(d_current), Ys - t + 1)
  }
  
  # compare predictions
  {
    qx_star_long <-
      expand.grid(x = 1:X, t = Yd + 1:Ys) %>%
      as.data.frame

    q_true_xt <- apply(qx_star_long, 1, function(x){
      x_current <- x[1]
      t_current <- x[2]
      q_true[x_current,t_current]
    })
    
    q_CI_LC <- apply(qx_star_long, 1, function(x){
        x_current <- x[1]
        t_current <- x[2]
        quantile(qx_star_LC[,x_current,t_current - Yd], probs = c(0.025, 0.975))
      })
    
    q_CI_M5 <- apply(qx_star_long, 1, function(x){
        x_current <- x[1]
        t_current <- x[2]
        quantile(qx_star_M5[,x_current,t_current - Yd], probs = c(0.025, 0.975))
      })
    
    qx_star_long %>%
      mutate(True = q_true_xt,
             l_CI = q_CI_M5[1,],
             u_CI = q_CI_M5[2,]) %>%
      filter(x %in% c(1,10)) %>%
      mutate(x = factor(x)) %>%
      ggplot(aes(x = t,
                 y = True,
                 ymin = l_CI,
                 ymax = u_CI,
                 group = x,
                 fill = x)) +
      geom_ribbon(alpha = .2) +
      geom_point() + theme_bw()
  }
  
  

  
}