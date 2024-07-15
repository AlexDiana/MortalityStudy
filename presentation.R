library(ggplot2)

X <- 1:40
ages <- 40 + X

abc <- c(-3, 0, 0.00125)

X_cov <- cbind(1, X, X^2)

y_true <- X_cov %*% abc

y_mort <- y_true + rnorm(length(X), sd = .1)

qplot(X, y_mort) + 
  theme_bw() + 
  xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
  ylab("Log-Mortality")


# additional models

# model 1
{
  X_cov <- cbind(1, X, X^2)
  abc <- c(-3.25, 1/20, 0)
  y_mort1 <- X_cov %*% abc
  abc_err <- abc + rnorm(length(abc), sd = c(.1, .005, exp(-50))) 
  y_mort1_err <- X_cov %*% abc_err
  
  predictions1_CI1 <- y_mort1 - .1
  predictions1_CI2 <- y_mort1 + .1
  predictions1_err_CI1 <- y_mort1_err - .1
  predictions1_err_CI2 <- y_mort1_err + .1  
}

# model 2
{
  set.seed(1)
  abc <- c(-3, -139/6000, 179 / 60000, -11 / 300000)
  X_cov <- cbind(1, X, X^2, X^3)
  y_mort2 <- X_cov %*% abc
  
  abc_err <- abc + rnorm(length(abc), sd = c(.1, exp(-4), exp(-15), exp(-15))) 
  y_mort2_err <- X_cov %*% abc_err
  
  predictions2_CI1 <- y_mort2 - .1
  predictions2_CI2 <- y_mort2 + .1
  predictions2_err_CI1 <- y_mort2_err - .1
  predictions2_err_CI2 <- y_mort2_err + .1  
}

# model 3
{
  set.seed(20)
  abc <- c(-3, 0, 0.00125)
  X_cov <- cbind(1, X, X^2)
  y_mort3 <- X_cov %*% abc
  abc_err <- abc + rnorm(length(abc), sd = c(.1, exp(-4), exp(-15))) 
  y_mort3_err <- X_cov %*% abc_err
  
  predictions3_CI1 <- y_mort3 - .1
  predictions3_CI2 <- y_mort3 + .1
  predictions3_err_CI1 <- y_mort3_err - .1
  predictions3_err_CI2 <- y_mort3_err + .1  
}

# model 4
{
  y_mort4 <- y_mort
  y_mort4_err <- y_mort + rnorm(length(y_mort), sd = .2)
  
  predictions4_CI1 <- y_mort4 - .1
  predictions4_CI2 <- y_mort4 + .1
  predictions4_err_CI1 <- y_mort4_err - .1
  predictions4_err_CI2 <- y_mort4_err + .1  
}

# prediction error
{
  
  predictions_CI1 <- y_true - .02
  predictions_CI2 <- y_true + .02
  
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions_CI1,
                    ymax = predictions_CI2), color = "grey", alpha = .2) + 
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11)) 
  
  predictions_CI1 <- y_true - .25
  predictions_CI2 <- y_true + .25
  
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions_CI1,
                    ymax = predictions_CI2), color = "grey", alpha = .2) + 
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality")  + 
    theme(axis.text = element_text(size = 11))
}

# model error
{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions1_CI1,
                    ymax = predictions1_CI2), fill = "red", alpha = .3) +
    # geom_ribbon(data = NULL, 
    #             aes(x = X, 
    #                 ymin = predictions2_CI1,
    #                 ymax = predictions2_CI2), fill = "blue", alpha = .2) +
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions3_CI1,
                    ymax = predictions3_CI2), fill = "black", alpha = .2) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))
  
  ggsave("plot1.jpeg")
}

# rjmcmc 1

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions1_CI1,
                    ymax = predictions1_CI2), fill = "red", alpha = .3) +
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions2_CI1,
                    ymax = predictions2_CI2), fill = "blue", alpha = .2) +
    # geom_text(data = NULL, aes(label = "Model 1", x = 30, y = -1.6), size = 5, color = "red") + 
    # geom_text(data = NULL, aes(label = "Model 2", x = 31, y = -2), size = 5, color = "blue") + 
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 2

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions2_err_CI1,
                    ymax = predictions2_err_CI2), fill = "blue", alpha = .4) +
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions2_CI1,
                    ymax = predictions2_CI2), fill = "blue", alpha = .2) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 3

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions2_err_CI1,
                    ymax = predictions2_err_CI2), fill = "blue", alpha = .4) +
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions3_CI1,
                    ymax = predictions3_CI2), fill = "green", alpha = .4) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 4

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions3_CI1,
                    ymax = predictions3_CI2), fill = "green", alpha = .4) +
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions3_err_CI1,
                    ymax = predictions3_err_CI2), fill = "green", alpha = .2) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 5

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions3_err_CI1,
                    ymax = predictions3_err_CI2), fill = "green", alpha = .4) +
    # geom_ribbon(data = NULL, 
    #             aes(x = X, 
    #                 ymin = predictions3_err_CI1,
    #                 ymax = predictions3_err_CI2), fill = "green", alpha = .2) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 6

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions3_err_CI1,
                    ymax = predictions3_err_CI2), fill = "green", alpha = .4) +
    geom_ribbon(data = NULL,
                aes(x = X,
                    ymin = predictions4_CI1,
                    ymax = predictions4_CI2), fill = "yellow", alpha = .5) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 7

{
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions4_err_CI1,
                    ymax = predictions4_err_CI2), fill = "yellow", alpha = .2) +
    geom_ribbon(data = NULL,
                aes(x = X,
                    ymin = predictions4_CI1,
                    ymax = predictions4_CI2), fill = "yellow", alpha = .5) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}

# rjmcmc 8

{
  ggplot() + 
    geom_ribbon(data = NULL, 
                aes(x = X, 
                    ymin = predictions4_CI1,
                    ymax = predictions4_CI2), fill = "yellow", alpha = .8) +
    geom_point(data = NULL, 
               aes(x = X, 
                   y = y_mort)) + 
    # geom_ribbon(data = NULL,
    #             aes(x = X,
    #                 ymin = predictions4_CI1,
    #                 ymax = predictions4_CI2), fill = "blue", alpha = .5) +
    theme_bw() + 
    xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
    ylab("Log-Mortality") + 
    theme(axis.text = element_text(size = 11))  
  
  ggsave("plot1.jpeg")
}


# model and precition error

abc <- c(-3.25, 1/20, 0)

abc_err <- abc + rnorm(length(abc), sd ) 

X_cov <- cbind(1, X, X^2)

y_mort1 <- X_cov %*% abc

abc <- c(-3, -139/6000, 179 / 60000, -11 / 300000)

X_cov <- cbind(1, X, X^2, X^3)

y_mort2 <- X_cov %*% abc

predictions1_CI1 <- y_mort1 - .1
predictions1_CI2 <- y_mort1 + .1
predictions2_CI1 <- y_mort2 - .1
predictions2_CI2 <- y_mort2 + .1

ggplot() + 
  geom_point(data = NULL, 
             aes(x = X, 
                 y = y_mort)) + 
  geom_ribbon(data = NULL, 
              aes(x = X, 
                  ymin = predictions1_CI1,
                  ymax = predictions1_CI2), fill = "red", alpha = .3) +
  geom_ribbon(data = NULL, 
              aes(x = X, 
                  ymin = predictions2_CI1,
                  ymax = predictions2_CI2), fill = "blue", alpha = .2) +
  theme_bw() + 
  xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
  ylab("Log-Mortality") + 
  theme(axis.text = element_text(size = 11))

# model error

abc <- c(-3.25, 1/20, 0)

X_cov <- cbind(1, X, X^2)

y_mort1 <- X_cov %*% abc

abc <- c(-3, -139/6000, 179 / 60000, -11 / 300000)

X_cov <- cbind(1, X, X^2, X^3)

y_mort2 <- X_cov %*% abc

predictions1_CI1 <- y_mort1 - .1
predictions1_CI2 <- y_mort1 + .1
predictions2_CI1 <- y_mort2 - .1
predictions2_CI2 <- y_mort2 + .1

ggplot() + 
  geom_point(data = NULL, 
             aes(x = X, 
                 y = y_mort)) + 
  geom_ribbon(data = NULL, 
              aes(x = X, 
                  ymin = predictions1_CI1,
                  ymax = predictions1_CI2), fill = "red", alpha = .3) +
  geom_ribbon(data = NULL, 
              aes(x = X, 
                  ymin = predictions2_CI1,
                  ymax = predictions2_CI2), fill = "blue", alpha = .2) +
  theme_bw() + 
  xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
  ylab("Log-Mortality") + 
  theme(axis.text = element_text(size = 11))

# model error with true

predictions3_CI1 <- y_true - .1
predictions3_CI2 <- y_true + .1

predictions4_CI1 <- y_mort - .1
predictions4_CI2 <- y_mort + .1

ggplot() + 
  # geom_ribbon(data = NULL,
  #             aes(x = X,
  #                 ymin = predictions3_CI1,
  #                 ymax = predictions3_CI2), fill = "grey", alpha = .6) +
  geom_ribbon(data = NULL,
              aes(x = X,
                  ymin = predictions4_CI1,
                  ymax = predictions4_CI2), fill = "grey", alpha = .6) +
  # geom_ribbon(data = NULL, 
              # aes(x = X, 
                  # ymin = predictions1_CI1,
                  # ymax = predictions1_CI2), fill = "red", alpha = .3) +
  # geom_ribbon(data = NULL, 
  #             aes(x = X, 
  #                 ymin = predictions2_CI1,
  #                 ymax = predictions2_CI2), fill = "blue", alpha = .2) +
  geom_point(data = NULL, 
             aes(x = X, 
                 y = y_mort)) + 
  theme_bw() + 
  xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
  ylab("Log-Mortality") + 
  theme(axis.text = element_text(size = 11))

# rjmcmc

predCurrent1 <- predictions3_CI1
predCurrent2 <- predictions3_CI2

ggplot() + 
  # geom_ribbon(data = NULL,
  #             aes(x = X,
  #                 ymin = predictions3_CI1,
  #                 ymax = predictions3_CI2), fill = "grey", alpha = .6) +
  geom_ribbon(data = NULL,
              aes(x = X,
                  ymin = predCurrent1,
                  ymax = predCurrent2), fill = "grey", alpha = .5) +
 # geom_ribbon(data = NULL, 
              # aes(x = X, 
                  # ymin = predictions1_CI1,
                  # ymax = predictions1_CI2), fill = "red", alpha = .3) +
  # geom_ribbon(data = NULL, 
  #             aes(x = X, 
  #                 ymin = predictions2_CI1,
  #                 ymax = predictions2_CI2), fill = "blue", alpha = .2) +
  
  geom_point(data = NULL, 
             aes(x = X, 
                 y = y_mort)) + 
  ylim(c(-3.4, -.7)) + 
  theme_bw() + 
  xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
  ylab("Log-Mortality") + 
  theme(axis.text = element_text(size = 11))

ggsave("plot1.jpeg")

# predictions

subExcl <- 30:33

predCurrent1 <- predictions2_CI1
predCurrent2 <- predictions2_CI2

ggplot() + 
  # geom_ribbon(data = NULL,
  #             aes(x = X,
  #                 ymin = predictions3_CI1,
  #                 ymax = predictions3_CI2), fill = "grey", alpha = .6) +
  geom_ribbon(data = NULL,
              aes(x = X[-subExcl],
                  ymin = predictions1_CI1[-subExcl],
                  ymax = predictions1_CI2[-subExcl]), fill = "blue", alpha = .3) +
  geom_ribbon(data = NULL,
              aes(x = X[-subExcl],
                  ymin = predictions2_CI1[-subExcl],
                  ymax = predictions2_CI2[-subExcl]), fill = "red", alpha = .3) +
  geom_ribbon(data = NULL,
              aes(x = X[-subExcl],
                  ymin = predictions3_CI1[-subExcl],
                  ymax = predictions3_CI2[-subExcl]), fill = "yellow", alpha = .3) +
  geom_ribbon(data = NULL,
              aes(x = X[-subExcl],
                  ymin = predictions4_CI1[-subExcl],
                  ymax = predictions4_CI2[-subExcl]), fill = "green", alpha = .3) +
  geom_rect(aes(xmin=30, xmax=33, ymin=-Inf, ymax=Inf), fill = "black",
            fill = "black", alpha = .2) +
  # geom_ribbon(data = NULL, 
              # aes(x = X, 
                  # ymin = predictions1_CI1,
                  # ymax = predictions1_CI2), fill = "red", alpha = .3) +
  # geom_ribbon(data = NULL, 
  #             aes(x = X, 
  #                 ymin = predictions2_CI1,
  #                 ymax = predictions2_CI2), fill = "blue", alpha = .2) +
  geom_point(data = NULL, 
             aes(x = X[-subExcl], 
                 y = y_mort[-subExcl])) + 
  theme_bw() + 
  xlab("Age") + scale_x_continuous(breaks = X, labels = ages) + 
  ylab("Log-Mortality") + 
  theme(axis.text = element_text(size = 11))

ggsave("plot1.jpeg")
