
library(ggplot2)

X <- 20
Y <- 5
ages <- seq(41, 60, length.out = X)
years <- 2000 + 1:Y

# Choice 1 (Base age effect)

cx <- (ages - mean(ages)) / (mean(ages) - ages[1])

{
  # choice1 <- 1 # unconstrainted
  choice1 <- 2 # line

  if(choice1 == 2){
    term1 <- -3 + 1 * cx
  } else {
    ax <- seq(-5, -0, length.out = X) + rnorm(X, sd = .15)
    term1 <- ax
  }
}

k1t <- seq(2, -2, length.out = Y)
k1t <- k1t - mean(k1t)

# Choice 2

{
  # choice2 <- 1 # 0
  choice2 <- 2 # 1
  # choice2 <- 3 # k^2t

  if(choice2 == 1){
    term2 <- rep(0, Y)
  } else if(choice2 == 2){
    # term2 <- rep(1, Y)
    term2 <- rep(.25, Y)
  } else if(choice2 == 3){
    k2t <- seq(0, 1, length.out = Y) + rnorm(Y, sd = .2)
    k2t_tilde <- k2t * k1t
    if(choice1 == 1 | choice2 == 1){
      k2t_tilde <- k2t_tilde - mean(k2t_tilde)
    }
    k2t <- k2t_tilde / k1t
    k2t[is.nan(k2t)] <- 0
    term2 <- k2t
  }
  
}

# Choice 3 (Varying year-age effect)

{
  # choice3 <- 0 # 0
  choice3 <- 1 # cx
  # choice3 <- 2 # bx

  if(choice3 == 0){
    term3 <- rep(0, X)
  } else if(choice3 == 1){
    b_bar <- .5
    a_bar <- 1
    term3 <- (a_bar + b_bar * cx)
  } else if(choice3 == 2){
    bx <- seq(0.2, 1, length.out = X) + rnorm(X, sd = .05)
    bx <- bx - mean(bx)
    # bx <- bx / (-bx[1])
    term3 <- bx
  }

}

# Choice 4 (Cohort effect)

{
  choice4 <- 1 # no cohort effect
  # choice4 <- 2 # cohort effect

  if(choice4 == 2){
    # gtx <- seq(0, 1, length.out = X + Y) + rnorm(X + Y, sd = .1)
    gtx <- rnorm(X + Y, sd = .075)
    # gtx <- rnorm(X + Y, sd = .15)
    term4 <- gtx
  } else {
    term4 <- rep(0, X + Y)
  }

}

data <- expand.grid(x = 1:X, t = 1:Y)

mxtp <- apply(data, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  # term1[x]
  # term1[x] + k1t[t]
  term1[x] + k1t[t] * (1 + term2[t] * term3[x]) + term4[t - x + X + 1]
  # term1[x] + k1t[t] * (1 + term2[t] * term3[x]) + term4[t - x + X + 1]
})

data$x <- ages[data$x]
data$t <- years[data$t]

df <- data.frame(Age = data$x,
                 Year = factor(data$t),
                 m = mxtp)

ggplot(data = df, aes(x = Age,
                      y = m,
                      color = Year,
                      group = Year)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:X, labels = ages) +
  # ylim(c(-10,0))+
  theme_bw() + xlab("Age") + ylab("Log mortality")

{
  
  df$m2 <- df$m + rnorm(nrow(df), sd = .1)
  
  ggplot(data = df, aes(x = Age,
                        y = m2,
                        color = Year,
                        group = Year)) +
    # geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:X, labels = ages) +
    # ylim(c(-10,0))+
    theme_bw() + xlab("Age") + ylab("Log mortality")
  
  k1t0 <- c(2.4, 1.1, 0, -1.1, -2.4)
  
  data <- expand.grid(x = 1:X, t = 1:Y)
  df$m0 <- apply(data, 1, function(dat){
    x <- dat[1]
    t <- dat[2]
    # term1[x]
    term1[x] + k1t0[t]
    # term1[x] + k1t[t] * (1 + term2[t] * term3[x]) + term4[t - x + X + 1]
  })
  
  df$CI1_m1 <- df$m0 - .23
  df$CI2_m1 <- df$m0 + .23
  
  df$CI1_m2 <- df$m - .12
  df$CI2_m2 <- df$m + .12
  
  ggplot(data = df, aes(x = Age,
                        y = m2,
                        ymin = CI1_m1,
                        ymax = CI2_m1,
                        # ymin = CI1_m2,
                        # ymax = CI2_m2,
                        color = Year,
                        group = Year,
                        fill = Year)) +
    geom_ribbon(alpha = .5) +
    # geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:X, labels = ages) +
    ylim(c(-7,1))+
    theme_bw() + xlab("Age") + ylab("Log mortality")
  
  ggsave("plot2.jpeg")
  
  }



# data.frame(Age = factor(data$x),
#                  Year = data$t,
#                  m = mxtp) %>% 
#   filter(Age %in% 40:45) %>% 
# ggplot(aes(x = Year,
#                       y = m,
#                       color = Age,
#                       group = Age)) +
#   geom_line() +
#   geom_point(size = 2) +
#   scale_x_continuous(breaks = 1:Y, labels = years) +
#   # ylim(c(-10,0))+
#   theme_bw() + xlab("Years") + ylab("Log mortality")
