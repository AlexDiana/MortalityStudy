
library(ggplot2)

X <- 20
Y <- 5
ages <- seq(41, 60, length.out = X)
years <- 2000 + 1:Y

# baseline 
mu0 <- -3

# Choice 1 (Base age effect)

{
  # choice1 <- "U" # unconstrainted
  choice1 <- "L" # line
  
  ax <- seq(-2, 1, length.out = X) + rnorm(X, sd = .25)
  cx <- -3 + .1 * (ages - mean(ages))
  if(choice1 == "L"){
    term1 <- cx
  } else {
    term1 <- ax
  }
}

# Choice 2 (Additive year effect)

{
  choice2 <- "A" # additive
  # choice2 <- "NA" # non additive
  
  k1t <- seq(-1.5, 1.5, length.out = Y)
  if(choice2 == "A"){
    term2 <- k1t
  } else {
    term2 <- rep(0, Y)
  }
  
}

# Choice 3 (Varying year-age effect)

{
  choice3 <- 3 # 1, varying line, 2 lee carter, 3 no effect
 
  bx <- seq(.5, 0.0, length.out = X) + rnorm(X, sd = .05)
  c2x <- -.04 * (ages - mean(ages))
    
  k2t <- seq(-11, -6, length.out = Y)
  k2t - mean(k2t)
  
  if(choice3 == 1){
    term3 <- matrix(c2x, X, Y, byrow = F) * matrix(k2t, X, Y, byrow = T)
  } else if(choice3 == 2){
    term3 <- matrix(bx, X, Y, byrow = F) * matrix(k2t, X, Y, byrow = T)
  } else if(choice3 == 3){
    term3 <- matrix(0, X, Y)  
  }
  
}

# Choice 4 (Cohort effect)

{
  # choice4 <- "C" # cohort effect
  choice4 <- "NC" # no cohort effect
  
  gtx <- rnorm(X + Y, sd = .1) 
  # gtx <- seq(0, 2, length.out = X + Y) 
  
  if(choice4 == "C"){
    term4 <- gtx
  } else {
    term4 <- rep(0, X + Y)
  }
  
}

data <- expand.grid(x = 1:X, t = 1:Y)

mxtp <- apply(data, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  # mu0 + term1[x] #+ term2[t] + term3[x,t] + term4[t - x + X + 1]
  # mu0 + term1[x] + term2[t] #+ term3[x,t] + term4[t - x + X + 1]
  mu0 + term1[x] + term2[t] + term3[x,t] + term4[t - x + X + 1]
})

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


# df <- data.frame(Age = factor(data$x),
#                  Year = data$t,
#                  m = mxtp)
# 
# ggplot(data = df, aes(x = Year,
#                       y = m,
#                       color = Age,
#                       group = Age)) + 
#   geom_line() +
#   geom_point(size = 2) + 
#   scale_x_continuous(breaks = 1:Y, labels = years) +
#   # ylim(c(-10,0))+
#   theme_bw() + xlab("Years") + ylab("Log mortality")
