X <- 5
TT <- 20
years <- seq(41, 60, length.out = TT) 
P <- 2

ax <- seq(-5, -2, length.out = X)
c1p <- rnorm(P, sd = .5)
c2p <- rnorm(P, sd = .1)
bx <- seq(.04, 0.02, length.out = X)
# bx <- seq(.2, .04, length.out = X)
kt <- seq(30, -40, length.out = TT)
c3p <- rnorm(P, sd = .015)
c3p <- c3p - mean(c3p)
c4p <- rnorm(P, sd = 10) 
c4p <- c4p - mean(c4p)

data <- expand.grid(x = 1:X, t = 1:TT, p = 1:P)

mxtp <- apply(data, 1, function(dat){
  x <- dat[1]
  t <- dat[2]
  p <- dat[3]
  ax[x] + bx[x] * kt[t]
  ax[x] + bx[x] * kt[t] + c1p[p] # (1A)
  # ax[x] * (1 + c2p[p]) + bx[x] * kt[t] # (1M)
  # ax[x] + (bx[x] + c3p[p]) * kt[t] # (2.1A)
  ax[x] + bx[x] * (kt[t] + c4p[p]) # (2.2A)
})

data$p <- paste0("Product ", 1:P)[data$p]

df <- data.frame(Year = factor(data$x),
                 t = data$t,
                 p = data$p,
                 m = mxtp)

library(ggplot2)

ggplot(data = df, aes(x = t,
                      y = m,
                      color = Year,
                      group = Year)) + 
  geom_line() + 
  facet_grid(cols = vars(p)) + 
  scale_x_continuous(breaks = 1:TT, labels = years) +
  theme_bw() + xlab("Age")


# df[df$x %in% 1,3] - df[df$x %in% 2,3]
# 
# (ax[1] - ax[2]) + kt * (bx[1] - bx[2])


# 
