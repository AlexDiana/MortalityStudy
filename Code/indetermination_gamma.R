X <- 7
Y <- 3

ax <- rnorm(X)
kt <- rnorm(Y)
gtx <- rnorm(X + Y - 1)

term1 <- matrix(ax, X, Y, byrow = F)
term2 <- matrix(kt, X, Y, byrow = T)
term4 <- createTerm4(gtx, X, Y)
term4 <- term4 - mean(term4)

mxt <- term1 + term2 + term4

cx <- rnorm(X)
cx <- cx - mean(cx)
cx <- 1:X - mean(1:X)

dt <- rnorm(Y)
dt <- dt - mean(dt)
dt <- 1:Y - mean(1:Y)

ax_tilde <- ax + cx
kt_tilde <- kt - dt

term4_tilde <- term4
for (x in 1:X) {
  for (t in 1:Y) {
    term4_tilde[x,t] <- 
      term4_tilde[x,t] +
      dt[t] - cx[x]
  }
}

term1_tilde <- matrix(ax_tilde, X, Y, byrow = F)
term2_tilde <- matrix(kt_tilde, X, Y, byrow = T)

mxt_tilde <- term1_tilde + term2_tilde + term4_tilde


gtx_diff <- c((term4_tilde - term4)[X:1,1], (term4_tilde - term4)[1,2:Y])

qplot(1:(X+Y-1), gtx_diff)

#  START

convertTerm4ToGtx <- function(term4){
  c(term4[X:1,1], term4[1,2:Y])
  
}

X <- 7
Y <- 3

ax <- rnorm(X)
kt <- rnorm(Y)
gtx <- rnorm(X + Y - 1)

# renormalise kt and gtx
kt <- kt - mean(kt)
term4 <- createTerm4(gtx, X, Y)
term4 <- term4 - mean(term4)
gtx <- c(term4[X:1,1], term4[1,2:Y])

# find kg

t0 <- mean(1:Y)
x0 <- mean(1:X)

idx_0 <- 1
# idx_0 <- floor(t0) - floor(x0) + X - 1
kg <- - gtx[idx_0] / (idx_0 - X - t0 + x0)

# 

gtx_diff <- rep(NA, X + Y - 1)

for (i in 1:(X+Y-1)) {
  gtx_diff[i] <- kg * ((i - X) - t0 + x0) 
}

term4_diff <- createTerm4(gtx_diff, X, Y)

term4_tilde <- term4_diff + term4

gtx_new <- convertTerm4ToGtx(term4_tilde)
kt_new <- kt - kg * (1:Y - t0)
ax_new <- ax + kg * (1:X - x0)

mxt_transformed <- matrix(ax_new, X, Y, byrow = F) + 
  matrix(kt_new, X, Y, byrow = T) + 
  createTerm4(gtx_new, X, Y)

mxt <- matrix(ax, X, Y, byrow = F) + 
  matrix(kt, X, Y, byrow = T) + 
  createTerm4(gtx, X, Y)

max(abs(mxt - mxt_transformed))
