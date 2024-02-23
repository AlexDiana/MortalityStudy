library(nimble)

# FULL OCCUPANCY ----

nimbleCodeCMI <- nimbleCode({
  
  for(i in 1:n){
    
    # logit(q[i]) <- a + b * k;
    # logit(q[i]) <- gamma1 * a[age[i]] + gamma2 * b[age[i]] * k[year[i]];
    # d[i] ~ dbinom(prob = q[i], size = E[i])
    d[i] ~ dbinom(prob = q, size = E[i])
    
  }
  
  # q ~ dbeta(1,1)
  a ~ dnorm(0,1)
  b ~ dnorm(0,1)
  k ~ dnorm(0,1)
  # for(i in 1:A){
  #   a[i] ~ dnorm(0, 1)
  #   b[i] ~ dnorm(0, 1)
  # }
  # 
  # for(i in 1:Y){
  #   k[i] ~ dnorm(0, 1)
  # }
  
  # gamma1 ~ dbern(p1)
  # gamma2 ~ dbern(p2)
  
})

# TRUE DATA -----

load(here("data_summarised.rda"))

n <- nrow(data_summarised)

E <- data_summarised$Exposure
d <- data_summarised$Claim

# refactor years
{
  year <- data_summarised$Year
  year <- factor(year)
  originalYears <- levels(year)
  year <- as.numeric(year)
  Y <- length(originalYears)
}


# refactor ages
{
  age <- data_summarised$Age
  age <- factor(age)
  originalAges <- levels(age)
  age <- as.numeric(age)
  A <- length(originalAges)
}

# refactor products

{
  product <- data_summarised$Product
  product <- factor(product)
  originalProducts <- levels(product)
  product <- as.numeric(product)
  P <- length(originalProducts)
}

cmi_dat <- list(n = n,
                A = A,
                Y = Y,
                P = P,
                age = age,
                year = year,
                product = product,
                E = E,
                d = d)

# SIMULATED DATA ----------

n <- 1000
A <- 5
Y <- 10

E <- rpois(n, lambda = 1000)

a <- rnorm(A)
b <- rnorm(A)
k <- rnorm(Y)

b <- b - mean(b) + 1 / A
k <- k - mean(k)

year <- sample(1:Y, n, replace = T)
age <- sample(1:A, n, replace = T)

logistic <- function(x){
  1 / (1 + exp(-x))
}

q <- logistic(a[age] + b[age] * k[year])

d <- sapply(1:n, function(i){
  rbinom(1, E[i], q[i])
})

# FITTING -------

n <- 1000
A <- 5
Y <- 10

year <- sample(1:Y, n, replace = T)
age <- sample(1:A, n, replace = T)

E <- rpois(n, lambda = 1000)

constants <- list(
  n = n,
  A = A,
  Y = Y,
  age = age,
  year = year,
  E = E
)

cmiModel <- nimbleModel(nimbleCodeCMI, constants)

# a <- rnorm(A)
# b <- rnorm(A)
# k <- rnorm(Y)
# E <- rpois(n, lambda = 1000)

a <- 1
b <- 0
k <- 0

q <- .5

initVars <- list(
  q = q
  # a = a,
                 # b = b,
                 # k = k
                 # gamma1 = 1,
                 # gamma2 = 1
                 )

cmiModel$setInits(initVars)

{
  
  logistic <- function(x){
    1 / (1 + exp(-x))
  }
  
  q <- logistic(a[age] + b[age] * k[year])
  
  d <- sapply(1:n, function(i){
    # rbinom(1, E[i], q[i])
    rbinom(1, E[i], .2)
  })
  
}

cmiModel$calculate()

# nodesToSimulate <- cmiModel$getDependencies(c("a","b",
                                              # "k","gamma1","gamma2","q"), self = F)
# cmiModel$simulate(nodesToSimulate)
# newData <- cmiModel$d

cmiModel$setData(d = d)

cmiModel$getLogProb()

cmiModel$d
cmiModel$E
cmiModel$q

mcmcConf <- configureMCMC(cmiModel)

mcmcConf$printSamplers()

Rmcmc <- buildMCMC(mcmcConf)
RmodelComp <- compileNimble(cmiModel)

Cmcmc <- compileNimble(Rmcmc, project = RmodelComp)

output <- runMCMC(Cmcmc)

str(output)

output[,3]
