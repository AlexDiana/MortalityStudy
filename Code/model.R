
library(here)
library(ggplot2)
library("rstan") 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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
P <- 2

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

cmi_dat <- list(n = n,
                A = A,
                A2 = A,
                Y = Y,
                P = P,
                age = age,
                year = year,
                E = E,
                d = d)

# MODEL -----

model <- stan_model(file = here('Code','code.stan'), verbose = T)

# mcmc
{
  mcmc_fit <- sampling(model, data = cmi_dat, 
                       chains = 1, iter = 10000,
                       pars = c("a","b","k","p"))
  
  mcmc_matrix_of_draws <- as.matrix(mcmc_fit)
}

# vb

{
  niter <- 2500
  
  vb_fit <- vb(model, data = cmi_dat, 
               algorithm = "meanfield", 
               pars = c("a","b","k","p"),
               include = T,
               iter = 5000,
               elbo_samples = 500,
               tol_rel_obj = 0.00001,
               output_samples = niter)
  
  vb_matrix_of_draws <- as.matrix(vb_fit)
  
}

# OUTPUT ------

param_samples <- vb_matrix_of_draws[,-ncol(vb_matrix_of_draws)]

colnames(param_samples)

param_qtl <- apply(param_samples, 2, function(x){
  quantile(x, probs = c(0.025, 0.975))
})

param_mean <- apply(param_samples, 2, mean)


subsetParams <- 1:Y
# trueParams <- k

ggplot(data = NULL, aes(x = 1:length(subsetParams),
                        ymin = param_qtl[1,subsetParams],
                        ymax = param_qtl[2,subsetParams]#,
                        # y = trueParams
                        )
) + 
  geom_errorbar(size = .7) + theme_bw() + 
  # geom_point(color = "red") +
  # ylab("Occupancy probability") + 
  theme(axis.text = element_text(angle = 90, size = 12),
        axis.title = element_text(size = 15))

# create lee carter trajectory

K <- 1

a_mean <- param_mean[1:A]
b_mean <- param_mean[A + 1:A]
k_mean <- param_mean[1]
d_mean <- param_mean[2 * A + Y + 1:P]

qx_1 <- logistic(a_mean + b_mean * k_mean + d_mean[1])
qx_2 <- logistic(a_mean + b_mean * k_mean + d_mean[2])
qx_3 <- logistic(a_mean + b_mean * k_mean + d_mean[3])
qx_4 <- logistic(a_mean + b_mean * k_mean + d_mean[4])

originalProducts

ggplot() + 
  geom_line(data = NULL, aes(x = 1:A, 
                                      y = qx_1)) + 
  geom_line(data = NULL, aes(x = 1:A, 
                                      y = qx_2), color = "red") + 
  geom_line(data = NULL, aes(x = 1:A, 
                                      y = qx_3), color = "blue") + 
  geom_line(data = NULL, aes(x = 1:A, 
                                      y = qx_4), color = "green")  + 
  coord_cartesian(xlim = c(60,80))


# PRODUCTS --

originalProducts

subsetParams <- 2 * A + Y + 2:P

ggplot(data = NULL, aes(x = 1:length(subsetParams),
                        ymin = param_qtl[1,subsetParams],
                        ymax = param_qtl[2,subsetParams]#,
                        # y = trueParams
)
) + 
  geom_errorbar(size = .7) + theme_bw() + 
  # geom_point(color = "red") +
  ylab("Occupancy probability") + 
  theme(axis.text = element_text(angle = 90, size = 12),
        axis.title = element_text(size = 15))
