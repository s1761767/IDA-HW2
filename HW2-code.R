# University of Edinburgh
# Incomplete Data Analysis
# Assignment 2 (R code)

# Daniela Olejnikova
# s1761767





##### Q2(b) #####


load("dataex2.Rdata")

# Define the log likelihood function as described in Q2(a)
loglik2 <- function(y, mu) {
  x <- y[,1]
  r <- y[,2]
  sum(r*dnorm(x, mean=mu, sd=1.5, log=T) + (1-r)*pnorm(x, mean=mu, sd=1.5, log=T))
}

# Maximize the log likelihood with respect to the mean mu
muopt <- optim(par=10, fn=loglik2, y=dataex2,
               method="Brent", lower=-10, upper=10,
               control=list("fnscale"=-1), hessian=T)

muopt





##### Q4 #####


load("dataex4.Rdata")

# Expected log likelihood function for the E-step
# Arguments:  - x: vector of covariates (fully observed)
#             - y: vector of the observed responses only (no NA's)
#             - b: beta0 and beta1
#             - bhat: beta0(t) and beta1(t) at step t
ell4 <- function(x, y, b, bhat) {
  n <- length(x); m <- length(y)
  xobs <- x[1:m]; xmis <- x[(m+1):n]
  b0 <- b[1]; b1 <- b[2]
  bhat0 <- bhat[1]; bhat1 <- bhat[2]
  
  # E(ymis | yobs, bhat)
  pbhat <- exp(bhat0 + bhat1*xmis)/(1 + exp(bhat0 + bhat1*xmis))
  
  # Corresponding expected log likelihood with respect to Ymis
  sum(y*(b0 + b1*xobs)) + sum(pbhat*(b0 + b1*xmis)) - sum(log(1 + exp(b0 + b1*x)))
}

# M-step function outputting the MLE estimates for beta0 and beta1
# Arguments:  - data: pre-processed source data
#             - obs: number of responses that are observed, i.e. dim(Yobs)
#             - bt0: initial values for beta, i.e. beta0(0), beta1(0)
#             - eps: convergence criterion parameter
maxlikEM4 <- function(data, obs, bt0, eps) {
  x <- data[,1]; y <- data[1:obs,2]
  bt <- bt0
  diff <- 1
  
  while(diff > eps) {
    # Use optim to update beta at each step t
    bt.new <- optim(par=bt, fn=ell4, x=x, y=y, bhat=bt, control=list("fnscale"=-1))
    
    diff <- sum(abs(bt.new$par - bt))
    bt <- bt.new$par
  }
  return(bt)
}

# Pre-processing the data - reorderring, so that the fully observed responses
# and corresponding covariates come first
Xorg <- dataex4$X
Yorg <- dataex4$Y
firstm <- which(!is.na(Yorg))

n <- length(Xorg)
m <- length(firstm)

X <- rep(0, n)
X[1:m] <- Xorg[firstm]
X[(m+1):n] <- Xorg[-firstm]

Y <- rep(0, n)
Y[1:m] <- Yorg[firstm]
Y[(m+1):n] <- Yorg[-firstm]

dataex4ord <- data.frame(X, Y)


res4 <- maxlikEM4(data=dataex4ord, obs=m, bt0=c(1,1), eps=0.00001); res4


# Since we assume ignorability, we can check our results from the EM algorithm
# by simply maximizing the observed log likelihood

loglikobs4 <- function(data, beta) {
  b0 <- beta[1]; b1 <- beta[2]
  x <- data[,1]; y <- data[,2]
  sum(y*(b0 + b1*x) - log(1 + exp(b0 + b1*x)))
}

# Extract the responses and corresponding covariates that are fully observed
Xobs <- Xorg[firstm]
Yobs <- Yorg[firstm]
data4check <- data.frame(Xobs, Yobs)

res4check <- optim(par=c(1,1), fn=loglikobs4, data=data4check, control=list("fnscale"=-1))
res4check





##### Q5(b) #####


# M-step function outputting the MLE estimates for the parameters of interest
# Arguments:  - y: data vector
#             - theta0: initial values of the parameters
#             - eps: convergence criterion parameter
maxlikmix5 <- function(y, theta0, eps) {
  theta <- theta0
  p <- theta[1]; mu <- theta[2]; sig <- theta[3]; lam <- theta[4]
  
  diff <- 1
  while(diff > eps) {
    theta.old <- theta   # store current theta before updating
    
    # E-step - calculate E(Z | y, theta(t))
    ptil1 <- p*dlnorm(y, meanlog=mu, sdlog=sig)
    ptil2 <- (1-p)*dexp(y, rate=lam)
    ptil <- ptil1/(ptil1 + ptil2)
    
    # M-step
    p <- mean(ptil)
    mu <- sum(log(y)*ptil)/sum(ptil)
    sig <- sqrt(sum(((log(y) - mu)^2)*ptil)/sum(ptil))
    lam <- sum(1 - ptil)/sum(y*(1 - ptil))
    
    theta <- c(p, mu, sig, lam)
    diff <- sum(abs(theta - theta.old))
  }
  return(theta)
}

# Prepare the data
load("dataex5.Rdata")
y <- dataex5
theta.inits <- c(0.1, 1, 0.5, 2)

# Obtain and extract the results
res5 <- maxlikmix5(y=y, theta0=theta.inits, eps=0.00001); res5
p <- res5[1]
mu <- res5[2]
sig <- res5[3]
lam <- res5[4]


# Plot the histogram of the data with the estimated density superimposed
hist(y, freq=F, main="Q5(b): Data Histogram + Estimated Mixture pdf")
curve(p*dlnorm(x, mu, sig) + (1-p)*dexp(x, lam), add=T, from=0, to=120)