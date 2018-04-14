########################################################################
#### Homework 4 ################################
#### By Blain Morin ################
#########################

require(readr) == T || install.packages("readr", repos="https://cloud.r-project.org/")

library(readr)

authors <- function() {
  c("Blain Morin")
}


#######################################################################
#### Question 1 ################################
####################################


### Import data

x = read_csv("X.csv", col_names = FALSE)

y = read_csv("y.csv", col_names = FALSE)


### ADMM function from class

softThresh <- function(x, lambda) {
  sign(x)*pmax(0, abs(x) - lambda)
}


admmLasso <- function(X, y, tau, maxit = 1000, tol=1e-4) {
  XX <- t(X) %*% X
  Xy <- t(X) %*% y
  
  p <- ncol(X)
  lambda <- rep(0, p)
  maxRho <- 5
  rho <- 4
  
  z0 <- z <- beta0 <- beta <- rep(0, p)
  Sinv <- solve(XX + rho*diag(rep(1, p)) )
  
  for (it in 1:maxit) {
    ## update beta
    ## beta <- solve(XX + rho*diag(rep(1, p)) ) %*% (Xy + rho * z - lambda)
    beta <- Sinv %*% (Xy + rho * z - lambda)
    
    ## update z
    z <- softThresh(beta + lambda/rho, tau/rho)
    ## update lambda
    lambda <- lambda + rho* (beta - z ) 
    ## increase rho
    ## rho <- min(maxRho, rho*1.1)
    
    change <- max(  c( base::norm(beta - beta0, "F"),
                       base::norm(z - z0, "F") ) )
    if (change < tol || it > maxit) {
      break
    }
    beta0 <-  beta
    z0 <-  z
    
  }
  z
}


## Using the solution from Zou and Hastie Paper
## Put lambda1 and lambda2 in terms of the lamba and alpha formula from class
## Basically, the function transforms the x and y matrix and then solves them like lasso

nothing_but_net2 = function(x, y, lambda, alpha) {
  
  xstandard = as.matrix(x)
  ystandard = as.matrix(y)
  lamb1 = lambda * alpha 
  lamb2 = .5*(1-alpha)*lambda 
  xstar = ((1 +lamb2)**-.5) * rbind(xstandard, (lamb2^.5) * diag(nrow = ncol(x)))
  ystar = rbind(ystandard, as.matrix(rep(0, ncol(x))))
  betastar = admmLasso(xstar, ystar, tau = nrow(x))
  netbeta = ((1 + lamb2)^.5) * betastar
  return(netbeta)
  
}

## Run net with lambda = .01, .1, 1, and 10

b1 = NULL

lambdas = c(.01, .1, 1, 10)


for(i in 1:length(lambdas)){
  
  b1 = cbind(b1, nothing_but_net2(x, y, lambdas[i], .95))
  
}


## Prepare for export
b1 = as.data.frame(b1)
names(b1) = c(".01", ".1", "1", "10")

## Write csv
write.csv(b1, file = "b1.csv", row.names = FALSE)

###########################################################################
#### Question 2 ##########################################
#############################################

b2 = NULL

for(i in 1:length(lambdas)){
  
  b2 = cbind(b2, nothing_but_net2(x, y, lambdas[i], 1))
  
}

## Prepare for export
b2 = as.data.frame(b2)
names(b2) = c(".01", ".1", "1", "10")

## Write csv
write.csv(b2, file = "b2.csv", row.names = FALSE)
