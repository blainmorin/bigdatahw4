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


soft_thresholding <- function(x,a){
  ## This could be done more efficiently using vector multiplication
  ## See the forumula in slides
  ##  sign(x)*pmax(abs(x) - a, 0)
  result <- numeric(length(x))
  result[which(x > a)] <- x[which(x > a)] - a
  result[which(x < -a)] <- x[which(x < -a)] + a
  return(result)
}



lasso_kkt_check <- function(X,y,beta,lambda, tol=1e-3){
  ## check convergence 
  beta <- as.matrix(beta); X <- as.matrix(X)
  ## Assuming no intercepts 
  G <- t(X)%*%(y-X%*%beta)/length(y)
  ix <- which(beta == 0 )
  iy <- which(beta != 0)
  if (any(abs(G[ix]) > (lambda + tol) )) { return(pass=0) }
  if (any(abs( G[iy] - lambda*sign(beta[iy] ) ) > tol)) { return(pass=0) }  
  return(pass=1)
}



lasso.cd <- function(X,y,beta,lambda,tol=1e-6,maxiter=1000,quiet=FALSE){
  # note that the LS part  in this function is the one in slides divided by length(y) = n 
  ## Or equivalently  lambda here = n * lambda in class
  beta <- as.matrix(beta); X <- as.matrix(X)
  obj <- numeric(length=(maxiter+1))
  betalist <- list(length=(maxiter+1))
  betalist[[1]] <- beta
  
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r <- y - X[,-k]%*%beta[-k]
      beta[k] <- (1/norm(as.matrix(X[,k]),"F")^2)*soft_thresholding(t(r)%*%X[,k],length(y)*lambda)
    }
    betalist[[(j+1)]] <- beta
    obj[j] <- (1/2)*(1/length(y))*norm(y - X%*%beta,"F")^2 + lambda*sum(abs(beta))
    if (norm(betalist[[j]] - beta,"F") < tol) { break }
  } 
  check <- lasso_kkt_check(X,y,beta,lambda) 
  
  if (quiet==FALSE){
    if (check==1) {
      cat(noquote("Minimum obtained.\n"))
    }
    else { cat(noquote("Minimum not obtained.\n")) } 
  }
  return(list(obj=obj[1:j],beta=beta)) 
}


## Using the solution from Zou and Hastie Paper
## Put lambda1 and lambda2 in terms of the lamba and alpha formula from class
## Basically, the function transforms the x and y matrix and then solves them like lasso


## This one centers y and standardizes x
nothing_but_net = function(x, y, lambda, alpha) {
  
  xstandard = scale(x)
  ystandard = scale(y, scale = FALSE)
  lamb1 = lambda * alpha 
  lamb2 = .5*(1-alpha)*lambda 
  xstar = ((1 +lamb2)**-.5) * rbind(xstandard, (lamb2^.5) * diag(nrow = ncol(x)))
  ystar = rbind(ystandard, as.matrix(rep(0, ncol(x))))
  betastar = lasso.cd(xstar, ystar, lambda = lambda, beta = rep(0, ncol(x)))
  netbeta = ((1 + lamb2)^.5) * betastar$beta
  return(netbeta)
  
}


## This one does not scale x or y
nothing_but_net2 = function(x, y, lambda, alpha) {
  
  xstandard = as.matrix(x)
  ystandard = as.matrix(y)
  lamb1 = lambda * alpha 
  lamb2 = .5*(1-alpha)*lambda 
  xstar = ((1 +lamb2)**-.5) * rbind(xstandard, (lamb2^.5) * diag(nrow = ncol(x)))
  ystar = rbind(ystandard, as.matrix(rep(0, ncol(x))))
  betastar = lasso.cd(xstar, ystar, lambda = lambda, beta = rep(0, ncol(x)))
  netbeta = ((1 + lamb2)^.5) * betastar$beta
  return(netbeta)
  
}

## Run net with lambda = .01, .1, 1, and 10

b1 = NULL

lambdas = c(.01, .1, 1, 10)

for(i in 1:length(lambdas)){
  
  b1 = cbind(b1, nothing_but_net2(x, y, lambdas[i], .95))
  
}

# glmnet check
#for(i in 1:length(lambdas)){
#  fit = glmnet(as.matrix(x), as.matrix(y), lambda = lambdas[i], alpha = 1, standardize = FALSE, intercept = FALSE)
#  b1 = cbind(b1, as.numeric( coef(fit) )[-1])
#  
#}

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
