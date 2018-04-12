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


nothing_but_net = function(x, y, lambda, alpha) {
  
  xstandard = scale(x)
  ystandard = scale(y)
  lamb1 = lambda * alpha
  lamb2 = .5*(1-alpha)*lambda
  xstar = rbind(xstandard, lamb2 * diag(nrow = ncol(x)))
  ystar = rbind(ystandard, as.matrix(rep(0, ncol(x))))
  betas = lasso.cd(xstar, ystar, beta = rep(0, ncol(x)), lambda = lambda)
  netbeta = ((1 + lamb2)^.5) * betas$beta
  return(netbeta)
  
}

test0 = nothing_but_net(x, y, 10, .95)
test1 = glmnet(as.matrix(x), as.matrix(y), lambda = 10, alpha = .95)
re = data.frame(me = test0, theirs = as.vector(test1$beta))
View(re)

### Import data

x = read_csv("X.csv", col_names = FALSE)

y = read_csv("y.csv", col_names = FALSE)

