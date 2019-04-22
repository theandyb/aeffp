#' Fit l1-penalized logistic regression model using ADMM
#' 
#' Use an ADMM approach to find the parameters for a l1-penalized logistic regression model.
#' Finds solution to argmin_beta sum(log(1+-yX beta)) + lambda sum(|beta|)
#' 
#' @param X Covariate matrix (no column for intercept)
#' @param y Vector of observations (coded in -1/1)
#' @param u Current value of u vector (ADMM optimization)
#' @param z Current value of z vector (ADMM optimization)
#' @param rho Tuning parameter for ADMM optimization
#' @return Vector containing updated estimate of beta vector
#' @export
admmlasso_log=function(X, y, lam, rho=1e-3, maxit=1000, tol=1e-3){
  n <- nrow(X)
  p <- ncol(X)
  
  # Parameters
  z=rep(0,p+1)
  w=rep(0,p+1)
  betas=rep(0,p+1)

  objVal <- -Inf
  
  iterations <- 0
  for(i in 1:maxit){
    iterations <- iterations + 1
    betas=b_update(X,y,w,z,rho)
    
    b0 <- betas[1]
    beta <- betas[2:(p+1)]
    
    zold <- z
    z <- betas+w
    z[2:length(z)] <- shrinkage(z[2:length(z)], n*lam/rho)
    w <- w + (betas - z)
    
    # Convergence check
    old_objVal <- objVal
    objVal <- sum(log(1+exp(-y*X%*%beta - y*b0)))+n*lam*sum(abs(z))
    if(abs(objVal-old_objVal) < tol){break}
    
  }
  return(list("beta"=z,"total.iterations"=iterations))
}

#' Update beta estimates using Newton-Raphson algorithm
#' 
#' The beta-update step requires optimizing a convex function. This version of the update function
#' uses a Newton-Raphson approach to minimizing the objective function.
#' 
#' @param X Covariate matrix (no column for intercept)
#' @param y Vector of observations (coded in -1/1)
#' @param u Current value of u vector (ADMM optimization)
#' @param z Current value of z vector (ADMM optimization)
#' @param rho Tuning parameter for ADMM optimization
#' @return Vector containing updated estimate of beta vector
b_update <- function(X, y, u, z, rho, maxiter=50, toler=1e-5, b = 0.5, alpha = 0.1){
  m <- nrow(X)
  n <- ncol(X)
  
  x = rep(0,n+1)
  
  I <- diag(1,n+1)
  C <- cbind(-y, -X)
  
  f <- function(w) {sum(log(1 + exp(C %*% w))) + (rho/2) * sum((w - z + u)^2)}
  
  for(i in 1:maxiter){
    fx <- f(x)
    
    g <- t(C) %*% (exp(C %*% x) / (1 + exp(C%*%x))) + rho*(x - z + u)
    
    H = t(C) %*% diag(as.numeric(exp(C%*%x)/(1 + exp(C%*%x))^2)) %*% C + rho*I
    dx = -solve(H) %*% g
    dfx = t(g)%*%dx
    
    if(abs(dfx) < toler) break
    
    t = 1;
    while (f(x + t*dx) > fx + alpha*t*dfx){
      t <- b*t;
    }
    x = x + t*dx
  }
  return(x)
}

#' Update beta estimates using L-BFGS-B
#' 
#' The beta-update step requires optimizing a convex function. This version of the update function
#' uses R's built-it optim function to do this for us
#' 
#' @importFrom stats optim
#' @param X Covariate matrix (no column for intercept)
#' @param y Vector of observations (coded in -1/1)
#' @param u Current value of u vector (ADMM optimization)
#' @param z Current value of z vector (ADMM optimization)
#' @param rho Tuning parameter for ADMM optimization
#' @return Vector containing updated estimate of beta vector
b_update2 <- function(X, y, u, z, rho){
  # function to minimize
  C <- cbind(-y, -X)
  m <- nrow(X)
  n <- ncol(X)
  f <- function(w) (sum(log(1 + exp(C %*% w))) + (rho/2) * sqrt(sum(w - z + u)^2))
  
  optim(rep(0.1,n+1),f,method="L-BFGS-B",lower=rep(-100,n+1),upper=rep(100,n+1))$par
}

#' Shrinkage function
#' 
#' Computes soft-thresholded z-update
#' @param z z vector from ADMM problem
#' @param kappa value to be (soft-thresholdedly) subtracted from each element of z
#' @return result of soft thresholding operation
shrinkage <- function(a, kappa){
  pmax(0, a-kappa) - pmax(0, -a-kappa);
}