logisticLik <- function(X,y,betas){
  b0 <- betas[1]
  beta <- matrix(betas[-1],ncol=1)
  sum(-y*log(1 + exp(-(b0 + X%*%beta))) - (1-y)*log(1+exp(b0 + X%*%beta)))
}

penLogLik <- function(X,y,betas,lambda=0.01){
  b0 <- betas[1]
  beta <- matrix(betas[-1],ncol=1)
  sum(-y*log(1 + exp(-(b0 + X%*%beta))) - (1-y)*log(1+exp(b0 + X%*%beta))) - lambda*sum(abs(beta))
}

library(tidyverse)
library(glmnet)
library(Rcpp)
library("RColorBrewer")
library(factoextra)

myocarde = read.table("http://freakonometrics.free.fr/myocarde.csv",head=TRUE, sep=";")
myocarde$PRONO = (myocarde$PRONO=="SURVIE")*1
y = myocarde$PRONO
X = myocarde[,1:7]
for(j in 1:7) X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
X = as.matrix(X)
y <- ifelse(y==1, 1, -1)

sourceCpp("aeffp/src/test.cpp")

# PennegLogLik = function(bbeta,lambda=0){
#   b0=bbeta[1]
#   beta=bbeta[-1]
#   -sum(-y*log(1 + exp(-(b0+X%*%beta))) -
#          (1-y)*log(1 + exp(b0+X%*%beta)))+lambda*sum(abs(beta))
# }

opt_lasso = function(lambda){
  beta_init = lm(PRONO~.,data=myocarde)$coefficients
  logistic_opt = optim(par = beta_init*0, 
                       function(x) -PenLogLikeC(X,y,x,lambda), 
                       hessian=TRUE, 
                       method = "BFGS",
                       control=list(abstol=1e-9))
  logistic_opt$par[-1]
}
v_lambda=c(exp(seq(-4,2,length=61)))
est_lasso=Vectorize(opt_lasso)(v_lambda)
colrs=brewer.pal(7,"Set1")
plot(v_lambda,est_lasso[1,],col=colrs[1],type="l", ylim=c(-5,5))
for(i in 2:7) lines(v_lambda,est_lasso[i,],col=colrs[i],lwd=2)

glm_lasso = glmnet(X, y, alpha=1)
plot(glm_lasso,xvar="lambda",col=colrs,lwd=2)
glmnet(X, y, alpha=1,lambda=exp(-4))$beta

# Orthogonal Covariates
pca <- princomp(X)
pca_X <- get_pca_ind(pca)$coord
glm_lasso <- glmnet(pca_X, y, alpha=1)
plot(glm_lasso, xvar = "lambda", col = colrs)
plot(glm_lasso, col = colrs)

# Interior Point Approach

##  Standard Lasso

soft_thresholding <- function(x, a){
  sign(x)*max(0,abs(x)-a)
}

lasso_coord_desc <- function(X,y,beta,lambda,tol=1e-6,maxiter=1000){
  beta = as.matrix(beta)
  X = as.matrix(X)
  omega = rep(1/length(y),length(y))
  obj = numeric(length=(maxiter+1))
  betalist = list(length(maxiter+1))
  betalist[[1]] = beta
  beta0list = numeric(length(maxiter+1))
  beta0 = sum(y-X%*%beta)/(length(y))
  beta0list[1] = beta0
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r = y - X[,-k]%*%beta[-k] - beta0*rep(1,length(y))
      beta[k] = (1/sum(omega*X[,k]^2))*soft_thresholding(t(omega*r)%*%X[,k],length(y)*lambda)
    }
    beta0 = sum(y-X%*%beta)/(length(y))
    beta0list[j+1] = beta0
    betalist[[j+1]] = beta
    obj[j] = (1/2)*(1/length(y))*norm(omega*(y - X%*%beta - 
                                               beta0*rep(1,length(y))),'F')^2 + lambda*sum(abs(beta))
    if (norm(rbind(beta0list[j],betalist[[j]]) - rbind(beta0,beta),'F') < tol) { break } 
  } 
  return(list(obj=obj[1:j],beta=beta,intercept=beta0)) 
}

# Lasso logistic

soft_thresholding = function(x,a){
  result = numeric(length(x))
  result[which(x > a)]  a)] - a
  result[which(x < -a)] <- x[which(x < -a)] + a
return(result)
}

lasso_coord_desc <- function(X,y,beta,lambda,tol=1e-6,maxiter=1000){
  beta = as.matrix(beta)
  X = as.matrix(X)
  obj = numeric(length=(maxiter+1))
  betalist = list(length(maxiter+1))
  betalist[[1]] = beta
  beta0 = sum(y-X%*%beta)/(length(y))
  p = exp(beta0*rep(1,length(y)) + X%*%beta)/(1+exp(beta0*rep(1,length(y)) + X%*%beta))
  z = beta0*rep(1,length(y)) + X%*%beta + (y-p)/(p*(1-p))
  omega = p*(1-p)/(sum((p*(1-p))))
  beta0list = numeric(length(maxiter+1))
  beta0 = sum(y-X%*%beta)/(length(y))
  beta0list[1] = beta0
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r = z - X[,-k]%*%beta[-k] - beta0*rep(1,length(y))
      beta[k] = (1/sum(omega*X[,k]^2))*soft(t(omega*r)%*%X[,k],length(y)*lambda)
    }
    beta0 = sum(y-X%*%beta)/(length(y))
    beta0list[j+1] = beta0
    betalist[[j+1]] = beta
    obj[j] = (1/2)*(1/length(y))*norm(omega*(z - X%*%beta - 
                                               beta0*rep(1,length(y))),'F')^2 + lambda*sum(abs(beta))
    p = exp(beta0*rep(1,length(y)) + X%*%beta)/(1+exp(beta0*rep(1,length(y)) + X%*%beta))
    z = beta0*rep(1,length(y)) + X%*%beta + (y-p)/(p*(1-p))
    omega = p*(1-p)/(sum((p*(1-p))))
    if (norm(rbind(beta0list[j],betalist[[j]]) - 
             rbind(beta0,beta),'F') < tol) { break } 
  } 
  return(list(obj=obj[1:j],beta=beta,intercept=beta0)) 
}

# ADMM lasso logistic

admm_lasso_logistic <- function(X, y, mu, rho, alpha, 
                                max.iter = 1000, abstol = 1e-4, reltol = 1e-2){
  m <- nrow(X)
  n <- ncol(X)
  
  # init
  x <- rep(0, n+1)
  z <- x
  u <- x
  
  if(0 %in% unique(y)) y <- ifelse(y == 0, -1, 1)
  
  for(i in 1:max.iter){
    print(i)
    # x update
    x = update_x(X, y, x, u, z, rho)
    print(x)
    
    # z update
    z_old <- z
    x_hat <- alpha*x + (1 - alpha)*z_old
    z <- x_hat + u
    z[2:(n+1)] <- shrinkage(z[2:(n+1)], (m*mu)/rho)
    print(z)
    
    # u update
    u <- u + (x_hat - z)
    print(u)
    
    # diagnostics, reporting, termination checks
    
    r_norm  = sqrt(sum((x - z)^2))
    print(r_norm)
    # s_norm  = sqrt(sum((rho*(z - z_old))^2))
    # 
    # eps_pri = sqrt(n)*abstol + reltol*max(sqrt(sum(x^2)), sqrt(sum(z^2)))
    # eps_dual= sqrt(n)*abstol + reltol*sqrt(sum((rho%*%u)^2))
    
    # print(r_norm)
    # print(eps_pri)
    # print(s_norm)
    # print(eps_dual)
    # if (r_norm < eps_pri & s_norm < eps_dual) break
  }
}

objective <- function(A, b, mu, x, z){
  m <- ncol(A)
  sum(log(1 + exp(-A%*%x[2:(m+1)] - b*x[1]))) + m*mu*sqrt(sum((z-1)^2))
}

sub_objective <- function(X, y, b, z, u, rho){
  logistic_loss(X,y,b) + (rho/2)*sum((b - z + u)^2)
}

logistic_loss <- function(X,y,b){
  beta <- b[2:length(b)]
  b0 <- b[1]
  p <- log(1/(1+exp(-(b0+X%*%beta))))
  
  sum((1-y)*p + y * (1 - p))
}

update_x <- function(X, y, b, u, z, rho, x0 = NULL){
# solve the x update
# minimize [ -logistic(x_i) + (rho/2)||x_i - z^k + u^k||^2 ]
# via Newton's method; for a single subsystem only.
  optim(b, function(x) sub_objective(X, y, x, u, z, rho), method = "BFGS")$par
}



shrinkage <- function(a, kappa){
    pmax(0, a-kappa) - pmax(0, -a-kappa)
}


n <- 50
m <- 200

w <- rsparsematrix(n, 1, 0.1)
v <- rnorm(1)

X <- rsparsematrix(m, n, 10/n)
btrue <- sign(X%*%w + v)

b <- sign(X%*%w + v + sqrt(0.1)*rnorm(m))

A <- diag(as.numeric(b)) %*% X
