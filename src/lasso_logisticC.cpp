#include "RcppArmadillo.h"

arma::vec predC(arma::mat X, arma::vec w);
arma::colvec softT(arma::colvec a, double kappa);
double obj1(arma::mat C, arma::colvec w, double rho, arma::colvec z, arma::colvec u);

//' Compute loglikelihood for a logistic regression model
//' 
//' @param X numeric matrix of observed predictors (not including a column of 1s for the intercept term)
//' @param y numeric vector of responses (0 or 1)
//' @param betas numeric vector of coefficients (including intercept term in first position)
// [[Rcpp::export]]
double logLikeC(arma::mat X, arma::vec y, arma::vec betas){
  int p = (int) X.n_cols;
  int n = (int) X.n_rows;
  double b0 = betas(0);
  arma::vec beta = betas.subvec(1,p);
  double ret = sum((-y % log(1 + exp(-1*(b0 + X * beta)))) - ((1-y) % log(1 + exp(b0 + X * beta)) ));
  return(ret);
}

//' Compute penalized loglikelihood for a logistic regression model
//' 
//' @param X numeric matrix of observed predictors (not including a column of 1s for the intercept term)
//' @param y numeric vector of responses (0 or 1)
//' @param betas numeric vector of coefficients (including intercept term in first position)
//' @param lambda penalty tuning parameter
//' @return penalized loglikelihood value
// [[Rcpp::export]]
double PenLogLikeC(arma::mat X, arma::vec y, arma::vec betas, double lambda = 0.01){
  int p = (int) X.n_cols;
  int n = (int) X.n_rows;
  double b0 = betas(0);
  arma::vec beta = betas.subvec(1,p);
  double ret = sum((-y % log(1 + exp(-1*(b0 + X * beta)))) - ((1-y) % log(1 + exp(b0 + X * beta)) ));
  ret -= lambda*sum(abs(beta));
  return(ret);
}

// [[Rcpp::export]]
double costC(arma::mat X, arma::vec y, arma::vec betas){
  int p = (int) X.n_cols;
  int n = (int) X.n_rows;
  arma::vec pred = predC(X, betas);
  double ret = sum((-y % log(pred)) - ((1 - y) % log(1 - pred))) / n; 
  return ret;
}

// [[Rcpp::export]]
arma::vec predC(arma::mat X, arma::vec w){
  return 1 / (1 + exp(-(X * w)));
}

// [[Rcpp::export]]
double obj1(arma::mat C, arma::colvec w, double rho, arma::colvec z, arma::colvec u){
  return(sum( arma::log(1 + arma::exp(C * w)) ) + (rho / 2) * sum(arma::pow((w - z + u),2) ));
}

//' Update beta estimates using Newton-Raphson algorithm
//' 
//' The beta-update step requires optimizing a convex function. This version of the update function
//' uses a Newton-Raphson approach to minimizing the objective function.
//' 
//' @param X Covariate matrix (no column for intercept)
//' @param y Vector of observations (coded in -1/1)
//' @param u Current value of u vector (ADMM optimization)
//' @param z Current value of z vector (ADMM optimization)
//' @param rho Tuning parameter for ADMM optimization
//' @return Vector containing updated estimate of beta vector
// [[Rcpp::export]]
arma::colvec b_updateC(arma::mat X, arma::colvec y, arma::colvec u, arma::colvec z,
                       double rho, unsigned int maxiter=50, double toler=1e-5, 
                       double b = 0.5, double alpha = 0.1){
  int m = X.n_cols;
  int n = X.n_rows;
  double fx, dfx, t;
  arma::colvec g;
  arma::colvec dx;
  arma::mat H;
  
  arma::colvec x = arma::vec(m+1, arma::fill::zeros);
  arma::mat I = arma::mat(m+1, m+1, arma::fill::eye);
  arma::mat C = arma::join_rows(-y, -X);
  
  for(int i = 0;  i < maxiter; ++i){
    fx = obj1(C, x, rho, z, u);

    g = C.t() * (arma::exp(C * x) / (1 + arma::exp(C * x))) + rho * (x - z + u);

    H = C.t() * arma::diagmat(arma::exp(C * x) / arma::pow(1 + arma::exp(C * x), 2)) * C + rho*I;
    dx = -H.i() * g;
    dfx = arma::dot(g, dx);

    if(abs(dfx) < toler) {break;};

    t = 1;
    while (obj1(C, x + t*dx, rho, z, u) > fx + alpha*t*dfx){
      t = b * t;
    }
    x += t * dx;
  }
  return(x);
}

//' Fit l1-penalized logistic regression model using ADMM
//' 
//' Use an ADMM approach to find the parameters for a l1-penalized logistic regression model.
//' Finds solution to argmin_beta sum(log(1+-yX beta)) + lambda sum(|beta|)
//' 
//' @param X Covariate matrix (no column for intercept)
//' @param y Vector of observations (coded in -1/1)
//' @param u Current value of u vector (ADMM optimization)
//' @param z Current value of z vector (ADMM optimization)
//' @param rho Tuning parameter for ADMM optimization
//' @return Vector containing updated estimate of beta vector
// [[Rcpp::export]]
arma::colvec admmlasso_logC(arma::mat X, arma::colvec y,double lam,
                           double rho=1e-3, unsigned int maxit=1000,
                           double tol=1e-3){
  int n = X.n_rows;
  int m = X.n_cols;
  
  // Parameters
  arma::colvec z = arma::vec(m+1, arma::fill::zeros); 
  arma::colvec zold = arma::vec(m+1, arma::fill::zeros); 
  arma::colvec w = arma::vec(m+1, arma::fill::zeros); 
  arma::colvec betas = arma::vec(m+1, arma::fill::zeros); 
  arma::colvec beta = arma::vec(m, arma::fill::zeros);
  double b0;
    
  double objVal = -999999;
  double old_objVal;
    
  int iterations = 0;
  for(int i =0;  i < maxit; ++i){
    iterations++;
    
    betas=b_updateC(X,y,w,z,rho);
    
    b0 = betas(0);
    beta = betas.subvec(1,m);
    
    zold = z;
    z = betas + w;
    
    z.subvec(1,m) = softT(z.subvec(1, m), n*lam/rho);
    w = w + (betas - z);
    
    // Convergence check
    old_objVal = objVal;
    objVal = sum(arma::log(1 + arma::exp(-y % (X * beta) - y * b0)))+ n * lam *sum(arma::abs(z));
    if(abs(objVal-old_objVal) < tol){break;}
      
  }
  return(z);
}

arma::colvec softT(arma::colvec a, double kappa){
  arma::colvec comp = arma::vec(a.n_elem, arma::fill::zeros);
  return arma::max(comp, a - kappa) - arma::max(comp, -a-kappa);
}