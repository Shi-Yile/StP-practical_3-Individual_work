## Yile Shi (s2168022)

## Linear Regression Model:
linmod <- function(formula, dat){

  y <- model.frame(formula, dat)[[1]]   ## the vector of response variable
  yname <- all.vars(formula)[1]   ## the name of response variable
  
  xname <- all.vars(formula)[-1]  ## the name(s) of explanatory variable(s)
  x <- model.matrix(~. , dat[xname])    ## model matrix
  
  p <- ncol(x)    ## number of parameters
  n <- nrow(x)    ## number of observations
  qrx <- qr(x)    ## get QR decomposition of model matrix
  beta <- backsolve(qr.R(qrx), qr.qty(qrx, y)[1 : p])   ## R^{-1} Q^T y
  names(beta) <- colnames(x)    ## name beta to identify the model components
  
  y_fitted <- beta %*% t(x)    ## the vector of 'fitted values'
  mu <- mean(y_fitted)    ## expected value of 'fitted values'
  names(mu) <- 'fitted'   ## name to identify
  
  e <- y - y_fitted     ## residuals
  var_e <- as.numeric(e %*% t(e) / (n - p))     ## the variance of residuals
  inv_R <- backsolve(qr.R(qrx), diag(p))    ## the inverse of R
  V <- (inv_R %*% t(inv_R)) * var_e    ## the estimated covariance matrix of least square estimators
  colnames(V) <- rownames(V) <- colnames(x)   ## name the rows and columns of estimated covariance matrix
  
  sigma <- c(sqrt(var_e))   ## the estimated standard deviation of response and residuals
  names(sigma) <- 'residuals'     ## name sigma to identify response and residuals
  
  # flev <- 
  
  re <- list(beta = beta, V = V, mu = mu, y = y, yname = yname, formula = formula, sigma = sigma) 
  class(re) <- 'linmod'
  return(re)
}

## Print Function
print.linmod <- function(x){
  
  est_bind <- cbind(x$beta, sqrt(diag(x$V)))    ## calculate standard deviations and bind with parameters estimates 
  colnames(est_bind) <- c('Estimate', 's.e.')   ## name the columns of binded matrix    
  print(x$formula)    ## print formula 
  cat("\n")     ## print blank line
  print(est_bind)     ## print binded matrix
}

## Plot Function
plot.linmod <- function(x){
  
  
}