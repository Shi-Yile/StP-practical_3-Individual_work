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
  
  y_fitted <- beta %*% t(x)   ## the vector of 'fitted values'
  mu <- mean(y_fitted)    ## expected value of 'fitted values'
  names(mu) <- 'fitted'   ## name to identify
  
  resid <- y - y_fitted     ## residuals
  var_resid <- as.numeric(resid %*% t(resid) / (n - p))     ## the variance of residuals
  inv_R <- backsolve(qr.R(qrx), diag(p))    ## the inverse of R
  V <- (inv_R %*% t(inv_R)) * var_resid    ## the estimated covariance matrix of least square estimators
  colnames(V) <- rownames(V) <- colnames(x)   ## name the rows and columns of estimated covariance matrix
  
  sigma <- c(sqrt(var_resid))   ## the estimated standard deviation of response and residuals
  names(sigma) <- 'residuals'     ## name sigma to identify response and residuals
  
  # flev <- 
  
  re <- list(beta = beta, V = V, mu = mu, y = y, yname = yname, fitted_values = drop(y_fitted), residuals = drop(resid), formula = formula, sigma = sigma) 
  class(re) <- 'linmod'
  return(re)
}

## Print Function
print.linmod <- function(x){
  
  est_bind <- cbind(x$beta, sqrt(diag(x$V)))    ## calculate standard deviations and bind with parameters estimates 
  colnames(est_bind) <- c('Estimate', 's.e.')   ## name the columns of bound matrix    
  print(x$formula)    ## print formula 
  cat("\n")     ## print blank line
  print(est_bind)     ## print bound matrix
}

## Plot Function
plot.linmod <- function(x){
  
  plot(x = x$fitted_values, y = x$residuals, xlab = 'fitted values', ylab = 'residuals', xlim = )    ## plot model residuals against fitted values
  abline(h = 0, col = 'blue', lty = 'dashed')  ## dashed horizontal line
}

## Predict Function
predict.linmod <- function(x, newdata){
   
  x_new <- model.matrix(~ . , newdata)    ## model matrix from newdata
  y_predict <- as.vector(x$beta %*% t(x_new))   ## predictions of response variable
  return(y_predict)
}