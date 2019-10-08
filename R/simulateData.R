#' Simulate Data For CAR Function
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme
#' @param mean.s Mean of the short term response
#' @param mean.t Mean long-term response
#' @param sigma0 sigma0 in the bivariate normal distribution
#' @param sigma sigma in bivariate normal distribution
#' @param rho rho in bivariate normal distribution
#' @param tau1 Tau1 is the weight on the first covariate
#' @param tau2 Tau2 is the weight on the second covariate
#' @param treat List containing the treatment assignments and the base matrix
#' @param covValues covValues is the dummy variable matrix from the covariate values
#' @param data Is the data from the previous analysis if applicable. The default is NULL
#' @export
simulatedata.car  <- function(mean.x,mean.y, sigma0, sigma, rho, tau1, tau2, beta0, beta1,treat,covValues)
  # simulate multiple treatment groups
{

  trts = unique(treat)
  trts = trts[order(trts)]
  # if (is.null(data)) keeps = rep(TRUE, nrow(covValues)) else keeps = data$treat %in% trts
  # if (is.null(data)) covValues = covValues else covValues = covValues[(length(data$treat)+1):nrow(covValues),]
  n.trt = length(unique(treat))-1
  n.x = table(treat)
  n.y = ceiling(n.x/2)
  diff = n.x-n.y
  mean.full.x = rep(rep(0,n.trt+1),n.x)
  mean.full.y = rep(rep(0,n.trt+1),n.y)
  treat = tail(treat,sum(n.x))
  for (trt in trts)
  {
    mean.full.y[treat==trt][1:n.x[which(trts==trt)]] = c(rep(mean.y[which(trts==trt)], n.y[which(trts==trt)]),rep(NA,diff[which(trts==trt)]))
    mean.full.x[treat==trt] = rep(mean.x[which(trts==trt)], n.x[which(trts==trt)])
  }

  mean.full = cbind(mean.full.x,mean.full.y)

  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  # covmatrix = matrix(c(tau1*covValues[,1],tau2*covValues[,2]),ncol=2)
  full.data  = full.error + mean.full
  new.data = data.frame(x=full.data[,1],y=full.data[,2],treat=treat)
}

