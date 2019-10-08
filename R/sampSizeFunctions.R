#' Sample Size Function
#'
#' This function is used to calculate the sample size for a given z.alpha and z.beta
#' @param z.alpha z alpha is the set Type I error level
#' @param z.beta Z beta is the set power level
#' @param rho Is the allocation ratio of available long-term endpoints to short-term endpoints
#' @param beta1 Beta1 is the beta from the relationship between the X and Y variable
#' @param beta0 Beta 0 is the beta from the relationship between the X and Y variable
#' @param sigma Sigma is either one sigma or two sigmas if two treatments
#' @param tau Tau is used in the calculation of r
#' @param mu.delta Mu.delta is the set level between the treatment and control, or the two treatments
#' @param n.trt N.trt is the number of treatments.
#' @export
get.samplesize = function(z.alpha, z.beta, rho, beta1, beta0, sigma, tau, mu.delta, n.trt, gam = NULL) {

  # m = rho*n
  # N = (1+rho)*n
  #Sample size calculations
  # var.mu.gd = (1/(n/(beta.1^2*tau^2) + (m/sigma^2)))*(1+2*omega*(1-omega)*((1/(n-1))+(1/(m-1))) + O*((1/(n^2))+(1/(m^2))))
  if (n.trt == 1) N = sampleSizeOne(z.alpha, z.beta, rho, beta1, beta0, sigma, tau, mu.delta) else N = sampleSizetwo(z.alpha,z.beta, rho, beta1, beta0, sigma, tau, mu.delta, gam)

  N
}
#' samplesizeone calculates the sample size for a one treatment case
#'
#' This function is used to calculate the sample size for a one treatment case
#' @param z.alpha z alpha is the set Type I error level
#' @param z.beta Z beta is the set power level
#' @param rho Is the allocation ratio of available long-term endpoints to short-term endpoints
#' @param beta1 Beta1 is the beta from the relationship between the X and Y variable
#' @param beta0 Beta 0 is the beta from the relationship between the X and Y variable
#' @param sigma Sigma is either one sigma or two sigmas if two treatments
#' @param tau Tau is used in the calculation of r
#' @param mu.delta Mu.delta is the set level between the treatment and control, or the two treatments
#' @param gam Gamma is the ratio between the two treatment groups and is not used in this case
sampleSizeOne = function(z.alpha, z.beta, rho, beta1, beta0, sigma, tau, mu.delta, gam = NULL) {
  M.fixed = ((qnorm(z.alpha/2)+qnorm(z.beta))^2*sigma^2)/(mu.delta)^2
  r = beta1^2*tau^2*(1/sigma^2)
  n = ((M.fixed/(2*(rho+r^-1)))*(1+sqrt(1+(8*(rho+1))/((1+rho*r)*M.fixed))))
  N = (1+rho)*ceiling(n)
  N
}

#' samplesizetwo calculates the sample size for a two treatment case
#'
#' This function is used to calculate the sample size for a two treatment case
#' @param z.alpha z alpha is the set Type I error level
#' @param z.beta Z beta is the set power level
#' @param rho Is the allocation ratio of available long-term endpoints to short-term endpoints
#' @param beta1 Beta1 is the beta from the relationship between the X and Y variable
#' @param beta0 Beta 0 is the beta from the relationship between the X and Y variable
#' @param sigma Sigma is either one sigma or two sigmas if two treatments
#' @param tau Tau is used in the calculation of r
#' @param mu.delta Mu.delta is the set level between the treatment and control, or the two treatments
#' @param gam Gam is the ratio between the two treatment groups
sampleSizetwo = function(z.alpha, z.beta, rho, beta1, beta0, sigma, tau, mu.delta, gam) {

  sigma1 = sigma[1]
  sigma2 = sigma[2]

  r = beta1^2*tau^2/sigma^2
  A = ((qnorm(z.alpha/2)+qnorm(z.beta))^2)/(mu.delta^2)
  B = sigma1^2/(rho+r[1]^-1) + sigma2^2/(gam*(rho+r[2]^-1))
  C = B^-2*((sigma1^2/(r[1]*(rho+r[1]^-1)^3))+sigma2^2/(gam^2*r[2]*(rho+r[2]^-1)^3))
  n1 = .5*A*B*(1+sqrt(1+ 8*(1+rho)*(A^-1)*C))

  N = (1+rho)*(1+gam)*n1
  N
}
