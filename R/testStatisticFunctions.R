#' get.t
#'
#' This function gets the test statistic for comparing a treatment against a mean
#' @param data Data is the dataset containing short-term and long-term endpoints
#' @param beta0 Beta0 is the beta_0 coefficient from the known relationship between the Y and X variables
#' @param beta1 Beta1 is the beta_1 coefficient from the
#' @param tau Tau is the known tau from the weight relationship
#' @param sigma Sigma is the known variance
#' @param alpha Alpha is the Type I error rate level
#' @export
get.t = function(beta0, beta1, mu.null  = 0, sigma = NULL, tau = NULL, z.alpha = 0.05, z.beta = 0.2) {

  mu.gd.vec = rep(0,n.trt)
  var.mu.gd.vec = rep(0,n.trt)
  z.alpha = 0.05
  z.beta = .2

  for (i in 1:n.trt) {

    n = length(data[data$treat==i,]$x)
    m = length(na.omit(data[data$treat==i,]$y))

    data$yhat = data$x*beta1 + beta0

    y.hat.bar = mean((data[data$treat==i,]$yhat))
    s.1.sq = var(data[data$treat==i,]$yhat)
    y.bar = mean(na.omit(data[data$treat==i,]$y))
    s.2.sq = var(na.omit(data[data$treat==i,]$y))
    # mu.hat = omega*y.hat.bar + (1-omega)*y.bar

    # omega = (n/(beta.1^2*tau^2))/((n/(beta.1^2*tau^2)) + m/sigma^2)

    if (is.null(tau)) omega.hat = (n/(s.1.sq))/((n/s.1.sq)+(m/(s.2.sq))) else omega.hat = (n/((beta0^2)*(tau^2)))/((n/((beta0^2)*(tau^2)))+(m/(sigma^2)))

    mu.gd = omega.hat*y.hat.bar + (1-omega.hat)*y.bar

    if (is.null(tau)) var.mu.gd = (1/((n/s.1.sq)+(m/s.2.sq)))*(1+4*omega.hat*(1-omega.hat)*((1/(n-1))+(1/(m-1)))) else var.mu.gd = 1/((n/((beta0^2)*(tau^2)))+(m/(sigma^2)))
    mu.gd.vec[i] = mu.gd
    var.mu.gd.vec[i] = var.mu.gd
    }
  if (n.trt == 1) t.ci = get.t.one(mu.gd.vec, var.mu.gd.vec, z.alpha, mu.null) else t.ci = get.t.two(mu.gd.vec, var.mu.gd.vec, z.alpha)
  T.stat = t.ci[1]
  ci = t.ci[2:3]

}

#' Two Treatment Test Statistic
#'
#' This function finds the test statistic for a two-treatment case of Chow & Tse (2007)
#' @param var.mu.gd.vec This is the vector of the variance of weighted means
#' @param mu.gd.vec This is the vector of the means of weighted means
#' @param z.alpha This is the set alpha level.
get.t.two = function(mu.gd.vec, var.mu.gd.vec,z.alpha, mu.null ) {
  T.stat = (mu.gd.vec[2] - mu.gd.vec[1]-mu.null)/(sqrt(sum(var.mu.gd.vec)))
  ci.upper = mu.gd.vec[1]-mu.gd.vec[2] +qnorm(z.alpha/2)*sqrt(sum(var.mu.gd.vec))
  ci.lower = mu.gd.vec[1]-mu.gd.vec[2] -qnorm(z.alpha/2)*sqrt(sum(var.mu.gd.vec))
  c(T.stat, ci.upper, ci.lower)
}

#' One Treatment Test Statistic
#'
#' This function finds the test statistic for a two-treatment case of Chow & Tse (2007)
#' @param var.mu.gd.vec This is the vector of the variance of weighted means
#' @param mu.gd.vec This is the vector of the means of weighted means
#' @param z.alpha This is the set alpha level.
#' @param mu.null Mu.null is the mu under the null distribution
get.t.one = function(mu.gd.vec, var.mu.gd.vec,z.alpha, mu.null ) {
  T.stat = (mu.gd.vec[1]-mu.null)/(sqrt(var.mu.gd.vec))
  ci.upper = mu.gd.vec[1]-mu.null +qnorm(z.alpha/2)*sqrt(var.mu.gd.vec)
  ci.lower = mu.gd.vec[1]-mu.null -qnorm(z.alpha/2)*sqrt(var.mu.gd.vec)
  c(T.stat, ci.upper, ci.lower)

}
