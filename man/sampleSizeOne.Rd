% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampSizeFunctions.R
\name{sampleSizeOne}
\alias{sampleSizeOne}
\title{samplesizeone calculates the sample size for a one treatment case}
\usage{
sampleSizeOne(z.alpha, z.beta, rho, beta1, beta0, sigma, tau, mu.delta,
  gam = NULL)
}
\arguments{
\item{z.alpha}{z alpha is the set Type I error level}

\item{z.beta}{Z beta is the set power level}

\item{rho}{Is the allocation ratio of available long-term endpoints to short-term endpoints}

\item{beta1}{Beta1 is the beta from the relationship between the X and Y variable}

\item{beta0}{Beta 0 is the beta from the relationship between the X and Y variable}

\item{sigma}{Sigma is either one sigma or two sigmas if two treatments}

\item{tau}{Tau is used in the calculation of r}

\item{mu.delta}{Mu.delta is the set level between the treatment and control, or the two treatments}

\item{gam}{Gamma is the ratio between the two treatment groups and is not used in this case}
}
\description{
This function is used to calculate the sample size for a one treatment case
}
