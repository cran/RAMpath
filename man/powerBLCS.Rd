\name{powerBLCS}
\alias{powerBLCS}
\title{
Power analysis for bivariate latent change score models
}
\description{
Calculate power for bivariate latent change score models based on Monte Carlo simulation.
}
\usage{
powerBLCS(N=100, T=5, R=1000, betay=0, my0=0, mys=0, varey=1, 
vary0=1, varys=1, vary0ys=0, alpha=0.05, betax=0, mx0=0, 
mxs=0, varex=1, varx0=1, varxs=1, varx0xs=0, varx0y0=0,  
varx0ys=0, vary0xs=0, varxsys=0, gammax=0, gammay=0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Sample size, can be a scalar or a vector. For better performance, make sure N is at least two times of T}
  \item{T}{Number of times, occasions or waves of measurements, can be a scalar or a vector}
  \item{R}{Number of replications to run in Monte Carlo simulation. Recommended 1000 or more}
  \item{betay}{Population parameter values}
  \item{my0}{Population parameter values}
  \item{mys}{Population parameter values}
  \item{varey}{Population parameter values}
  \item{vary0}{Population parameter values}
  \item{varys}{Population parameter values}
  \item{vary0ys}{Population parameter values}
  \item{betax}{Population parameter values}
  \item{mx0}{Population parameter values}
  \item{mxs}{Population parameter values}
  \item{varex}{Population parameter values}
  \item{varx0}{Population parameter values}
  \item{varxs}{Population parameter values}
  \item{varx0xs}{Population parameter values}
  \item{gammax}{Population parameter values}
  \item{gammay}{Population parameter values}
  \item{varx0y0}{Population parameter values}
  \item{varx0ys}{Population parameter values}
  \item{vary0xs}{Population parameter values}
  \item{varxsys}{Population parameter values}
  \item{alpha}{Significance level}
  \item{\dots}{Options can be used for lavaan}
}

\value{
A matrix with power for each parameter.
}
\references{
Zhang, Z., & Liu, H. (2016). Sample Size Planning for Latent Change Score Models through Monte Carlo Simulation.
}

\examples{
\dontrun{
powerBLCS(R=1000)
}
}