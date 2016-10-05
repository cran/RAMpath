\name{powerLCS}
\alias{powerLCS}
\title{
Power analysis for univariate latent change score models
}
\description{
Calculate power for univariate latent change score models based on Monte Carlo simulation.
}
\usage{
powerLCS(N=100, T=5, R=1000, betay=0, my0=0, mys=0, 
varey=1, vary0=1, varys=1, vary0ys=0, alpha=0.05, ...)
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
  \item{alpha}{Significance level}
  \item{\dots}{Options can be used for lavaan}
}
\value{
\item{model}{The lavaan model specification of the bivariate latent change score model}
\item{lavaan}{The lavaan output}
\item{ram}{Output in terms of RAM matrices}
}
\references{
Zhang, Z., & Liu, H. (2016). Sample Size Planning for Latent Change Score Models through Monte Carlo Simulation.
}

\examples{
\dontrun{
powerLCS(R=1000)
}
}