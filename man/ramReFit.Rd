\name{ramReFit}
\alias{ramReFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Refit a model with additional paths}
\description{
Generate a vector field plot based on the bivariate lcsm
}
\usage{
ramReFit(object, add, ram.out=FALSE, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Output from any data analysis
}
  \item{add}{
Additional paths to be added, e.g., add='X1~~X2'.
}
\item{ram.out}{Whether to print the RAM matrices}
  \item{\dots}{
Options for plot and arrows function.
}
}
\examples{
data(ex3)
gcm.l<-ramLCM(ex3, 1:6, model='linear', ram.out=TRUE)
## Add correlated errors
ramReFit(gcm.l, add='X1~~X2')
}