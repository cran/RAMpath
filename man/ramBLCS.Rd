\name{ramBLCS}
\alias{ramBLCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conduct bivariate latent change score analysis
}
\description{
Conduct bivariate latent change score analysis}
\usage{
ramBLCS(data, y, x, timey, timex, ram.out = FALSE, betax, 
betay, gammax, gammay, mx0, mxs, my0, mys, varex, varey, 
varx0, vary0, varxs, varys, varx0y0, varx0xs, vary0ys, 
varx0ys, vary0xs, varxsys, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data
}
  \item{y}{
Indices for y variables
}
  \item{x}{
Indices for x variables
}
  \item{timey}{
Time for y variables
}
  \item{timex}{
Time for x variables}
  \item{ram.out}{
whether print ram matrices
}
  \item{betax}{
Starting value
}
  \item{betay}{
Starting value
}
  \item{gammax}{
Starting value
}
  \item{gammay}{
Starting value
}
  \item{mx0}{
Starting value
}
  \item{mxs}{
Starting value
}
  \item{my0}{
Starting value
}
  \item{mys}{
Starting value
}
  \item{varex}{
Starting value
}
  \item{varey}{
Starting value
}
  \item{varx0}{
Starting value
}
  \item{vary0}{
Starting value
}
  \item{varxs}{
Starting value
}
  \item{varys}{
Starting value
}
  \item{varx0y0}{
Starting value
}
  \item{varx0xs}{
Starting value
}
  \item{vary0ys}{
Starting value
}
  \item{varx0ys}{
Starting value
}
  \item{vary0xs}{
Starting value
}
  \item{varxsys}{
Starting value
}
  \item{\dots}{
Options can be used for \code{\link[lavaan]{lavaan}}
}
}

\value{
\item{model}{The lavaan model specification of the bivariate latent change score model}
\item{lavaan}{The lavaan output}
\item{ram}{Output in terms of RAM matrices}
}
\references{
Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2013). Using R Package RAMpath for Tracing SEM Path Diagrams and Conducting Complex Longitudinal Data Analysis. Structural Equation Modeling.
}

\examples{
data(ex3)
## Test the bivariate latent change score model ramBLCS

test.blcs<-ramBLCS(ex3, 7:12, 1:6, ram.out=TRUE)
summary(test.blcs$lavaan, fit=TRUE)

bridge<-ramPathBridge(test.blcs$ram, allbridge=FALSE,indirect=FALSE)
plot(bridge, 'blcs')


## Test the vector field plot
## test.blcs is the output of the ramBLCS function.
ramVF(test.blcs, c(0,80),c(0,80), length=.05, xlab='X', ylab='Y',scale=.5, ninterval=9)
}