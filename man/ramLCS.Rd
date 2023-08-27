\name{ramLCS}
\alias{ramLCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate latent change score model
}
\description{
Univariate latent change score model}
\usage{
ramLCS(data, y, timey, ram.out = FALSE, betay, my0, mys, 
varey, vary0, varys, vary0ys, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data
}
  \item{y}{
y data
}
  \item{timey}{
time of y
}
  \item{ram.out}{
Whether print ram matrices
}
  \item{betay}{
Starting value
}
  \item{my0}{
Starting value
}
  \item{mys}{
Starting value
}
  \item{varey}{
Starting value
}
  \item{vary0}{
Starting value
}
  \item{varys}{
Starting value
}
  \item{vary0ys}{
Starting value
}
  \item{\dots}{
Options can be used for lavaan
}
}
\value{
\item{model}{The lavaan model specification of the bivariate latent change score model}
\item{lavaan}{The lavaan output}
\item{ram}{Output in terms of RAM matrices}
}
\references{
Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2015). Using R package RAMpath for tracing SEM path diagrams and conducting complex longitudinal data analysis. Structural Equation Modeling, 22(1), 132-147. https://doi.org/10.1080/10705511.2014.935257
}

\examples{
data(ex3)
test.lcs<-ramLCS(ex3, 7:12)
summary(test.lcs$lavaan, fit=TRUE)

bridge<-ramPathBridge(test.lcs$ram, allbridge=FALSE, indirect=FALSE)
## uncomment to plot
## plot(bridge, 'lcs')
}