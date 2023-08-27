\name{ramVF}
\alias{ramVF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a vector field plot based on the bivariate lcsm
}
\description{
Generate a vector field plot based on the bivariate lcsm
}
\usage{
ramVF(ramout, ylim, xlim, ninterval=10, scale=.1, length=.25, 
scatter=TRUE, n=20, alpha=.95, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ramout}{
Output from the ramBLCS function
}
  \item{ylim}{
Range of y data, for example, c(0,80) from 0 to 80
}
  \item{xlim}{
Range of x data, for example, c(0,80) from 0 to 80
}
  \item{ninterval}{
Number of intervals for plotting. The default is 10.
}
  \item{scale}{
Time interval to calculate vector fields.
}
  \item{length}{
The length of arrows to plot
}
\item{scatter}{Whether to plot the data points}
\item{n}{The number of data points to be plotted}
\item{alpha}{The confidence level to calculate the ellipse}
  \item{\dots}{
Options for plot and arrows function.
}
}

\references{
Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2015). Using R package RAMpath for tracing SEM path diagrams and conducting complex longitudinal data analysis. Structural Equation Modeling, 22(1), 132-147. https://doi.org/10.1080/10705511.2014.935257
}

\examples{
data(ex3)
test.blcs<-ramBLCS(ex3, 1:6, 7:12, ram.out=TRUE)
ramVF(test.blcs, c(0,80),c(0,80), length=.05, xlab='X', ylab='Y',scale=.5, ninterval=9)
}