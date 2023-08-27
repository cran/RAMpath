\name{plot.RAMpath}
\alias{plot.RAMpath}
\alias{plot}
\alias{plot.blcs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the path diagram according to RAM path and bridges or Plot the vector field for the bivariate latent change score model
}
\description{
Plot the path diagram according to RAM path and bridges or Plot the vector field for the bivariate latent change score model}
\usage{
\method{plot}{RAMpath}(x, file, from, to, type = c("path", "bridge"), 
size = c(8, 8), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10),
rank.direction = c("LR", "TB"), digits = 2, output.type = c("graphics", "dot"),
graphics.fmt = "pdf", dot.options = NULL, ...)

\method{plot}{blcs}(x, ylim, xlim, ninterval=10, scale=.1, length=.25, 
scatter=TRUE, n=20, alpha=.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Output from the \code{\link{ramPathBridge}} function or Output from the \code{\link{ramBLCS}} function
}
  \item{file}{
File name for the generated figures
}
  \item{from}{
from variable: path starts from this variable
}
  \item{to}{
to variable: path ends on this variable
}
  \item{type}{
\code{path}: to plot the effect path; 

\code{bridge}: to plot the bridges path
}
  \item{size}{
The size of the plot in inches
}
  \item{node.font}{
The size of the text for the variables
}
  \item{edge.font}{
The size of the text on the pahts
}
  \item{rank.direction}{
LR: from left to right; TB: from top to bottom.
}
  \item{digits}{
Digits of numbers to plot
}
  \item{output.type}{
If "graphics", the default, both a ".dot" file and a graphics file will be created.
}
  \item{graphics.fmt}{
a graphics format recognized by the dot program; the default is "pdf"; graphics.fmt is also used for the extension of the graphics file that is created.
}
  \item{dot.options}{
options to be passed to the dot program, given as a character string.
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
plot(test.blcs, c(0,80),c(0,80), length=.05, xlab='X', ylab='Y',scale=.5, ninterval=9)
}