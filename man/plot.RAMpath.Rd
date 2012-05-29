\name{plot.RAMpath}
\alias{plot.RAMpath}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the path diagram according to RAM path and bridges
}
\description{
Plot the path diagram according to RAM path and bridges}
\usage{
\method{plot}{RAMpath}(x, file, from, to, type = c("path", "bridge"), size = c(8, 8), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = c("LR", "TB"), digits = 2, output.type = c("graphics", "dot"), graphics.fmt = "pdf", dot.options = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Output from the \code{\link{ramPathBridge}} function
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
  \item{\dots}{
Other options for plot
}
}
