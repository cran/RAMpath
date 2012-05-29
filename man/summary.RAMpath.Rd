\name{summary.RAMpath}
\alias{summary.RAMpath}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the total and individual contribution for each path and bridge
}
\description{
Calculate the total and individual contribution for each path and bridge}
\usage{
\method{summary}{RAMpath}(object, from, to, type = c("path", "bridge"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Output from the ramPathBridge function
}
  \item{from}{
from variable: starting from this variable
}
  \item{to}{
to variable: end on this variable
}
  \item{type}{
path: to calculate the effect; bridge: to calculate the bridges
}
  \item{...}{
Other options
}
}
