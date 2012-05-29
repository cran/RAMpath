\name{makePathList}
\alias{makePathList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a list of effects
}
\description{
Make a list of effects}
\usage{
makePathList(AMatrix, Ase, indirect = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AMatrix}{
A matrix from the ram matrices
}
 \item{Ase}{
Standard error matrix for A matrix from the ram matrices
}
  \item{indirect}{
Whether to generate all indirect effects}
}
\references{
Boker, S. M., McArdle, J. J. & Neale, M. C. (2002) An algorithm for the hierarchical organization of path diagrams and calculation of components of covariance between variables. Structural Equation Modeling, 9(2), 174-194

Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2013). Using R Package RAMpath for Tracing SEM Path Diagrams and Conducting Complex Longitudinal Data Analysis. Structural Equation Modeling.
}
