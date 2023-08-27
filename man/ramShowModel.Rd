\name{ramShowModel}
\alias{ramShowModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show the model using Lavvan model syntax}
\description{
Show the model using Lavvan model syntax
}
\usage{
ramShowModel(object)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Output from any data analysis
}
}
\examples{
data(ex3)
gcm.l<-ramLCM(ex3, 1:6, model='linear', ram.out=TRUE)
## Add correlated errors
ramShowModel(gcm.l)
}

\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. URL http://www.jstatsoft.org/v48/i02/.

Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2015). Using R package RAMpath for tracing SEM path diagrams and conducting complex longitudinal data analysis. Structural Equation Modeling, 22(1), 132-147. https://doi.org/10.1080/10705511.2014.935257
}