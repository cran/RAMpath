\name{ramLCM}
\alias{ramLCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conduct growth curve analysis
}
\description{
Conduct growth curve analysis}
\usage{
ramLCM(data, outcome, model = c("all", "no", "linear", "quadratic", "latent"), 
basis = 0:(length(outcome) - 1), predictor, equal.var = TRUE, digits = 3, 
ram.out = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data
}
  \item{outcome}{
Outcome variable indices
}
  \item{model}{
Models to fit
}
  \item{basis}{
Basis coefficients
}
  \item{predictor}{
Covariates as predictors
}
  \item{equal.var}{
Set residual variances to be equal
}
  \item{digits}{
Print digits
}
  \item{ram.out}{
Print ram matrices
}
  \item{\dots}{
Options can be used for lavaan
}
}
\value{
\item{model}{The lavaan model specification of the bivariate latent change score model}
\item{lavaan}{The lavaan output}
\item{ram}{Output in terms of RAM matrices}
\item{fit}{Model fit}
}
\references{
Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2015). Using R package RAMpath for tracing SEM path diagrams and conducting complex longitudinal data analysis. Structural Equation Modeling, 22(1), 132-147. https://doi.org/10.1080/10705511.2014.935257
}

\examples{
data(ex3)
## Example 3. Growth curve models
gcm.all<-ramLCM(ex3, 1:6, ram.out=TRUE)
## plot the path diagram
bridge<-ramPathBridge(gcm.all$ram$latent, FALSE, FALSE)
## uncomment to plot
## plot(bridge, 'latent')

##unequal variance
gcm.all<-ramLCM(ex3, 1:6, ram.out=TRUE, equal.var=FALSE)

## missing data
gcm.all<-ramLCM(ex3, c(1,2,4,6), basis=c(1,2,4,6), ram.out=TRUE)

gcm.l<-ramLCM(ex3, 1:6, model='linear', ram.out=TRUE)

## with a predictor
gcm.pred<-ramLCM(ex3, c(1,2,4,6), model='linear', basis=c(1,2,4,6), 
                 predictor=c(3,5), ram.out=TRUE)
bridge3<-ramPathBridge(gcm.pred$ram$linear)
## uncomment to plot
## plot(bridge3, 'gcmlinear')
}