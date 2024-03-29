\name{ramFit}
\alias{ramFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a model using lavaan based on ram input
}
\description{
Fit a model using lavaan based on ram input}
\usage{
ramFit(ramModel, data, type=c('ram','lavaan'), digits = 3, zero.print = "0", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ramModel}{
An ram model
}
  \item{data}{
data
}
\item{type}{
  ram: specify a ram model; lavaan: specify a lavaan model
}
  \item{digits}{
Digits for print
}
  \item{zero.print}{
Format of zeros
}
  \item{\dots}{
Options for lavaan
}
}

\value{
\item{A and Ase}{A matrix and its standard error}
\item{S and Sse}{S matrix and its standard error}
\item{lavaan}{Original lavaan output}
\item{fit}{Model fit statistics and indices}
}
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. URL http://www.jstatsoft.org/v48/i02/.

Zhang, Z., Hamagami, F., Grimm, K. J., & McArdle, J. J. (2015). Using R package RAMpath for tracing SEM path diagrams and conducting complex longitudinal data analysis. Structural Equation Modeling, 22(1), 132-147. https://doi.org/10.1080/10705511.2014.935257
}
\examples{
## Example 1. A path model
data(ex1)
m1<-'
manifest=3
label=age,hvlt,ept
arrow(2,1)=?
arrow(3,1)=?
arrow(3,2)=?
sling(1,1)=?
sling(2,2)=?
sling(3,3)=?
'

## Fit the model
res1<-ramFit(m1, ex1)

## More output from Lavaan
summary(res1$lavaan, fit=TRUE)

## Effects and variance decomposition
bridge<-ramPathBridge(res1, allbridge=TRUE, indirect=TRUE)
summary(bridge)
summary(bridge, type='bridge')

## plot the path diagram
## uncomment to plot
## plot(bridge, 'ex1')

## plot the effects from age to ept
## uncomment to plot
## plot(bridge, 'ex1effect', 'age','ept')

## plot the bridges for ept
## uncomment to plot
## plot(bridge, 'ex1bridge', 'ept','hvlt', type='bridge')

## summarize
summary(bridge)
summary(bridge, type='bridge')


## Example 2: An SEM model (MIMIC model)
data(ex2)
## Using lavaan directly for model estimation and specification
mimic<-'
R =~ ws1 + ls1 + lt1
R ~ edu + gender
'

mimic.res<-sem(mimic, data=ex2)

mimic.ram<-lavaan2ram(mimic.res)

## plot the path diagram
bridge<-ramPathBridge(mimic.ram, allbridge=FALSE, indirect=FALSE)
## uncomment to plot
## plot(bridge, 'mimic')

}