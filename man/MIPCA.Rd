\name{MIPCA}
\alias{MIPCA}
\title{Multiple Imputation with PCA}
\description{
MIPCA performs Multiple Imputation with a PCA model. Can be used as a preliminary step to perform Multiple Imputation in PCA.}
\usage{
MIPCA(X, ncp = 2, scale = TRUE, method=c("Regularized","EM"), threshold = 1e-04, 
    nboot = 100,  method.mi="Boot", Lstart=1000, L=100, verbose=FALSE)
}
\arguments{
  \item{X}{a data.frame with continuous variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the PCA reconstruction formulae}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable}
  \item{method}{"Regularized" by default or "EM"}
  \item{threshold}{the threshold for the criterion convergence}
  \item{nboot}{the number of imputed datasets}
  \item{method.mi}{a string. If "Bayes", the uncertainty on the parameters of the imputation model is taken into account using a Bayesian treatment of PCA. By default "Boot" leading to a MI which reflect uncertainty a bootstrap procedure. See details.}
  \item{Lstart}{number of iterations for the burn-in period (only used if method.mi="Bayes")}
  \item{L}{number of skipped iterations to keep one imputed data set after the burn-in period (only used if method.mi="Bayes")}
  \item{verbose}{use verbose=TRUE for screen printing of iteration numbers}
  
}
\details{
MIPCA generates nboot imputed datasets from a PCA model. The observed values are the same from one dataset to the others whereas the imputed values change. The variation among the imputed values reflects the variability with which missing values can be predicted.
The multiple imputation is proper in the sense of Little and Rubin (2002) since it takes into account the variability of the parameters. Two versions are available: multiple imputation using a parametric bootstrap (Josse, J., Husson, F. (2010)) and multiple imputation using a Bayesian treatment of the PCA model (Audigier et al 2015). The methods differ by the way in which the variability due to missing values is reflected. The method used is controlled by the \code{method.mi} argument. By default, MIPCA uses the parametric bootstrap \code{method.mi="Boot"}. This bootstrap method is more recommended to evaluate uncertainty in PCA (through confidence ellipses). Otherwise, the Bayesian method can be used by specifying the argument \code{method.mi="Bayes"}. It is based on an iterative algorithm which alternates imputation of the data set and draw of the PCA parameters in a posterior distribution. These steps are repeated \code{Lstart} times to reach a convergence. Then, one imputed data set is kept each \code{L} iterations to ensure independence between imputed values from a data set to another. The Bayesian method is more recommanded to apply a statistical method on an incomplete data set.
}
\value{
\item{res.imputePCA}{A matrix corresponding to the imputed dataset obtained with the function imputePCA (the completed dataset)}
\item{res.MI}{A list of data frames corresponding to the nboot imputed data sets}
\item{call}{the matched call}
}

\references{Josse, J., Husson, F. (2011). Multiple Imputation in PCA. Advances in Data Analysis and Classification.

Audigier, V. Josse, J., Husson, F. (2015). Multiple imputation for continuous variables using a Bayesian principal component analysis. Journal of Statistical Computation and Simulation.

Little R.J.A., Rubin D.B. (2002) Statistical Analysis with Missing Data. Wiley series in probability and statistics, New-York.
}

\author{Francois Husson  \email{francois.husson@institut-agro.fr}, Julie Josse  \email{julie.josse@polytechnique.edu} and Vincent Audigier}

\seealso{\code{\link{imputePCA}},\code{\link{plot.MIPCA}},\code{\link{Overimpute}},\code{\link{MIMCA}},\code{\link{with.mids}},\code{\link{pool}},\code{\link{summary.mira}}}
\examples{
\dontrun{
#########################################################
## Multiple Imputation for visualization on the PCA map
#########################################################

data(orange)
## First the number of components has to be chosen 
##   (for the reconstruction step)
nb <- estim_ncpPCA(orange,ncp.max=4)

## Multiple Imputation
resMI <- MIPCA(orange,ncp=2)

## Visualization on the PCA map
plot(resMI)

#########################################################
## Multiple Imputation for applying statistical methods
(Bayesian method)
#########################################################
data(ozone)

## First the number of components has to be chosen 
nb <- estim_ncpPCA(ozone[,1:11])

## Multiple Imputation with Bayesian method
res.BayesMIPCA<-MIPCA(ozone[,1:11],ncp=2,method.mi="Bayes",verbose=TRUE)

## Regression on the multiply imputed data set and pooling with mice
require(mice)
imp<-prelim(res.mi=res.BayesMIPCA,X=ozone[,1:11])#creating a mids object
fit <- with(data=imp,exp=lm(maxO3~T9+T12+T15+Ne9+Ne12+Ne15+Vx9+Vx12+Vx15+maxO3v))#analysis
res.pool<-pool(fit);summary(res.pool)#pooling

## Diagnostics
res.over<-Overimpute(res.BayesMIPCA)
}
}
\keyword{multivariate}
\keyword{imputation}
