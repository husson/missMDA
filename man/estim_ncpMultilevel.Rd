\name{estim_ncpMultilevel}
\alias{estim_ncpMultilevel}
\title{Estimate the number of dimensions for the Multilevel PCA, multlevel MCA or Multilevel FAMD by cross-validation}
\description{Estimate the number of dimensions for Multilevel principal component (multilevel PCA, multilevel MCA ormultilevel Factorial Analysis of Mixed Data) by cross-validation}
\usage{
estim_ncpMultilevel(X,  ifac=1, ncpW.min = 1, ncpW.max = 5, ncpB.min = 1, 
    ncpB.max = 5, scale = TRUE, nbsim=100, pNA=0.05, threshold=1e-4, 
	nb.cores = NULL, verbose = TRUE)
}
\arguments{
  \item{X}{a data.frame with categorical variables; with missing entries or not}
  \item{ifac}{index of the group variable}
  \item{ncpB.min}{integer corresponding to the minimum number of components to test for the between matrix}
  \item{ncpB.max}{integer corresponding to the maximum number of components to test for the between matrix}
  \item{ncpW.min}{integer corresponding to the minimum number of components to test for the within matrix}
  \item{ncpW.max}{integer corresponding to the maximum number of components to test for the within matrix}
  \item{scale}{if all the variables are continuous, should they be standardized? Yes if true.}
  \item{nbsim}{number of simulations}
  \item{pNA}{percentage of missing values added in the data set, useful only if method.cv="Kfold"}
  \item{threshold}{the threshold for assessing convergence}
  \item{nb.cores}{Integer, number of core used. By default, NULL and the number of cores used are the number of cores of your computer  minus 1}
  \item{verbose}{boolean. TRUE means that a progressbar is writtent}
}

\details{
pNA percentage of missing values is inserted at random in the data matrix and predicted with a multilevel model using ncpB.min to ncpB.max and ncpW.min to ncpW.max dimensions. 
This process is repeated nbsim times. The number of components which leads to the smallest MSEP is retained. 
More precisely, the missing entries are predicted using the imputeMultilevel function.
}
\value{
  \item{ncp}{the number of components retained for the FAMD}
  \item{criterion}{the criterion (the MSEP) calculated for each number of components}
}

\seealso{\code{\link{imputeMultilevel}}}
\examples{
\dontrun{
data(ozone)
result <- estim_ncpMultilevel(ozone, ifac=12)
}
}
\keyword{multivariate}
