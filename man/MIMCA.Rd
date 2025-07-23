\name{MIMCA}
\alias{MIMCA}
\title{Multiple Imputation with MCA}
\description{
MIMCA performs multiple imputations for categorical data using Multiple Correspondence Analysis.}
\usage{
MIMCA(X, nboot=100, ncp, coeff.ridge=1, threshold = 1e-06, maxiter = 1000, verbose=FALSE)
}
\arguments{
  \item{X}{a data.frame with categorical variables containing missing values}
  \item{nboot}{the number of imputed datasets}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the MCA reconstruction formulae}
  \item{coeff.ridge}{1 by default to perform the regularized imputeMCA algorithm. Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of an EM method) or more than 1 to regularized more (to get closer to the results of the proportion imputation)}
  \item{threshold}{the threshold for assessing convergence for the (regularized) iterative MCA algorithm}
  \item{maxiter}{integer, maximum number of iterations for the (regularized) iterative MCA algorithm}
  \item{verbose}{use verbose=TRUE for screen printing of iteration numbers}
}
\details{
MIMCA generates nboot imputed data sets from MCA. The observed values are the same from one dataset to the others whereas the imputed values change. First, nboot weightings are defined for the individuals. Then, the iterative regularized MCA algorithm (Josse, 2012) is applied according to each weighting, leading to nboot imputed tables. These imputed tables are scaled to verify the constraint that the sum is equal to one per variable and per individual. Lastly, missing categories are drawn from the probabilities given by the imputed tables. Thus, nboot imputed categorical data sets are obtained. The variation among the imputed values reflects the variability with which missing values can be predicted. The multiple imputation is proper in the sense of Little and Rubin (2002) since it takes into account the variability of the parameters using a non-parametric bootstrap approach.
}
\value{
\item{res.MI}{A list of data frames corresponding to the nboot imputed categorical data sets}
\item{res.imputeMCA}{A matrix corresponding to the single imputed disjunctive table obtained with the function imputeMCA}
\item{call}{The matched call}

}

\references{Audigier, V., Husson, F., Josse, J. (2015). MIMCA: Multiple imputation for categorical variables with multiple correspondence analysis.

Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis, Journal of Classification, 29 (1), pp. 91-116.

Little R.J.A., Rubin D.B. (2002) Statistical Analysis with Missing Data. Wiley series in probability and statistics, New-York.

}

\author{Vincent Audigier \email{vincent.audigier@cnam.fr}, Francois Husson  \email{francois.husson@institut-agro.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\seealso{\code{\link{imputeMCA}},\code{\link{MIPCA}},\code{\link{estim_ncpMCA}},\code{\link[mice]{with.mids}},\code{\link[mice]{pool}},\code{\link[mice]{summary.mira}}}
\examples{
\dontrun{
data(TitanicNA)

## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(TitanicNA) ## Time-consuming, nb = 5

## Multiple Imputation
res.mi <- MIMCA(TitanicNA, ncp=5, verbose=TRUE)

## First completed data matrix
res.mi$res.MI[[1]]
 
## Analysis and pooling with mice
require(mice)
imp<-prelim(res.mi,TitanicNA)
fit <- with(data=imp,exp=glm(SURV~CLASS+AGE+SEX,family = "binomial"))
res.pool<-pool(fit)
summary(res.pool)

}
}
\keyword{multivariate}
\keyword{imputation}
