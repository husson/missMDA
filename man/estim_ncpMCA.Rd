\name{estim_ncpMCA}
\alias{estim_ncpMCA}
\title{Estimate the number of dimensions for the Multiple Correspondence Analysis by cross-validation}
\description{Estimate the number of dimensions for the Multiple Correspondence Analysis by cross-validation}
\usage{
estim_ncpMCA(don, ncp.min=0, ncp.max=5,  method = c("Regularized","EM"), 
     method.cv = c("Kfold","loo"), nbsim=100, pNA=0.05, ind.sup=NULL,
	 quanti.sup=NULL, quali.sup=NULL, threshold=1e-4,verbose = TRUE)
}
\arguments{
  \item{don}{a data.frame with categorical variables; with missing entries or not}
  \item{ncp.min}{integer corresponding to the minimum number of components to test}
  \item{ncp.max}{integer corresponding to the maximum number of components to test}
  \item{method}{"Regularized" by default or "EM"}
  \item{method.cv}{"Kfold" for cross-validation or "loo" for leave-one-out}
  \item{nbsim}{number of simulations, useful only if method.cv="Kfold"}
  \item{pNA}{percentage of missing values added in the data set, useful only if method.cv="Kfold"}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{quanti.sup}{a vector indicating the indexes of the quantitative supplementary variables}
  \item{quali.sup}{a vector indicating the indexes of the categorical supplementary variables}
  \item{threshold}{the threshold for assessing convergence}
  \item{verbose}{boolean. TRUE means that a progressbar is writtent}
}

\details{
For leave-one-out cross-validation (method.cv="loo"), each cell of the data matrix is alternatively removed and predicted with a MCA model using ncp.min to ncp.max dimensions. The number of components which leads to the smallest mean square error of prediction (MSEP) is retained.
For the Kfold cross-validation (method.cv="Kfold"), pNA percentage of missing values is inserted at random in the data matrix and predicted with a MCA model using ncp.min to ncp.max dimensions. This process is repeated nbsim times.  The number of components which leads to the smallest MSEP is retained. More precisely, for both cross-validation methods, the missing entries are predicted using the imputeMCA function, it means using it means using the regularized iterative MCA algorithm (method="Regularized") or the iterative MCA algorithm (method="EM"). The regularized version is more appropriate to avoid overfitting issues.
}
\value{
  \item{ncp}{the number of components retained for the MCA}
  \item{criterion}{the criterion (the MSEP) calculated for each number of components}
}
\references{Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis, Journal of Clcassification, 29 (1), pp. 91-116.}

\author{Francois Husson  \email{francois.husson@institut-agro.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}
\seealso{\code{\link{imputeMCA}}}
\examples{
\dontrun{
data(vnf)
result <- estim_ncpMCA(vnf,ncp.min=0, ncp.max=5)
}
}
\keyword{multivariate}
