\name{estim_ncpPCA}
\alias{estim_ncpPCA}
\title{Estimate the number of dimensions for the Principal Component Analysis by cross-validation}
\description{Estimate the number of dimensions for the Principal Component Analysis by cross-validation}
\usage{
estim_ncpPCA(X, ncp.min = 0, ncp.max = 5, method = c("Regularized","EM"), 
       scale = TRUE, method.cv = c("gcv","loo","Kfold"), nbsim = 100, 
	   pNA = 0.05, threshold=1e-4, verbose = TRUE)
}
\arguments{
  \item{X}{a data.frame with continuous variables; with missing entries or not}
  \item{ncp.min}{integer corresponding to the minimum number of components to test}
  \item{ncp.max}{integer corresponding to the maximum number of components to test}
  \item{method}{"Regularized" by default or "EM"}
  \item{scale}{boolean. TRUE implies a same weight for each variable}
  \item{method.cv}{string with the values "gcv" for generalised cross-validation, "loo" for leave-one-out or "Kfold" cross-validation}
  \item{nbsim}{number of simulations, useful only if method.cv="Kfold"}
  \item{pNA}{percentage of missing values added in the data set, useful only if method.cv="Kfold"}
  \item{threshold}{the threshold for assessing convergence}
  \item{verbose}{boolean. TRUE means that a progressbar is writtent}
}

\details{

For leave-one-out (loo) cross-validation, each cell of the data matrix is alternatively removed and predicted with a PCA model using ncp.min to ncp.max dimensions. The number of components which leads to the smallest mean square error of prediction (MSEP) is retained.
For the Kfold cross-validation, pNA percentage of missing values is inserted and predicted with a PCA model using ncp.min to ncp.max dimensions. This process is repeated nbsim times. The number of components which leads to the smallest MSEP is retained. \cr
For both cross-validation methods, missing entries are predicted using the imputePCA function, it means using the regularized iterative PCA algorithm (method="Regularized") or the iterative PCA algorithm (method="EM"). The regularized version is more appropriate when there are already many missing values in the dataset to avoid overfitting issues.
\cr
Cross-validation (especially method.cv="loo") is time-consuming. The generalised cross-validation criterion (method.cv="gcv") can be seen as an approximation of the loo cross-validation criterion which provides a straightforward way to estimate the number of dimensions without resorting to a computationally intensive method. \cr

This argument scale has to be chosen in agreement with the PCA that will be performed. If one wants to perform a normed PCA (where the variables are centered and scaled, i.e. divided by their standard deviation), then the argument scale has to be set to the value TRUE.


}
\value{
  \item{ncp}{the number of components retained for the PCA}
  \item{criterion}{the criterion (the MSEP) calculated for each number of components}
}
\references{Bro, R., Kjeldahl, K. Smilde, A. K. and Kiers, H. A. L. (2008) Cross-validation of component models: A critical look at current methods. Analytical and Bioanalytical Chemistry, 5, 1241-1251.\cr

Josse, J. and Husson, F. (2011). Selecting the number of components in PCA using cross-validation approximations. Computational Statistics and Data Analysis. 56 (6), pp. 1869-1879.}


\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}
\seealso{\code{\link{imputePCA}}}
\examples{
\dontrun{
data(orange)
nb <- estim_ncpPCA(orange,ncp.min=0,ncp.max=4) 
}
}
\keyword{multivariate}
