\name{imputePCA}
\alias{imputePCA}
\title{Impute dataset with PCA}
\description{Impute the missing values of a dataset with the Principal Components Analysis model. Can be used as a preliminary step before performing a PCA on an completed dataset.}
\usage{imputePCA(X, ncp = 2, scale = TRUE, method = c("Regularized","EM"), 
       row.w = NULL, ind.sup=NULL,quanti.sup=NULL,quali.sup=NULL,
	   coeff.ridge = 1, threshold = 1e-06, seed = NULL, nb.init = 1,  
	   maxiter = 1000, ...)
}

\arguments{
  \item{X}{a data.frame with continuous variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to  to predict the missing entries}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable}
  \item{method}{"Regularized" by default or "EM"}
  \item{row.w}{row weights  (by default, a vector of 1 for uniform row weights)}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{quanti.sup}{a vector indicating the indexes of the quantitative supplementary variables}
  \item{quali.sup}{a vector indicating the indexes of the categorical supplementary variables}
  \item{coeff.ridge}{1 by default to perform the regularized imputePCA algorithm; useful only if method="Regularized". Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of the EM method) or more than 1 to regularized more (to get closer to the results of the mean imputation)}
  \item{threshold}{the threshold for assessing convergence}
  \item{seed}{integer, by default seed = NULL implies that missing values are initially imputed by the mean of each variable. Other values leads to a random initialization}
  \item{nb.init}{integer corresponding to the number of random initializations; the first initialization is the initialization with the mean imputation}
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Impute the missing entries of a mixed data using the iterative PCA algorithm (method="EM") or the regularised iterative PCA algorithm (method="Regularized"). The (regularized) iterative PCA algorithm first consists imputing missing values with initial values such as the mean of the variable. If the argument seed is set to a specific value, a random initialization is performed: the initial values are drawn from a gaussian distribution
with mean and standard deviation calculated from the observed values. nb.init different random initialization can be drawn. In such a situation, the solution giving the smallest objective function (the mean square error between the fitted matrix and the observed one) is kept.  The second step of the (regularized) iterative PCA algorithm is to perform PCA on the completed dataset. Then, it imputes the missing values with the (regularized) reconstruction formulae of order ncp (the fitted matrix computed with ncp components for the (regularized) scores and loadings). These steps of estimation of the parameters via PCA and imputation of the missing values using the (regularized) fitted matrix are iterate until convergence. The iterative PCA algorithm is also known as the EM-PCA algorithm since it corresponds to an EM algorithm of the fixed effect model where the data are generated as a fixed structure (with a low rank representation) corrupted by noise. The number of components used in the algorithm can be found using cross-validation criteria implemented in the function estim_ncpPCA.\cr

We advice to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. In the regularized algorithm, the singular values of the PCA are shrinked.\cr

The output of the algorithm can be used as an input of the PCA function of the FactoMineR package in order to perform PCA on an incomplete dataset.}
\value{
  \item{completeObs}{the imputed dataset; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones.}
  \item{fittedX}{the reconstructed data}
}

\references{Josse, J & Husson, F. (2013). Handling missing values in exploratory multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.\cr
Josse, J. and Husson, F. missMDA (2016). A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31 <doi:10.18637/jss.v070.i01>}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\seealso{\code{\link{estim_ncpPCA}}, \code{\link{MIPCA}},\cr
\href{https://www.youtube.com/watch?v=YDbx2pk9xNY&list=PLnZgp6epRBbQzxFnQrcxg09kRt-PA66T_&index=2}{Video showing how to perform PCA on an incomplete dataset}}
\examples{
\dontrun{
data(orange)
## First the number of components has to be chosen 
##   (for the imputation step)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

## Imputation
res.comp <- imputePCA(orange,ncp=2)

## A PCA can be performed on the imputed data 
res.pca <- PCA(res.comp$completeObs)
}
}
\keyword{models}
\keyword{multivariate}
\keyword{imputation}
