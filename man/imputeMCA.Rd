\name{imputeMCA}
\alias{imputeMCA}

\title{Impute categorical dataset}
\description{Impute the missing values of a categorical dataset using Multiple Correspondence Analysis (MCA). Can be used as a preliminary step before performing MCA on an incomplete dataset.}
\usage{
imputeMCA(don, ncp=2, method = c("Regularized","EM"), row.w=NULL, coeff.ridge=1, 
    threshold=1e-06, ind.sup = NULL, quanti.sup=NULL, quali.sup=NULL,
	seed=NULL, maxiter=1000)
}

\arguments{
  \item{don}{a data.frame with categorical variables containing missing values}
  \item{ncp}{integer corresponding to the number of dimensions used  to predict the missing entries}
  \item{method}{"Regularized" by default or "EM"}
  \item{row.w}{row weights  (by default, a vector of 1 for uniform row weights)}
  \item{coeff.ridge}{1 by default to perform the regularized imputeMCA algorithm; useful only if method="Regularized". Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of the EM method) or more than 1 to regularized more (to get closer to the results of the proportion imputation)}
  \item{threshold}{the threshold for assessing convergence}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{quanti.sup}{a vector indicating the indexes of the quantitative supplementary variables}
  \item{quali.sup}{a vector indicating the indexes of the categorical supplementary variables}
  \item{seed}{integer, by default seed = NULL implies that missing values are initially imputed by the proportion of the category for the categorical variables coded with indicator matrices of dummy variables. Other values leads to a random initialization}
  \item{maxiter}{integer, maximum number of iterations for the regularized iterative MCA algorithm}
}
\details{
Impute the missing entries of a categorical data using the iterative MCA algorithm (method="EM") or the regularised iterative MCA algorithm (method="Regularized"). The (regularized) iterative MCA algorithm first consists in coding the categorical variables using the indicator matrix
of dummy variables. Then, in the initialization step, missing values are imputed with initial values such as the proportion of the category for each category using the non-missing entries. This imputation corresponds also to using the algorithm with ncp=0 and is sometimes called in the literature the "missing fuzzy average method". If the argument seed is set to a specific value, a random initialization is performed: random values are drawn in such a way that the constraint that the sum of the entries corresponding to one individual and one variable is equal to one in the indicator matrix of dummy variables. 
The second step of the (regularized) iterative MCA algorithm consists in performing MCA on the completed dataset. Then, it imputes the missing values with the (regularized) reconstruction formulae of order ncp (the fitted matrix computed with ncp components for the (regularized) scores and loadings). These steps of estimation of the parameters via MCA and imputation of the missing values using the (regularized) fitted matrix are iterate until convergence. \cr
We advice to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. In the regularized algorithm, the singular values of the MCA are shrinked. \cr
The number of components ncp used in the algorithm can be selected using the function ncpMCA. A small number of components can also be seen as a way to regularize more and consequently may be advices to get more stable predictions.
\cr
The output of the algorithm can be used as an input of the MCA function of the FactoMineR package in order to perform MCA on an incomplete dataset.

}

\value{
\item{tab.disj}{The imputed indicator matrix; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones.  The imputed values are real numbers but they  but they met the constraint that the sum of the entries corresponding to one individual and one variable is equal to one. Consequently they can be seen as degree of membership to the corresponding category.
}
\item{completeObs}{The categorical imputed dataset; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones. Missing values are imputed with the most plausible categories according to the values in the tab.disj output}

}

\author{
Francois Husson  \email{francois.husson@institut-agro.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\references{Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis, Journal of Clcassification, 29 (1), pp. 91-116.\cr
Josse, J. and Husson, F. missMDA (2016). A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31 \doi{10.18637/jss.v070.i01}}

\seealso{\code{\link{estim_ncpMCA}},\cr
\href{https://www.youtube.com/watch?v=_Wa6R4PM9dY&list=PLnZgp6epRBbQzxFnQrcxg09kRt-PA66T_&index=1}{Video showing how to perform MCA on an incomplete dataset}}


\examples{
\dontrun{
data(vnf)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(vnf,ncp.max=5) ## Time-consuming, nb = 4

## Impute the indicator matrix and perform a MCA
res.impute <- imputeMCA(vnf, ncp=4)

## The imputed indicator matrix can be used as an input of the MCA function of the
## FactoMineR package to perform the MCA on the incomplete data vnf 
require(FactoMineR)
res.mca <- MCA(vnf,tab.disj=res.impute$tab.disj) 

## With supplementary variables (var 11 to 14), impute the active ones
res.impute <- imputeMCA(vnf[,1:10], ncp=4)
res.mca <- MCA(vnf,tab.disj=res.impute$tab.disj,quali.sup=11:14) 
}
}
\keyword{models}
\keyword{multivariate}
\keyword{imputation}
