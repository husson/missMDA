\name{imputeFAMD}
\alias{imputeFAMD}
\title{Impute mixed dataset}
\description{Impute the missing values of a mixed dataset (with continuous and categorical variables) using the principal component method "factorial analysis for mixed data" (FAMD). Can be used as a preliminary step before performing FAMD on an incomplete dataset.}

\usage{imputeFAMD(X, ncp = 2, method=c("Regularized","EM"), row.w = NULL,
      coeff.ridge=1,threshold = 1e-06, ind.sup=NULL, sup.var=NULL,
	  seed = NULL, maxiter = 1000,...)
	   }

\arguments{
  \item{X}{a data.frame with continuous and categorical variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to predict the missing entries}
  \item{method}{"Regularized" by default or "EM"}
  \item{row.w}{row weights (by default, uniform row weights)}
  \item{coeff.ridge}{1 by default to perform the regularized imputeFAMD algorithm; useful only if method="Regularized". Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of the EM method) or more than 1 to regularized more}
  \item{threshold}{the threshold for assessing convergence}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{sup.var}{a vector indicating the indexes of the supplementary variables (quantitative or categorical)}
  \item{seed}{integer, by default seed = NULL implies that missing values are initially imputed by the mean of each variable for the continuous variables and by the proportion of the category for the categorical variables coded with indicator matrices of dummy variables. Other values leads to a random initialization}
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Impute the missing entries of a mixed data using the iterative FAMD algorithm (method="EM") or the regularised iterative FAMD algorithm (method="Regularized"). The (regularized) iterative FAMD algorithm first consists in coding the categorical variables using the indicator matrix
of dummy variables. Then, in the initialization step, missing values are imputed with initial values such as the mean of the variable for the continuous variables and the proportion of the category for each category using the non-missing entries. If the argument seed is set to a specific value, a random initialization is performed: the initial values are drawn from a gaussian distribution
with mean and standard deviation calculated from the observed values for each continuous variable. The second step of the (regularized) iterative FAMD algorithm is to perform FAMD on the completed dataset. Then, it imputes the missing values with the (regularized) reconstruction formulae of order ncp (the fitted matrix computed with ncp components for the (regularized) scores and loadings). These steps of estimation of the parameters via FAMD and imputation of the missing values using the (regularized) fitted matrix are iterate until convergence. \cr

We advice to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. In the regularized algorithm, the singular values of the FAMD are shrinked.\cr

The output of the algorithm can be used as an input of the FAMD function of the FactoMineR package in order to perform FAMD on an incomplete dataset.
}
\value{
 \item{tab.disj}{the imputed matrix; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones. The categorical variables are coded with the indicator matrix of dummy variables. In this indicator matrix, the imputed values are real numbers but they met the constraint that the sum of the entries corresponding to one individual and one variable is equal to one. Consequently they can be seen as degree of membership to the corresponding category}
  \item{completeObs}{the mixed imputed dataset; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones. For the continuous variables, the values are the same as in the tab.disj output; for the categorical variables missing values are imputed with the most plausible categories according to the values in the tab.disj output}
 
  \item{call}{the matched call}
}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\references{Audigier, V., Husson, F. & Josse, J. (2013). A principal components method to impute mixed data.  Advances in Data Analysis and Classification, 10(1), 5-26. https://arxiv.org/abs/1301.4797}

\seealso{\code{\link{imputePCA}}}
\examples{
\dontrun{
data(ozone)
res.impute <- imputeFAMD(ozone, ncp=3) 
## The output can be used as an input of the FAMD function of the FactoMineR package 
##to perform the FAMD on the incomplete data ozone 
require(FactoMineR)
res.afdm <- FAMD(ozone,tab.disj=res.impute$tab.disj) 
}
}
\keyword{models}
\keyword{multivariate}
