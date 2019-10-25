\name{imputeMFA}
\alias{imputeMFA}
\title{Impute dataset with variables structured into groups of variables (groups of continuous or categorical variables)}
\description{Impute the missing values of a dataset with Multiple Factor Analysis (MFA). The variables are structured a priori into groups of variables. The variables can be continuous or categorical but within a group the nature of the variables is the same.  Can be used as a preliminary step before performing MFA on an incomplete dataset.}

\usage{imputeMFA(X, group, ncp = 2, type=rep("s",length(group)), method = c("Regularized","EM"), 
       row.w = NULL, coeff.ridge = 1,threshold = 1e-06, ind.sup = NULL,  
	   num.group.sup = NULL, seed = NULL, maxiter = 1000, ...)
	   }

\arguments{
  \item{X}{a data.frame with groups of continuous or categorical variables containing missing values}
  \item{group}{a vector indicating the number of variables in each group}
  \item{ncp}{integer corresponding to the number of components used to predict the missing entries}
  \item{type}{the type of variables in each group; three possibilities: "c" or "s" for continuous variables
  (for "c" the variables are centered and for "s" variables are scaled to unit variance), "n" for categorical variables}
  \item{method}{"Regularized" by default or "EM"}
  \item{row.w}{row weights (by default, a vector of 1 for uniform row weights)}
  \item{coeff.ridge}{1 by default to perform the regularized imputeMFA algorithm; useful only if method="Regularized". Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of the EM method) or more than 1 to regularized more}
  \item{threshold}{the threshold for assessing convergence}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{num.group.sup}{a vector indicating the group of variables that are supplementary}
  \item{seed}{integer, by default seed = NULL implies that missing values are initially imputed by the mean of each variable for the continuous variables and by the proportion of the category for the categorical variables coded with indicator matrices of dummy variables. Other values leads to a random initialization}
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{

Impute the missing entries of a data with groups of variables using the iterative MFA algorithm (method="EM") or the regularised iterative MFA algorithm (method="Regularized"). The (regularized) iterative MFA algorithm first consists in coding the categorical variables using the indicator matrix
of dummy variables. Then, in the initialization step, missing values are imputed with initial values such as the mean of the variable for the continuous variables and the proportion of the category for each category using the non-missing entries. If the argument seed is set to a specific value, a random initialization is performed: the initial values are drawn from a gaussian distribution
with mean and standard deviation calculated from the observed values for each continuous variable. The second step of the (regularized) iterative MFA algorithm is to perform MFA on the completed dataset. Then, it imputes the missing values with the (regularized) reconstruction formulae of order ncp (the fitted matrix computed with ncp components for the (regularized) scores and loadings). These steps of estimation of the parameters via MFA and imputation of the missing values using the (regularized) fitted matrix are iterate until convergence. \cr

We advice to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. In the regularized algorithm, the singular values of the MFA are shrinked.\cr

The output of the algorithm can be used as an input of the MFA function of the FactoMineR package in order to perform the MFA on an incomplete dataset.}
\value{
  \item{tab.disj}{the imputed matrix; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones. The categorical variables are coded with the indicator matrix of dummy variables. In this indicator matrix, the imputed values are real numbers but they met the constraint that the sum of the entries corresponding to one individual and one variable is equal to one. Consequently they can be seen as degree of membership to the corresponding category}
  \item{completeObs}{the imputed dataset; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones. For the continuous variables, the values are the same as in the tab.disj output; for the categorical variables missing values are imputed with the most plausible categories according to the values in the tab.disj output}

  \item{call}{the matched call}
}

\references{F. Husson, J. Josse (2013) Handling missing values in multiple factor analysis. Food Quality and Preferences, 30 (2), 77-85.}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}


\seealso{\code{\link{imputePCA}}}
\examples{
\dontrun{
data(orange)
## Impute the data and perform a MFA
## with groups of continuous variables only
res.impute <- imputeMFA(orange, group=c(5,3), type=rep("s",2),ncp=2) 
res.mfa <- MFA(res.impute$completeObs,group=c(5,3),type=rep("s",2)) 
}
\dontrun{
data(vnf)
## Impute the indicator matrix and perform a MFA 
## with groups of categorical variables only
res.comp <- imputeMFA(vnf,group=c(6,5,3),type=c("n","n","n"),ncp=2)
res.mfa <- MFA(vnf,group=c(6,5,3),type=c("n","n","n"),tab.comp=res.comp)
}
}
\keyword{models}
\keyword{multivariate}
