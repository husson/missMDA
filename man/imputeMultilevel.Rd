\name{imputeMultilevel}
\alias{imputeMultilevel}
\title{Impute a multilevel mixed dataset}
\description{Impute the missing values of a multilevel mixed dataset (with a variable that groups the individuals, and with continuous and categorical variables) using the principal component method "multilevel factorial analysis for mixed data".}

\usage{imputeMultilevel(X, ifac = 1, ncpB = 2, ncpW=2, method=c("Regularized","EM"), 
    scale=TRUE, row.w = NULL, threshold = 1e-04, maxiter = 1000,...)
}

\arguments{
  \item{X}{a data.frame with continuous and categorical variables containing missing values}
  \item{ifac}{integer corresponding to the index of the group variable}
  \item{ncpB}{integer corresponding to the number of components used for the between group}
  \item{ncpW}{integer corresponding to the number of components used for the within group}
  \item{method}{"Regularized" by default or "EM"}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable. This is useful only when all the variables are continuous.}
  \item{row.w}{row weights (by default, uniform row weights)}
  \item{threshold}{the threshold for assessing convergence}
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Impute the missing entries of a multilevel mixed data using the iterative multilevel FAMD algorithm (method="EM") or the regularised iterative multilevel FAMD algorithm (method="Regularized"). 

We advice to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. In the regularized algorithm, the singular values of the FAMD are shrinked.\cr
}
\value{
  \item{completeObs}{the mixed imputed dataset; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones. For the continuous variables, the values are the same as in the tab.disj output; for the categorical variables missing values are imputed with the most plausible categories according to the values in the tab.disj output}
 }

\author{Francois Husson  \email{francois.husson@institut-agro.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\references{F. Husson, J. Josse, B. Narasimhan, G. Robin (2019). Imputation of mixed data with multilevel singular value decomposition. Journal of Computational and Graphical Statistics, 28 (3), pp. 552-566 <DOI:10.1080/10618600.2019.1585261>}

\seealso{\code{\link{imputePCA}},\code{\link{imputeFAMD}}}
\examples{
\dontrun{
## Example on artificial data
data(ozone)
res <- imputeMultilevel(ozone, ifac=12, ncpB=2, ncpW=2)
}
}
\keyword{models}
\keyword{multivariate}
\keyword{imputation}
