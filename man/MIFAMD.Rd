\name{MIFAMD}
\alias{MIFAMD}
\title{
Multiple Imputation with FAMD
}
\description{
MIFAMD performs multiple imputations for mixed data (continuous and categorical) using Factorial Analysis of Mixed Data.}
\usage{
MIFAMD(X, ncp = 2, method = c("Regularized", "EM"), coeff.ridge = 1, threshold = 1e-06, 
    seed = NULL, maxiter = 1000, nboot = 20, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data.frame with continuous AND categorical variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the FAMD reconstruction formulae}
  \item{method}{"Regularized" by default or "EM"}
  \item{coeff.ridge}{1 by default to perform the regularized imputeFAMD algorithm. Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of an EM method) or more than 1 to regularized more (to get closer to the results of the proportion imputation)}
  \item{threshold}{the threshold for the criterion convergence}
  \item{seed}{integer, by default seed = NULL implies that missing values are initially imputed by the mean of each variable for the continuous variables and by the proportion of the category for the categorical variables coded with indicator matrices of dummy variables. Other values leads to a random initialization}
  \item{maxiter}{integer, maximum number of iterations for the algorithm}
  \item{nboot}{the number of imputed datasets}
  \item{verbose}{use verbose=TRUE for screen printing of iteration numbers}
}
\details{MIFAMD generates nboot imputed data sets using FAMD. The observed values are the same from one dataset to the others, whereas the imputed values change. The algorithm is as follows: first, nboot weightings are defined for the individuals (equivalent to a non-parametric bootstrap). Then, the iterative regularized FAMD algorithm (Audigier et al., 2016) is applied according to each weighting, leading to nboot imputed tables. Dummy variables (coding for categorial variables) of these imputed tables are scaled to verify the constraint that the sum is equal to one per variable and per individual. Lastly, missing categories are drawn from the probabilities given by the imputed tables, and gaussian noise is added to the prediction of continuous variables. Thus, nboot imputed mixed data sets are obtained. The variation among the imputed values reflects the variability with which missing values can be predicted.
}
\value{
\item{res.MI}{A list of data frames corresponding to the nboot imputed mixed data sets}
\item{res.imputeFAMD}{A list corresponding to the output obtained with the function imputeFAMD (single imputation)}
\item{call}{The matched call}
}
\references{
Audigier, V., Husson, F. & Josse, J. (2015). A principal components method to impute mixed data. Advances in Data Analysis and Classification, 10(1), 5-26. <doi:10.1007/s11634-014-0195-1>

Audigier, V., Husson, F., Josse, J. (2017). MIMCA: Multiple imputation for categorical variables with multiple correspondence analysis. <doi:10.1007/s11222-016-9635-4>

Little R.J.A., Rubin D.B. (2002) Statistical Analysis with Missing Data. Wiley series in probability and statistics, New-York
}
\author{
Vincent Audigier \email{vincent.audigier@lecnam.net}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{imputeFAMD}},\code{\link{MIPCA}},\code{\link{MIMCA}},\code{\link{estim_ncpFAMD}},\code{\link[mice]{with.mids}},\code{\link[mice]{pool}},\code{\link[mice]{summary.mira}}}

\examples{
\dontrun{
data(ozone)

## First the number of components has to be chosen 
##   (for the reconstruction step)
 nb <- estim_ncpFAMD(ozone) ## Time-consuming, nb = 2


## Multiple Imputation
res.mi<-MIFAMD(ozone,ncp = 2,nboot=50)


## First completed data matrix
head(res.mi$res.MI[[1]])

## Analysis and pooling with mice
require(mice)
imp<-prelim(res.mi,ozone)
fit <- with(data=imp,exp=lm(maxO3~T9+T12+T15+Ne9+Ne12+Ne15+Vx9+Vx12+Vx15+maxO3v+vent+pluie))
res.pool<-pool(fit)
summary(res.pool)
}
}

\keyword{multivariate}
\keyword{imputation}
