\name{imputeCA}
\alias{imputeCA}

\title{Impute contingency table}
\description{Impute the missing entries of a contingency table using Correspondence Analysis (CA). Can be used as a preliminary step before performing CA on an incomplete dataset.}
\usage{
imputeCA(X, ncp = 2, threshold = 1e-08, maxiter = 1000, row.sup=NULL, 
     col.sup=NULL, quanti.sup=NULL, quali.sup=NULL)
}

\arguments{
  \item{X}{a data.frame that is a contingency table containing missing values}
  \item{ncp}{integer corresponding to the number of dimensions used  to predict the missing entries}
  \item{threshold}{the threshold for assessing convergence}
  \item{maxiter}{integer, maximum number of iterations for the regularized iterative CA algorithm}
  \item{row.sup}{a vector indicating the indexes of the supplementary rows}
  \item{col.sup}{a vector indicating the indexes of the supplementary columns}
  \item{quanti.sup}{a vector indicating the indexes of the quantitative supplementary variables}
  \item{quali.sup}{a vector indicating the indexes of the categorical supplementary variables}
}
\details{
Impute the missing entries of a contingency table using a regularized CA algorithm. 
The (regularized) iterative CA algorithm first consists in initializing missing values with random initial values. 
The second step of the (regularized) iterative CA algorithm consists in performing CA on the completed dataset. Then, it imputes the missing values with the (regularized) reconstruction formulae of order ncp (the fitted matrix computed with ncp components for the (regularized) scores and loadings). These steps of estimation of the parameters via CA and imputation of the missing values using the (regularized) fitted matrix are iterate until convergence. \cr
In this regularized algorithm, the singular values of the CA are shrinked. \cr
The number of components ncp used in the algorithm should be small. 
A small number of components can also be seen as a way to regularize more and consequently may be advices to get more stable predictions.
\cr
The output of the algorithm can be used as an input of the CA function of the FactoMineR package in order to perform CA on an incomplete dataset.
}

\value{
The imputed contingency table; the observed values are kept for the non-missing entries and the missing values are replaced by the predicted ones.
}

\author{
Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\examples{
\dontrun{
data(children)

## Impute the indicator matrix and perform a CA
res.impute <- imputeCA(children, ncp=2)
res.ca <- CA(res.impute) 
}
}
\keyword{models}
\keyword{multivariate}
\keyword{imputation}
