\name{vnf}
\alias{vnf}
\docType{data}
\title{
Questionnaire done by 1232 individuals who answered 14 questions
}
\description{
A user satisfaction survey of pleasure craft operators on the ``Canal des Deux Mers'', located in South of France, was carried out by the public corporation ``Voies Navigables de France'' responsible for managing and developing the largest network of navigable waterways in Europe. Some values are missing.
}
\usage{data(vnf)}

\format{
A data frame with 1232 observations on 14 categorical variables.
}

\references{Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis, Journal of Classification, 29 (1), pp. 91-116.}

\examples{
\dontrun{
data(vnf)
tab.disj.impute <- imputeMCA(vnf, ncp=4)$tab.disj
require(FactoMineR)
res.mca <- MCA(vnf,tab.disj=tab.disj.impute)
}
}
\keyword{datasets}
