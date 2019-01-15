\name{orange}
\alias{orange}
\docType{data}
\title{
Sensory description of 12 orange juices by 8 attributes.
}
\description{
Sensory description of 12 orange juices by 8 attributes. Some values are missing.
}

\usage{data(orange)}

\format{
A data frame with 12 rows and 8 columns. Rows represent the different orange juices,
columns represent the attributes.
}
\details{
A sensory data frame.
}
\source{
Francois Husson, Agrocampus Rennes
}

\examples{
\dontrun{
data(orange)
nb <- estim_ncpPCA(orange,ncp.min=0,ncp.max=5,method.cv="Kfold",nbsim=20,pNA=0.05)
res.comp <- imputePCA(orange,ncp=nb$ncp)
require(FactoMineR)
res.pca <- PCA(res.comp$completeObs)
resMI <- MIPCA(orange,ncp=nb$ncp)
plot(resMI)
}
}
\keyword{datasets}
