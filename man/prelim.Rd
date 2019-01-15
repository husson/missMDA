\name{prelim}
\alias{prelim}
\title{Converts a dataset imputed by MIMCA, MIPCA or MIFAMD into a mids object}
\description{
This function performs grouping and sorting operations on a multiply imputed dataset. It creates a mids object that is needed for input to with.mids, which allows analyse of the multiply imputed data set. The original incomplete data set needs to be available so that we know where the missing data are.}
\usage{
prelim(res.mi,X)
}
\arguments{
  \item{res.mi}{an output of the functions MIPCA, MIMCA or MIFAMD}
  \item{X}{the original incomplete data set corresponding to the res.mi argument}
}

\value{
\item{imp.mids}{An object of type mids}
}

\author{Vincent Audigier \email{vincent.audigier@lecnam.net}, Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{josse@agrocampus-ouest.fr}
}

\seealso{\code{\link{MIFAMD}},\code{\link{MIPCA}},\code{\link{MIMCA}},\code{\link{with.mids}},\code{\link{pool}},\code{\link{summary.mira}}}
\examples{
\dontrun{
data(TitanicNA)

## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(TitanicNA,ncp.max=5) ## Time-consuming, nb = 5

## Multiple Imputation
res.mi <- MIMCA(TitanicNA, ncp=5, verbose=T)

#Analysis
imp<-prelim(res.mi,TitanicNA)
fit <- with(data=imp,exp=glm(SURV~CLASS+AGE+SEX,family = "binomial"))

#Pooling
res.pool<-pool(fit)
summary(res.pool)
}
}

\keyword{multivariate,imputation,categorical,nominal}
