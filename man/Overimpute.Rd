\name{Overimpute}
\alias{Overimpute}

\title{Overimputation diagnostic plot}
\description{
Assess the fit of the predictive distribution after performing multiple imputation with the function \code{MIPCA}.
}
\usage{
Overimpute(output, plotvars)

}


\arguments{
  \item{output}{output from the function \code{MIPCA}.}
  \item{plotvars}{column number of the variable to overimpute.}
}

\details{This function imputes each observed values from each of the parameters of the imputation model obtained from the \code{MIPCA} procedure. The comparison between the ``overimputed" values and the observed values is made by constructing a confidence interval for each observed value using the quantiles of the overimputed values (Blackwell et al. (2015)). Note that confidence intervals constructed with quantiles require a large number of imputations. If the model fits well the data, then the 90\% confidence interval should contain the observed value in 90\% of the cases. The function \code{Overimpute} takes as an input the output of the MIPCA function (\code{output}) and the indices of the variables that are plotted (\code{plotvars}).
}

\value{A list of 6-column matrix that contains (1) the row in the original data,
(2) the observed value of that observation, (3) the mean of the
overimputations, (4) the lower bound of the 90\% confidence interval of
the overimputations, (5) the upper bound of the 90\% confidence interval
of the overimputations, and (6) the proportion of the other variables that were
missing for that observation in the original data.}

\references{
Blackwell, M., Honaker, J. and King. G. 2015. A Unified Approach to Measurement Error and Missing Data: Overview and Applications. Sociological Methods and Research, 1-39.
}

\seealso{
  \code{\link{MIPCA}}}
  
\examples{
\dontrun{
require(Zelig)
data(ozone)

# First the number of components has to be chosen 
nb <- estim_ncpPCA(ozone[,1:11])

# Multiple Imputation with Bayesian method
res.BayesMIPCA<-MIPCA(ozone[,1:11],ncp=2,method.mi="Bayes",verbose=T)

# Regression on the multiply imputed data set and pooling
z.out <- zelig(maxO3~., model = "ls", data = res.BayesMIPCA$res.MI,cite=F)
summary(z.out,digits=5)

# Diagnostics
res.over<-Overimpute(res.BayesMIPCA)

}
}

