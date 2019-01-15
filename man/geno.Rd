\name{geno}
\alias{geno}
\docType{data}
\title{
Genotype-environment data set with missing values
}
\description{
A data-frame with 16 rows corresponding to genotypes (triticale lines) and 10 columns corresponding to different environments where the genotypes were sown. 
Each cell of the data-frame corresponds to the grain yield (kilograms per hectare) for one genotype in an environment.  
The first six genotypes correspond to the so-called ``complete'' type, while the next eight were of the ``substituted'' type and two 
check genotypes were included. Such data sets are often incomplete. 
Indeed, it frequently happens that all the varieties are not assessed on all the environments.  
}

\format{
A data-frame with 16 rows and 10 columns}

\source{
Royo C, Rodriguez A, Romagosa I (1993). Differential adaptation of complete and substitute
triticale. Plant Breeding, 111, 113- 119.}

\examples{
\dontrun{
data(geno)
ncomp <- estim_ncpPCA(geno)
res.imp <- imputePCA(geno, ncp= ncomp$ncp)
res.pca <- PCA(res.imp$completeObs)
}
}
\keyword{datasets}
