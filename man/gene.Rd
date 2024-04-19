\name{gene}
\alias{gene}
\docType{data}
\title{
Gene expression
}
\description{
A data frame with 53 brain tumors of 4 different types defined by the standard World Health 
Organization (WHO) classification (O, oligodendrogliomas; A, astrocytomas; OA, mixed oligo-astrocytomas 
and GBM, glioblastomas)  are described by information at the transcriptome level with expression data 
(356 continuous variables for microarrays) and at the genome level (76 continuous variables for CGH data) as illustrated.
10 rows are missing for the expression data.
}

\format{
A data-frame with 53 rows and 432 continuous variables}

\source{
de Tayrac M, Le S, Aubry M, Mosser J, Husson F. (2009). Simultaneous analysis of distinct Omics data sets with integration of biological knowledge: Multiple Factor Analysis approach. BMC Genomics, 10.}

\examples{
\dontrun{
data(gene)
res.impute <- imputeMFA(gene[,-1], group = c(76,356), 
    type = rep("s",2), ncp = 2) 
res.mfa <- MFA(cbind.data.frame(gene[,1], res.impute$completeObs), 
      group = c(1,76,356), type=c("n",rep("s",2)), 
	  name.group = c("WHO","CGH","expr"), num.group.sup = 1)
plot(res.mfa, habillage = 1, lab.ind = FALSE)
plot(res.mfa, habillage = "group", invisible = "ind", partial = "all")
plot(res.mfa, habillage = "group", lab.ind = FALSE, partial = "all")
plot(res.mfa, choix = "var", habillage = "group", lab.var = FALSE)
plot(res.mfa, choix = "group", habillage = "group")
}
}
\keyword{datasets}
