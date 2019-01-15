\name{snorena}
\alias{snorena}
\docType{data}
\title{
Characterization of people who snore
}
\description{
This dataset contains 100 individuals and 7 variables (age, weight, size, alcohol, sex, snore, tobacco). There are 4 continuous variables and 3 categorical variables with 2 levels. Some values are missing.}
\usage{data(snorena)}

\format{
A data frame with 100 observations on 7 variables.
}

\source{
Cornillon, P.-A., Guyader, A., Husson, F., Jegou, N., Josse, J., Kloareg, M., Matzner-Lober, E., Rouviere, L., (2012). R for Statistics. Chapman & Hall/CRC Computer Science & Data Analysis, Rennes.
}

\examples{
\dontrun{
data(snorena)
res.comp <- imputeFAMD(snorena, ncp=3)
require(FactoMineR)
res.afdm <- FAMD(snorena, tab.comp = res.comp$tab.disj)
}
}
\keyword{datasets}
