\name{ozone}
\alias{ozone}
\docType{data}
\title{
Daily measurements of meteorological variables and ozone concentration
}
\description{
This dataset contains 112 daily measurements of meteorological variables (wind speed, temperature, rainfall, etc.) and ozone concentration recorded in Rennes (France) during the summer 2001. There are 11 continuous variables and 2 categorical variables with 2 or 4 levels. Some values are missing.}
\usage{data(ozone)}

\format{
A data frame with 112 observations on 13 variables.
}

\source{
Cornillon, P.-A., Guyader, A., Husson, F., Jegou, N., Josse, J., Kloareg, M., Matzner-Lober, E., Rouviere, L., (2012). R for Statistics. Chapman & Hall/CRC Computer Science & Data Analysis, Rennes.
}

\examples{
\dontrun{
data(ozone)
res.comp <- imputeFAMD(ozone, ncp=3)
require(FactoMineR)
res.afdm <- FAMD(ozone,tab.comp=res.comp$tab.disj)
}
}
\keyword{datasets}
