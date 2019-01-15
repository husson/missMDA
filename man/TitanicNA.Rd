\name{TitanicNA}
\alias{TitanicNA}
\docType{data}
\title{
Categorical data set with missing values: Survival of passengers on the Titanic
}
\description{
This data set provides information on the fate of passengers on the fatal maiden voyage of the ocean liner Titanic, summarized according to economic status (class), sex, age and survival. Twenty percent of values are missing completely at random on each variable.
}
\usage{data(TitanicNA)}

\format{
A data frame with 2201 observations on the following 4 variables:

CLASS
0 = crew, 1 = first, 2 = second, 3 = third, denoting the economic status of the subject

AGE
1 = adult, 0 = child, denoting if the subject is an adult or a child

SEX
1 = male, 0 = female, denoting the sex of the subject

SURV
1 = yes, 0 = no, denoting if the passenger lived through the fatal maiden voyage of the ocean liner Titanic
}

\source{
British Board of Trade (1990), Report on the Loss of the Titanic (S.S.). British Board of Trade Inquiry Report (reprint). Gloucester, UK: Allan Sutton Publishing.
}

\examples{
data(TitanicNA)
}
\keyword{datasets}
