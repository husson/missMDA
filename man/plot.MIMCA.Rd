\name{plot.MIMCA}
\alias{plot.MIMCA}

\title{Plot the graphs for the Multiple Imputation in MCA}

\description{From the multiple imputed datasets, the function plots graphs for the individuals, categories and dimensions for the Multiple Correspondance Analysis (MCA)}

\usage{\method{plot}{MIMCA}(x, choice = "all", axes = c(1, 2), new.plot = TRUE, 
   main = NULL, level.conf = 0.95, \dots)
}

\arguments{
  \item{x}{an object of class MIMCA}
  \item{choice}{the graph(s) to plot. By default "all" the graphs are plotted. "ind.proc" the procrustean representation of the individuals, "dim" the representation of the dimensions of the MCA, "ind.supp" the projection of the individuals as supplementary individuals, "mod.supp" the projection of the categories}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{main}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{level.conf}{confidence level used to construct the ellipses. By default, 0.95}
  \item{\dots}{further arguments passed to or from other methods}
}

\details{
Plots the multiple imputed datasets obtained by the function MIMCA.
The idea is to represent the multiple imputed dataset on a reference configuration (the map obtained from the MCA on the incomplete dataset). Different ways are available to take into account and visualize the supplement variability due to missing values.
}

\value{
Four graphs can be drawn:
  \item{ind.supp}{The individuals of the imputed datasets are projected as supplementary individuals onto the reference MCA map; then confidence ellipses are drawn}
  \item{mod.supp}{The individuals of the imputed datasets are projected as supplementary individuals onto the reference MCA map, but only categories are plotted; then confidence ellipses are drawn}
  \item{ind.proc}{A PCA is performed on each imputed dataset and each configuration of scores is rotated onto the reference MCA map with procrustes rotation; then confidence ellipses are drawn}
  \item{dim}{The dimensions of each imputed dataset are projected as supplementary variables onto the dimensions of the reference MCA dimensions}
}
\references{Audigier, V., Husson, F., Josse, J. (2016). MIMCA: Multiple imputation for categorical variables with multiple correspondence analysis}

\author{Audigier Vincent \email{vincent.audigier@cnam.fr}, Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\seealso{\code{\link{MIMCA}},\code{\link{imputeMCA}}}
\examples{
\dontrun{
data(TitanicNA)

## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(TitanicNA) ## Time-consuming, nb = 5

## Multiple Imputation
res.mi <- MIMCA(TitanicNA, ncp=5, verbose=TRUE)

## Plot the graphs
plot(res.mi)
}
}
\keyword{dplot}
