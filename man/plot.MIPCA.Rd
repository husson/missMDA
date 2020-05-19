\name{plot.MIPCA}
\alias{plot.MIPCA}

\title{Plot the graphs for the Multiple Imputation in PCA}

\description{From the multiple imputed datasets, the function plots graphs for the individuals, variables and dimensions for the Principal Component Analysis (PCA)}

\usage{\method{plot}{MIPCA}(x, choice = "all", axes = c(1, 2), new.plot = TRUE, 
   main = NULL, level.conf = 0.95, graph.type = c("ggplot","classic"), \dots)
}

\arguments{
  \item{x}{an object of class MIPCA}
  \item{choice}{the graph(s) to plot. By default "all" the graphs are plotted. "ind.proc" the procrustean representation of the individuals, "dim" the representation of the dimensions of the PCA, "ind.supp" the projection of the individuals as supplementary individuals, "var" the projection of the variables as supplementary variables}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{main}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{level.conf}{confidence level used to construct the ellipses. By default, 0.95}
  \item{graph.type}{a character that gives the type of graph used: "ggplot" or "classic"}
  \item{\dots}{further arguments passed to or from other methods}
}

\details{
Plots the multiple imputed datasets obtained by the function MIPCA.
The idea is to represent the multiple imputed dataset on a reference configuration (the map obtained from the PCA on the incomplete dataset). Different ways are available to take into account and visualize the supplement variability due to missing values.
}

\value{
Four graphs can be drawn:
  \item{ind.supp}{The individuals of the imputed datasets are projected as supplementary individuals onto the reference PCA map; then confidence ellipses are drawn}
  \item{var}{The variables of the imputed datasets are projected as supplementary variables onto the reference PCA map}
  \item{ind.proc}{A PCA is performed on each imputed dataset and each configuration of scores is rotated onto the reference PCA map with procrustes rotation; then confidence ellipses are drawn}
  \item{dim}{The dimensions of each imputed dataset are projected as supplementary variables onto the dimensions of the reference PCA dimensions}
  \item{graph}{if the \code{graph.type=ggplot} an object with a list with the ggplot graphs that can be modified and plotted}
}
\references{Josse, J., Husson, F. (2010). Multiple Imputation in PCA}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{julie.josse@polytechnique.edu}
}

\seealso{\code{\link{MIPCA}},\code{\link{imputePCA}}}
\examples{
\dontrun{
data(orange)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2
resMI <- MIPCA(orange,ncp=2)
plot(resMI)
}
}
\keyword{dplot}
