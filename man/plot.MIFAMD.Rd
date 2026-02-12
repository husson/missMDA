\name{plot.MIFAMD}
\alias{plot.MIFAMD}

\title{Plot the graphs for the Multiple Imputation in FAMD}

\description{From the multiple imputed datasets, the function plots graphs for the individuals, categories and dimensions for the Factor Analysis of Mixed Data (FAMD)}

\usage{\method{plot}{MIFAMD}(x, choice = "all", axes = c(1, 2), new.plot = TRUE, 
   main = NULL, level.conf = 0.95, graph.type = c("ggplot","classic"), \dots)
}
					  
\arguments{
  \item{x}{an object of class MIFAMD}
  \item{choice}{the graph(s) to plot. By default "all" the graphs are plotted. "ind.proc" the procrustean representation of the individuals, "dim" the representation of the dimensions of the FAMD, "ind.supp" the projection of the individuals as supplementary individuals, "mod.supp" the projection of the categories, "var" the projection of the quantitative variables as supplementary variables}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{main}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{level.conf}{confidence level used to construct the ellipses. By default, 0.95}
  \item{graph.type}{a character that gives the type of graph used: "ggplot" or "classic"}
  \item{\dots}{further arguments passed to or from other methods}
}

\details{
Plots the multiple imputed datasets obtained by the function MIFAMD.
The idea is to represent the multiple imputed dataset on a reference configuration (the map obtained from the FAMD on the incomplete dataset). Different ways are available to take into account and visualize the supplement variability due to missing values.
}

\value{
Four graphs can be drawn:
  \item{ind.supp}{The individuals of the imputed datasets are projected as supplementary individuals onto the reference FAMD map; then confidence ellipses are drawn}
  \item{mod.supp}{The individuals of the imputed datasets are projected as supplementary individuals onto the reference FAMD map, but only categories are plotted; then confidence ellipses are drawn}
  \item{ind.proc}{A FAMD is performed on each imputed dataset and each configuration of scores is rotated onto the reference FAMD map with procrustes rotation; then confidence ellipses are drawn}
  \item{var}{The quantitative variables of the imputed datasets are projected as supplementary variables onto the reference FAMD map}
  \item{dim}{The dimensions of each imputed dataset are projected as supplementary variables onto the dimensions of the reference FAMD dimensions}
}

\references{
Audigier, V., Husson, F. & Josse, J. (2015). A principal components method to impute mixed data. Advances in Data Analysis and Classification, 10(1), 5-26. <doi:10.1007/s11634-014-0195-1>
}
\author{
Francois Husson \email{francois.husson@institut-agro.fr}
}

\seealso{\code{\link{MIFAMD}},\code{\link{imputeFAMD}}}
\examples{
\dontrun{
data(ozone)

## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpFAMD(ozone) ## Time-consuming, nb = 2


## Multiple Imputation
res.mi<-MIFAMD(ozone,ncp = 2,nboot=50)

## Plot the graphs
plot(res.mi)
}
}
\keyword{dplot}
