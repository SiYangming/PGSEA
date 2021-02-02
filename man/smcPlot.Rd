\name{smcPlot}
\alias{smcPlot}
\title{ Plot PGSEA results}
\description{
  This basic function will plot results from PGSEA with easy altering of margins, colors, and text.
}
\usage{
smcPlot(m, ff = NULL, skip = "NO", scale = c(-3, 3), na.color = par("bg"), margins = NULL, r.cex = NULL, c.cex = NULL, show.grid = F, cnames = TRUE, rnames = TRUE, grid.lty = 3, clust = FALSE, ...)
}

\arguments{
  \item{m}{ matrix - your results from PGSEA (or any other numeric matrix of data)}
  \item{ff}{ factor - this factor corresponds to the subtypes of your samples and will control the column names}
  \item{skip}{ character - which subtype(s) to skip from "ff"}
  \item{scale}{ vector, length 2 - this vector sets the minimum and maximum values for the graph scale (at bottom of plot)}
  \item{na.color}{ character - color to display in the result of an NA}
  \item{margins}{ vector, length 4 - this vector gives the expansion values for the margins}
  \item{r.cex}{ numeric - number giving the amount by which row names should be scaled relative to the default}
  \item{c.cex}{ numeric - number giving the amount by which column names should be scaled relative to the default}
  \item{show.grid}{ boolean - show grid outlines within plot?}
  \item{cnames}{ boolean or character - vector of alternative column names}
  \item{rnames}{ boolean or character - vector of alternative row names}
  \item{grid.lty}{ numeric - line type of the grid lines}
  \item{clust}{ boolean - want to cluster?}

  \item{\dots}{ additional graphical parameters passed along to the plotting function }
}

\author{ Karl Dykema <karl.dykema@vai.org>}

\examples{

library(PGSEA)
datadir <- system.file("extdata", package = "PGSEA")
sample <- readGmt(file.path(datadir, "sample.gmt"))
data(nbEset)

pg <- PGSEA(nbEset,cl=sample,ref=1:5)
sub <- factor(c(rep(NA,5),rep("NeuroB",5),rep("NeuroB_MYC+",5)))

smcPlot(pg,sub,scale=c(-10,10),col=.rwb,margins=c(1,1,8,13))

}
\keyword{ hplot }
