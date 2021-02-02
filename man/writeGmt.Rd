\name{writeGmt}
\alias{writeGmt}
\title{ writeGmt }
\description{
  This function writes out SMC objects into .gmt file format
}
\usage{
writeGmt(fname, cl)
}

\arguments{
  \item{fname}{ name of the file to be written out }
  \item{cl}{ list of SMC objects }
}

\details{

The .gmt file format is a tab delimited file format used to store gene lists. These gene lists are stored row by row. The first column is the gene set name. The second column is a brief description, and every entry after that is a gene within that gene set. 

}

\references{ \url{http://www.broad.mit.edu/gsea/doc/data_formats.html#gmt }}

\author{ Kyle Furge <kyle.furge@vai.org> and Karl Dykema <karl.dykema@vai.org> }

\seealso{ \code{\link{readGmt}} }

\examples{

	datadir <- system.file("extdata", package = "PGSEA")
	sample <- readGmt(file.path(datadir, "sample.gmt"))
	str(sample)

\dontrun{
	writeGmt(paste(datadir,"/output.gmt",sep=""),sample)
}


}
\keyword{ data }
