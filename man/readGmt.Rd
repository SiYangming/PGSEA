\name{readGmt}
\alias{readGmt}

\title{ readGmt }
\description{
  This function will read a "gmt" file into R, returning results as a list of SMC objects.
}
\usage{
readGmt(fname)
}
\arguments{
  \item{fname}{ File name of concepts in .gmt format }
}

\value{
	A list of SMC objects
}

\details{

The .gmt file format is a tab delimited file format used to store gene lists. These gene lists are stored row by row. The first column is the gene set name. The second column is a brief description, and every entry after that is a gene within that gene set. 

}

\references{ \url{http://www.broad.mit.edu/gsea/doc/data_formats.html#gmt } }

\author{ Karl Dykema <karl.dykema@vai.org> }

\seealso{ \code{\link{writeGmt}} }

\examples{

	datadir <- system.file("extdata", package = "PGSEA")
	sample <- readGmt(file.path(datadir, "sample.gmt"))
	str(sample)


}
\keyword{data}
