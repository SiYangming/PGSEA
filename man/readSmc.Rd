\name{readSmc}
\alias{readSmc}
\alias{smc-class}

\title{ Read SMC files}
\description{
  This function reads in SMCs (simple molecular concepts) from individual text files.
}

\usage{
readSmc(files)
}

\arguments{
  \item{files}{ a character vector of file names }
}

\value{
	A list of SMC objects

}

\references{ ?? }

\author{ Kyle Furge <kyle.furge@vai.org> and Karl Dykema <karl.dykema@vai.org> }

\seealso{ \code{\link{writeSmc}} }

\examples{


datadir <- system.file("extdata", package = "PGSEA")
sample <- readGmt(file.path(datadir, "sample.gmt"))
str(sample)
sample[[1]]@reference <- "fileName"

\dontrun{
	writeSmc(sample[[1]])
	smc <- readSmc("fileName-0.txt")
}

}
\keyword{ data}

