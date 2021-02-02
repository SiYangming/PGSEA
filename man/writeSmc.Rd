\name{writeSmc}
\alias{writeSmc}
\title{ writeSmc }
\description{
  This function will write out SMC objects to individual text files
}
\usage{
writeSmc(x)
}

\arguments{
  \item{x}{ an object of class SMC  }
}
\details{
  The file name is determined by the reference slot of the SMC object.
}


\author{ Kyle Furge <kyle.furge@vai.org> and Karl Dykema <karl.dykema@vai.org>}

\seealso{ \code{\link{writeSmc}} }

\examples{

datadir <- system.file("extdata", package = "PGSEA")
sample <- readGmt(file.path(datadir, "sample.gmt"))
str(sample)
sample[[1]]@reference <- "fileName"

\dontrun{
	writeSmc(sample[[1]])
}



}
\keyword{ data }
