\name{go2smc}
\alias{go2smc}
\title{ Gene Ontology 2 "smc"}
\description{
  This function creates "smc" objects from the "GO" Bioconductor library.
}
\usage{
go2smc(min = 50, max = 200,organism="human")
}

\arguments{
  \item{min}{ numeric - minimum length of ids to be included }
  \item{max}{ numeric - maximum length of ids to be included  }
  \item{organism}{ character - organism }
}

\value{
a list of "smc" objects
}
\author{ Karl Dykema <karl.dykema@vai.org>}

\examples{
	if(require(GO.db)){
		mcs <- go2smc()[1:2]
		str(mcs)
	}
}
\keyword{ manip }

