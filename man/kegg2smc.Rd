\name{kegg2smc}
\alias{kegg2smc}

\title{ KEGG pathway to "smc"}
\description{
  This function creates "smc" objects from the "KEGG" Bioconductor library.
}
\usage{
 kegg2smc(min = 1, max = 284,organism="human")
}

\arguments{
  \item{min}{ numeric - minimum length of ids to be included }
  \item{max}{ numeric - maximum length of ids to be included  }
  \item{organism}{ character - organism }
}

\value{
a list of "smc" objects
}

\author{ Karl Dykema <karl.dykema@vai.org> and Richard Birnie <richard.birnie@pro-curetherapeutics.com>}

\examples{
	if(require(KEGG.db)){
		mcs <- kegg2smc(min=20,max=284)
		length(mcs)
		str(mcs[[1]])
	}
}
\keyword{ manip }

