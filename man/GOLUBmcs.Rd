\name{GOLUBmcs}
\alias{GOLUBmcs}
\docType{data}
\title{ Molecular Concepts prepared at VAI from data created by Golub et al.}
\description{
  386 molecular concepts generated at VAI. The data these concepts were generate from is available from \url{http://www.broad.mit.edu/cmap/}.
}
\usage{data(GOLUBmcs)}
\format{
	a list of "smc" objects
}

\details{
	These concepts were generated using the limma BioConductor package. The code used for generation of these concepts is available upon request.

}
\source{
 \url{http://www.broad.mit.edu/cmap/}
}
\examples{
data(GOLUBmcs)
str(GOLUBmcs[1:4]) 
}
\keyword{datasets}
