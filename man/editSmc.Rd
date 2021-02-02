\name{editSmc}
\alias{editSmc}

\title{ Edit "smc" objects}
\description{
  This function will edit a single or list of "smc" objects.
}
\usage{
editSmc(smcList, attName = "creator", newAtt = "changed!!")
}
\arguments{
  \item{smcList}{ a list of "smc" objects }
  \item{attName}{ character - which slot to change}
  \item{newAtt}{ character - what to change the slot to }
}
\value{
a list of edited "smc" objects 
}
\author{ Karl Dykema <karl.dykema@vai.org}

\examples{

	datadir <- system.file("extdata", package = "PGSEA")
	sample <- readGmt(file.path(datadir, "sample.gmt"))
	str(sample[1:2])
	
	temp <- editSmc(sample[1:2],"creator","Joe Smith")
	
	str(temp)


}
\keyword{ manip }


