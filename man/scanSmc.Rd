\name{scanSmc}
\alias{scanSmc}

\title{ Scan through smc objects}
\description{
  This function scans through smc objects and returns those with specified attributes.
}
\usage{
scanSmc(smcList, scanSlot = "private", scanFor = "no")
}
\arguments{
  \item{smcList}{ list of "smc" objects }
  \item{scanSlot}{ character - which smc slot to investigate }
  \item{scanFor}{ character - what character string to look for }
}
\value{
	a list of "smc" objects with the desired attribute
}

\author{ Karl Dykema <karl.dykema@vai.org>}

\examples{

datadir <- system.file("extdata", package = "PGSEA")
sample <- readGmt(file.path(datadir, "sample.gmt"))
sample[1:2] <- editSmc(sample[1:2],"creator","Joe Smith")

scanned <- scanSmc(sample,"creator","Joe Smith")

str(scanned)


}
\keyword{ manip }

