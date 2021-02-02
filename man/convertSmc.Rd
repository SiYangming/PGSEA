\name{convertSmc}
\alias{convertSmc}

\title{ Convert Entrez ID based "smc" object }
\description{
  This function will convert the Entrez IDs of an smc object to the corresponding Entrez IDs from a different species.
  Data from the homologene project is downloaded and used within this function. 
}
\usage{
convertSmc(mcs, fromSpecies = "h", toSpecies = "r",hgX="./homologene.data")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcs}{ a list of "smc" objects}
  \item{fromSpecies}{ character - a single letter describing the species to convert from ie, h=human, r= rat, etc..}
  \item{toSpecies}{ character - a single letter describing the species to convert to ie, h=human, r= rat, etc.. }
  \item{hgX}{ character - file name of homologene data file }
}
\details{
This function will not work if you have not downloaded the homologene data file. Please use this command to do so:
download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")
}
\value{
  a list of converted "smc" objects
}
\author{ Karl Dykema <karl.dykema@vai.org>}
\examples{

\dontrun{
	download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")
	
	datadir <- system.file("extdata", package = "PGSEA")
	sample <- readGmt(file.path(datadir, "sample.gmt"))
	
	converted <- convertSmc(sample[1:2],"h","r")
	str(converted)
}

}
\keyword{ manip }
