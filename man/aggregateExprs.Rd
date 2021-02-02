\name{aggregateExprs}
\alias{aggregateExprs}
\title{Aggregate expression data}
\description{
  This function removes duplicates row names from an expression set, summarizing them with a function of the users choice.
  The "absMax" function located in package "reb" we have found to be useful.
}
\usage{
aggregateExprs(x, package = "hgu133plus2", using = "ENTREZID", FUN, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ expression data - matrix, eSet, or ExpressionSet }
  \item{package}{ annotation package of expression data }
  \item{using}{ format type that gene IDs are converted to  }
  \item{FUN}{ function by which to summarize duplicated values }
  \item{\dots}{ extra parameters passed on to FUN }
}
\value{
  A matrix of expression data with the rows aggregated to a unique format chosen by the user.
  The new identifiers of the returned matrix are those specified with the "using" argument.
  To see possible values, use the ls() command illustrated below in the examples. 
}
\author{ Kyle Furge <kyle.furge@vai.org> and Karl Dykema <karl.dykema@vai.org> }
\examples{
	if (require(hgu95av2.db) & require(annaffy)) {
		library(annaffy)
		data(aafExpr)
		class(exprs(aafExpr))
		exprs(aafExpr)[1:4, 1:4]

		#list possible values for the "using" argument
		ls(pos=which(search()=="package:hgu95av2.db"))

		convert <- aggregateExprs(exprs(aafExpr),"hgu95av2.db",FUN=mean,na.rm=TRUE)
		convert[1:4,1:4]
	}
}
\keyword{datasets}
