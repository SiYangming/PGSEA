\name{PGSEA}
\alias{PGSEA}
\title{Parametric Gene Set Enrichment Analysis}
\description{
    This package contains functions for an exploratory parametric analysis of gene
    expression data. This type of analysis can assist in determining of
    lists of genes, such as those deregulated in defined experimental
    systems, are similarly deregulated in other data sets. 

    This function subsets the data based on lists of genes, computes a
    summary statistic for each gene list, and returns the results in a
    convenient form.
}
\usage{
PGSEA(exprs, cl, range = c(25, 500), ref = NULL, center = TRUE, p.value = 0.005, weighted = TRUE, enforceRange=TRUE, \dots)
}
\arguments{
  \item{exprs}{ matrix expression data, a numeric matrix, eSet, or ExpressionSet }
  \item{cl}{ gene set list - "GeneSetCollection" or list of "SMC" objects  }
  \item{range}{ a 2 element vector describing the min and max length of concepts to analyze  }
  \item{enforceRange}{ boolean - if TRUE, the expression matrix must contain data for the proper number of genes
  as set by the range argument to return a significant result. (this argument is used for data that contains NA's...) }
  \item{ref}{a vector containing the index of reference samples from
  which to make comparisons. Defaults to NULL (internally referenced samples) }
  \item{center}{ boolean - median center gene expression matrix columns
  prior to analysis. Can be helpful if 'ref' is used }
  \item{p.value}{ numeric p.value threshold or NA to return all data or TRUE to return a matrix of p.values }
  \item{weighted}{ boolean - weight results by the size of each gene list }
  \item{\dots}{ extra arguments passed along to FUN }
}
\details{
  Gene expression values are separated into subsets based on the lists
  of genes contained in the cl argument. This can be a "GeneSetCollection" or a list of "SMC" (Simple Molecular Concept) objects.
  For example, \code{readGmt} can be used to produce a 'smc'
  object list from a simple tab-delimited text file. The gene
  expression values from each of these gene lists is extracted and a
  summary statistic is computed for each subset (or region in the case of chromosomal bands/arms).

  The expression data must have the same identifiers as the list of
  genes being tested. If they are not, the expression data can be
  converted using the \code{aggregateExprs} function, that can use a
  current annotation environment to convert and condense the gene
  expression data.

  By default the method set out by Kim and Volsky \url{http://www.biomedcentral.com/1471-2105/6/144} is applied to the gene set.
  If weighted==FALSE than the default \code{t.test} function is used.

  The function is set up to perform the analysis on individual
  samples. For convenient method to analyze groups of samples, see the
  "Limma User's Guide" for more information on how to see up a contrast
  matrix and perform a linear model fit.  The coefficients of the fit
  can then be used a input into the \code{PGSEA} function.
  
  This package has not been extensively tested beyond a set of well defined curated pathways using the 
  Affymetrix platform and significance values represent approximations. Any results should be confirmed 
  by additional gene set testing methodologies.
  
}
\value{
	If p.value is set to a number, a matrix of results that pass at that significance is returned,
	of size <number of samples> x <number of molecular concepts>.
	
	If p.value is set to NA, all results are returned.
	
	If p.value is set to TRUE, then a list is returned that consists of the PGSEA results as well as their p.values.
}
\references{PGSEA: Parametric Analysis of Gene Set Enrichment}
\author{ Kim SY, Volsky DJ., kyle.furge@vai.org and karl.dykema@vai.org }
\note{\url{http://www.biomedcentral.com/1471-2105/6/144}}
\examples{
	datadir <- system.file("extdata", package = "PGSEA")
	sample <- readGmt(file.path(datadir, "sample.gmt"))
	data(nbEset)
	pg <- PGSEA(nbEset,cl=sample,ref=1:5)
	
	print(pg[,-c(1:5)])
}
\keyword{htest}
