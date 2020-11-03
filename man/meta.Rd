\name{meta}
\alias{meta}
\alias{Dsetnames}
\alias{genesets}
\alias{groups}
\docType{data}
\title{ Auxiliary data objects required by the tutorial }
\description{
  Three objects dedicated for gene sets, sample groupping, and data set names, respectively.
}
\usage{data(meta)}
\format{
	Three R objects are loaded under the superset name meta.
\itemize{
  \item{genesets}{ A list of three components, each component comprising gene members of one pathway}
  \item{groups}{ A list of 11 components, each component indicating phenotype labels of one TCGA dataset}
  \item{Dsetnames}{ Eleven dataset names, e.g., BRCA}
}
}
\details{
	*genesets* includes the gene lists for three signaling pathways. One gene set can be supplied as a primary input to MetaGSCA(), and two or three gene sets can be supplied altogether as a primary input to MetaGSCA_Multi_Geneset().

	*groups* stores the binary grouping information for all samples of all datasets. The names of 11 components correspond to the dataset names, which coincide to the root names of the dataset RData files.

	*Dsetnames* is a character vector of 11 items. It contains the individual dataset names, or the cancer names. Examples include BRCA and KIRP.
}
\author{ Yan Guo, Hui Yu }
\examples{
data(meta)
sapply(genesets,length)
sapply(groups,table)
Dsetnames
}
