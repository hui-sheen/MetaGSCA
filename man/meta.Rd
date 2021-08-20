\name{meta}
\alias{meta}
\alias{datasets}
\alias{genesets}
\alias{groups}
\docType{data}
\title{ Minimum example inputs to MetaGSCA function }
\description{
  Three objects dedicated to gene sets, datasets, and sample groupping, respectively.
}
\usage{data(meta)}
\format{
	Three R objects are loaded under the superset name meta.
\itemize{
  \item{genesets}{ A list of three components, each component comprising gene members of one pathway}
  \item{groups}{ A list of 11 components, each component indicating phenotype labels of one TCGA dataset}
  \item{datasets}{ A list of 11 case-vs-control cancer expression datasets, e.g., BRCA}
}
}
\details{
	*genesets* is a list of three gene sets (signaling pathways). One gene set or multiple gene sets from this list can be supplied as a primary input to MetaGSCA().

	*groups* stores the binary grouping information for all samples of all datasets. The names of 11 components correspond to the dataset names, which coincide to the root names of the dataset RData files.

	*datasets* is a list of 11 cancer-specific gene expression datasets. Each dataset has roughly half samples for cancer tissues and the other half for normal controls, with the designation specified in the companion data object (groups). Each component dataset is named with TCGA acronyms (e.g., BRCA and KIRP).
}
\author{ Yan Guo, Hui Yu }
\examples{
data(meta)
sapply(genesets,length)
sapply(groups,table)
names(datasets)
}
