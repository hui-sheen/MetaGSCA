\name{Individual TCGA transcriptome datasets}
\alias{BRCA}
\alias{COAD}
\alias{HNSC}
\alias{KIRC}
\alias{KIRP}
\alias{LIHC}
\alias{LUAD}
\alias{LUSC}
\alias{PRAD}
\alias{STAD}
\alias{THCA}
\docType{data}
\title{ RNA-Seq transcriptome data for paired tumor and adjacent samples}
\description{
  Eleven  cancer cohorts each had a minimum of 30 adjacent normal samples, and these normal samples were coupled with the corresponding tumor samples to yield 11 paired  expression matrices. 
}
\usage{data(BRCA)}
\format{
	A numeric data matrix with around 3800 gene expression vectors across tens to hundreds of variables. Rows correspond to genes and Columns correspond to samples (tumor or normal).
}
\details{
	Raw RNA-Seq data were downloaded from The Cancer Genome Atlas (TCGA). The minimum sample size is with KIRP and STAD, 64; the largest sample size is with BRCA, 226. The 11 data matrices share the same format, and almost all the covered genes (~3800). These 11 data matrices are provided as parallel datasets for the user to try out the meta-analysis offered by MetaGSCA, though the user can always select fewer (>1) datasets to try out the functions.  
}
\author{ Yan Guo, Hui Yu }
\examples{
data(BRCA)
head(BRCA[,1:5])
head(row.names(BRCA)) # gene names
data(KIRP)
dim(KIRP)
}
\keyword{datasets}
