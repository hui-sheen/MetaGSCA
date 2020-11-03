\name{input2PWcTalk}
\alias{input2PWcTalk}
\docType{data}
\title{ Example input data to pathway crosstalk analysis }
\description{
  Data frame outputted from prior steps of multi-pathway meta-analysis in our pan-cancer case study
}
\usage{data(input2PWcTalk)}
\format{
A data frame with 330 observations across 14 variables. Each observation corresponds to a pathway.
\itemize{
  \item{Number.of.Genes}{ Number of Genes included in the pathway }
  \item{meta.p}{ Meta-analysis summarized permutation p-value across 11 datasets}
  \item{bootstrap.p}{ Median of Meta-analysis summarized bootstrap p-values across 11 datasets}
  \item{BRCA}{ permutation p-value in Breast invasive carcinoma}
  \item{COAD}{ permutation p-value in Colon adenocarcinoma}
  \item{HNSC}{ permutation p-value in Head and Neck squamous cell carcinoma}                                   
  \item{KIRC}{ permutation p-value in Kidney renal clear cell carcinoma}
  \item{KIRP}{ permutation p-value in Kidney renal papillary cell carcinoma}
  \item{LIHC}{ permutation p-value in Liver hepatocellular carcinoma}
  \item{LUAD}{ permutation p-value in Lung adenocarcinoma}
  \item{LUSC}{ permutation p-value in Lung squamous cell carcinoma}
  \item{PRAD}{ permutation p-value in Prostate adenocarcinoma}
  \item{STAD}{ permutation p-value in Stomach adenocarcinoma}
  \item{THCA}{ permutation p-value in Thyroid carcinoma}
}
}
\details{
We conducted a pan-cancer case study with 330 pathways containing at least 10 genes. Eleven cancer types from TCGA that had  a minimum of 30 normal samples were included in the meta-analysis. This dataset is a data frame outputted from function MetaGSCA_Multi_Geneset.
}
\author{ Yan Guo, Hui Yu }
\examples{
data(input2PWcTalk)
head(input2PWcTalk)
row.names(input2PWcTalk) # Names of 330 pathways
}
\keyword{datasets}
