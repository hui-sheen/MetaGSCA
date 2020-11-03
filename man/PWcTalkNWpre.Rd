\name{PWcTalkNWpre}
\alias{PWcTalkNWpre}
\title{
Cascaded steps in pathway crosstalk analyis prior to network drawing.
}
\description{
Identify crosstalk pathway pairs from pathway-wise statistical significance values (p-values) across multiple datasets.
}
\usage{
PWcTalkNWpre(input2PWcTalk,test=c('binary','binary')[1],pTh.dataset,pTh.pwPair,pTh.pw)
}
\arguments{
  \item{input2PWcTalk}{ Input CSV file name or a data frame object. The matrix within the file or the data frame must contain gene set analysis results (p-values) across multiple datasets, as well as a bootstrap.p column which represents the meta-analysis result. }
  \item{test}{ The test method used to quantify pathway similarity between two binary vectors. Default is binary.}
  \item{pTh.dataset}{ P-values in the pathway-by-dataset matrix less than pTh.dataset are converted to 1 and those greater than the threshold are converted to 0.  }
  \item{pTh.pwPair}{ Pathway pairs with similarity value (p-value or an analogy) less than pTh.pwPair are retained as edges of the graph of pathway crosstalk. }
  \item{pTh.pw}{ Pathways with meta-analysis p-value less than pTh.pw are retained as vertices of the graph of pathway crosstalk. }
}
\details{
PWcTalkNWpre first infers pathway pairwise similarity from a dichotomized pathway-by-dataset p-value matrix, then apply an edge filter and a node filter to obtain a network of discrete pathway connections.

Firstly, the algorithm converts the pathway-dataset p-value matrix to a binary matrix on account of pTh.dataset. Then, it quantifies similarity between all possible pathway pairs using either the asymmetric binary similarity (the binary method of R function dist) or Pearson's phi (enabled through R function crosstable_statistics from R package sjstats). If test is set as phi, the pairwise distance value is inverted to a similarity metric through 1-dist operation and converted to percentages in a high-to-low sorted list. Finally, a graph (network) is defined by retaining pathway connections with similarity measure lower than pTh.pwPair, provided that the involved pathways have their meta-analysis p-value less than pTh.pw. If pTh.pw is set to 1, the vertex filter is waived.
}
\value{
  a list object with two components:
  \item{g}{A graph object as defined in R package igraph.}
  \item{tkid}{The identifier of the current tkplot window. This should be supplied to function tk_coords once layout adjustment is completed.}
}
\seealso{
\code{\link{PWcTalkNW}}
}
\author{
Hui Yu \email{schena.yu@gmail.com}
}
\examples{
data(input2PWcTalk)
## One code block to execute pathway crosstalk analysis. More code lines, but enabling interactive layout tuning.
preNW <- PWcTalkNWpre(input2PWcTalk,test='binary',
  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01)
g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p)
##### PAUSE here: adjust the network layout on the pop-out window to reach a satisfaction #####
coords <- tk_coords(g_tkid$tkid)
g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p,layout=coords,pdfW=14,pdfH=10,figname='PWcTalk',asp=0.5) 
}




