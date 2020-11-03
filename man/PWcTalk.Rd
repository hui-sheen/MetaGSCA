\name{PWcTalk}
\alias{PWcTalk}
\title{
Pathway crosstalk analysis on pathway-dataset p-value matrix
}
\description{
Perform pathway crosstalk analysis based on setwise statistical significance (p-values) across multiple datasets, and draw pathway crosstalk network.
}
\usage{
PWcTalk(input2PWcTalk,test=c('binary','phi')[1],pTh.dataset,pTh.pwPair,pTh.pw,
  figname,pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2)
}
\arguments{
  \item{input2PWcTalk}{ Input CSV file name or a data frame object. The matrix within the file or the data frame must contain gene set analysis results (p-values) across multiple datasets, as well as a bootstrap.p column which represents the meta-analysis result. }
  \item{test}{ The test method used to quantify pathway similarity between two binary vectors. Default is binary.}
  \item{pTh.dataset}{ P-values in the pathway-by-dataset matrix less than pTh.dataset are converted to 1 and those greater than the threshold are converted to 0.  }
  \item{pTh.pwPair}{ Pathway pairs with similarity value (p-value or an analogy) less than pTh.pwPair are retained as edges of the graph of pathway crosstalk. }
  \item{pTh.pw}{ Pathways with meta-analysis p-value less than pTh.pw are retained as vertices of the graph of pathway crosstalk. }
  \item{figname}{ Name of the PDF file to draw pathway crosstalk network. }
  \item{pdfW}{ Width of the PDF file in inch. It is an argument of R function pdf. }
  \item{pdfH}{ Height of the PDF file in inch. It is an argument of R function pdf. }
  \item{asp}{ Relative ratio between width and height of the graph to be plotted. It is an argument of R function plot.igraph. }
  \item{vbase}{ Base size of vertices in the network graph. vbase will be multiplied by a coefficient proportional to the meta-analysis statistical significance to determine the vertex size in the pathway crosstalk network. }
  \item{ebase}{ Base width of edges in the network graph. ebase will be multiplied by a coefficient proportional to pathway pair's statistical significance to determine the edge width in the pathway crosstalk network. }
  \item{vlbase}{ Base character expansion ratio (cex) for the vertex labels. vlbase will be multiplied by a coefficient proportional to vertex degree to determine the vertex label size in the pathway crosstalk network. }
  \item{power}{ Paramter applied on (converted) meta-analysis p-values, (converted) gene pair p-values, and vertex degree values to attenuate extreme disparity. The original value x changes to x^power before being multiplied to vbase, ebase, or vlbase. }
}
\details{
PWcTalk first infers pathway pairwise similarity from a dichotomized pathway-by-dataset p-value matrix, then apply an edge filter and a node filter to obtain a network of discrete pathway connections. The resultant pathway crosstalk network is shown in a pop-up tkplot window and also in a PDF file.

Firstly, the algorithm converts the pathway-dataset p-value matrix to a binary matrix on account of pTh.dataset. Then, it quantifies similarity between all possible pathway pairs using either the asymmetric binary similarity (the binary method of R function dist) or Pearson's phi (enabled through R function crosstable_statistics from R package sjstats). If test is set as phi, the pairwise distance value is inverted to a similarity metric through 1-dist operation and converted to percentages in a high-to-low sorted list. Finally, a graph (network) is defined by retaining pathway connections with similarity measure lower than pTh.pwPair, provided that the involved pathways have their meta-analysis p-value less than pTh.pw. If pTh.pw is set to 1, the vertex filter is waived.

Note that PWcTalk is a wrapper function combining two sequential functions (PWcTalkNWpre and PWcTalkNW). PWcTalk function allows less flexibility than executing PWcTalkNWpre and PWcTalkNW sequentially. Users are generally recommended to execute PWcTalkNWpre once and for all, then repeat the network drawing step (through function PWcTalkNW) multiple times to refine the layout/appearance of the network.

}
\value{
	a graph object defined in R package igraph.
}
\seealso{
\code{\link{PWcTalkNW}}, \code{\link{PWcTalkNWpre}}
}
\author{
Hui Yu \email{schena.yu@gmail.com}
}

\examples{
data(input2PWcTalk)
## One code line to execute pathway crosstalk analysis. Simple, but no flexibility for layout tuning. 
PWcTalk(input2PWcTalk,test='binary',
  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01,figname='PWcTalk',
  pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2)
}

