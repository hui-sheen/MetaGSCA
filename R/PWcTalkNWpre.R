#' Cascaded steps in pathway crosstalk analyis prior to network drawing
#'
#' Identify crosstalk pathway pairs from pathway-wise statistical significance values (p-values) across multiple datasets.
#' 
#' PWcTalkNWpre first infers pathway pairwise similarity from a dichotomized pathway-by-dataset p-value matrix, then apply an edge filter and a node filter to obtain a network of discrete pathway connections.
#' Firstly, the algorithm converts the pathway-dataset p-value matrix to a binary matrix on account of pTh.dataset. Then, it quantifies similarity between all possible pathway pairs using either the asymmetric binary similarity (the binary method of R function dist) or Pearson's phi (enabled through R function crosstable_statistics from R package sjstats). If test is set as phi, the pairwise distance value is inverted to a similarity metric through 1-dist operation and converted to percentages in a high-to-low sorted list. Finally, a graph (network) is defined by retaining pathway connections with similarity measure lower than pTh.pwPair, provided that the involved pathways have their meta-analysis p-value less than pTh.pw. If pTh.pw is set to 1, the vertex filter is waived.
#'
#' @export
#' @family aggregations
#' @seealso [PWcTalkNW()] for steps post this function, [PWcTalk()] for overall compacted pathway crosstalk analysis module.
#' @examples
#' data(input2PWcTalk)
#' # One code block to execute pathway crosstalk analysis, enabling interactive layout tuning.
#' preNW <- PWcTalkNWpre(input2PWcTalk,test='binary',
#'  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01)
#' #Code requires XMing support (x11 server), thus being turned off.
#' #g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p)
#' ## PAUSE here: adjust the network layout on the pop-out window to reach a satisfaction ##
#' #coords <- tk_coords(g_tkid$tkid)
#' #g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p,layout=coords,
#' #pdfW=14,pdfH=10,figname='PWcTalk',asp=0.5) 
#'
#' @param input2PWcTalk Input CSV file name or a data frame object. The matrix within the file or the data frame must contain gene set analysis results (p-values) across multiple datasets, as well as a bootstrap.p column which represents the meta-analysis result.
#' @param test The test method used to quantify pathway similarity between two binary vectors. Default is binary.
#' @param pTh.dataset P-values in the pathway-by-dataset matrix less than pTh.dataset are converted to 1 and those greater than the threshold are converted to 0.
#' @param pTh.pwPair Pathway pairs with similarity value (p-value or an analogy) less than pTh.pwPair are retained as edges of the graph of pathway crosstalk.
#' @param pTh.pw Pathways with meta-analysis p-value less than pTh.pw are retained as vertices of the graph of pathway crosstalk.

PWcTalkNWpre <- function(input2PWcTalk,test='binary',
  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01)
{
  dichotGSAR.res <- dichotGSAR(input2PWcTalk,pTh.dataset=pTh.dataset)
  PWsim.res <- PWsim(dichotGSAR.res$dichotP,test=test) #crossTab
  edges <- PWsimDichot(PWsim.res,dichotGSAR.res$metaP,pTh.pwPair=pTh.pwPair,pTh.pw=pTh.pw)
  PWpair <- subset(edges,connect==1)[,seq_len(3)]
  if (nrow(PWpair)>0) {
    pws <- unique(as.character(unlist(PWpair[,seq_len(2)]))) 
    PWp <- data.frame(pws,pwP=dichotGSAR.res$metaP[pws])
  } else {
    warning('Zero edge survives specified thresholds! Impossible to draw pathway crosstalk network!')
  }
  list(PW.pair=PWpair,PW.p=PWp)
}
