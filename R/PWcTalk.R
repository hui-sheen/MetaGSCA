#' Streamlined functions for Pathway Crosstalk Network analysis
#'
#' Overall main function packed as a cascade of major sub-functions
#'
#' Nevertheless, it is recommended to collapse the function and run the sub-functions in a script instead.  
#'
#' @export
#' @family aggregations
#' @seealso [PWcTalkNWpre()] for the first module of prior steps, [PWcTalkNW()] for the second module of post steps.
#' @examples
#' data(input2PWcTalk)
#' #Execute pathway crosstalk analysis in one command. 
#' #Code requires XMing support (x11 server), thus being turned off.
#' #PWcTalk(input2PWcTalk,test='binary',
#' # pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01,figname='PWcTalk',
#' # pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2)
#' @param input2PWcTalk A data frame with 330 observations across 14 variables. Each observation corresponds to a pathway
#' @param test Similarity measure over two binary vectors, either asymmetric Binary similarity or Pearson's Phi index 
#' @param pTh.dataset Threshold for permutation p-statistic of a gene set co-expression in a dataset
#' @param pTh.pwPair Threshold for the similarity measure between two pathways
#' @param pTh.pw Threshold for the meta-analysis p-statistic of a pathway’s co-expression across multiple datasets
#' @param figname Stem file name for the resultant pathway crosstalk network figure
#' @param pdfW Width of PDF figure file 
#' @param pdfH Height of PDF figure file
#' @param asp Aspect ratio (between width and height). Asp parameter of function plot.igraph
#' @param vbase A base value regarding vertex size. Vertex size depends on vbase and a soft thresholding Power parameter 
#' @param ebase A base value regarding edge width. Edge width depends on ebase and a soft thresholding Power parameter
#' @param vlbase A base value regarding verbex label. Vertex label font depends on vlbase and a soft thresholding Power parameter
#' @param power A power paramter involved in soft-thresholding treatment. Power is implicated in vertex size, edge width, and vertex label font
#' @return An igraph graph (network) object 
PWcTalk <- function(input2PWcTalk,test='binary',
  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01,figname='PWcTalk',
  pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2) 
{
  preNW <- PWcTalkNWpre(input2PWcTalk,test,pTh.dataset=pTh.dataset,pTh.pwPair=pTh.pwPair,pTh.pw=pTh.pw)
  g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p)
  coords <- tk_coords(g_tkid$tkid)
  g <- PWcTalkNW(preNW$PW.pair,preNW$PW.p,layout=coords,pdfW=pdfW,pdfH=pdfH,figname=figname,asp=asp)
}
