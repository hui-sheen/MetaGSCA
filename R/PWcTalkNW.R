#' Draw pathway crosstalk network given a finite set of pathway connections
#'
#' This function draws the pathway crosstalk network either in a pop-up tkplot window or on a disk file (PDF format).
#' 
#' PWcTalkNW relies heavily on R package igraph to enable plotting of pathway crosstalk network on an interactive GUI or on a disk PDF file. It is recommended to execute PWcTalkNW twice sequentially. At the first time, set argument layout as NULL, and a tkplot interactive window will pop up with a draft of the network. Users can adjust the layout of network to improve the appearance. If the layout has reached a level of satisfaction, users can call tk_coords to record the layout coordinates, and then execute PWcTalkNW one more time while setting layout to the recorded coordinates. In this way, the PDF file will use the same layout as adjusted by the user on the interactive tkplot window. However, depending on the specific setting of other graphing arguments (pdfW, pdfH, asp, vbase, ebase, vlbase, and power), the network in PDF file may deviate the layout tuned on tkplot window. If a severe distortion occurs, users need to adjust graphing arguments in a trial and error way and may need to execute PWcTalkNW function multiple times.
#' The function decides the graphing destination by setting of layout. When layout is set as NULL, it plots the graph on a tkplot window; when it is not NULL, it plots the graph on a PDF file. Edge information and node information will always be written as CSV files with the figname keyword in the file names.
#'
#' @export
#' @family aggregations
#' @seealso [PWcTalkNWpre()] for prior steps, [PWcTalk()] for overall compacted pathway crosstalk analysis module.
#' @examples
#' data(input2PWcTalk)
#' ## One code block to execute pathway crosstalk analysis, enabling interactive layout tuning.
#' preNW <- PWcTalkNWpre(input2PWcTalk,test='binary',
#' pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01)
#' g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p)
#' ##### PAUSE here: adjust the network layout on the pop-out window to reach a satisfaction #####
#' coords <- tk_coords(g_tkid$tkid)
#' g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p,layout=coords,
#' pdfW=14,pdfH=10,figname='PWcTalk',asp=0.5) 
#'
#' @param PWpair A data frame of three columns. First two columns define the pairs of pathways to be connected in the resultant network. Third column quantifies pathway similarity, either p-value (out of Pearson's phi) or percentage value (out of asymmetric binary distance).
#' @param PWp A data frame of two columns. First column has pathway names, and Second column contains meta-analysis p-values of individual pathways.
#' @param layout  Coordinate information of a tkplot object. Typically obtained via call of function tk_coords.
#' @param figname Name of the PDF file to draw pathway crosstalk network.
#' @param pdfW Width of the PDF file in inch. It is an argument of R function pdf.
#' @param pdfH Height of the PDF file in inch. It is an argument of R function pdf.
#' @param asp Relative ratio between width and height of the graph to be plotted. It is an argument of R function plot.igraph.
#' @param vbase Base size of vertices in the network graph. vbase will be multiplied by a coefficient proportional to the meta-analysis statistical significance to determine the vertex size in the pathway crosstalk network.
#' @param ebase Base width of edges in the network graph. ebase will be multiplied by a coefficient proportional to pathway pair's statistical significance to determine the edge width in the pathway crosstalk network.
#' @param vlbase Base character expansion ratio (cex) for the vertex labels. vlbase will be multiplied by a coefficient proportional to vertex degree to determine the vertex label size in the pathway crosstalk network.
#' @param power Paramter applied on (converted) meta-analysis p-values, (converted) gene pair p-values, and vertex degree values to attenuate extreme disparity. The original value x changes to x^power before being multiplied to vbase, ebase, or vlbase.
#' @return a list object with two components:
#' `g` A graph object defined in R package igraph.
#' `tkid` The identifier of the current tkplot window. This should be supplied to function tk_coords once layout adjustment is completed.   
PWcTalkNW <- function(PWpair,PWp,layout=NULL,figname='PWcTalk',
  pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2)
{
  write.csv(PWpair,paste0('PWpair.',figname,'.csv'),row.names=FALSE,quote=FALSE)
  write.csv(PWp,paste0('PWp.',figname,'.csv'),row.names=FALSE,quote=FALSE)
	if (max(PWpair[,3])<=1)
		PWpair[,3] <- -log10(PWpair[,3])
	colnames(PWpair)[3] <- 'minus.logP'
	PW.w <- -log10(unlist(PWp[,2]))
  names(PW.w) <- PWp[,1]
	g <- igraph::graph_from_data_frame(PWpair,directed=FALSE)
	V(g)$weight <- PW.w[V(g)$name]
	E(g)$weight <- PWpair$minus.logP
 	E(g)$color <- 'grey'
	if (is.null(layout)) {
    cat('Graph of PWcTalkNW has',vcount(g),'vertices and',ecount(g),'undirectional edges. See two CSV files saved on disk.\n')
  	cat('node weight (-log10(p)), min & max:\t',min(signif(PW.w,2)), '\t',max(signif(PW.w,2)),'\n')
  	cat('edge weight (-log10(p)), min & max:\t',min(signif(PWpair$minus.logP,2)),'\t',max(signif(PWpair$minus.logP,2)),'\n')
		tkid=igraph::tkplot(g,layout=layout_with_fr,vertex.color='lightgray',vertex.frame.color=NA,
      vertex.size=vbase*(V(g)$weight)^(power),vertex.label.cex=vlbase*degree(g)^(power),edge.width=ebase*(E(g)$weight)^(power))
 	  return(list(g=g,tkid=tkid))
	} else {
		pdf(paste(figname,'pdf',sep='.'),width=pdfW,height=pdfH)
    par(mar=par()$mar+c(0,10,0,10))
		igraph::plot.igraph(g,vertex.label=V(g)$name,asp=asp,layout=layout,
      vertex.size=vbase*(V(g)$weight)^(power),vertex.label.cex=vlbase*degree(g)^(power),
			vertex.label.color='black',vertex.color='lightgray',vertex.frame.color=NA,
			edge.color='lightgray',vertex.label.font=2,
			edge.width=ebase*(E(g)$weight)^(power))
		dev.off()
    return(list(g=g,tkid=NULL))
	}
}
