#' Check gene-wise STD of each one matrix
#'
#' Check if data matrix per group has a minimum or higher per-gene STD
#'
#' @export
#' @examples
#' data(meta)
#' STAD <- datasets[['STAD']]
#' N <- ncol(STAD)
#' n1 <- N%/%2
#' objt1 <- t(STAD[,1:n1])
#' objt2 <- t(STAD[,(n1+1):N])
#' genes <- rownames(STAD)
#' check.res <- check_sd(objt1,objt2,genes,0.1)
#'
#' @param objt1 dataset for condition 1, genes in columns
#' @param objt2 dataset for condition 2, genes in columns
#' @param genes gene names or colnames of objt1/objt2
#' @param min.sd a valid data matrix per group must have at least this much per-feature STD
#' @return list of retained genes and to-be-removed genes for STD reason
check_sd <- function(objt1,objt2,genes,min.sd=0.001) {
	sd1 <- matrixStats::colSds(objt1,na.rm=TRUE) #sd1 <- apply(objt1, 2, 'sd', na.rm=TRUE)
	sd2 <- matrixStats::colSds(objt2,na.rm=TRUE)
	delcol <- (sd1 < min.sd) | (sd2 < min.sd)
	delcol[is.na(delcol)] <- TRUE
	genes.removed <- paste(genes[delcol], collapse = ',')
	if (sum(delcol,na.rm=T) == 0) genes.removed <- '-'
	if (sum(!delcol,na.rm=T) < 2) {
		print(delcol)
		stop(paste0('There must be at least 2 genes with standard deviation of gene expression data bigger than ', min.sd, ' in group 1 and group 2'))
	}
	genes.kept <- genes[!delcol]
	list(genes.kept=genes.kept,genes.removed=genes.removed)
}
