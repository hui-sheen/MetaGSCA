#' Test statistics of Gene Sets Net Correlation Analysis (GSNCA) 
#'
#' Derive distance statistics by solving eigenvector of two pairise correaltion matrix.
#'
#' @export
#' @seealso [gsnca_p()] for the external function that calls on the current function and returns permutation p-statistic.
#' @examples
#' data(meta)
#' BRCA <- datasets[['BRCA']]
#' N <- ncol(BRCA)
#' n1 <- floor(N/2)
#' objt1 <- t(BRCA[1:100,1:n1])
#' objt2 <- t(BRCA[1:100,(n1+1):N])
#' gStat.res <- gsnca_stat(objt1,objt2)
#'
#' @param objt1 dataset for condition 1, genes in columns
#' @param objt2 dataset for condition 2, genes in columns
#' @param cor.method correlation coefficient method, the same as in function cor
#' @return L1-norm distance between two weight vectors
gsnca_stat <- function(objt1,objt2,cor.method='pearson') {
	cormat1 <- abs(cor(objt1, method = cor.method))
	cormat2 <- abs(cor(objt2, method = cor.method))
	e1 <- eigen(cormat1)        
	e2 <- eigen(cormat2)
	p1 <- abs(e1$vectors[, 1])
	p2 <- abs(e2$vectors[, 1])
	D_obs <- sum(abs(p1 * norm(matrix(p1)) - p2 * norm(matrix(p2))))
	D_obs
}
