#' Gene Sets Net Correlation Analysis (GSNCA) including permutation
#'
#' @export
#' @seealso [gsnca.stat()] for the internal algorithm that returns the coexpression distance statistics.
#' @examples
#' data(BRCA)
#' smpCode <- substr(colnames(BRCA),14,15)
#' grp1 <- which(smpCode=='01')
#' grp2 <- which(smpCode=='11')
#' object <- BRCA[1:25,c(grp1,grp2)]
#' group <- c(rep(1,length(grp1)),rep(2,length(grp1)))
#' perm.list <- vector('list',500)
#' for (i in seq_len(500)) {perm.list[[i]] <- sample(ncol(object))}
#' gsnca.p.res <- gsnca.p(object,group,perm.list)
#'
#' @param object Gene expression matrix. Genes in rows.
#' @param group original groupping of samples, vector of 1's and 2's.
#' @param perm.list list of permutation specs. Each component gives permutated sample indices
#' @param max.skip maximum number of repeated permutation/bootstrap times to avoid zero STD
#' @param min.sd a valid data matrix per group must have at least this much per-feature STD
#' @param cor.method correlation method
#' @return list of p-value and d-statistics of GSNCA
gsnca.p <- function(object,group,perm.list,cor.method='pearson',max.skip=50,min.sd=0.001) {
    objt1 <- t(object[,group==1])
    objt2 <- t(object[,group==2])
    D_obs <- gsnca.stat(objt1,objt2,cor.method=cor.method)
    nv <- ncol(object)
    nperm <- length(perm.list)
    D_perm <- numeric(length=nperm)
    for (itr in seq_len(nperm)) {
			perm <- perm.list[[itr]]
			object <- object[,perm]
			objt1 <- t(object[,group==1])
			objt2 <- t(object[,group==2])
			tmp <- check_sd(objt1,objt2,rownames(object),min.sd=0.001)
			skip.counter=0        
			while ( (length(tmp$genes.kept)<ncol(objt1)) && (skip.counter <= max.skip) ) {
				if (skip.counter == max.skip){
					stop('number of skipped bootstraps or permutations exceeded max.skip')
				}
				skip.counter=skip.counter+1
				perm <- sample(ncol(object))
				objt1 <- t(object[,perm][,group==1])
				objt2 <- t(object[,perm][,group==2])
				tmp <- check_sd(objt1,objt2,rownames(object),min.sd=0.001)      
			}        
			D_perm[itr] <- gsnca.stat(objt1,objt2,cor.method=cor.method)
    }
    pvalue <- (sum(D_perm >= D_obs) + 1)/(length(D_perm) + 1)
    res <- list(p=pvalue,stat=D_obs)
}
