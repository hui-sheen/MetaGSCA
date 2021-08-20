#' Run GSNCA over multiple gene sets
#'
#' For each gene-set loop, genes are overlaid to expression matrix and too small gene sets are declined for GSNCA run.
#'
#' Due to too small gene set size (with consideration of intersection with expression data), certain gene sets have NA as p and stat results.
#'
#' @export
#' @seealso [gsnca_p()] for the GSNCA algorithm, which further calls on [gsnca_stat()] for coexpression distance statistics.
#' @examples
#' data(meta)
#' BRCA <- datasets[['BRCA']]
#' smpCode <- substr(colnames(BRCA),14,15)
#' grp1 <- which(smpCode=='01')
#' grp2 <- which(smpCode=='11')
#' object <- BRCA[1:100,c(grp1,grp2)]
#' group <- c(rep(1,length(grp1)),rep(2,length(grp1)))
#' perm.list <- vector('list',500)
#' for (i in seq_len(500)) {perm.list[[i]] <- sample(ncol(object))}
#' gsets <- split(rownames(object),rep(1:4,each=25))
#' res <- gsnca_gsets(gsets,object,group,perm.list)
#'
#' @param object gene expression matrix covering two groups. Row names are gene symbols.
#' @param gsets list for multiple gene sets.
#' @param group original groupping of samples, vector of 1's and 2's.
#' @param perm.list list of permutation specs. Each component gives permutated sample indices
#' @param minGsize considered gene set must have this minimum size after overlaying with gene expression matrix.
#' @param max.skip maximum number of repeated permutation/bootstrap times to avoid zero STD
#' @param min.sd a valid data matrix per group must have at least this much per-feature STD
#' @param cor.method correlation method
#' @return geneset-wise GSNCA results, each consisting of p-value and statistics out of GSNCA.
gsnca_gsets <- function(gsets,object,group,perm.list,cor.method='pearson',max.skip=50,min.sd=0.001,minGsize=3) {
  nPerm <- length(perm.list)
  nGset <- length(gsets)
  gsetNames <- names(gsets)
  gset.lst <- vector('list',nGset)
  names(gset.lst) <- gsetNames
  for (gsetName in gsetNames) {
    gset <- intersect(rownames(object),gsets[[gsetName]])
    if (length(gset)>=minGsize) {
      gem <- object[gset,]
      gset.lst[[gsetName]] <- gsnca_p(gem,group,perm.list,cor.method=cor.method,max.skip=max.skip,min.sd=min.sd)
    } else {
      gset.lst[[gsetName]] <- list(p=NA,stat=NA)
    }
  }  
  gset.lst
}
