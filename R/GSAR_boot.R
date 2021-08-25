#' Perform differential coexpression analysis on original and bootstrap datasets
#' 
#' Differential coexpression analysis is conducted on one original dataset and multiple bootstrap datasets.
#'   
#' 1. Genes/rows in original data failing STD check are simply dropped. Genes/rows in bootstrap data failing STD check are simply dropped, we did not bother trying a second bootstrap for presence of certain minimal STD genes. 2. When original dataset does not contain sufficient genes for running GSNCA on one particular gene set, warning will be shown as GSET: Results NOT Captured!!! The output list simply does not contain component for that GSET.
#'
#' @param object gene expression matrix covering two groups. Row names are gene symbols.
#' @param gsets list of multiple gene sets.
#' @param group original groupping of samples, vector of 1's and 2's.
#' @param R number of bootstrap times.
#' @param nperm times of sample indix permutation, necessitated by GSNCA
#' @param minGsize considered gene set must have this minimum size after overlaying with gene expression matrix.
#' @param max.skip maximum number of repeated permutation/bootstrap times to avoid zero STD
#' @param min.sd a valid data matrix per group must have at least this much per-feature STD
#' @param cor.method correlation method
#' @return list of GSARboot_format() results for multiple gene sets. Component has multiple elements plus a dsetRow which cover all ultimate output 
GSAR_boot <- function(R,gsets,object,group,nperm=100,cor.method='pearson',max.skip=50,min.sd=0.001,minGsize=3) {
  call <- match.call()
  object1 <- object[, c(which(group == 1))]
  object2 <- object[, c(which(group == 2))]
  tmp <- check_sd(t(object1),t(object2),rownames(object),min.sd=min.sd)
  genes.removed <- tmp$genes.removed
  object <- object[tmp$genes.kept,]
  nGset <- length(gsets)
  perm.list <- vector('list',nperm)
  for (i in seq_len(nperm)) {
    perm.list[[i]] <- sample(ncol(object)) 
  }
  res <- gsnca_gsets(gsets,object,group,perm.list,max.skip=max.skip,min.sd=min.sd,minGsize=minGsize)
  D_obs.vec <- sapply(res,function(x) x$stat)
  pvalue.vec <- sapply(res,function(x) x$p)
  object1 <- object[, c(which(group == 1))] # Re-derive group 1/2 matrices after STD check and gene dropping.
  object2 <- object[, c(which(group == 2))]  
  n1 <- ncol(object1) 
	n2 <- ncol(object2)
  D_obs_boot.mat <- matrix(NA,nrow=nGset,ncol=R,dimnames=list(names(gsets),1:R))
  pvalue_boot.mat <- matrix(NA,nrow=nGset,ncol=R,dimnames=list(names(gsets),1:R))
  for (j in seq_len(R)) {
    btrp.idx1 <- sample.int(n1, n1, replace=TRUE) #Bootstrap, group 1
    btrp.idx2 <- sample.int(n2, n2, replace=TRUE) #Bootstrap, group 2
    object_boot1 <- object1[, btrp.idx1] 
    object_boot2 <- object2[, btrp.idx2] 
    tmp <- check_sd(t(object_boot1),t(object_boot2),rownames(object),min.sd=min.sd)
    object_boot <- cbind(object_boot1[tmp$genes.kept,],object_boot2[tmp$genes.kept,]) # Simply dropping STD-failed genes in bootstrap data matrix.    
    res.j <- gsnca_gsets(gsets,object_boot,c(rep(1,n1),rep(2,n2)),perm.list,max.skip=max.skip,min.sd=min.sd,minGsize=minGsize)
    D_obs_boot.mat[,j] <- sapply(res.j,function(x) x$stat)
    pvalue_boot.mat[,j] <- sapply(res.j,function(x) x$p)
  }
  cat('\n')
  # Format result for individual gset, like output of GSAR_boot():
  gset.lst <- vector('list',length=nGset)# list results with geneset-specific components
  names(gset.lst) <- names(gsets)
  for (i in seq_len(nGset)) {
    D_obs <- round(D_obs.vec[i], 6)
    pvalue <- round(pvalue.vec[i], 6)
    gset <- gsets[[i]]
    genes.not.found <- paste(gset[!(gset %in% rownames(object))], collapse = ",")
    if(sum(!(gset %in% rownames(object)) == 0)) genes.not.found <- "-"
    if (!is.na(pvalue)) { # Result is not NA due to overly small gene set size
      cat('\t',names(gsets)[i],': Results Captured\n')
      gset.lst[[i]] <- GSARboot_format(D_obs,pvalue,D_obs_boot.mat[i,],pvalue_boot.mat[i,],
        genes.not.found,genes.removed,nperm=nperm,R=R,call=call)
    } else {
      cat('\t','!!!',names(gsets)[i],': Results NOT Captured!!!\n')
    }
  }
  gset.lst <- gset.lst[!sapply(gset.lst,is.null)]
  return(gset.lst)
}
