#' Calculate similarity measure for every pair of PWs
#'
#' @param dichotP MetaGSAR p-value matrix for n PWs across k datasets.
#' @param test Choice of similarity measure over two binary vectors. For binary, it uses 1-dist(...,method='binary')
#' @return A data frame for all possible pathway paris [n*(n-1)/2]. 1st & 2nd cols: the two vertices comprising an edge. 3rd col: edge connection strength. It is full connection still.
PWsim <- function(dichotP,test=c('phi','binary')[1]) {
  pws <- rownames(dichotP)
  nPW <- length(pws)
  nComb <- nPW*(nPW-1)/2
  edge <- matrix(nr=nComb,nc=2)
  edgeVal <- numeric(nComb)
  cnt <- 0
  for (i in seq_len(nPW-1)) {
    for (j in (i+1):nPW) {
      cnt <- cnt+1
      edge[cnt,] <- c(pws[i],pws[j])
      tab <- table(factor(dichotP[i,],levels=c(0,1)),factor(dichotP[j,],levels=c(0,1)))
      edgeVal[cnt] <- switch(tolower(test),
        phi=sjstats::crosstable_statistics(tab,statistics='phi')$p.value,
        binary=1-dist(rbind(dichotP[i,],dichotP[j,]),method='binary') # inverse distance
      )           
    }
  }
  if (test=='binary') edgeVal <- rank(-edgeVal,ties.method='average')/length(edgeVal)
  res <- data.frame(edge,edgeVal)
  colnames(res) <- c('PW1','PW2','p_pairSim')
  res
}
