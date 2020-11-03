# dichotGSAR(): dichotomize p-value matrix to 0-1 matrix for each PW across datasets£»output metaP and metaFDR as well.
# INPUT input2PWcTalk: file name of the CSV output file from prior metaGSAR analysis.
### includes k p-values for each pathway represented in rows (k: number of repetitive datasets).
### Several columns of meta-P values follow the dataset-specific p-value columns.
dichotGSAR <- function(input2PWcTalk,pTh.dataset=0.01) { #,meta=c('glmm','inverse')[1]
  if (is.null(nrow(input2PWcTalk))) {# input is a file name, rather the data frame
    input2PWcTalk <- read.csv(input2PWcTalk,as.is=TRUE,row.names=1)
  }
  metaGSAR <- input2PWcTalk
  GSAR <- metaGSAR[,-seq_along(c(1:3))]
  dichotP <- GSAR<=pTh.dataset
  dichotP[is.na(dichotP)] <- 0
  metaP <- metaGSAR[,'bootstrap.p']
  names(metaP) <- rownames(metaGSAR)
  res <- list(dichotP=dichotP,metaP=metaP)
}
