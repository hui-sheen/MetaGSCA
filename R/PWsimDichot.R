#' Dichotomize p-values of PW pairs according to pair p threshold (say 0.01) [and vertex metaP threshold (say 0.01)]
#'
#' @param pwPairP output from function PWsim
#' @param metaP list component $metaP of function dichotGSAR output
#' @param pTh.pwPair threshold value applied on the overall meta-analysis p-value (p-statistic) 
#' @param pTh.pw threshold value applied on the p-value (p-statistics) of individual datasets
PWsimDichot <- function(pwPairP,metaP=NULL,pTh.pwPair=0.01,pTh.pw=0.01) {
  pwScope <- names(metaP)[metaP<=pTh.pw]
  dichotP0 <- pwPairP[,3]<=pTh.pwPair
  dichotP0[is.na(dichotP0)] <- 0
  cat(sum(dichotP0), 'pathway pairs survived pair threshold',round(100*sum(dichotP0)/length(dichotP0),1),'%\n')
  if (pTh.pw<1 & !is.null(metaP)) {
    dichotP <- dichotP0 & (pwPairP[,1]%in%pwScope) & (pwPairP[,2]%in%pwScope)
    cat(sum(dichotP),'pathway pairs survived further pathway threshold',round(100*sum(dichotP)/length(dichotP),1),'%\n')  
  } else {
    dichotP <- dichotP0
  }
  res <- data.frame(pwPairP,connect=dichotP)
  res
}
