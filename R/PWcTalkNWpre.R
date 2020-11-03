# PWcTalkNWpre(): Cascaded steps in pathway crosstalk analyis prior to network drawing.
PWcTalkNWpre <- function(input2PWcTalk,test='binary',
  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01) #,figname='PWcTalk',
  #pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2 
{
  if (!require('sjstats') | !require('igraph')) stop('Please install libraries sjstats and igraph before invoking PWcTalk() function.')
  dichotGSAR.res <- dichotGSAR(input2PWcTalk,pTh.dataset=pTh.dataset)
  PWsim.res <- PWsim(dichotGSAR.res$dichotP,test=test) #crossTab
  edges <- PWsimDichot(PWsim.res,dichotGSAR.res$metaP,pTh.pwPair=pTh.pwPair,pTh.pw=pTh.pw)
  PWpair <- subset(edges,connect==1)[,seq_len(3)]
  if (nrow(PWpair)>0) {
    pws <- unique(as.character(unlist(PWpair[,seq_len(2)]))) 
    PWp <- data.frame(pws,pwP=dichotGSAR.res$metaP[pws])
  } else {
    warning('Zero edge survives specified thresholds! Impossible to draw pathway crosstalk network!')
  }
  list(PW.pair=PWpair,PW.p=PWp)
}
