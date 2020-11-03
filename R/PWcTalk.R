# PWcTalk(): the overall main function comprising a cascade of major sub-functions.
### however it is recommended to collapse the function and run the sub-functions in a script instead.  
PWcTalk <- function(input2PWcTalk,test='binary',
  pTh.dataset=0.01,pTh.pwPair=0.01,pTh.pw=0.01,figname='PWcTalk',
  pdfW=10,pdfH=10,asp=0.7,vbase=15,ebase=2,vlbase=1,power=1/2) 
{
  preNW <- PWcTalkNWpre(input2PWcTalk,test,pTh.dataset=pTh.dataset,pTh.pwPair=pTh.pwPair,pTh.pw=pTh.pw)
  g_tkid <- PWcTalkNW(preNW$PW.pair,preNW$PW.p)
  coords <- tk_coords(g_tkid$tkid)
  g <- PWcTalkNW(preNW$PW.pair,preNW$PW.p,layout=coords,pdfW=pdfW,pdfH=pdfH,figname=figname,asp=asp)
}