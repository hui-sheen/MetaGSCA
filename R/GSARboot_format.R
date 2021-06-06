#' Assemble GSNACA results on original dataset and bootstrap datasets
#'
#' Collate statistics and p-value of GSNCA on original dataset and counterpart statistics or median and CIs from bootstrap datasets
#'
#' Only one gene set is considered in the current function
#'
#' @param D_obs distance statistics from original dataset, by GSNCA
#' @param pvalue p-value from original dataset, by GSNCA
#' @param D_obs_boot distance statistics from multiple bootstrap datasets
#' @param pvalue_boot p-values from multiple bootstrap datasets
#' @param nperm number of permutations passed on to GSNCA
#' @param R number of bootstrapping
#' @param level confidence interval of this level. Defaults to 0.05 (95% CI)
#' @param call external ultimate command exertion statement for GSAR_boot(). A character string
#' @param genes.not.found string of excluded genes for absence in expression matrix 
#' @param genes.removed string of genes that were removed for STD disqualification
#' @return list of multiple objects of statistics value, where component dsetRow consists of 14 elements to appear as a row in output file
GSARboot_format <- function(D_obs,pvalue,D_obs_boot,pvalue_boot,genes.not.found='-',genes.removed='-',nperm=100,R=200,level=0.05,call='Command used') {
      D_obs_boot <- round(pvalue_boot, 6)
      pvalue_boot <- round(pvalue_boot, 6)
      ## 1, Test statistics/ Mean
      ts.mean <- round(mean(pvalue_boot), 6)
      ts.bias <- round((mean(pvalue_boot)-D_obs), 6)
      ts.std <- round(sd(pvalue_boot), 6)
      ts.norm.ci <- paste((1-level)*100, "% CI: [", round(mean(pvalue_boot)-qnorm(1-level/2)*sd(pvalue_boot), 6), ", ", round(mean(pvalue_boot)+qnorm(1-level/2)*sd(pvalue_boot), 6), "]", sep = "")
      ## 2, Test statistics/ Median
      ts.percentile <- round(quantile(pvalue_boot, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)), 6)
      ts.pt.ci <- paste((1-level)*100, "% CI: [", round(as.numeric(quantile(pvalue_boot,  probs = level/2)), 6), ", ", round(as.numeric(quantile(pvalue_boot,  probs = (1-level/2))), 6), "]",  sep = "")
      ## 3, P value/ Mean
      p.mean <- round(mean(pvalue_boot), 6)
      p.bias <- round((mean(pvalue_boot)-pvalue), 6)
      p.std <- round(sd(pvalue_boot), 6)
      p.norm.ci <- paste((1-level)*100, "% CI: [", round(mean(pvalue_boot)-qnorm(1-level/2)*sd(pvalue_boot), 6), ", ", round(mean(pvalue_boot)+qnorm(1-level/2)*sd(pvalue_boot), 6), "]", sep = "")
      ## 4, P value/ Meidan
      p.percentile <- round(quantile(pvalue_boot, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)), 6)
      p.pt.ci <- paste((1-level)*100, "% CI: [", round(as.numeric(quantile(pvalue_boot,  probs = level/2)), 6), ", ", round(as.numeric(quantile(pvalue_boot,  probs = (1-level/2))), 6), "]",  sep = "")
      ## dsetRow denotes a row to be rbinded for multiple datasets. 
      dsetRow <- data.frame(Original.TS = D_obs,
                                 Original.Pvalue = pvalue,
                                 Boot.TS.Mean = ts.mean,
                                 Boot.TS.Bias = ts.bias,
                                 Boot.TS.Std = ts.std,
                                 Boot.TS.Norm.CI = ts.norm.ci,
                                 Boot.TS.PT.CI = ts.pt.ci,
                                 Boot.Pvalue.Mean = p.mean,
                                 Boot.Pvalue.Bias = p.bias,
                                 Boot.Pvalue.Std = p.std,
                                 Boot.Pvalue.Norm.CI = p.norm.ci,
                                 Boot.Pvalue.PT.CI = p.pt.ci,
                                 Genes.not.found = genes.not.found,
                                 Genes.removed = genes.removed)
      output<-list(Original.TS = D_obs,
                   Original.Pvalue = pvalue,
                   Boot.TS = D_obs_boot,
                   Boot.Pvalue = pvalue_boot,
                   Boot.TS.Mean = ts.mean,
                   Boot.TS.Bias = ts.bias,
                   Boot.TS.Std = ts.std,
                   Boot.TS.Norm.CI = ts.norm.ci,
                   Boot.TS.Percentile = ts.percentile,
                   Boot.TS.PT.CI = ts.pt.ci,
                   Boot.Pvalue.Mean = p.mean,
                   Boot.Pvalue.Bias = p.bias,
                   Boot.Pvalue.Std = p.std,
                   Boot.Pvalue.Norm.CI = p.norm.ci,
                   Boot.Pvalue.Percentile = p.percentile,
                   Boot.Pvalue.PT.CI = p.pt.ci,
                   call = call,
                   nperm = nperm,
                   R = R,
                   dsetRow = dsetRow)  
      output
}
