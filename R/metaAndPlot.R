#' Perform GLMM or inverse variance meta-analysis
#'
#' Perform GLMM or inverse variance meta-analysis on original and bootstrap analysis results
#'
#' 1.GLMM may reject to apply if all individual datasets have identical proportion, so there is always special treatment on GLMM track considering this marginal condition. For bootstrap sample meta-analysis, used tryCatch() to ignore occasional marginal situation on individual bootstrap time. 2. Forest plot is generated in a PDF file. 3. method designates GLMM or Inverse. If inappropriate, change method to method.choice throughout. 4. Use single generic term to designate either fixed or random series.
#'
#' @param lst list of GSAR_boot() results for multiple datasets regarding one gene set
#' @param name.geneset name of the concerned single gene set
#' @param method method for meta analysis, either GLMM or Inverse (inverse)
#' @param effect choice of effects model, either fixed or random
#' @param nperm times of sample indix permutation, necessitated by GSNCA
#' @param nboot number of bootstrap times
#' @return metaprop function result object together with p-values for multiple bootstrap meta-analyses.
metaAndPlot <- function(lst, ## list of GSAR_boot() results for multiple datasets.
                     name.geneset, ##  the names of a single gene sets, used for output file name
                     #names.dataset,  ## a list of dataset names corresponding to list.dataset, used for forest plot
                     nperm = 500,
                     nboot = 200,
										 method = c('GLMM','Inverse')[1],
										 effect = c('random','fixed')[1]
) {
	names.dataset <- names(lst)
	effect=tolower(effect)
	Event <- numeric(length(lst))
	for (j in seq_len(length(lst))){
			Event[j] <- ceiling(lst[[j]]$Original.Pvalue * (nperm + 1) - 1)
	}
	N <- rep(nperm, length(lst))	
  ## STEP 1, Meta Analysis for original p-value
  if (length(unique(Event))==1 & method=='GLMM') {
		meta <- Event[1]/nperm
		cat("\n---------------------- GLMM Warning--------------------\n")
		cat("If all input datasets have identical estimates of proportions, the GLMM approach cannot fit ML model.\n")
  } else {
    meta <- metaprop(
			event = Event,
			n = N,
			studlab = names.dataset,
			comb.fixed = ifelse(effect=='fixed',TRUE,FALSE),
			comb.random = ifelse(effect=='random',TRUE,FALSE),
			prediction = TRUE,
			method = method,
			method.tau = "ML"
    )
  }
  ## STEP 2, Meta Analysis for bootstrap p-value
  meta.p <- numeric(nboot)
  for (i in seq_len(nboot)) {
		event.boot <- numeric(length(lst))
		for(j in seq_len(length(lst))) {
				event.boot[j] <- ceiling(lst[[j]]$Boot.Pvalue[i] * (nperm + 1) - 1)
		}
		N.boot <- rep(nperm, length(lst))
		m.boot <- tryCatch(metaprop(
			event = event.boot,
			n = N.boot,
			studlab = names.dataset,
			comb.fixed = ifelse(effect=='fixed',TRUE,FALSE),
			comb.random = ifelse(effect=='random',TRUE,FALSE),
			prediction = TRUE,
			method = method,
			method.tau = "ML"),
			error = function(e) return(NA)
		)
		if(is.na(m.boot)){
			meta.p[i] <- NA
		}else{
			meta.p[i] <-ifelse(effect=='fixed', exp(m.boot$TE.fixed),exp(m.boot$TE.random))
		}
  }
	distinct.label=paste(method,' ',effect,' effect model bootstrapping result (nperm=', nperm, ', nboot=', nboot, ') ', sep='')
  if (is.numeric(meta) & method=='GLMM') {
		print(paste('GLMM bootstrap result (permutation time=', nperm, ', bootstrap time=', nboot, ') ',
								round(median(meta.p,na.rm=TRUE),4), ' [', round(quantile(meta.p, probs=0.025,na.rm=TRUE),4),
								', ', round(quantile(meta.p, probs=0.975,na.rm=TRUE),4), ']',
								sep=''))  
  } else {
		png(paste0(name.geneset, '_',method,'_',effect,'_forestPlot.png'), width=960, height=(3+length(lst)*0.25)*60)
    #pdf(file = paste0(name.geneset, '_',method,'_',effect,'_forestPlot.pdf'), width=10, height=3+length(lst)*0.25);
    MetaGSCA_plot(meta,meta.p,distinct.label)
    dev.off()
  }
  list(meta=meta,meta.p=meta.p)
}
