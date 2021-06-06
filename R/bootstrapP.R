#' Append final two bootstrap statistics
#'
#' meta.p & bootstrap.p are appended to series bootstrap statistics.
#'
#' Either Fixed Effects model or Random Effects model can be designated. Between the two options, result structures are identical but names are different.
#'
#' @param bootstrap R object flowing from prior content of MetaGSCA function
#' @param meta.res output from metaAndPlot function
#' @param fixed.effect logic variable for designation of Fixed Effect Model
#' @param random.effect logic variable for designation of Random Effect Model
#' @return Series of bootstrap result statistics with two new items appended: meta.p & bootstrap.p   

bootstrapP <- function(bootstrap,meta.res,fixed.effect=TRUE,random.effect=FALSE) {
	meta <- meta.res$meta
	meta.p <- meta.res$meta.p

	if(fixed.effect == TRUE){
      if (is.numeric(meta)) {
        bootstrap$meta.p.fixed = meta
        bootstrap$bootstrap.p.fixed = median(meta.p, na.rm = TRUE)
      } else {
  			fixed <- c(meta$TE.fixed, meta$lower.fixed, meta$upper.fixed)
  			bootstrap$meta.p.fixed = meta:::backtransf(fixed, sm="PLOGIT")[1]
  			bootstrap$bootstrap.p.fixed = median(meta.p, na.rm = TRUE)
      }
	}
	if(random.effect == TRUE){
    if (is.numeric(meta)) {
      bootstrap$meta.p.random = meta
      bootstrap$bootstrap.p.random = median(meta.p, na.rm = TRUE)
    } else {
			random <- c(meta$TE.random, meta$lower.random, meta$upper.random)
			bootstrap$meta.p.random = meta:::backtransf(random, sm="PLOGIT")[1]
			bootstrap$bootstrap.p.random = median(meta.p, na.rm = TRUE)
    }
  }
	bootstrap
}
