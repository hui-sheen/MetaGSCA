#' Meta-analysis of differential coexpression analysis 
#'
#' Meta-analysis of differential coexpression analysis across multiple datasets on multiple gene sets
#'
#' `MetaGSCA` systematically assesses the coexpression disturbance of a gene set by pooling the results from individual studies. In the kernel, a nonparametric approach named GSNCA tests whether a gene set is differentially coexpressed between two comparative conditions, and it produces a set-wise coexpression distance statistics as well as a permutation-based p-statistics. A meta-analysis is then performed to combine individual study results with either of two options: a random-intercept logistic regression model and the inverse variance method. The current function can take as input both a single gene set and multiple gene sets.
#' In the meta-analysis of proportion method, we could perform fixed- and random- effects model. If we choose the fixed-effects model, we assume that the parameter of interest is identical across studies and any difference between the observed proportion is only due to sampling error. If we choose the random-effects model, we assume that the observed difference between the proportions cannot be entirely attributed to sampling error and other factors such as differences in study popula-tion, study designs, etc. To examine con-sistency of findings across studies, a statistical test of heterogeneity is reported in the meta-analysis. If heterogeneity is a concern, the random-effects model is recommended. In this case, each study estimates a different parameter, and the pooled estimate describes the mean proportion of the distribution. The variance parameter de-scribes the heterogeneity among the studies (when the variance is zero, this is equivalent to the fixed-effects model).
#' Logit transformation is applied to the binary outcome and regressed on the study variable in fixed effects or random effects logistic model. Define the logit(p) as effect size (ES), it can be estimated using two approaches: in-verse variance method and generalized linear mixed model (GLMM) with fixed intercept. The interval estimation is provided from these methods as well as from the bootstrap approach.
#'
#' @export
#' @examples
#' data(meta)
#' data3 <- c('BRCA','COAD','HNSC')
#' testSingle <- MetaGSCA(list.geneset = genesets[2],
#'         list.dataset = datasets[data3], 
#'         list.group = groups[data3],
#'         names.geneset = names(genesets)[2],
#'         names.dataset = data3,
#'         nperm = 100,
#'         nboot = 100,
#'         method = 'GLMM',
#'         effect = 'random')
#'
#' @param list.geneset a list of gene sets (one or multiple).
#' @param list.dataset a list of datasets, first column is gene name.
#' @param list.group a list of samples/patients subgroup or condition (e.g. (1,1,1,2,2,2)).
#' @param names.geneset gene set names, corresponding to list.geneset
#' @param names.dataset dataset names, corresponding to list.dataset
#' @param nperm number of permutations used to estimate the null distribution of the test statistic. If not given, a default value 500 is used.
#' @param nboot number of bootstraps used to estimate the point and interval estimate. If not given, a default value 200 is used.
#' @param method meta-analysis method. Must be either `GLMM` or `Inverse`.
#' @param effect statistical model applied in meta-analysis. Must be either `random` or `fixed`.
MetaGSCA <- function(list.geneset,  ## a pre-specified gene list from a gene set or pathway
                     list.dataset,  ## a list of datasets, first column is gene name
                     list.group,  ## a list of samples/patients subgroup or condition (e.g. (1,1,1,2,2,2))
                     names.geneset=names(list.geneset), ##  the names of gene sets, used for output file name
                     names.dataset=names(list.dataset),  ## a list of dataset names corresponding to list.dataset, used for forest plot
                     nperm = 500,
                     nboot = 200,
                     method = c('GLMM','Inverse')[1],
                     effect = c('random','fixed')[1]
) {
	if (is.null(names.dataset)) names.dataset <- paste0('Dset',1:length(list.dataset))
	if (is.null(names.geneset)) names.geneset <- paste0('Gset',1:length(list.geneset))
	method.Inverse <- ifelse(method=='Inverse',TRUE,FALSE)
	method.GLMM <- !method.Inverse
	fixed.effect <- ifelse(effect=='fixed',TRUE,FALSE)
	random.effect <- !fixed.effect
	if (length(list.group) != length(list.dataset))
			stop("The length of the 'list.dataset' and the length of the 'list.group' must be equal")
	if (length(names.dataset) != length(list.dataset))
			stop("The length of the 'list.dataset' and the length of the 'names.dataset' must be equal")
	names(list.geneset) <- names.geneset
	names(list.dataset) <- names.dataset
	# GSAR_boot Analysis
	GSAR_boot.res <- list()
	res <- list()
	for (j in seq_len(length(list.dataset))){
		current.time <- Sys.time()
		cat(paste0('\n',current.time, " Processing dataset: ", names.dataset[j], "...\n"))
		GSAR_boot.res[[j]] <- GSAR_boot(R=nboot,gsets=list.geneset,object=list.dataset[[j]],group=list.group[[j]],nperm=nperm)
		if (j==1) {
			gsets <- names(GSAR_boot.res[[j]])
		} else {
			gsets <- intersect(gsets,names(GSAR_boot.res[[j]]))
		}
	}
	list.bootstrap <- vector('list',length(gsets))
	names(list.bootstrap) <- gsets
	for (gset in gsets) { # Intersect genesets where all datasets have results.
		dsets.res <- lapply(GSAR_boot.res,function(x) x[[gset]])
		dsetMat <- do.call(rbind,lapply(dsets.res,function(x) x$dsetRow))
		rownames(dsetMat) <- names.dataset
		write.csv(dsetMat, paste0(gset, "_Individual Dataset GSAR Results with Bootstrapping.csv"))
		meta.res=metaAndPlot(dsets.res,gset,nperm,nboot,
									 method,effect)
		bootstrap <- list()
		bootstrap$Geneset.Name <- gset
		bootstrap$Number.of.Genes <- length(list.geneset[[gset]])
		bootstrap <- bootstrapP(bootstrap,meta.res,fixed.effect,random.effect)
		bootstrap[(length(bootstrap)+1):(length(list.dataset)+length(bootstrap))] <- dsetMat$Original.Pvalue
		names(bootstrap)[(length(bootstrap)-length(list.dataset)+1):length(bootstrap)] <- names.dataset
		list.bootstrap[[gset]] <- bootstrap
	}
	res <- do.call('rbind',list.bootstrap)
	if(!is.null(res))
		write.table(res, "_Meta Analysis bootstrap result.tsv", row.names = FALSE,sep='\t',quote=F)    
}
