\name{MetaGSCA}
\alias{MetaGSCA}
\alias{MetaGSCA_Multi_Geneset}
\title{
Meta Gene Set Co-expression Analysis
}
\description{
\code{MetaGSCA} is used to perform meta-analysis of single geneset, help to aggregate gene-expression data from multiple studies. The meta-analysis is implemented with the R package \code{meta} and \code{metafor}.
\code{MetaGSCA_Multi_Geneset} function to perform multiple geneset meta-analysis.
The difference is: In \code{MetaGSCA}, \code{list.geneset} is a prespecified gene list, could only report one geneset analysis in the end;
In \code{MetaGSCA_Multi_Geneset}, \code{list.geneset} is several prespecified gene lists, could report several geneset analysis in the same time.
}

\usage{
MetaGSCA(list.geneset, list.dataset, list.group, name.geneset = "Geneset",
         name.dataset, nperm = 500, nboot = 200, method.Inverse = FALSE,
         method.GLMM = TRUE, fixed.effect = FALSE, random.effect = TRUE)

MetaGSCA_Multi_Geneset(list.geneset, list.dataset, list.group,
                       name.geneset = "Geneset", name.dataset, nperm = 500,
                       nboot = 200, method.Inverse = FALSE, method.GLMM = TRUE,
                       fixed.effect = FALSE, random.effect = TRUE)
}

\arguments{
\item{list.geneset}{a/several pre-specified gene list(s) from a gene set or pathway.}

\item{list.dataset}{a list of datasets, first column is gene name.}

\item{list.group}{a list of samples/patients subgroup or condition (e.g. (1,1,1,2,2,2)).}

\item{name.geneset}{the name of the gene set(s), used for output file name.}

\item{name.dataset}{a list of dataset names corresponding to list.dataset, used for forest plot.}

\item{nperm}{number of permutations used to estimate the null distribution of the test statistic. If not given, a default value 500 is used.}

\item{nboot}{number of bootstraps used to estimate the point and interval estimate. If not given, a default value 200 is used.}

\item{method.Inverse}{logical. Indicating whether "Inverse" method is to be used for pooling of studies.}

\item{method.GLMM}{logical. Indicating whether "GLMM" method is to be used for pooling of studies.}

\item{fixed.effect}{logical. Indicating whether a fixed effect meta-analysis should be conducted.}

\item{random.effect}{logical. Indicating whether a random effects meta-analysis should be conducted.}

}
\details{
In the meta-analysis of proportion method, we could perform fixed- and random- effects model. If we choose the fixed-effects model, we assume that the parameter of interest is identical across studies and any difference between the observed proportion is only due to sampling error. If we choose the random-effects model, we assume that the observed difference between the proportions cannot be entirely attributed to sampling error and other factors such as differences in study popula-tion, study designs, etc. To examine con-sistency of findings across studies, a statistical test of heterogeneity is reported in the meta-analysis. If heterogeneity is a concern, the random-effects model is recommended. In this case, each study estimates a different parameter, and the pooled estimate describes the mean proportion of the distribution. The variance parameter de-scribes the heterogeneity among the studies (when the variance is zero, this is equivalent to the fixed-effects model).

Logit transformation is applied to the binary outcome and regressed on the study variable in fixed effects or random effects logistic model. Define the logit(p) as effect size (ES), it can be estimated using two approaches: in-verse variance method and generalized linear mixed model (GLMM) with fixed intercept. The interval estimation is provided from these methods as well as from the bootstrap approach.
}

\author{
Yan Guo \email{YaGuo@salud.unm.edu}, Fei Ye \email{fei.ye@vumc.org}
}

\references{
Balduzzi, S., Rucker, G. and Schwarzer, G. How to perform a meta-analysis with R: a practical tutorial. Evid-Based Ment Heal 2019;22(4):153-160.

Viechtbauer, W. Conducting Meta-Analyses in R with the metafor Package. J Stat Softw 2010;36(3):1-48.
}

\examples{

###########################
## step 1, read in pathways
mygeneset <- list.files(system.file("extdata/pathway29", package = "MetaGSCA"),
                        pattern = "*.lst", full.names=TRUE)
genesets <- list(); genesetnames <- list()
for (i in 1:length(mygeneset)) {
    g <- fread(mygeneset[i], header = FALSE, data.table = FALSE)
    genesets[[i]] <- as.character(g[, 1])
    genesetnames[[i]] <- tools::file_path_sans_ext(mygeneset[i])
}
genesetnames <- gsub(paste(system.file("extdata/pathway29", package = "MetaGSCA"),
                     "/", sep = ""), "", genesetnames)

##########################
## step 2,read in datasets
mydataset <- c(paste(system.file("extdata", package = "MetaGSCA"),
                     "/KIRP.rds", sep=""),
               paste(system.file("extdata", package = "MetaGSCA"),
                     "/LIHC.rds", sep=""))
datasets <- list(); groups <- list(); datasetnames <- list()
for (j in 1:length(mydataset)) {
    d <- readRDS(file = mydataset[j])
    datasets[[j]] <- d[-1,]
    groups[[j]] <- as.character(d[1, -1])
    datasetnames[[j]] <- tools::file_path_sans_ext(mydataset[j])
}
datasetnames <- gsub(paste(system.file("extdata", package = "MetaGSCA"),
                     "/", sep = ""), "", datasetnames)

## set your workdir
# workdir <- setwd()

## step 3.1,call MetaGSCA function to perform Meta-analysis
# dir.create(file.path(workdir, "result_single_geneset"), showWarnings = FALSE)
# setwd(file.path(workdir, "result_single_geneset"))
MetaGSCA(list.geneset = genesets[[2]],
         list.dataset = datasets,
         list.group = groups,
         name.geneset = genesetnames[[2]],
         name.dataset = datasetnames,
         nperm = 100,
         nboot = 100,
         method.Inverse = FALSE,
         method.GLMM = TRUE,
         fixed.effect = FALSE,
         random.effect = TRUE)

## step 3.2,call MetaGSCA_Multi_Geneset function to perform Multiple geneset Meta-analysis
# dir.create(file.path(workdir, "result_multiple_geneset"), showWarnings = FALSE)
# setwd(file.path(workdir, "result_multiple_geneset"))
MetaGSCA_Multi_Geneset(list.geneset = genesets[2:3],
                       list.dataset = datasets,
                       list.group = groups,
                       name.geneset = genesetnames[2:3],
                       name.dataset = datasetnames,
                       nperm = 100,
                       nboot = 100,
                       method.Inverse = FALSE,
                       method.GLMM = TRUE,
                       fixed.effect = FALSE,
                       random.effect = TRUE)

}
