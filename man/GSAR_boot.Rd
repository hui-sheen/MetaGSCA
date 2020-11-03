
\name{GSAR_boot}
\alias{GSAR_boot}
\title{Gene set analysis in R with bootstrap}

\description{
\code{GSAR_boot} is used to perform correlation analysis for gene set. It consists of two-sample nonparametric statistical methods, to test correlation structure and whether the gene set is differentially co-expressed between two conditions. To provide reproducible conclusions, we include bootstrap method and provide corresponding point and interval estimates.
}

\usage{
GSAR_boot(projectname, gematrix, group, genelist, nperm = 500,
          cor.method = "pearson", check.sd = TRUE, min.sd = 0.001,
          max.skip = 50, R = 200, level = 0.05)
}

\arguments{
\item{projectname}{a list of dataset/project names, which could be used for meta analysis forest plot.}

\item{gematrix}{a list of dataset/project matrix, contains gene set informations.}

\item{group}{patient subgroup or condition (values: 1 or 2).}

\item{genelist}{a pre-specified gene list or signature or defined by pathway.}

\item{nperm}{number of permutations used to estimate the null distribution of the test statistic. If not given, a default value 500 is used.}

\item{cor.method}{a character string indicating which correlation coefficient is to be computed. Possible values are "pearson" (default), "spearman" and "kendall".}

\item{check.sd}{logical. Should the standard deviations of features checked for small values before the intergene correlations are computed? Default is TRUE (recommended).}

\item{min.sd}{the minimum allowed standard deviation for any feature. If any feature has a standard deviation smaller than min.sd the execution stops and an error message is returned.}

\item{max.skip}{maximum number of skipped random permutations which yield any feature with a standard deviation less than min.sd.}

\item{R}{number of bootstraps used to estimate the point and interval estimate. If not given, a default value 200 is used.}

\item{level}{the significance level, also denoted as alpha or ??, is the probability of rejecting the null hypothesis when it is true.}
}

\details{
\code{GSNCA} is a nonparametric test that assesses multivariate changes in the gene coexpression network between two conditions. Gene sets should be defined a priori, often derived from preliminary analysis on the basis of biological relevance and/or statistical significance, a pathway of interest, or a signature reported in the literature.

Suppose we have \code{N} gene expression studies containing a given gene set with \code{p} genes, and within each study samples can be divided into two biological conditions. Denote the weight of a gene as \code{wi}, and let the weight vector \code{w} be the first eigenvector of p×p gene correlation matrix. The test statistic is defined as the L1 norm between the 2 weight vectors. The p-values for the test statistics are obtained by comparing the observed value of the test statistic to its null distribution, which is constructed using permutation. This permutation test p-value can be interpreted as the proportion of test statistic values from the permuted data that are at least as extreme as the test statistic from the original data. This estimated proportion is a point estimate without an accompanying measure of precision.

Then we use a bootstrap approach to supplement point estimates by confidence intervals. \code{Bootstrap} is a resampling procedure that constructs a sampling distribution of the observed data, which allows making inference about sample esti-mates. In this work, we used the bootstrap method by repeatedly drawing random samples from original dataset with replacement to construct confidence intervals of the sample estimates. With each bootstrap sample, we calculate a permutation p-value.
}

\author{
Yan Guo \email{YaGuo@salud.unm.edu}, Fei Ye \email{fei.ye@vumc.org}
}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net correlations analysis (GSNCA): a multivariate differential coexpression test for gene sets. Bioinformatics 30, 360-368.

Walters, S.J. and Campbell, M.J. The use of bootstrap methods for estimating sample size and analysing health-related quality of life outcomes. Stat Med 2005;24(7):1075-1102.
}



