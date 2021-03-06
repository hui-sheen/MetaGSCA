% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrapP.R
\name{bootstrapP}
\alias{bootstrapP}
\title{Append final two bootstrap statistics}
\usage{
bootstrapP(bootstrap, meta.res, fixed.effect = TRUE, random.effect = FALSE)
}
\arguments{
\item{bootstrap}{R object flowing from prior content of MetaGSCA function}

\item{meta.res}{output from metaAndPlot function}

\item{fixed.effect}{logic variable for designation of Fixed Effect Model}

\item{random.effect}{logic variable for designation of Random Effect Model}
}
\value{
Series of bootstrap result statistics with two new items appended: meta.p & bootstrap.p
}
\description{
meta.p & bootstrap.p are appended to series bootstrap statistics.
}
\details{
Either Fixed Effects model or Random Effects model can be designated. Between the two options, result structures are identical but names are different.
}
