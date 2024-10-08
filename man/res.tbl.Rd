% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{res.tbl}
\alias{res.tbl}
\title{Result output from full.mle()}
\format{
A dataframe with columns "mutant.start","nmito.start","generation","null.nlogL","alot.nlogL","beta3","beta2","beta1","beta0"
}
\usage{
res.tbl
}
\description{
The example of result table from full.mle() for some functions in this package
}
\details{
Model estimations result output from full.mle(). The dataset contains 14398 models (rows), the mitochondrial copy number of parent cell (nmito.start) is 312,
the initial mutant mitochondria read count for each single cell (mutant.start) is from 0 to 312, the replication is from 0 to 45 generation. For each model, the negative log-likelihood of
non-selection (null.nlogL) and selection (alt.nlog.L) is calculated. Also, the betas are estimated (beta0-3) for selection model cubic function for log.psi
}
\seealso{
\code{\link{full.mle}}
}
\keyword{datasets}
