% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plt.mtDNA.model}
\alias{plt.mtDNA.model}
\title{Plot the modeling of mtDNA mutation distribution}
\usage{
plt.mtDNA.model(
  res.tbl.best,
  read.data,
  plot.type = "CDF",
  read.clr = "black",
  brks = 0:25/25,
  clr.scheme = "rainbow"
)
}
\arguments{
\item{res.tbl.best}{A data frame from get.best.model() output with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional), each row represents a model}

\item{read.data}{A matrix with columns: mut.read and wt.read, representing mutant and wildtype mitochondria read count for each single cell}

\item{plot.type}{Options must be one of "CDF" or "OR"}

\item{read.clr}{Color for observed data in CDF plot}

\item{brks}{Breaks for bar plot of observed reads}

\item{clr.scheme}{Default color scheme for models if no clr specified in res.tbl.best}
}
\value{
Plot the selection and no selection modeling of mtDNA mutation distribution respect to CDF or log odds ratio
}
\description{
Plot the modeling of mtDNA mutation distribution respect to CDF or log odds ratio
}
\examples{
res.tbl.best = get.best.model(res.tbl = res.tbl,K = 10,show.best = TRUE)
plt.mtDNA.model(res.tbl.best = res.tbl.best,
                 read.data = reads,
                 plot.type = "CDF")
}
\seealso{
\code{\link{get.best.model}} to get the best selection model and best no selection model
}
