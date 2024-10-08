% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_function.R
\name{dnchg}
\alias{dnchg}
\title{Probability mass function of non-central hypergeometric distribution}
\usage{
dnchg(x, n1, n2, m1, psi)
}
\arguments{
\item{x}{The location to evaluate the probability}

\item{n1}{The size of group one}

\item{n2}{The size of group two}

\item{m1}{The size of both two groups}

\item{psi}{Odds ratio}
}
\value{
The probability at point x
}
\description{
Evaluate probability of a single point from non-central hypergeometric distribution
}
\examples{
prob = dnchg(1,2,38,20,1)
}
\references{
J. G. Liao and Ori Rosen. 2001. “Fast and Stable Algorithms for Computing and Sampling From the Noncentral Hypergeometric Distribution." The American Statistician. 55: 366-369.
}
