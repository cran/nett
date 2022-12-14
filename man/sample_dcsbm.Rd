% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{sample_dcsbm}
\alias{sample_dcsbm}
\title{Sample from a DCSBM}
\usage{
sample_dcsbm(z, B, theta = 1)
}
\arguments{
\item{z}{Node labels (\eqn{n * 1})}

\item{B}{Connectivity matrix (\eqn{K * K})}

\item{theta}{Node connectivity propensity vector (\eqn{n * 1})}
}
\value{
An adjacency matrix following DCSBM
}
\description{
Sample an adjacency matrix from a degree-corrected block model (DCSBM)
}
\examples{
B = pp_conn(n = 10^3, oir = 0.1, lambda = 7, pri = rep(1,3))$B
head(sample_dcsbm(sample(1:3, 10^3, replace = TRUE), B, theta = rexp(10^3)))
}
\seealso{
\link{sample_dcpp}, \link{fast_sbm}, \link{sample_tdcsbm}
}
