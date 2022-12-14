% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{sample_tdcsbm}
\alias{sample_tdcsbm}
\title{Sample truncated DCSBM (fast)}
\usage{
sample_tdcsbm(z, B, theta = 1)
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
Sample an adjacency matrix from a truncated degree-corrected block model
(DCSBM) using a fast algorithm.
}
\details{
The function  samples an adjacency matrix from a truncated DCSBM, with
entries having Bernoulli distributions with mean \deqn{ E[A_{ij} | z] =
B_{z_i, z_j} \min(1, \theta_i \theta_j).} The approach uses the masking idea
of Aiyou Chen, leading to fast sampling for sparse networks. The masking,
however, truncates \eqn{\theta_i \theta_j} to at most 1, hence
we refer to it as the truncated DCSBM.
}
\examples{
B = pp_conn(n = 10^4, oir = 0.1, lambda = 7, pri = rep(1,3))$B
head(sample_tdcsbm(sample(1:3, 10^4, replace = TRUE), B, theta = rexp(10^4)))
}
