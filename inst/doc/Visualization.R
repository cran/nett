## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(nett)
library(igraph)

## -----------------------------------------------------------------------------
n = 1500
Ktru = 4
lambda = 15 # expected average degree
oir = 0.1
pri = 1:Ktru

set.seed(1234)
theta <- EnvStats::rpareto(n, 2/3, 3)
B = pp_conn(n, oir, lambda, pri=pri, theta)$B
z = sample(Ktru, n, replace=T, prob=pri)

# sample the adjacency matrix
A = sample_dcsbm(z, B, theta)

## -----------------------------------------------------------------------------
original = par("mar")

gr = igraph::graph_from_adjacency_matrix(A, "undirected") # convert to igraph object 
par(mar = c(0,0,0,0))
out = nett::plot_net(gr, community = z)

par(mar = original)

## -----------------------------------------------------------------------------
nett::plot_deg_dist(gr)
summary(igraph::degree(out$gr))

## -----------------------------------------------------------------------------
d = Ktru
labels = sample(Ktru, n, replace = T, prob = pri)
labels = sort(labels)
mu = diag(Ktru)
x = 2*mu[labels, ] + 0.75*matrix(rnorm(n*d), n)

A = sample_dclvm(x, lambda, theta)

## -----------------------------------------------------------------------------
original = par("mar")

gr = igraph::graph_from_adjacency_matrix(A, "undirected") # convert to igraph object 
par(mar = c(0,0,0,0))
out = nett::plot_net(gr, community = labels)

par(mar = original)

## -----------------------------------------------------------------------------
nett::plot_deg_dist(gr)
summary(igraph::degree(out$gr))

## -----------------------------------------------------------------------------
original = par("mar")

par(mar = c(0,0,0,0))
out = nett::plot_net(polblogs, community = igraph::V(polblogs)$community)

par(mar = original)

## -----------------------------------------------------------------------------
nett::plot_deg_dist(polblogs)
summary(igraph::degree(polblogs))

