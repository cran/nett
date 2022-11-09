## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(nett)

## -----------------------------------------------------------------------------
n = 1500
Ktru = 4
lambda = 15 # expected average degree
oir = 0.1
pri = 1:Ktru

set.seed(1234)
theta <- EnvStats::rpareto(n, 2/3, 3)
B = pp_conn(n, oir, lambda, pri=pri, theta)$B
z = sample(Ktru, n, replace=T, prob=pri) # randomly smaple "true community labels" 
A = sample_dcsbm(z, B, theta) # sample the adjacency matrix

## -----------------------------------------------------------------------------
zh = spec_clust(A, K=4)

## -----------------------------------------------------------------------------
compute_mutual_info(z, zh)

## ----  message = FALSE, warning= FALSE----------------------------------------
set.seed(1234)
nrep = 20
nlam = 12
lamvec = 10^seq(log10(1), log10(50), length.out = nlam)  # the vector of logarithmically spaced lambda
runs = expand.grid(rep = 1:nrep, lambda = lamvec)

res = do.call(rbind, lapply(1:nrow(runs), function(j) {
 lambda = runs[j,"lambda"]
 B = pp_conn(n, oir, lambda, pri=pri, theta)$B
 A = sample_dcsbm(z, B, theta)
 zh = spec_clust(A, K = Ktru) # defaults to tau = 0.25 for the  regularization parameter
 zh_noreg = spec_clust(A, K = Ktru, tau = 0)
 data.frame(rep = runs[j,"rep"], lambda = lambda, 
            nmi = compute_mutual_info(z, zh), nmi_noreg = compute_mutual_info(z, zh_noreg))
}))

agg_nmi = aggregate(res, by = list(res$lambda), FUN = mean)

## ----fig.width=6, fig.height=5, fig.align="center"----------------------------
plot(agg_nmi$lambda, agg_nmi$nmi, log="x",
     type = "b", col = "blue", ylab = "NMI", xlab = "lambda", pch=19,
     main="Specral clustering performance")
lines(agg_nmi$lambda, agg_nmi$nmi_noreg, col="red", lty=2, pch=18, type = "b")
legend(1, 1, legend = c("0.25 regularization","No regularization"), 
       col = c("blue","red"), lty=1:2, box.lty=0)

## ----  message = FALSE, warning= FALSE----------------------------------------
set.seed(1234)
nrep = 20
nlam = 12
lamvec = 10^seq(log10(1), log10(200), length.out = nlam)  # the vector of logarithmically spaced lambda
runs = expand.grid(rep = 1:nrep, lambda = lamvec)

res = do.call(rbind, lapply(1:nrow(runs), function(j) {
 lambda = runs[j,"lambda"]
 B = gen_rand_conn(n, Ktru, lambda = lambda, gamma = 0.1, pri=pri)
 A = sample_dcsbm(z, B, theta)
 zh = spec_clust(A, K = Ktru) # defaults to tau = 0.25 for the  regularization parameter
 zh_noreg = spec_clust(A, K = Ktru, tau = 0)
 data.frame(rep = runs[j,"rep"], lambda = lambda, 
            nmi = compute_mutual_info(z, zh), nmi_noreg = compute_mutual_info(z, zh_noreg))
}))

agg_nmi = aggregate(res, by = list(res$lambda), FUN = mean)

## ----fig.width=6, fig.height=5, fig.align="center"----------------------------
plot(agg_nmi$lambda, agg_nmi$nmi, log="x",
     type = "b", col = "blue", ylab = "NMI", xlab = "lambda", pch=19,
     main="Specral clustering performance")
lines(agg_nmi$lambda, agg_nmi$nmi_noreg, col="red", lty=2, pch=18, type = "b")
legend(1, max(agg_nmi$nmi), legend = c("0.25 regularization","No regularization"), 
       col = c("blue","red"), lty=1:2, box.lty=0)

