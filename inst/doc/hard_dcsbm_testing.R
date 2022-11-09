## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
Ktru = 4  # number of true communities
oir = 0.15 # out-in-ratio
h = 1:Ktru 
h = h / sum(h) # class prior or `pri` in our notation

P_diag = rep(1, Ktru)
P = oir + diag(P_diag - oir)

r = nett::sinkhorn_knopp(P, sums=h, sym=T)$r
d = r / h
diag(d) %*% P %*% diag(d) %*% h # verify Sinkhorn's scaling


## -----------------------------------------------------------------------------
n = 2000 # number of nodes
lambda = 5 # expected average degree

# Generate from the degree-corrected Erdős–Rényi model
set.seed(1234)
theta <- EnvStats::rpareto(n, 2/3, 3)
alpha2 = (sum(theta)^2 - sum(theta^2))/n
theta =  theta * sqrt(lambda/alpha2)

A = nett::sample_dcer(theta)  # sample the adjacency matrix matrix 
mean(Matrix::rowSums(A))
hist(Matrix::rowSums(A), 15, main  = "Degree distribution", xlab = "degree")

## -----------------------------------------------------------------------------
z = sample(Ktru, n, replace=T, prob=h) # randomly smaple "true community labels" 
tht = d[z] * theta
A = nett::sample_dcsbm(z, P, tht) # sample the adjacency matrix
hist(Matrix::rowSums(A), 15, main = "Degree distribution", xlab = "degree")

## -----------------------------------------------------------------------------
nrep = 10
degs1 = degs2 = rep(0, n)
for (r in 1:nrep) {
  degs1 = degs1 + Matrix::rowSums(nett::sample_dcer(theta))
  z = sample(Ktru, n, replace=T, prob=h) # randomly smaple "true community labels" 
  degs2 = degs2 + Matrix::rowSums(nett::sample_dcsbm(z, P,  d[z] * theta))
}
degs1 = degs1 / nrep
degs2 = degs2 / nrep
# cbind(degs1, degs2)
mean(abs(degs1 - degs2)/(degs1 + degs2)/2)

## -----------------------------------------------------------------------------
K_H0 = 1
K_H1 = Ktru
apply_methods = function(A) {
  z0 = nett::spec_clust(A, K_H0)
  z0p = nett::spec_clust(A, K_H0+1)
  z1 = nett::spec_clust(A, K_H1)
  tibble::tribble(
    ~method, ~tstat, ~twosided,
    "NAC+", nett::nac_test(A, K_H0, z=z0, y=z0p)$stat, F,
    "SNAC+", nett::snac_test(A, K_H0, z=z0)$stat, F,
    "AS", nett::adj_spec_test(A, K_H0, z=z0), F,
    "LR", nett::eval_dcsbm_loglr(A, cbind(z0, z1), poi=T), F
  )
}

## -----------------------------------------------------------------------------
gen_null_data = function() {
  nett::sample_dcer(theta)
}

gen_alt_data = function() {
  z = sample(Ktru, n, replace=T, prob=h)
  nett::sample_dcsbm(z, P,  d[z] * theta)
}

## -----------------------------------------------------------------------------
# try reducing nruns for faster bu cruder estimation of the ROC curves
roc_res = nett::simulate_roc(apply_methods,
                       gen_null_data = gen_null_data,
                       gen_alt_data = gen_alt_data,
                       nruns = 100, core_count = 1) 
nett::printf('time = %3.3f', roc_res$elapsed_time)



## -----------------------------------------------------------------------------
plot_res <- roc_res$roc
plot_res$method <- factor(plot_res$method, levels = c("NAC+","SNAC+","AS", "LR"))
nett::plot_roc(plot_res)

