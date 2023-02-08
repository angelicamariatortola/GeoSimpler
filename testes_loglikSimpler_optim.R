
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 1000)
x2 <- seq(0,1, l = 1000)
coord <- cbind(x1,x2)

# Building the distance matrix:
U <- dist(coord)

cov.model = "matern" #c("exp", "matern")
cov.pars = list(c(5, 0.3, 0.9)) # list(c(0.5, 0.15), c(0.2, 0.08, 0.7)) # para modelos exp, gaussian, etc.
SigmaB =  NULL #matrix(c(1, 0.9, 0.9, 1), nc = 2);
nugget = 10 # c(1,1)
p = 1

Sigma2 <- CovSimpler(coords = coord, nugget = nugget,
                     cov.model = cov.model, cov.pars = cov.pars, p = p, 
                     SigmaB = SigmaB) 

set.seed(1234)
y <- matrix(mvnfast::rmvn(n = 1, mu = rep(0, p*nrow(coord)),
                          sigma = Sigma2$varcov), nc = p)

# data = y; coords = coord; dists.lowertri = NULL
# cov.model = cov.model;
# cov.pars = cov.pars; est_mean = T;
# SigmaB = SigmaB
# 
# 
# trend = list("cte", "1st")

system.time(loglik.GRF(coords = coord, data = y, cov.model = cov.model,
           cov.pars = unlist(cov.pars)[1:2], kappa = unlist(cov.pars)[3], nugget = nugget,
           trend = "1st"))

system.time(loglikSimpler(data = y, coords = coord, trend = list("1st"),
                         cov.model = cov.model, nugget = nugget,
                         cov.pars = cov.pars, est_mean = T,
                         SigmaB = SigmaB))
ll$loglik

 