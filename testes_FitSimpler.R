
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 20)
x2 <- seq(0,1, l = 20)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))

### Univariate
# Covariance matrix
cov.model = c("matern")
cov.pars1 = list(c(0.09, 0.2, 0.5))
nloc = nrow(U)

Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model, 
                     cov.pars = cov.pars1, SigmaB = NULL)

# Simulating data
set.seed(1234)
y <- as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, nloc), 
                              sigma = Sigma2))

loglikSimpler_optim_uni(par = unlist(cov.pars1), data = y, dist.matrix = U, 
                        cov.model = cov.model, 
                        nloc = nrow(U), logpars = FALSE)

est_uni <- FitSimpler(data = y, dist.matrix = U, 
           cov.model = cov.model, cov.pars = cov.pars1, p = 1,
           method = "Nelder-Mead", hessian = F, logpars = F)
loglik.GRF(coords = coord, data = y, cov.pars =c(0.09, 0.2),
       kappa = 0.5, nugget = 0, cov.model = cov.model)
