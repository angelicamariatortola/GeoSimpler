
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 500)
x2 <- seq(0,1, l = 500)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))

## Multivariate scenario, com o mesmo modelo de covariância para todas as respostas

cov.model = c("exp")
cov.pars1 = list(c(0.1, 0.2))
nloc = nrow(U)
p = 2
SigmaB = matrix(c(1,0.9,0.9,1),nc = 2)
geoR::practicalRange(cov.model = "exponential", phi = 0.2)


# Simulating data
Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model, 
                     cov.pars = cov.pars1, SigmaB = SigmaB, p = p)
set.seed(1234)
y <- as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, p*nloc), 
                              sigma = Sigma2))


est_multi1 <- FitSimpler(data = y, dist.matrix = U, SigmaB = SigmaB,
                      cov.model = cov.model, cov.pars = cov.pars1, p = p,
                      method = "Nelder-Mead", hessian = F, logpars = F)
est_multi1$est_Simpler$par
est_multi1$est_Simpler$value
est_multi1$time_Simpler


## Multivariate scenario, modelos distintos de covariância

cov.model = c("exp", "exp")
cov.pars1 = list(c(0.1, 0.2), c(0.1, 0.2))
nloc = nrow(U)
SigmaB = matrix(c(1,0.9,0.9,1), nc = 2)
p = 2
#geoR::practicalRange(cov.model = "matern", phi = 0.3, kappa = 0.3)

est_multi2 <- FitSimpler(data = y, dist.matrix = U, SigmaB = SigmaB,
                         cov.model = cov.model, cov.pars = cov.pars1, p = p,
                         method = "Nelder-Mead", hessian = F, logpars = F)
est_multi2$est_Simpler$par
est_multi2$est_Simpler$value
est_multi2$time_Simpler


## TRV
-2*(est_multi1$est_Simpler$value - est_multi2$est_Simpler$value)
gl <- length(est_multi2$est_Simpler$par)-length(est_multi1$est_Simpler$par)

qchisq(0.95, gl)

# ### Univariate 
# cov.model = c("matern")
# cov.pars1 = list(c(0.1, 0.3, 0.3))
# nloc = nrow(U)
# geoR::practicalRange(cov.model = "matern", phi = 0.3, kappa = 0.3)
# 
# # Simulating data
# Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model, 
#                      cov.pars = cov.pars1, SigmaB = NULL, p = 1)
# set.seed(1234)
# y <- as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, nloc), 
#                               sigma = Sigma2))
# 
# est_uni <- FitSimpler(data = y, dist.matrix = U, 
#            cov.model = cov.model, cov.pars = cov.pars1, p = 1,
#            method = "Nelder-Mead", hessian = F, logpars = F)
# est_uni$est_Simpler$par
