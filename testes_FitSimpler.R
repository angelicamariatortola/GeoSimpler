

rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 2)
x2 <- seq(0,1, l = 2)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))

library(geoR)

# p = 2
cov.model = c("exp")
cov.pars1 = list(c(0.5, 0.3))
nloc = nrow(U)
SigmaB =  NULL #matrix(c(1,0.9, 0.9, 1), nc = 2)
p = 1
# Simulating data
Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model,
                     cov.pars = cov.pars1, SigmaB = SigmaB, p = p)
set.seed(1234)
y <- matrix(as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, p*nloc),
                                     sigma = Sigma2)), nc = p)

# "gencauchy", "gneiting.matern" --> nparam = 4
# "matern", "cauchy", "powered.exponential"   --> nparam = 3
# “exponential”, “gaussian”, “spherical”, “circular”, “cubic”, 
# “wave”, “linear”, “power”, “stable”, “gneiting”, “pure.nugget”  --> nparam = 2

# data = y, dist.matrix = U, SigmaB = SigmaB,
# cov.model = cov.model, cov.pars = cov.pars1,
# method = "Nelder-Mead", hessian = F, logpars = F

est_multi1 <- FitSimpler(data = y, dist.matrix = U, SigmaB = SigmaB,
                         cov.model = cov.model, cov.pars = cov.pars1,
                         method = "Nelder-Mead", hessian = F, logpars = F)
print(est_multi1)
summary(est_multi1)









# cov.model = c("gencauchy","gencauchy","gencauchy")
# cov.pars1 = list(c(0.1, 0.2, 0.6, 0.5), c(0.1, 0.2, 0.6, 0.5), c(0.1, 0.2, 0.6, 0.5))
# nloc = nrow(U)
# SigmaB = matrix(c(1,0.9,0.8,
#                   0.9, 1, 0.5,
#                   0.8, 0.5, 1), nc = 3)
# p = 3
#geoR::practicalRange(cov.model = "matern", phi = 0.3, kappa = 0.3)

# data = y; dist.matrix = U; SigmaB = SigmaB;
# cov.model = c("matern","gencauchy"); cov.pars = list(c(0.1, 0.2, 0.6), c(0.1, 0.2, 0.6, 0.5));
# method = "Nelder-Mead"; hessian = F; logpars = F


est_multi1 <- FitSimpler(data = y, dist.matrix = U, SigmaB = SigmaB,
                      cov.model = "gencauchy", #c("gencauchy","gencauchy","gencauchy"), 
                      cov.pars = list(c(0.1, 0.2, 0.6, 0.5)), # , c(0.1, 0.2, 0.6, 0.5), c(0.1, 0.2, 0.6, 0.5)),
                      method = "Nelder-Mead", hessian = F, logpars = F)
print(est_multi1)
summary(est_multi1)



est_multi1$est_Simpler$par
est_multi1$est_Simpler$value
est_multi1$time_Simpler


# modelo completo
est_multi2 <- FitSimpler(data = y, dist.matrix = U, SigmaB = SigmaB,
                         cov.model = cov.model, cov.pars = cov.pars1, 
                         method = "Nelder-Mead", hessian = F, logpars = F)
est_multi2$est_Simpler$par
est_multi2$est_Simpler$value
est_multi2$time_Simpler


## TRV
# lrtest(fullmodel, reducedmodel)
# lrtest.FitSimpler

(x1 <- -2*(est_multi1$est_Simpler$value - est_multi2$est_Simpler$value))
gl <- length(est_multi2$est_Simpler$par)-length(est_multi1$est_Simpler$par)

pchisq(x1, gl, lower.tail = F)

lrtest.FitSimpler(fullmodel = c(est_multi2$est_Simpler$value,
                                length(est_multi2$est_Simpler$par)), 
                  reducedmodel = c(est_multi1$est_Simpler$value,
                                length(est_multi1$est_Simpler$par)))


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
