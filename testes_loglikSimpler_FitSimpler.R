
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 500)
x2 <- seq(0,1, l = 500)
coord <- cbind(x1,x2)

# Building the distance matrix:
U <- dist(coord)

# # Univariado
cov.model = "exp"
cov.pars = list(c(0.1, 0.2))
SigmaB =  NULL
nugget = 1
mean = 0
p = 1
trend = list("cte")

# # # Multivariado
# cov.model = c("exp", "exp")
# cov.pars = list(c(0.5, 0.15), c(5, 0.15)) # para modelos exp, gaussian, etc.
# SigmaB =  matrix(c(1, 0.9, 0.9, 1), nc = 2);
# nugget = c(1,10)
# mean = c(2,3)
# p = 2
# trend = list("cte", "cte")

Sigma2 <- CovSimpler(coords = coord, nugget = nugget,
                     cov.model = cov.model, cov.pars = cov.pars, p = p, 
                     SigmaB = SigmaB) 

set.seed(56789)
# y <- matrix(mvnfast::rmvn(n = 1, mu = rep(0, p*nrow(coord)),
#                           sigma = Sigma2$varcov), nc = p)

y <- grfSimpler(n = 2, coords = coord, mean = mean, nsim = 1,
  cov.model = cov.model, cov.pars = cov.pars, p = p, SigmaB = SigmaB)
# y$data$sim1

# ini.cov.pars = cov.pars; data = y$data$sim1; coords = coord;
# cov.model = cov.model; trend = trend
# nugget = nugget; fix.kappa = T; fix.nugget = T;
# SigmaB = SigmaB; limits = pars.limits()
# print.pars = F
# est_mean = F

# y2 <- grfSimpler(n = 25, nx = 5, ny = 5, mean = c(2,10), nsim = 1,
#            cov.model = cov.model, cov.pars = cov.pars, p = p, SigmaB = SigmaB)
# y2$data

ff = FitSimpler(ini.cov.pars = cov.pars, data = y$data$sim1, coords = coord,
           cov.model = cov.model, trend = trend, est_mean = T,
           nugget = nugget, fix.nugget = F,
           SigmaB = SigmaB, limits = pars.limits(),
           print.pars = F)
ff$parameters.summary
ff$loglik


likfit(data = y$data$sim1, coords = coord, trend = unlist(trend),
       ini.cov.pars = unlist(cov.pars)[1:2], kappa = unlist(cov.pars)[3],
       fix.nugget = T, nugget = nugget, fix.kappa = T,
       components = F, nospatial = F)


# likfit(data = y$data$sim1, coords = coord, trend = unlist(trend),
#        ini.cov.pars = unlist(cov.pars),
#        fix.nugget = T, nugget = nugget, fix.kappa = T,
#        components = F, nospatial = F)

# loglik.GRF(coords = coord, data = y$data$sim1, cov.model = cov.model,
#            cov.pars = unlist(cov.pars)[1:2], kappa = unlist(cov.pars)[3], nugget = nugget,
#            trend = unlist(trend))
# 
# loglikSimpler(data = y$data$sim1, coords = coord, trend = trend,
#                          cov.model = cov.model,
#                          cov.pars = cov.pars, est_mean = T,
#                          SigmaB = SigmaB)


# $betahat
# [1] 2.070584
# 
# $loglik
# [1] -19.7065

 