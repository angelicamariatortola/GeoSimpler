
rm(list = ls())

library(geoR)

# Generating the coordinates:
x1 <- seq(0,1, l = 3)
x2 <- seq(0,1, l = 3)
coord <- cbind(x1,x2)

# # Building the distance matrix:
# U <- dist(coord)

## Univariado
# cov.model = "exp"
# cov.pars = list(c(0.5, 0.1))
# SigmaB =  NULL
# nugget = 0
# mean = 0
# p = 1
# trend = list("cte")

## Multivariado p = 2 --> mesma estrutura de covariância para todas as respostas
# cov.model = "matern"
# cov.pars = list(c(0.2, 0.4, 0.7)) # para modelos exp, gaussian, etc.
# SigmaB =  matrix(c(1, 0.9, 0.9, 1), nc = 2);
# nugget = 2
# mean = 0
# p = 2
# trend = rep(list("cte"), p)
# list("cte", "1st")

## Multivariado p = 3 --> mesma estrutura de covariância para todas as respostas
# cov.model = "exp"
# cov.pars = list(c(0.5, 0.2)) # para modelos exp, gaussian, etc.
# SigmaB =  matrix(c(1, 0.4, 0.9,
#                    0.4, 1, 0.6,
#                    0.9, 0.6, 1), nc = 3);
# nugget = 0
# mean = 0
# p = nrow(SigmaB)
# n = nrow(coord)
# trend = rep(list("cte"), p)

## Multivariado p = 2 --> diferentes estruturas para as respostas
cov.model = c("exp", "exp")
cov.pars = list(c(0.5, 0.1), c(0.7, 0.2)) # para modelos exp, gaussian, etc.
SigmaB =  matrix(c(1, 0.9, 0.9, 1), nc = 2);
nugget = c(1,2)
mean = c(0,0)
p = 2
trend = list("cte", "cte")

Sigma2 <- CovSimpler(coords = coord, nugget = nugget,
                     cov.model = cov.model, cov.pars = cov.pars, p = p, 
                     SigmaB = SigmaB) 
v1 <- varcov.spatial2(coords = coord, nugget = nugget[1]/cov.pars[[1]][1],
                      cov.model = cov.model[1], cov.pars = c(1, cov.pars[[1]][2]))
v2 <- varcov.spatial2(coords = coord, nugget = nugget[2]/cov.pars[[2]][1],
                      cov.model = cov.model[2], cov.pars = c(1, cov.pars[[2]][2]))

Sigma2$varcov
s1 = sqrt(0.5)*matrix(1, nc = 3, nr = 3)
s2 = sqrt(0.7)*matrix(1, nc = 3, nr = 3)
bds = bdiag(s1,s2) # ou bds = kronecker(bdiag(sqrt(0.5), sqrt(0.7)), matrix(1, nc = 3, nr = 3))
bdvchol <- bdiag(chol(v1$varcov), chol(v2$varcov))

(bds*t(bdvchol))%*%kronecker(SigmaB, diag(3))%*%(bds*(bdvchol))

# det(as.matrix(bds)) = 0



set.seed(56789)
# y <- matrix(mvnfast::rmvn(n = 1, mu = rep(0, p*nrow(coord)),
#                           sigma = Sigma2$varcov), nc = p)
y <- grfSimpler(n = 2, coords = coord, mean = mean, nsim = 1, nugget = nugget,
  cov.model = cov.model, cov.pars = cov.pars, p = p, SigmaB = SigmaB)


ini.cov.pars = cov.pars; data = y$data$sim1; coords = coord;
cov.model = cov.model; trend = trend
nugget = nugget; fix.kappa = T; fix.nugget = T; fix.mean = F
SigmaB = SigmaB; limits = pars.limits2()
print.pars = F


ff <- FitSimpler(ini.cov.pars = cov.pars,  data = y$data$sim1, coords = coord, cov.model = cov.model,
                        fix.kappa = T, fix.nugget = T, fix.mean = F, 
                        nugget = nugget, SigmaB = SigmaB, trend = trend,
                        limits = pars.limits2(), print.pars = FALSE)
ff$parameters.summary





# system.time(
#   {f2 <- likfit(data = y$data$sim1, coords = coord, trend = unlist(trend),
#        ini.cov.pars = unlist(cov.pars),
#        fix.nugget = T, nugget = nugget, fix.kappa = TRUE,
#        components = F, nospatial = F)
#   })
# f2$parameters.summary



# > ff$parameters.summary
# status  values
# beta    estimated -0.0488
# tausq       fixed  0.0000
# phi     estimated  0.0685
# sigmasq estimated  0.3496
# kappa       fixed      NA

# likfit(data = y$data$sim1, coords = coord, trend = unlist(trend),
#        ini.cov.pars = unlist(cov.pars),
#        fix.nugget = F, nugget = nugget, fix.kappa = T,
#        components = F, nospatial = F)
# 
# likfit(data = y$data$sim1, coords = coord, trend = unlist(trend),
#        ini.cov.pars = unlist(cov.pars)[1:2], kappa = unlist(cov.pars)[3],
#        fix.nugget = F, nugget = nugget, fix.kappa = F,
#        components = F, nospatial = F)






# loglik.GRF(coords = coord, data = y$data$sim1, cov.model = cov.model,
#            cov.pars = unlist(cov.pars)[1:2], kappa = unlist(cov.pars)[3], nugget = nugget,
#            trend = unlist(trend))
# 
# loglikSimpler(data = y$data$sim1, coords = coord, trend = trend,
#                          cov.model = cov.model,
#                          cov.pars = cov.pars, est_mean = T,
#                          SigmaB = SigmaB)
