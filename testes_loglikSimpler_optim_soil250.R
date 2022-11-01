
rm(list = ls())

## Packages

# library(metaSEM)
# library(mvnfast)
# library(gdata) ## upperTriangle function

##  Loading data

data("soil250", package = "geoR")
dados12 <- data.frame(X=soil250$Linha/10, Y=soil250$Coluna/10, 
                      var1 = soil250$H, var2 = soil250$CTC)

## Transforming
dados <- transform(dados12,
                   var1 = resid(lm(var1 ~ 1, data = dados12)),
                   var2 = resid(lm(var2 ~ 1, data = dados12)))

## Initial value parameters
s12 <- round(c(sd(dados$var1), sd(dados$var2)), 2)
s12
rho12 <- round(cor(dados$var1, dados$var2), 2)
rho12

## Distance Matrix

U <- as.matrix(dist(as.matrix(dados[, 1:2]))) ##entra-se com long e lat
y <- as.numeric(as.matrix(dados[3:4]))

## Estimation
cov.pars <- list(c(s12[1]^2,1,0.3),c(s12[2]^2,2,0.5))
SigmaB <- matrix(c(1,rho12,rho12,1),nc = 2)
cov.model = c("matern","matern")
nparam <- nparam_covmodel(cov.model) 
p <- length(cov.model) ## número de variáveis
par <-c(unlist(cov.pars), gdata::upperTriangle(SigmaB))

par <- c(0.393129, 1.077, 0.543, 0.8464, 1.994, 0.513, 0.823)
loglikSimpler_optim(par = par, data = y, dist.matrix = U, n = nrow(U), p = p, nparam = nparam,
                    cov.model = c("matern","matern"), marg.error = NULL, logpars = F)
# 
# s2 <- CovSimpler(dist.matrix = U, cov.model = c("matern","matern"), cov.pars = cov.pars, SigmaB = SigmaB)
# mvnfast::dmvn(X = y, mu = rep(0,length(y)), sigma = s2, log = TRUE)

est_MatSimp_time <- system.time(est_MatSimp <- optim(par = par,
                     fn = loglikSimpler_optim, 
                     data = y, dist.matrix = U, p=p, n = nrow(U), nparam = nparam,
                     cov.model = cov.model, logpars = F, 
                     marg.error = NULL, #rep(list(c(5,5)), p),
                     control = list(fnscale = -1, maxit = 2000),
                     method="Nelder-Mead", hessian = F))
est_MatSimp$convergence
est_MatSimp$par
est_MatSimp$value
round(est_MatSimp_time[3],3)

