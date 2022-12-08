
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
y <- as.matrix(dados[3:4])

## Estimation
cov.pars <- list(c(s12[1]^2,1,0.3),c(s12[2]^2,2,0.5))
SigmaB <- matrix(c(1,rho12,rho12,1),nc = 2)
cov.model = c("matern","matern")
nparam <- nparam_covmodel(cov.model) 
p <- length(cov.model) ## número de variáveis
par <-c(unlist(cov.pars), gdata::upperTriangle(SigmaB))


est_multi1 <- FitSimpler(data = y, dist.matrix = U, SigmaB = SigmaB,
                         cov.model = cov.model, cov.pars = cov.pars,
                         method = "Nelder-Mead", hessian = F, logpars = F)
est_multi1$est_Simpler$par
est_multi1$est_Simpler$value
est_multi1$time_Simpler
