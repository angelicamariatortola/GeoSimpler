
rm(list = ls())

library(mvnfast)

setwd("/media/angelica/HD/Pesquisa/Projetos/Modelagem de dados espaciais multivariados/GeoSimpler package/GeoSimpler/R")

cov.model = c("matern")#, "matern")
nparam <- nparam_covmodel(cov.model) ## vetor com o número de parâmetros de cada modelo, e o num dos parâmetros de correlação, se houverem
## a soma do vetor nparam retorna o número total de parâmetros no modelo
p <- length(cov.model) ## número de variáveis
SigmaB <- matrix(c(1,0.75,0.75,1), nc = 2)
cov.pars1 = list(c(0.09, 0.2, 0.5))#, c(0.09, 0.2, 0.5))
par <- c(unlist(cov.pars1), 0.75)#, gdata::upperTriangle(SigmaB))





# Generating the coordinates:
ncord <- c(10, 15, 20) #, 25, 30)
# time_i <- c()
# est_i <- c()
# value_i <- c()

for (i in 1:length(ncord))
{
  # Coordinates
  x1 <- seq(0,1, l = ncord[i])
  x2 <- seq(0,1, l = ncord[i])
  coord <- expand.grid(x1,x2)
  
  # Distance matrix
  U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))
  
  # Covariance matrix
  Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model, cov.pars = cov.pars1,
                       SigmaB = SigmaB, marg.error = rep(list(c(0.1,0.1)), p), id_aux1 = 0)
  
  # Simulating data
  set.seed(1234)
  y <- as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, p*(ncord[i]^2)), 
                       sigma = Sigma2))
  
  # Estimating
  id_aux1 <- 0
  time_nvar <- system.time(est_nvar <- try(optim(par = ini_par,
                 fn = loglikSimpler_optim, data = y, dist.matrix = U, p=p, 
                 n = nrow(U), nparam = nparam, marg.error = rep(list(c(0.01,0.01)), p),
                 cov.model = cov.model, logpars = F, 
                 control = list(fnscale = -1, maxit = 2000),
                 method="Nelder-Mead", hessian = F)))
  
  
  
  est_pars <- c(est_nvar$par[c(2,5,3,6)], sqrt(est_nvar$par[c(1,4)]), est_nvar$par[7])
  
  if(class(est_nvar) != 'try-error')
  {
    time_i <- rbind(time_i, c(ncord[i], as.numeric(time_nvar[3])))
    est_i <- rbind(est_i, est_pars)
    value_i <- rbind(value_i, est_nvar$value)
  }
  results <- list(time_i = time_i, est_i = est_i, value_i = value_i)
  print(results)
}






# par = par; data = y1; dist.matrix = U; n = nrow(U); p = p; nparam = nparam;
# cov.model = cov.model; marg.error = NULL; logpars = FALSE

# loglikSimpler_optim(par = par,
#                     #   c(log(par[1:sum(nparam[1:p])]),
#                     # atanh(par[(sum(nparam[1:p])+1):length(par)])), 
#                     data = y1, dist.matrix = U, n = nrow(U), p = p, nparam = nparam,
#                     cov.model = cov.model, marg.error = NULL, logpars = F)

loglikSimpler_optim(par = par,
                    data = y1, dist.matrix = U, n = nrow(U), p = p, nparam = nparam,
                    cov.model = cov.model, marg.error = NULL, logpars = T)

# Sigma <- CovSimpler(dist.matrix = U, cov.model = c("matern", "exp"), SigmaB = SigmaB,
#                           cov.pars = list(c(0.5, 0.1, 0.8),c(2, 0.6)))
# out <- mvtnorm::dmvnorm(x = y1, mean = rep(0,length(y1)), sigma = Sigma, log = TRUE)
# out


#library(geoR)
# loglik.GRF(coords = coord, data = y1,
#            cov.model = "matern", cov.pars = c(0.5, 0.1, 0.8))


