
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 2)
x2 <- seq(0,1, l = 2)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))


#### Testes função CovSimpler

# ## Covariancia univariada:
cov_marg(dist.matrix = U, cov.model = "matern",
         cov.pars.uni = c(2,0.2,0.7))

CovSimpler(dist.matrix = U, cov.model = "matern",
           cov.pars = list(c(2,0.2,0.7)))

# is.positive.definite(CovSimpler(dist.matrix = U, cov.model = "gaussian",
#            cov.pars = list(c(1,0.25))))
 
# ## Covariancia Bivariada:
cov.model = c("matern", "matern")
CovSimpler(dist.matrix = U, cov.model = cov.model,
           cov.pars = list(c(1,0.2,0.7),c(3,0.2,0.701)), 
           SigmaB = matrix(c(1,0.9,
                             0.9,1),nc = 2))


## Covariancia Tri-variada - funções de correlação diferentes:       
cov.model = c("matern", "matern", "matern")
cov.pars = list(c(1,0.2,0.7), c(2,0.21,0.7),c(3,0.22,0.7))

system.time(cov_tri1 <- CovSimpler(dist.matrix = U, cov.model = cov.model,
                                   cov.pars = cov.pars,  SigmaB = matrix(c(1,0.5,0.7,
                                                                           0.5,1,0.8,
                                                                           0.7,0.8,1),nc = 3)))

system.time(cov_tri2 <- CovSimpler(dist.matrix = U, cov.model = cov.model,
                                   cov.pars = cov.pars,  SigmaB = matrix(c(1,0.5,0.7,
                                                                           0.5,1,0.8,
                                                                           0.7,0.8,1),nc = 3)))


## Covariancia Tri-variada - funções de correlação iguais:                                                            
CovSimpler(dist.matrix = U, cov.model = c("matern", "matern", "matern"),
           cov.pars = list(c(1.5,0.8,0.9), c(1,0.2,0.7),c(3,0.7,1)),  SigmaB = matrix(c(1,0.5,0.7,
                                                                                        0.5,1,0.8,
                                                                                        0.7,0.8,1),nc = 3))
## Covariancia Tri-variada --> Dados Independentes
CovSimpler(dist.matrix = U, cov.model = c("exp", "matern", "gaussian"), 
           cov.pars = list(c(1,0.2), c(1,0.2,0.5), c(1,0.5)), 
           SigmaB = diag(3))

## Covariancia Bi-variada --> modelo para dados separáveis 
## (mesma estrutura de correlação: phi1=phi2, v1=v2 (matern))
CovSimpler(dist.matrix = U, cov.model = c("matern", "matern"), 
           cov.pars = list(c(1,0.2,0.5), c(2,0.2,0.5)), 
           SigmaB = matrix(c(1,0.9,
                             0.9,1),nc = 2))

