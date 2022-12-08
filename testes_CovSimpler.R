
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 5)
x2 <- seq(0,1, l = 5)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))


#### Testes função CovSimpler

# ## Covariancia univariada: length(cov.model)==1 && (p == 1)
# cov_marg(dist.matrix = U, cov.model = "matern",
#          cov.pars.uni = c(2,0.2,0.7))

# CovSimpler(dist.matrix = U, cov.model = "matern",
#            cov.pars = list(c(2,0.2,0.7)), p=1, nugget = 2)

 
# ## Covariancia Bivariada: mesma correlação para variáveis diferentes
# cov.model = c("matern")
# CovSimpler(dist.matrix = U, cov.model = cov.model,
#            cov.pars = list(c(1,0.2,0.7)),
#            SigmaB = matrix(c(1,0.9,
#                              0.9,1),nc = 2), p = 2)


# ## Covariancia Bivariada: correlações diferentes
cov.model = c("matern", "matern")

# dist.matrix = U; cov.model = cov.model;
# cov.pars = list(c(1,0.2,0.7),c(3,0.2,0.701));
# SigmaB = matrix(c(1,0.9,
#                   0.9,1),nc = 2); p =2

CovSimpler(dist.matrix = U, cov.model = cov.model,
           cov.pars = list(c(1,0.2,0.7),c(3,0.2,0.701)), 
           SigmaB = matrix(c(1,0.9,
                             0.9,1),nc = 2), p =2, nugget = c(1,2))


## Covariancia Tri-variada: mesma correlação para variáveis diferentes    
cov.model = c("matern")
cov.pars = list(c(1,0.2,0.7))

CovSimpler(dist.matrix = U, cov.model = cov.model, p = 3,
           cov.pars = cov.pars,  SigmaB = matrix(c(1,0.5,0.7,
                                                   0.5,1,0.8,
                                                   0.7,0.8,1),nc = 3))


## Covariancia Tri-variada: correlações diferentes
cov.model = c("gaussian", "matern", "cauchy")
cov.pars = list(c(1,0.7), c(2,0.21,0.7),c(3,0.22,0.7))

CovSimpler(dist.matrix = U, cov.model = cov.model, p = 3,
           cov.pars = cov.pars,  SigmaB = matrix(c(1,0.5,0.7,
                                                   0.5,1,0.8,
                                                   0.7,0.8,1),nc = 3))
