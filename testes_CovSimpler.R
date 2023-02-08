
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 2)
x2 <- seq(0,1, l = 2)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- dist(coord)
U_mat = dist(coord, diag = T, upper = T)
# dists.lowertri = U
# coords = NULL
# 
p=1
cov.model = c("matern") #, "matern"); 
cov.pars = list(c(3,0.3,0.6)) #, c(3,0.3,0.6)); nugget = c(0,0);
SigmaB <- NULL
  # matrix(c(1,0.5,
  #                  0.5,1),nc=2)
# SigmaB <- matrix(c(1,0.3,0.5,
#                    0.3,1,0.2,
#                    0.5,0.2,1),nc=3)

# param_check_break(cov.pars = cov.pars, cov.model = cov.model, p = p) 

# s1 = CovSimpler(dists.lowertri = U, cov.model = cov.model, #nugget = nugget,
#            cov.pars = cov.pars, p = 1, SigmaB = SigmaB) 
CovSimpler(dists.lowertri = U, cov.model = cov.model, nugget = nugget,
           cov.pars = cov.pars, p = p, SigmaB = SigmaB) 






#### Testes função CovSimpler

# ## Covariancia Bivariada: mesma correlação para variáveis diferentes
cov.model = c("matern")
CovSimpler(dist.matrix = U, cov.model = cov.model,
           cov.pars = list(c(1,0.7,0.7)),
           SigmaB = matrix(c(1,0.9,
                             0.9,1),nc = 2), nugget = 2, p = 2)


# ## Covariancia Bivariada: correlações diferentes
cov.model = c("matern", "matern")

CovSimpler(dist.matrix = U, cov.model = cov.model,
           cov.pars = list(c(1,0,0.7),c(3,0.7,0.701)), 
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
