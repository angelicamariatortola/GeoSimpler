
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 5)
x2 <- seq(0,1, l = 5)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))

# p = 2
cov.model = "exp"
nugget = 2
cov.pars1 = c(0.5, 0.3)
nloc = nrow(U)
SigmaB =  NULL #matrix(c(1,0.9, 0.9, 1), nc = 2)
p = 1
# Simulating data
Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model, nugget = nugget,
                     cov.pars = cov.pars1, SigmaB = SigmaB, p = p)
set.seed(1234)
y <- matrix(as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, p*nloc),
                                     sigma = Sigma2)), nc = p)

ini.par = c(nugget, cov.pars1); data = y; mean = 2
dist.matrix = U; cov.model = "exp"
nloc = nloc; v.nugget = T; logpars = F

loglikSimpler_optim_uni(ini.par = c( cov.pars1), data = y, mean = 2,
                        dist.matrix = U, cov.model = "exp",
                        nloc = nloc, v.nugget = F, logpars = F)
