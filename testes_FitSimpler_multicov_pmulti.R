
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 100)
x2 <- seq(0,1, l = 100)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))
dist.matrix = U

### Univariate
# Covariance matrix
cov.model = c("matern", "matern")
nparam <- nparam_covmodel(cov.model)
cov.pars = list(c(0.09, 0.2, 0.5), c(0.09, 0.2, 0.5))
nloc = nrow(U)
ncov.model <- length(cov.model)
SigmaB <- matrix(c(1,0.9,0.9,1),nc = 2)
p = 2
# matrix(c(1,0.9,0.7,
#          0.9,1,0.8,
#          0.7,0.8,1),nc=3)

## Covariance matrix
Sigma2 <- CovSimpler(dist.matrix = U, cov.model = cov.model, 
                     cov.pars = cov.pars, SigmaB = SigmaB)
  
# Simulating data
set.seed(1234)
y <- as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, p*nloc), 
                              sigma = Sigma2))

ini_par <- c(unlist(cov.pars),gdata::upperTriangle(SigmaB, byrow = T))

# if(logpars)
# {
#   ini_par <- c(log(unlist(cov.pars)), atanh(gdata::upperTriangle(SigmaB, byrow = T)))
# }

loglikSimpler_optim_multicov_pmulti(par = ini_par, data = y, nparam = nparam, p=p,
                                    dist.matrix = U, cov.model = cov.model, nloc = nrow(U), 
                                    logpars = F)

ini_par_t <- c(log(unlist(cov.pars)), atanh(gdata::upperTriangle(SigmaB, byrow = T)))
loglikSimpler_optim_multicov_pmulti(par = ini_par_t, data = y, nparam = nparam, p=p,
                                    dist.matrix = U, cov.model = cov.model, nloc = nrow(U), 
                                    logpars = T)


fit1 <- FitSimpler(data = y, dist.matrix = U, cov.model = cov.model, 
           cov.pars = cov.pars, p = p,
           SigmaB = SigmaB, method = "Nelder-Mead", hessian = F, logpars = F)
  
SigmaB <- atanh(SigmaB)
diag(SigmaB) <- 1
fit2 <- FitSimpler(data = y, dist.matrix = U, cov.model = cov.model, 
           cov.pars = lapply(1:p, function(x){log(cov.pars[[x]])}), p = p,
           SigmaB = SigmaB, method = "Nelder-Mead", hessian = F, logpars = T)

fit1$est_Simpler$par
fit2$est_Simpler$par

fit1$est_Simpler$value
fit2$est_Simpler$value


sol1 <- c(0.02929125731, 0.10424420821, 1.32036569974, 0.01891575414, 0.13973918722, 0.69675448117, 0.87854006405)
sol2 <- c(0.02938140619, 0.10429398214, 2.30034375161, 0.01892987718, 0.13973390761, 3.21236957066, 0.87877093028)
loglikSimpler_optim_multicov_pmulti(par = sol1, data = y, nparam = nparam, p=p,
                                    dist.matrix = U, cov.model = cov.model, nloc = nrow(U), 
                                    logpars = F)
loglikSimpler_optim_multicov_pmulti(par = sol2, data = y, nparam = nparam, p=p,
                                    dist.matrix = U, cov.model = cov.model, nloc = nrow(U), 
                                    logpars = F)
