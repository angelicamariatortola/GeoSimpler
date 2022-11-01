
rm(list = ls())

# Generating the coordinates:
x1 <- seq(0,1, l = 100)
x2 <- seq(0,1, l = 100)
coord <- cbind(x1,x2)
nrow(coord)

# Building the distance matrix:
U <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))

### Univariate
# Covariance matrix
cov.model = c("matern")
nparam <- nparam_covmodel("matern")
cov.pars = list(c(0.09, 0.2, 0.5))
nloc = nrow(U)
ncov.model <- length(cov.model)
SigmaB <- matrix(c(1,0.9,0.9,1),nc = 2)
p = 2

  # matrix(c(1,0.9,0.7,
  #          0.9,1,0.8,
  #          0.7,0.8,1),nc=3)




if(ncov.model==1 && (p > 1))
{
  ## calcula a cholesky das matrizes de covariancia marginais para as p variáveis
  cov_uni <- Matrix::chol(cov_marg(dist.matrix = U, 
                        cov.model = cov.model, 
                        cov.pars.uni = unlist(cov.pars)))
  
  ct <- lapply(1:p, function(x) {cov_uni})
  
  cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(nloc)))
  Sigma <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
  Sigma2 <- as.matrix(Sigma)
}


# Simulating data
set.seed(1234)
y <- as.numeric(mvnfast::rmvn(n = 1, mu = rep(0, p*nloc), 
                              sigma = Sigma2))

ini_par <- c(unlist(cov.pars),gdata::upperTriangle(SigmaB, byrow = T))
loglikSimpler_optim_unicov_pmulti(par = ini_par, data = y, nparam = nparam, p=p,
                                    dist.matrix = U, cov.model = "matern", nloc = nrow(U), 
                                    logpars = F)

ini_par_trans <- c(log(unlist(cov.pars)), atanh(gdata::upperTriangle(SigmaB, byrow = T)))
loglikSimpler_optim_unicov_pmulti(par = ini_par_trans, data = y,  nparam = nparam, p=p,
                                    dist.matrix = U, cov.model = "matern", nloc = nrow(U), 
                                    logpars = T)


fit1 <- FitSimpler(data = y, dist.matrix = U, cov.model = "matern", 
           cov.pars = cov.pars, p = p,
           SigmaB = SigmaB, method = "Nelder-Mead", hessian = F, logpars = F)
  
fit2 <- FitSimpler(data = y, dist.matrix = U, cov.model = "matern", 
           cov.pars = cov.pars, p = p,
           SigmaB = SigmaB, method = "BFGS", hessian = F, logpars = T)
