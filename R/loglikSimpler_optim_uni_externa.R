#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler_optim_uni
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom mvtnorm dmvnorm
loglikSimpler_optim_uni <- function(data, coords, dist.matrix, cov.pars, nloc,
                                    cov.model, nugget = 0, compute.dists = TRUE)
{
  ## função de uso interno
  ## Univariate scenario --> ncov.model==1 && (ncol(data) == 1)
  
  ## ini.par --> vetor com todos os parametros.
  # Ex: modelo Exp --> se nugget = T --> par = c(nugget, var, phi),
  # senão --> par = c(var, phi)
  ## cov.model --> vetor de string contendo o modelo da variável
  ## nloc --> número de locações amostrais
  ## data --> vetor de valores da resposta y, que deve ter dimensão length = nloc (univariado)
  ## coords --> coordenadas dos dados: matrix nloc x 2
  ## dist.matrix --> matriz de distâncias euclidianas
  ## entrar com coords ou com dist.matrix
  ## compute.dists --> se as distancias devem ser calculadas ou não. 
  ## Se T o argumento coords deve ser dado
  # nloc <- nrow(data) --> calcular antes de entrar na função
  
  realisations <- as.factor(rep(1, length(data)))
  
  if(compute.dists)
  {
    dist.matrix <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  }
   
  Sigma <- as.matrix(CovSimpler(dist.matrix = dist.matrix, cov.model = cov.model,
                      cov.pars = cov.pars, p = 1, nugget = nugget))

  # ou
  D <- as.vector(rep(1, length = nrow(data))) ## 
  R <- chol(Sigma) 
  # tal que R'R = sigma, onde R = chol(sigma) --> triangular superior
  # nos cálculos a seguir será feita a substituição R'R = sigma, 
  # para calcular a estimativa de beta e a verossimilhança
  invRD <- backsolve(R, D, transpose = TRUE) # = t(D)%*%solve(R)
  invRy <- backsolve(R, data, transpose = TRUE) # = t(t(y)%*%solve(R))
  if(mean == 0)
  {
    bhat <- as.matrix(mean)
  } else
  {
    bhat <- solve(crossprod(invRD), # = invRD%*%t(invRD)
                  crossprod(invRD, invRy)) # = invRD%*%t(invRy)
    # = solve(crossprod(invRD))%*%crossprod(invRD, invRy)
  }
  invRe <- invRy - invRD %*% bhat
  out <- -drop(length(data) * log(2 * pi)/2 + sum(log(diag(R))) +
                 crossprod(invRe)/2)
  
  # sum(log(diag(R))) = 1/2 log|R'R| = 1/2 log|R^2| = log|R|
  # ou se mean = bhat -->
  # (out <-  mvtnorm::dmvnorm(data, mean=rep(bhat, nloc), sigma=Sigma, log=TRUE))
  return(out)
}









  # ## Cross-covariance
  # 
  # Sigma_chol <- lapply(1:p, function(x) 
  # {
  #   Matrix::chol(cov_marg(dist.matrix = dist.matrix, cov.model = cov.model[x],
  #                         cov.pars = par[npos[[x]]]))
  # })
  # 
  # 
  # ## Inverse of each Cholesky matrix
  # Sigma_inv_chol <- lapply(1:p, function(x){Matrix::solve(Sigma_chol[[x]])}) 
  # 
  # ## Log sum of cholesky diagonal elements
  # Sigma_log_sum_chol_diag <- sum(log(sapply(1:p, 
  #                                           function(x){Matrix::diag(Sigma_chol[[x]])})))
  # 
  # ## If p=1 (Univariate case)
  # if(all(is.na(SigmaB)))
  # {
  #   ## Transpose product of y and bdiag
  #   prod_y_bdiag <- Matrix::crossprod(data,Matrix::bdiag(Sigma_inv_chol))
  #   prod_y_bdiag_kron <- Matrix::tcrossprod(prod_y_bdiag, kronecker(1, diag(n)))
  #   nll_uni <- -0.5*drop(n_all*log(2*pi) + 2*Sigma_log_sum_chol_diag + # n*log(det(SigmaB)) +
  #                          + Matrix::tcrossprod(prod_y_bdiag_kron, prod_y_bdiag))
  #   return(nll_uni)
  # }
  # 
  # ## Transpose product of y and bdiag
  # prod_y_bdiag <- crossprod(data,bdiag(Sigma_inv_chol))
  # prod_y_bdiag_kron <- tcrossprod(prod_y_bdiag, kronecker(solve(SigmaB), diag(n)))
  # nll <- -drop(n_all*log(2*pi) + 2*Sigma_log_sum_chol_diag + n*log(det(SigmaB)) +
  #                tcrossprod(prod_y_bdiag_kron, prod_y_bdiag))/2
  # 

