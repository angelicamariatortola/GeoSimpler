#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom Matrix chol solve diag crossprod tcrossprod bdiag
#' 
#' @export
#' 
loglikSimpler <- function(dist.matrix, cov.model, cov.pars, data, SigmaB = NA)
{
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros de cada modelo
  ## SigmaB --> (opcional) matriz de correlação de dimensão pxp, necessária se p=length(cov.model) > 1
  
  
  ##-- Time
  time_start <- Sys.time()
  
  ##-- Tests
  
  n <- nrow(dist.matrix) ## número de locações amostrais
  p <- length(cov.model) ## número de variáveis
  n_all <- n*p
  
  ## Cross-covariance
  
  Sigma_chol <- lapply(1:p, function(x) 
  {
    Matrix::chol(cov_marg(dist.matrix = dist.matrix, cov.model = cov.model[x],
                          cov.pars = cov.pars[[x]]))
  })
  
  ## Inverse of each Cholesky matrix
  Sigma_inv_chol <- lapply(1:p, function(x){Matrix::solve(Sigma_chol[[x]])}) 
  
  ## Log sum of cholesky diagonal elements
  Sigma_log_sum_chol_diag <- sum(log(sapply(1:p, 
                                            function(x){Matrix::diag(Sigma_chol[[x]])})))
  
  ## If p=1 (Univariate case)
  if(all(is.na(SigmaB)))
  {
    ## Transpose product of y and bdiag
    prod_y_bdiag <- Matrix::crossprod(data,Matrix::bdiag(Sigma_inv_chol))
    prod_y_bdiag_kron <- Matrix::tcrossprod(prod_y_bdiag, kronecker(1, diag(n)))
    nll_uni <- as.numeric(-drop(n_all*log(2*pi) + 2*Sigma_log_sum_chol_diag + # n*log(det(SigmaB)) +
                                  Matrix::tcrossprod(prod_y_bdiag_kron, prod_y_bdiag))/2)
    return(nll_uni)
  }
  
  ## If p>1 (Multivariate case)
  ## Transpose product of y and bdiag
  prod_y_bdiag <- Matrix::crossprod(data,Matrix::bdiag(Sigma_inv_chol))
  prod_y_bdiag_kron <- Matrix::tcrossprod(prod_y_bdiag, kronecker(Matrix::solve(SigmaB), diag(n)))
  nll <- as.numeric(-drop(n_all*log(2*pi) + 2*Sigma_log_sum_chol_diag + n*log(det(SigmaB)) +
                            Matrix::tcrossprod(prod_y_bdiag_kron, prod_y_bdiag))/2)
  
  return(nll)
}
