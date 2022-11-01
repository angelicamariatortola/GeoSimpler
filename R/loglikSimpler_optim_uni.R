#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler_optim_uni
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom mvtnorm dmvnorm
loglikSimpler_optim_uni <- function(par, data, dist.matrix, cov.model, nloc, 
                                 logpars = FALSE)
{
  ## função de uso interno
  ## Univariate scenario --> ncov.model==1 && (ncol(data) == 1)
  
  ## par --> vetor com todos os parametros.
  # Ex: c(param_var1)
  ## cov.model --> vetor de string contendo o modelo da variável
  ## nloc --> número de locações amostrais
  ## data --> valores da resposta y, que deve ter dimensão nloc x 1
  
    if(logpars)
    {
      par <- c(exp(par[1:length(par)]))
    }
    
    Sigma <- cov_marg(dist.matrix = dist.matrix, cov.model = cov.model,
                      cov.pars.uni = par)
    out <- mvtnorm::dmvnorm(x = data, mean = rep(0,nloc), 
                            sigma = Sigma, log = TRUE)
  
    ## ou usar a função loglik.GRF {geoR}	??
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

