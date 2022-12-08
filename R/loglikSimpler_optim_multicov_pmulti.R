#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler_optim_multicov_pmulti
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom mvtnorm dmvnorm
loglikSimpler_optim_multicov_pmulti <- function(par, data, dist.matrix, p,
                                                cov.model, nloc, nparam,
                                 logpars = logpars)
{
  ## função de uso interno
  ## Multivariate scenario, com modelos distintos de covariância para as respostas
  ## --> ncov.model > 1 && (ncol(data) > 1)
  ## SigmaB != NULL
  
  ## par --> vetor com todos os parametros.
  # Ex: c(param_var1, param_va2,..., rho12, rho13, ...)
  ## cov.model --> vetor de string contendo os modelos das variáveis
  ## nloc --> número de locações amostrais
  ## data --> valores da resposta y, que deve ter dimensão nloc x p, p = número de respostas

  if(logpars)
  {
    par <- c(exp(par[1:sum(nparam[1:p])]), 
             tanh(par[(sum(nparam[1:p])+1):length(par)]))
  }
  
  # Initial parameters
  cov.pars_all <- split(par, rep(1:length(nparam), nparam))
  ## Correlation matrix
  SigmaB <- metaSEM::vec2symMat(cov.pars_all[[length(nparam)]], diag = F)
  
  ## Cholesky of each univariate matrix
  chol_cov <- lapply(1:p, function(x){Matrix::chol(cov_marg(dist.matrix = dist.matrix, 
                                        cov.model = cov.model[x], 
                                        cov.pars.uni = cov.pars_all[[x]]))})
  
  ## Inverse of each Cholesky matrix
  Sigma_inv_chol_cov <- lapply(1:p, function(x){spam::solve(chol_cov[[x]])})
  
  ## Log sum of cholesky diagonal elements (soma a diagonal de uma matriz e multiplica por p)
  Sigma_log_sum_chol_diag <- sum(log(sapply(1:p, function(x){spam::diag(chol_cov[[x]])})))

  ## Transpose product of y and bdiag
  y <- as.numeric(data)
  prod_y_bdiag <- crossprod(y, bdiag(Sigma_inv_chol_cov))
  ## = t(y)%*%bdiag(Sigma_inv_chol),
  # where bdiag(Sigma_inv_chol) is upper triangular
  
  prod_y_bdiag_kron <- tcrossprod(prod_y_bdiag, kronecker(solve(SigmaB), diag(nloc)))
  #tcrossprod(x,y) = x%*%t(y)
  ## pois kronecker(solve(SigmaB), diag(n_uni)) = t(kronecker(solve(SigmaB), diag(n_uni)))
  nll <- as.numeric(-drop((p*nloc)*log(2*pi) + 2*Sigma_log_sum_chol_diag + nloc*log(det(SigmaB)) +
                            tcrossprod(prod_y_bdiag_kron, prod_y_bdiag))/2)
  
  return(nll)
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

