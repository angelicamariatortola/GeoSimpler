#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler_optim_unicov_pmulti
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom mvtnorm dmvnorm
loglikSimpler_optim_unicov_pmulti <- function(par, data, dist.matrix, p,
                                              cov.model, nloc, nparam,
                                 logpars = logpars)
{
  ## função de uso interno
  ## Multivariate scenario, com o mesmo modelo de covariância para todas as respostas
  ## --> ncov.model==1 && (ncol(data) > 1)
  
  ## par --> vetor com todos os parametros.
  # Ex: c(param_var, rho12, rho13, ...)
  ## cov.model --> vetor de string contendo o modelo das variáveis
  ## nloc --> número de locações amostrais
  ## data --> valores da resposta y, que deve ter dimensão nloc x p, p = número de respostas
  
    if(logpars)
    {
      par <- c(exp(par[1:nparam]), tanh(par[(nparam+1):length(par)]))
    }
    
    # Initial parameters
    cov.pars <- par[1:nparam]
    ## Correlation matrix
    SigmaB <- metaSEM::vec2symMat(par[(nparam+1):length(par)], diag = F)
    
    ## Cholesky of the univariate covariance matrix
    chol_cov_uni <- Matrix::chol(cov_marg(dist.matrix = dist.matrix, 
                                     cov.model = cov.model, 
                                     cov.pars.uni = unlist(cov.pars)))
    
    ## Inverse of Cholesky matrix
    inv_chol_cov_uni <- spam::solve(chol_cov_uni)
    
    ## Atribuindo a mesma inv_chol_cov_uni para todas as respostas
    Sigma_inv_chol_cov <- lapply(1:p, function(x) {inv_chol_cov_uni}) 
  
    ## Log sum of cholesky diagonal elements (soma a diagonal de uma matriz e multiplica por p)
    Sigma_log_sum_chol_diag <- p*(sum(log(spam::diag(chol_cov_uni))))
    
    ## Transpose product of y and bdiag
    prod_y_bdiag <- crossprod(data, bdiag(Sigma_inv_chol_cov))
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

