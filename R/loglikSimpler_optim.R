#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler_optim
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom Matrix chol solve diag crossprod tcrossprod bdiag
loglikSimpler_optim <- function(par, data, dist.matrix, cov.model, nloc, 
                                ncov.model, nparam, logpars = FALSE)
{
  ## função de uso interno
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais

  ## par --> vetor com todos os parametros.
  # Ex: c(param_var1, param_var2, ..., rho12, rho13,...)
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## ncov.model = length(cov.model) --> número de modelos de covariancia
  ## nloc --> número de locações amostrais
  
  ## Marginal-covariance (univariate)
  if(ncov.model==1 && (ncol(data) == 1))
  {
    if(logpars)
    {
      par <- c(exp(par[1:length(par)]))
    }
    
    Sigma <- cov_marg(dist.matrix = dist.matrix, cov.model = cov.model,
                      cov.pars.uni = par)
    out <- mvtnorm::dmvnorm(x = data, mean = rep(0,length(data)), 
                            sigma = Sigma, log = TRUE)
    return(out)
  }
  
  ## Marginal-covariance (ncol(data) > 1), ncov.model==1 && (ncov.pars == 1))
  if(ncov.model==1 && (ncol(data) > 1))
  {
    if(logpars)
    {
      par <- c(exp(par[1:length(par)]))
    }
    

    return(out)
  }
  
  
  
  
  if(logpars)
  {
    par <- c(exp(par[1:sum(nparam[1:ncov.model])]),
             tanh(par[(sum(nparam[1:ncov.model])+1):length(par)]))
  }
  
  cov.pars_all <- split(par, rep(1:length(nparam), nparam))
  # quebra o vetor de parâmetros em uma lista, separando 
  # os parâmetros de cada variável e do parâmetro de correlação.
  cov.pars <- lapply(1:ncov.model, function(x) cov.pars_all[[x]])
  # separa apenas os parâmetros das covariâncias marginais
  SigmaB <- metaSEM::vec2symMat(cov.pars_all[[(ncov.model+1)]], diag = F)
  # separa apenas os parâmetros da função de correlação SigmaB
  
  #corr.pars <- lapply(1:length(cov.pars), function(x) cov.pars[[x]][-1]) 
  #lista dos parâmetros de correlação apenas

  ## Cholesky of each corr_marg
  ct <- lapply(1:ncov.model, function(x) {Matrix::chol(corr_marg(dist.matrix = dist.matrix, 
        cov.model = cov.model[x], corr.pars = cov.pars[[x]][2:length(cov.pars[[x]])]))})
  
  ## Inverse of each Cholesky univariate covariance matrix, 
  # pois inv(chol(sig2*R)) = (1/sqrt(sig2))*inv(chol(R))
  Sigma_inv_chol <- lapply(1:ncov.model, function(x){(1/sqrt(cov.pars[[x]][1]))*spam::solve(ct[[x]])}) 
  
  ## Log sum of cholesky diagonal elements: 
  # pois chol(sig2*R) = sqrt(sig2)*chol(R), onde ct = chol(R)
  Sigma_log_sum_chol_diag <- sum(log(sapply(1:ncov.model, function(x){
                             spam::diag(sqrt(cov.pars[[x]][1])*ct[[x]])})))
  
  prod_y_bdiag <- crossprod(data,bdiag(Sigma_inv_chol))
  
  prod_y_bdiag_kron <- tcrossprod(prod_y_bdiag, kronecker(solve(SigmaB), diag(n)))
  
  nll <- as.numeric(-drop((length(data))*log(2*pi) + 2*Sigma_log_sum_chol_diag + n*log(det(SigmaB)) +
                            tcrossprod(prod_y_bdiag_kron, prod_y_bdiag))/2)
  
  # chol_uni_models <- chol_uni_corr_dif_models(dist.matrix = dist.matrix, cov.model = cov.model,
  #                                             corr.pars = corr.pars, marg.error = marg.error, id_aux1 = id_aux1)
  # ct <- list_chol_uni_corr(chol.uni = chol_uni_models$modelos, cov.model = cov.model, corr.pars = corr.pars)
  # 
  # id_aux1 <- chol_uni_models$id_aux1
  # #print(id_aux1)

  #print(nll)
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

