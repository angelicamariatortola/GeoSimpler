#' @title Calculates the log-likelihood value for the Simpler Model
#' @name FitSimpler
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom Matrix chol solve diag crossprod tcrossprod bdiag
FitSimpler <- function(data, dist.matrix, cov.model, cov.pars, p,
                       SigmaB = NULL, method, hessian, logpars)
{
  
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  
  nloc <- nrow(dist.matrix) ## número de locações amostrais
  ncov.model <- length(cov.model) ## número de variáveis
  nparam <- nparam_covmodel(cov.model) #número de parâmetros apenas de 1 modelo, pois p>1, mas o modelo é o mesmo
  
  ## p deve ser igual a ncol(as.matrix(data))
  
  if(ncov.model == 1 && (p == 1))
  {
    ini_par <- unlist(cov.pars)
    
    time_Simpler <- system.time(est_Simpler <- try(optim(par = ini_par,
                     fn = loglikSimpler_optim_uni, data = data, 
                     dist.matrix = dist.matrix, nloc = nloc,
                     cov.model = cov.model, logpars = logpars, 
                     control = list(fnscale = -1),
                     method=method, hessian = hessian)))
    return(list(time_Simpler = time_Simpler, est_Simpler = est_Simpler))
  }
  
  if(ncov.model == 1 && (p > 1)) ## ncol(data) > 1
  {
    ini_par <- c(unlist(cov.pars), gdata::upperTriangle(SigmaB, byrow = T))
  
    time_Simpler <- system.time(est_Simpler <- try(optim(par = ini_par,
                     fn = loglikSimpler_optim_unicov_pmulti, data = data, 
                     dist.matrix = dist.matrix, nloc = nloc, nparam = nparam,
                     cov.model = cov.model, logpars = logpars, p = p,
                     control = list(fnscale = -1),
                     method=method, hessian = hessian)))
    return(list(time_Simpler = time_Simpler, est_Simpler = est_Simpler))
  }
  
  if((ncov.model > 1) && (p > 1) && (ncov.model == p)) ## ncol(data) > 1
  {
    ini_par <- c(unlist(cov.pars), gdata::upperTriangle(SigmaB, byrow = T))

    time_Simpler <- system.time(est_Simpler <- try(optim(par = ini_par,
                     fn = loglikSimpler_optim_multicov_pmulti, data = data, 
                     dist.matrix = dist.matrix, nloc = nloc, nparam = nparam,
                     cov.model = cov.model, logpars = logpars, p = p,
                     control = list(fnscale = -1),
                     method=method, hessian = hessian)))
    
    if(logpars)
    {
      est_Simpler$par <- c(exp(est_Simpler$par[1:sum(nparam[1:p])]), 
                           tanh(est_Simpler$par[(sum(nparam[1:p])+1):length(est_Simpler$par)]))
    }
    return(list(time_Simpler = time_Simpler, est_Simpler = est_Simpler))
  }
  return("error: wrong model specification")
}
