#' @title Calculates the log-likelihood value for the Simpler Model
#' @name FitSimpler
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom Matrix chol solve diag crossprod tcrossprod bdiag
FitSimpler <- function(data, dist.matrix, cov.model, ini.cov.pars, 
                       SigmaB = NULL, method, hessian, logpars)
{
  
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## ini.cov.pars --> lista de vetores iniciais, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  ## data --> dados (vetor, matriz ou data.frame). Se vetor p = 1, se matriz ou data.frame, p = ncol(data)
  
  if(is.vector(data) || (ncol(data) == 1)){p <- 1}
  if(is.matrix(data) || is.data.frame(data)){p <- ncol(data)}
  
  nloc <- nrow(dist.matrix) ## número de locações amostrais
  ncov.model <- length(cov.model) ## número de variáveis
  nparam <- nparam_covmodel(cov.model) 
  
  ## p deve ser igual a ncol(as.matrix(data))
  
  for (i in 1:length(cov.model))
  {
    if (any(cov.model[i] == c("power", "gneiting.matern", "gencauchy"))) 
    {
      stop(paste("parameter estimation for", cov.model[i], "is not yet implemented"))
    }
    
  }


  
  if(ncov.model == 1)
  {
    if(p == 1)
    {
      # Univariado
      ini_par <- unlist(ini.cov.pars)
      
      time_Simpler <- system.time(est_Simpler <- try(optim(par = ini_par,
                                 fn = loglikSimpler_optim_uni, data = data, 
                                 dist.matrix = dist.matrix, nloc = nloc,
                                 cov.model = cov.model, logpars = logpars, 
                                 control = list(fnscale = -1),
                                 method=method, hessian = hessian)))
    }else
    {
      # else (p > 1) --> multivariado: com o mesmo modelo de covariância para todas as respostas
      #neste caso o nparam é número de parâmetros apenas de 1 modelo, pois p>1, 
      # mas o modelo é o mesmo
      
      ini_par <- c(unlist(ini.cov.pars), gdata::upperTriangle(SigmaB, byrow = T))
      
      time_Simpler <- system.time(est_Simpler <- try(optim(par = ini_par,
                                 fn = loglikSimpler_optim_unicov_pmulti, data = data, 
                                 dist.matrix = dist.matrix, nloc = nloc, nparam = nparam,
                                 cov.model = cov.model, logpars = logpars, p = p,
                                 control = list(fnscale = -1),
                                 method=method, hessian = hessian)))
    }
  }else
  {# else (ncov.model > 1) --> multivariado --> tamanho do cov.model == p
    if(ncov.model == p)
    {
        ini_par <- c(unlist(ini.cov.pars), gdata::upperTriangle(SigmaB, byrow = T))

        time_Simpler <- system.time(est_Simpler <- try(optim(par = ini_par,
                       fn = loglikSimpler_optim_multicov_pmulti, data = data,
                       dist.matrix = dist.matrix, nloc = nloc, nparam = nparam,
                       cov.model = cov.model, logpars = logpars, p = p,
                       control = list(fnscale = -1),
                       method=method, hessian = hessian)))
        
    }else{stop("error: wrong model specification")}
  }

  fits_pars <- est_FitSimpler(est_pars = est_Simpler, p = p, ncov.model = ncov.model)
  
  out <- list(nloc = nloc, p = p,
              elapstime = as.numeric(time_Simpler[3]),
              cov.model = cov.model,
              method = method,
              hessian = hessian, 
              logpars = logpars,
              #sd = sd,
              est_pars = fits_pars)
  
  class(out) <- "FitSimpler"
  return(out)
}
