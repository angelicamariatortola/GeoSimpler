#' @title Calculates o número de parâmetros no modelo
#' @name nparam_covmodel
#'
#' @description Retorna um vetor com o número de parâmetros de cada modelo
#' da função cov.model e dos parâmetros da função de correlação SigmaB
#' 
nparam_covmodel <- function(cov.model)
{
  #### Função para uso interno
  
  modelos_param <- list(modelo = c("matern", "exp", "gaussian", "cauchy"),
                        np = c(3, 2, 2, 3))
  
  np_corr <- NULL
  if(length(cov.model) > 1)
  {
    np_corr <- choose(length(cov.model), 2)
  }
  
  np_marg <- c()
  for(i in 1:length(cov.model))
    {
      np_marg[i] <- modelos_param$np[which(cov.model[i] == modelos_param$modelo)]
  }
  
  nparam <- c(np_marg, np_corr)
  return(nparam)
}
