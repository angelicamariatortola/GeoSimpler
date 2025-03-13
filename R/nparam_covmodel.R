#' @title Calculates o número de parâmetros no modelo
#' @name nparam_covmodel
#'
#' @description Retorna um vetor com o número de parâmetros de cada modelo
#' da função cov.model e dos parâmetros da função de correlação SigmaB
#' 
#' @example 
nparam_covmodel <- function(cov.model)
{
  #### Função para uso interno
  
  np_corr <- NULL
  if(length(cov.model) > 1)
  {
    np_corr <- choose(length(cov.model), 2)
  }
  
  np_marg <- c()
  for(i in 1:length(cov.model))
    {
      if(cov.model[i] %in% c("gencauchy", "gneiting.matern"))
      {
        np_marg[i] <- 4
      }
      if(cov.model[i] %in% c("matern", "cauchy", "powered.exponential"))
      {
        np_marg[i] <- 3
      }
      if(!cov.model[i] %in% c("matern", "cauchy", "powered.exponential",
                              "gencauchy", "gneiting.matern"))
      {
        np_marg[i] <- 2
      }
    }
  
  nparam <- c(np_marg, np_corr)
  return(nparam)
}
