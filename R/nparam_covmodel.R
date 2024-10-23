#' @title Calculates o número de parâmetros do modelo CovSimpler
#' @name nparam_covmodel
#'
#' @description Retorna um vetor com o número de parâmetros de cada modelo marginal e da função de correlação entre as respostas SigmaB, dados os modelos de correlação marginais das respostas
#' 
#' @usage nparam_covmodel(model)
#' 
#' @param model vetor de modelos de função de correlação para cada resposta do modelo
#' 
#' @details
#' Retorna uma lista contendo 
#' 
#' @examples
#' library(geoR)
#' # Número de parâmetros do modelo de covariância CovSimpler considerando uma resposta com distribuição Matern
#' nparam_covmodel("matern")
#' 
#' # Número de parâmetros do modelo de covariância CovSimpler considerando duas resposta com distribuições Exponential e Matern
#' nparam_covmodel(c("exp","matern"))
#'   
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
