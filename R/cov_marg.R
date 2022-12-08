#' @title Calculates Value of the Marginal-Correlation Function
#' @name cov_marg
#'
#' @description Computes a marginal-correlation function,
#' given the vector of parameters and the distance of their locations.
#' Different correlation functions are available.
#'
#' @importFrom geoR matern
#' 
cov_marg <- function(dist.matrix, cov.model, cov.pars.uni)
{
  #### Função para uso interno
  ## calcula matrix de covariância univariada
  ## cov.model --> modelos de correlação da função geoR::cov.spatial
  ## cov.pars.uni --> vetor com os parâmetros da covariancia univariada: c(sig2, phi, kappa) - modelo Matern
  
  if(cov.model %in% c("matern", "powered.exponential", 
                      "cauchy", "gencauchy", "gneiting.matern"))
  {
   cov_fun <- geoR::cov.spatial(obj = dist.matrix, cov.model = cov.model,
                cov.pars = cov.pars.uni[1:2],
                kappa = cov.pars.uni[3:length(cov.pars.uni)])
  }else
  {
    cov_fun <- geoR::cov.spatial(obj = dist.matrix, cov.model = cov.model,
                                 cov.pars = cov.pars.uni)
  }

  return(cov_fun)
}
