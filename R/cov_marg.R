#' @title Calculates Value of the Marginal-Correlation Function
#' @name cov_marg
#'
#' @description Computes a marginal-correlation function,
#' given the vector of parameters and the distance of their locations.
#' Different correlation functions are available.
#'
#' @importFrom geoR matern
#' 
cov_marg <- function(coords = NULL, dists.lowertri = NULL, nugget,
                     cov.model, cov.pars.uni, log.det.chol = FALSE,
                     inv.chol = FALSE)
{
  #### Função para uso interno
  ## calcula matrix de covariância univariada
  ## cov.model --> modelos de correlação da função geoR::cov.spatial
  ## cov.pars.uni --> vetor com os parâmetros da covariancia univariada: c(sig2, phi, kappa) - modelo Matern
  ## coords ou dists.lowertri devem ser dados --> igual em varcov.spatial
  ## Se log.det.chol = T calcula o log do determinante da decomposição de cholesky da varcov
  ## Se inv.chol = T calcula a inversa da decomposição de cholesky da varcov
  
  if (is.null(coords) & is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) & !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  # if (!is.null(coords)) 
  #   n <- nrow(coords)
  # if (!is.null(dists.lowertri)) 
  #   n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  
  if(cov.model %in% c("matern", "powered.exponential", 
                      "cauchy", "gencauchy", "gneiting.matern"))
  {
    cov.pars <- as.vector(cov.pars.uni[1:2])
    kappa <- as.vector(cov.pars.uni[3:length(cov.pars.uni)])
    
    varcov_results <- varcov.spatial(cov.model = cov.model, coords = coords,
                             dists.lowertri = dists.lowertri,
                             cov.pars = cov.pars, kappa = kappa, nugget = nugget,
                             func.inv = "cholesky", sqrt.inv = inv.chol,
                             try.another.decomposition = FALSE,
                             det = log.det.chol)
  }else
  {
    cov.pars <- as.vector(cov.pars.uni)
    varcov_results <- varcov.spatial(cov.model = cov.model, coords = coords,
                             dists.lowertri = dists.lowertri,
                             cov.pars = cov.pars, nugget = nugget,
                             func.inv = "cholesky", sqrt.inv = inv.chol,
                             try.another.decomposition = FALSE,
                             det = log.det.chol)
  }

  return(varcov_results)
}
