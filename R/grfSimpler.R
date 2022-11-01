#' @title Simulation of Gaussian random fields with a Simpler Covariance Function
#' @name grfSimpler
#' 
#' @description Simulation of Gaussian random fields with a Simpler Covariance Function
#' 
#' @importFrom matrixcalc is.positive.definite

#' 
#' @export
#' 
grfSimpler <- function(n, nsim, cov.model, cov.pars, SigmaB = NULL)
{
  #### Função para uso externo
  ### n --> numero de locações espaciais
  ## nsim --> número de simulações
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  ## SigmaB --> (opcional) matriz de correlação de dimensão pxp, 
  #  necessária se p=length(cov.model) > 1
  
  p <- length(cov.model) ## número de variáveis
  
  # rseed <- set.seed(...)
  
  ## Marginal-covariance (univariate)
  if(p==1)
  {
    pars_uni <- unlist(cov.pars)
    sig2 <- pars_uni[1]
    Sigma <- sig2*corr_marg(dist.matrix = dist.matrix, cov.model = cov.model[1],
                            corr.pars = pars_uni[2:length(pars_uni)])
  }
  

  return(Sigma)
}
