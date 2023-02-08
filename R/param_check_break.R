param_check_break <- function (cov.pars, cov.model, p) 
{
  #### Função para uso interno
  ## verifica se o número de parâmetros está correto para o respectivo modelo
  # e separa os respectivos parâmetros: sigmap, phip, kappap,...
  # cov.pars --> lista
  # falta verificar se a dimensão de cada vetor cov.pars 
  # está de acordo com o repectivo modelo de cov.model
  
  ncov.model <- length(cov.model)
  ncov.pars <- length(cov.pars)
  
  if(is.list(cov.pars))
  {
    if((ncov.model == 1) && (ncov.pars > 1))
      stop("cov.pars must be a list with one vector")
    
    if((ncov.model > 1) && (ncov.model != p))
      stop("p must be equal to the size of the cov.model vetor. 
           Ex. If cov.model = \"matern\" then p=1.
               If cov.model = c(\"matern\", \"gaussian\") then p=2, etc.")
      
    if((ncov.model > 1) && (ncov.model == p) && (ncov.model != ncov.pars))
      stop("cov.pars must be a list of p vector of parameters, 
           one for each model of the cov.model vector")
  } else
  {
    stop("cov.pars must be a list of p vectors")
  }
  
  sigmap <- sapply(1:ncov.model, function(x) {cov.pars[[x]][1]})
  phip   <- sapply(1:ncov.model, function(x) {cov.pars[[x]][2]})
  kappap <- list()
  
  for(i in 1:ncov.model)
  {
    if(cov.model[i] %in% c("matern", "powered.exponential", 
                           "cauchy", "gencauchy", "gneiting.matern"))
    {
      kappap[[i]] <- cov.pars[[i]][3:length(cov.pars[[i]])]
    } else
    {
      kappap[[i]] <- NA
    }
  }
  
  return(list(sigmap = sigmap, phip = phip, kappap = kappap))
      
}
  
