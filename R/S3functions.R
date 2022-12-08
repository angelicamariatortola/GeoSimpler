#' @export
print.FitSimpler <- function(x, ...)
{
  
  var <- as.list(match.call())
  if(!("digits" %in% names(var))) 
  {digits <- 4}else{
      digits <- var$digits}
  
  est_round <- round(x$est_pars, digits)
  
  cat("number of locations:", x$nloc,"\n")
  cat("number of response variables:", x$p,"\n")
  cat("cov.model:", x$cov.model,"\n")
  
  ncov.model <- length(x$cov.model) ## número de variáveis
  nparam <- nparam_covmodel(x$cov.model)
  
  if(ncov.model == 1)
  {
    cat("var: ", est_round[1],"\n")
    cat("phi: ", est_round[2],"\n")
    
    # univariado --> p = 1
    if(x$cov.model %in% c("matern", "cauchy", "powered.exponential"))
    {
      cat("kappa: ", est_round[3],"\n")
    }
    if(x$cov.model %in% c("gencauchy", "gneiting.matern"))
    {
      cat("kappa1: ", est_round[3],"\n")
      cat("kappa2: ", est_round[4],"\n")
    }
    
    if(x$p > 1)
    { # multivariado --> p > 1 : considera o mesmo modelo para todas as respostas
      ncomb <- paste(combn(x$p,2)[1,], combn(x$p,2)[2,], sep = "")
      est_rho <- est_round[(nparam+1):length(est_round)]
      
      for (i in 1:length(ncomb)) 
      {
        cat(paste("rho", ncomb[i], ": ", est_rho, sep = ""),"\n")
      }
    }
  } else # ncov.model > 1 e ncov.model = p --> modelos diferentes para as respostas
  {
    est_list <- split(est_round, rep(1:length(nparam), nparam))
    
    for (i in 1:x$p) 
    {
      if(x$cov.model[i] %in% c("matern", "cauchy", "powered.exponential"))
      {
        cat(paste("var", i, ": ", est_list[[i]][1], sep = ""),"\n")
        cat(paste("phi", i, ": ", est_list[[i]][2], sep = ""),"\n")
        cat(paste("kappa", i, ": ", est_list[[i]][3], sep = ""),"\n")
      }
      if(x$cov.model[i] %in% c("gencauchy", "gneiting.matern"))
      {
        cat(paste("var", i, ": ", est_list[[i]][1], sep = ""),"\n")
        cat(paste("phi", i, ": ", est_list[[i]][2], sep = ""),"\n")
        cat(paste("kappa1", i, ": ", est_list[[i]][3], sep = ""),"\n")
        cat(paste("kappa2", i, ": ", est_list[[i]][4], sep = ""),"\n")
      }
      if(!x$cov.model[i] %in% c("matern", "cauchy", "powered.exponential",
                                "gencauchy", "gneiting.matern"))
      {
        cat(paste("var", i, ": ", est_list[[i]][1], sep = ""),"\n")
        cat(paste("phi", i, ": ", est_list[[i]][2], sep = ""),"\n")
      }
    }
    
    ncomb <- paste(combn(x$p,2)[1,], combn(x$p,2)[2,], sep = "")
    
    for (i in 1:length(ncomb)) 
    {
      cat(paste("rho", ncomb[i], ": ", 
                est_list[[length(est_list)]][i], sep = ""),"\n")
    }
  }
}


#' @export
summary.FitSimpler <- function(object, ...)
{
  var <- as.list(match.call())
  if(!("digits" %in% names(var))) 
  {digits <- 4}else{
    digits <- var$digits}
  
  cat('\n')
  cat('---------------------------------------------------\n')
  cat('  Likelihood Simpler Model Estimates \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("number of locations:", (object$nloc),"\n")
  cat("number of response variables:", object$p,"\n")
  cat("elapsed time:", object$elapstime,"\n")
  cat("cov.model:", object$cov.model,"\n")
  cat("method:", object$method,"\n")
  cat('------------------\n')
  cat('Paramter Estimates\n')
  cat('------------------\n')
  cat('\n')

  est_round <- round(object$est_pars, digits)
  ncov.model <- length(object$cov.model) 
  lab <- c()
  
  if(ncov.model == 1)
  {
    # univariado --> p = 1
    lab <- c("var", "phi") # labels para os parâmetros
    
    if(object$cov.model %in% c("matern", "cauchy", "powered.exponential"))
    {  lab <- append(lab, "kappa")}
    
    if(object$cov.model %in% c("gencauchy", "gneiting.matern"))
    {  lab <- append(lab, c("kappa1", "kappa2"))}
    
    nparam <- nparam_covmodel(object$cov.model)
    
    if(object$p > 1)
    { # multivariado --> p > 1 : considera o mesmo modelo para todas as respostas
      ncomb <- paste(combn(object$p,2)[1,], combn(object$p,2)[2,], sep = "")
      #est_rho <- est_round[(nparam+1):length(est_round)]
      for (i in 1:length(ncomb)) 
      {
        lab_rho <- paste("rho", ncomb[i], sep = "")
      }
      lab <- append(lab, lab_rho)
    }
  } else # ncov.model > 1 e ncov.model = p --> modelos diferentes para as respostas
  {
    for (i in 1:object$p) 
    {
      if(object$cov.model[i] %in% c("matern", "cauchy", "powered.exponential"))
      {
        lab <- append(lab, c(paste("var", i, sep = ""), 
                             paste("phi", i, sep = ""), 
                             paste("kappa", i, sep = ""))) # labels para os parâmetros
      }
      if(object$cov.model[i] %in% c("gencauchy", "gneiting.matern"))
      {
        lab <- append(lab, c(paste("var", i, sep = ""), 
                             paste("phi", i, sep = ""), 
                             paste("kappa1", i, sep = ""),
                             paste("kappa2", i, sep = ""))) # labels para os parâmetros
      }
      if(!object$cov.model[i] %in% c("matern", "cauchy", "powered.exponential",
                                     "gencauchy", "gneiting.matern"))
      {
        lab <- append(lab, c(paste("var", i, sep = ""), 
                             paste("phi", i, sep = "")))
      }
    }
    ncomb <- paste(combn(object$p,2)[1,], combn(object$p,2)[2,], sep = "")
    #est_rho <- est_round[(nparam+1):length(est_round)]
    for (i in 1:length(ncomb)) 
    {
      lab_rho <- paste("rho", ncomb[i], sep = "")
    }
    lab <- append(lab, lab_rho)
  }
  print(lab)
}

 

#                              
#   # tab <- matrix(NA, ncol=l, nrow=2)
#   # tab[1,] = ifelse(object$logpars, 
#   #                  rep(NA, l), round(object$sd, digits))
#   # t(round(est,digits))
#   # tab[2,] = ifelse(object$hessian, 
#   #                  rep(NA, l), round(object$sd, digits))
#   # colnames(tab)=t(lab)
#   # rownames(tab)=c("Estimate","Standard Error")

