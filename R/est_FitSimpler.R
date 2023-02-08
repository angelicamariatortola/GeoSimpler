est_FitSimpler <- function(est, p, ncov.model, nparam, logpars, hessian){
  if(logpars && hessian)
  {
    var_par <- diag(-solve(est$hessian))
    deriv_param <- exp(est$par[sum(nparam[1:(length(nparam)-1)])])
    # sum(nparam[1:(length(nparam)-1)]) --> soma o número de parâmetros dos modelos marginais
    
    # chol2inv(chol(x)) = solve(x) 
    
    if(p > 1)
    {
      
    }
   
    fits_pars <- rbind(c(exp(est$par[1:(3*nvar)]), 
                         tanh(est$par[(3*nvar+1):length(est$par)])),
                       sqrt((deriv_param^2)*var_par))
    
    
    if(ncov.model == 1)
    {
      if(p == 1)
      {
        
        
      }else
      {
        
        
        
      }
    }
  }

  
  if(x$logpars)
  {
    est <- c(exp(x$est_pars[1:(3*nvar)]), 
             tanh(est_MatSimpler$par[(3*nvar+1):length(est_MatSimpler$par)]))
  }

  }