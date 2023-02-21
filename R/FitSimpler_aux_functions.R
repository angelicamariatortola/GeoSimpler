"FitSimpler_aux" <-   function(pars, fp, ip, temp.list)
  ## função de uso interno!! 
  
  ### pars : values for the parameters to be estimated
  ## sequence is c(phi, tausq, kappa, sigmasq)
  ### fp = fixed pars: parameters considered fixed
  ### ip = ind.pars : list indicating which are fixed and which are to be estimated
  ##
  ## Warning:
  ##  if fix.nugget = TRUE and nugget > 0 ,
  ## sigmasq should be passed and fp$nugget is the value of the nugget
  ## otherwise the RELATIVE nugget should be passed
{
  # pbeta <- temp.list$beta.size
  # log.jacobian <- temp.list$log.jacobian
  
  ## Univariado
  ## Obligatory parameter:
  phi <- pars[1]
  
  ## Others
  if(ip$f.tausq){
    if(fp$tausq > 0){
      npars.min <- length(pars)
      sigmasq <- pars[npars.min]
    }  else sigmasq <- 1 # se tau é fixo e nulo, ele atribui 1 para o sigma (fixa o sigma também para depois estimá-lo usando as demais estimativas?)
  } else sigmasq <- 1 # se tau não é fixo, ele atribui 1 para o sigma 
  
  if(ip$f.tausq & ip$f.kappa){ # tau e kappa fixos
    tausq <- fp$tausq
    kappa <- fp$kappa
  }
  
  if(ip$f.tausq & !ip$f.kappa){ # tau fixo e kappa vai ser estimado
    tausq <- fp$tausq
    kappa <- pars[2]
  }
  
  if(!ip$f.tausq & ip$f.kappa){  # tau vai ser estimado e kappa fixo
    tausq <- pars[2]
    kappa <- fp$kappa
  }

  if(!ip$f.tausq & !ip$f.kappa){#ambos serão estimados
    tausq <- pars[2]
    kappa <- pars[3]
  }
  
  ##
  if(temp.list$print.pars){
    running.pars <- c(phi = phi, tausq = tausq, kappa =kappa)
    if(ip$f.tausq && fp$tausq > 0)
      running.pars <- c(sigmasq=sigmasq, running.pars)
    print(running.pars)
  }
  ##
  ## Absurd values
  ##
  if(!is.na(kappa))
  {
    if(kappa < 1e-04 | (tausq+sigmasq) < (.Machine$double.eps^0.5) |
       any(c(phi, tausq, sigmasq, kappa) < 0))
      return(.Machine$double.xmax^0.5)
  }

  ## Computing likelihood
  ## NOTE: Likelihood for Independent observations 
  ##       arbitrary criteria used here:
  ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
  ##
  
  n <- temp.list$n
  z <- temp.list$z
  est_mean <- temp.list$est_mean
  
    if((phi < 1e-16) | (sigmasq < 1e-16)){
      if(ip$f.tausq)
        v <- list(varcov = diag(x=(tausq+sigmasq), n), # a matrix de covariância fica nula fora da diagonal
                  log.det.to.half = (n/2) * log(tausq+sigmasq)) 
      # a cholesky de uma matriz diagonal é a raiz quadrada da diagonal
      else  # o determinante de uma matriz diagonal é o produto da diagonal = (tau2+sig2)^(n/2)
        # o log disso fica a espressão de log.det.to.half
        v <- list(varcov = diag(x=(1+tausq), n),
                  log.det.to.half = (n/2) * log(1+tausq))
    }  else
      v <- varcov.spatial2(coords = temp.list$coords, 
                           cov.model = temp.list$cov.model, kappa=kappa,
                           nugget = tausq, cov.pars = c(sigmasq, phi),
                           det = TRUE)
    
    #if(!is.null(v$crash.parms)) return(.Machine$double.xmax^0.5)
    
    if (est_mean)
    {
      xmat <- temp.list$xmat
      ivx <- solve(v$varcov,xmat)
      xivx <- crossprod(ivx,xmat)
      betahat <- try(.solve.geoR(xivx,crossprod(ivx,z)), silent=TRUE)
      
      if(inherits(betahat, "try-error"))
      { # se der erro no cálculo do beta
        t.ei <- eigen(xivx, symmetric = TRUE)
        betahat <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% crossprod(ivx,z), silent=TRUE)
      }
      
      if(inherits(betahat, "try-error")) # se der erro de novo no cálculo do beta
        stop("Covariates have very different orders of magnitude. Try to multiply and/or divide them to bring them to similar orders of magnitude") 
      
      res <- z - xmat %*% betahat
    } else  res <- z
    

    ssres <- drop(crossprod(res,solve(v$varcov,res))) # e^t inv(Sigma) e
    
    # if(temp.list$method.lik == "ML")
    # {
      if(ip$f.tausq & (tausq > 0))
        negloglik <- v$log.det.to.half +  0.5 * ssres  else
        negloglik <- (n/2) * log(ssres) +  v$log.det.to.half ## essa expressão vem da verossimilhança concentrada 
                                                             ## pois neste caso, o tau deve ser estimado (model based geo - p. 112)
    # }                                                      ## -(1/2)log|V| = -(1/2)log|Sigma|, quando sigma^2 = 1,
                                                             ## neste caso, -(1/2)log|V| = -(1/2)log|Sigma| = -v$log.det.to.half
    
    if(temp.list$method.lik == "RML")
    {
      if(length(as.vector(xivx)) == 1) {
        choldet <- 0.5 * log(xivx)
      }
      else {
        chol.xivx <- chol(xivx)
        choldet <- sum(log(diag(chol.xivx)))
      }
      
      if(ip$f.tausq & (tausq > 0))
        negloglik <- v$log.det.to.half +  0.5 * ssres + choldet
      else
          negloglik <- ((n-p)/2) * log(ssres) +  v$log.det.to.half + choldet
    }  
    
    negloglik <- negloglik - temp.list$loglik.cte
  
  if(negloglik > (.Machine$double.xmax^0.5) | negloglik == Inf | negloglik == -Inf)
    negloglik <- .Machine$double.xmax^0.5
  if(temp.list$print.pars)
    cat(paste("log-likelihood = ", -negloglik, "\n"))
  
  return(negloglik) 
}

vecdist <- function(x){as.vector(dist(x))}

