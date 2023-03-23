"FitSimpler_aux" <- function(pars, fp, ip, temp.list)
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
  
  if(temp.list$ncov.model == 1)
  {     ## Cenario univariado
    if(temp.list$p == 1)
    {
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
        if(kappa < 1e-04 || (tausq+sigmasq) < (.Machine$double.eps^0.5) ||
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
      fix.mean <- temp.list$fix.mean
      
        if((phi < 1e-16) || (sigmasq < 1e-16)){
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
        
        if (!fix.mean)
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
          if(length(as.vector(xivx)) == 1) 
            {
              choldet <- 0.5 * log(xivx)
            } else 
              {
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
        
    } else
    {
      ## Multivariado
      ## Obligatory parameter:
      phi <- pars[1]
      nrhos <- temp.list$nrhos
      
      ## Others
      if(ip$f.tausq){
        if(fp$tausq > 0){
          npars.min <- length(pars)
          sigmasq <- pars[npars.min - nrhos] #ultima posição que antecede os rhos
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
      
      rhos <- pars[(length(pars)-nrhos+1):length(pars)]
      
      # phi = 0.201; sigmasq = 1; tausq = 0; rhos <- c(0.4,0.9,0.6)
      # print(c(phi, sigmasq, tausq, rhos))
      
      m1 <- matrix(1, nc = p, nr = p)
      m1[lower.tri(m1)] <- rhos
      SigmaBq <- t(m1)
      SigmaBq[lower.tri(SigmaBq)] <- rhos
      
      if(matrixcalc::is.positive.definite(SigmaBq))
      {
        ##
        if(temp.list$print.pars){
          running.pars <- c(phi = phi, tausq = tausq, kappa =kappa, rhos = rhos)
          if(ip$f.tausq && fp$tausq > 0)
            running.pars <- c(sigmasq=sigmasq, running.pars)
          print(running.pars)
        }
        ##
        ## Absurd values
        ##
        if(!is.na(kappa))
        {
          if(kappa < 1e-04 || (tausq+sigmasq) < (.Machine$double.eps^0.5) ||
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
        fix.mean <- temp.list$fix.mean
        
        ## até aqui igual o univariado
        
        ###########
        
        if((phi < 1e-16) || (sigmasq < 1e-16))
        {
          if(ip$f.tausq)
          {
            v <- list(varcov = diag(x=(tausq+sigmasq), n),
                      # a matrix de covariância fica nula fora da diagonal
                      sqrt.varcov = sqrt(tausq+sigmasq)*diag(n),
                      sqrt.inverse = (1/sqrt(tausq+sigmasq))*diag(n),
                      log.det.to.half = (n/2) * log(tausq+sigmasq)) 
            # Considerando que, se phi<1e-16, a varcov é formada por rhoij*(tausq+sigmasq) nas diagonais das suas submatrizes
            # e zero no restante. Deste modo, det(chol(varcov)) = det(chol(kronecker((tausq+sigmasq)*SigmaB, diag(n))))
            # = det(kronecker(chol((tausq+sigmasq)*SigmaB), diag(n))), pois chol(kronecker(A, I)) = kronecker(chol(A), I)
            # = det(kronecker(sqrt(tausq+sigmasq)*chol(SigmaB), diag(n))), pois chol(a²A) = a*chol(A)
            # = det((sqrt(tausq+sigmasq)*chol(SigmaB))^n), pois det(A ⊗ B) = det(A)^(rank(B))*det(B)^(rank(A))
            # = det(sqrt(tausq+sigmasq)*chol(SigmaB))^n, pois det(A^n) = det(A)^n
            # aplicando o log, obtemos: log(det(chol(varcov))) = n*log(det(sqrt(tausq+sigmasq)*chol(SigmaB)))
            # = n*log(prod(sqrt(tausq+sigmasq)*diag(chol(SigmaB))))
            # = n*sum(log(sqrt(tausq+sigmasq)*diag(chol(SigmaB))))
            # o determinante de uma matriz diagonal é o produto da diagonal.
            # o log disso fica a espressão de log.det.to.half
          } else  
          {
            # Se tau não é fixo (ou seja, será estimado), ele faz a a divisão por sigmasq, assim tausq = nugget/sigmasq (confirmar !!!)
            v <- list(varcov = diag(x=(1+tausq), n),
                      sqrt.varcov = sqrt(tausq+1)*diag(n),
                      sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
                      log.det.to.half = (n/2) * log(1+tausq))
          }
        } else
        {
          # covariancia univariada
          # Neste caso, se tau não for fixo (será estimado) ou for igual a zero, já foi considerado como valores iniciais:
          # tau = nugget/params$sigmap e sigma2 = 1. Neste caso, Sigma11 = sigmasq*v$varcov (v$varcov é como se dividisse a cov univariada pela variancia)
          # Se tau é fixo ele entra com o valor do tau fornecido no argumento nugget
          # Já está fazendo os 2 casos
          
          #print(c(phi, sigmasq, rhos))
          v <- varcov.spatial2(coords = temp.list$coords, 
                               cov.model =  temp.list$cov.model, kappa = kappa, 
                               nugget = tausq, cov.pars = c(sigmasq, phi), det = T)
        }

        kroninv <- kronecker(solve(SigmaBq), solve(v$varcov))
        
        # 0.5*kronecker(SigmaBq, v$varcov) =  CovSimpler(coords = coord, cov.model = cov.model, cov.pars = cov.pars, 
                                              # nugget = nugget, p = 2, SigmaB = SigmaB)
        
        #loglikSimpler(data = data, coords = coord, cov.model = cov.model, cov.pars = cov.pars,  SigmaB = SigmaB)
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, nugget = 0,  
        #                   cov.pars = cov.pars, p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(0, length(z)), sigma = s22$varcov, log = T) #  --> fix.mean = T
        
        # print(c(phi, sigmasq, tausq, rhos))
        
        if (!fix.mean)
        {
          D <- temp.list$D
          kroninv_x <- crossprod(kroninv, D)
          tx_kroninv_x <- as.matrix(crossprod(kroninv_x, D)) # t(x)%*%kronecker(inv(Sigmab), inv(V))%*%x
          betahat <- try(.solve.geoR(tx_kroninv_x, as.matrix(crossprod(kroninv_x, z))), silent=TRUE)
          
          if(inherits(betahat, "try-error"))
          { # se der erro no cálculo do beta
            t.ei <- eigen(tx_kroninv_x, symmetric = TRUE)
            betahat <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% crossprod(kroninv_x,z), silent=TRUE)
          }
          if(inherits(betahat, "try-error") || any(is.nan(as.numeric(betahat)))) # se der erro de novo no cálculo do beta
          {
            stop("Covariates have very different orders of magnitude. Try to multiply and/or divide them to bring them to similar orders of magnitude")
          } else{ res <- z - D %*% betahat} # se não der erro no betahat ele calcula o res
        } else  res <- z
        
        ssres <- as.numeric(drop(crossprod(res, kroninv%*%res))) # e^t %*% kronecker(inv(Sigmab), inv(V)) %*% e
        # = t(res)%*%kroninv%*%(res)
        
        # if(temp.list$method.lik == "ML")
        # {
        if(ip$f.tausq & (tausq > 0))
          negloglik <- p*(v$log.det.to.half) + (n/2)*log(det(SigmaBq)) + 0.5 * ssres  else
          # N = n*p
          # negloglik = -0.5*(N*log(2*pi)+N*log(sig²)+n*log(det(SigmaB)) + p*log(det(v$varcov)) + ssres/sig²)
            
            negloglik <- 0.5*(n*p*log(ssres) +  n*log(det(SigmaBq))) + p*(v$log.det.to.half)
            # Aqui estima sig² por ssres/N:
            # negloglik = -0.5*(N*log(2*pi)+N*log(ssres/N)+n*log(det(SigmaB)) + p*log(det(v$varcov)) + N)
        
        ## essa expressão vem da verossimilhança concentrada 
        ## pois neste caso, o tau deve ser estimado (model based geo - p. 112)                                                      
        ## (1/2)p log|V| = (1/2)p log|R'R| = (1/2)p log|R|^2  = 2(1/2)p log|R| = p log|R| = p * v$log.det.to.half
        
        negloglik <- as.numeric(negloglik - temp.list$loglik.cte)
        
        # Se tau é fixo e igual a zero - fix.mean = T: betahat = 0
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, cov.pars = list(c(ssres/(n*p), phi)), p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(0, length(z)), sigma = s22$varcov, log = T) #  --> fix.mean = T
        
        # Se tau é fixo e igual a zero - fix.mean = F: estima betahat
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, cov.pars = list(c(ssres/(n*p), phi)), p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(betahat, each = nrow(coord)), sigma = s22$varcov, log = T)
        
        # Se tau é fixo e maior do que zero - fix.mean = T: betahat = 0
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, cov.pars = list(c(sigmasq, phi)), # estima o sigmasq no optim
        #                   nugget = nugget, p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(0, length(z)), sigma = s22$varcov, log = T)
        
        # Se tau é fixo e maior do que zero - fix.mean = F: estima betahat
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, cov.pars = list(c(sigmasq, phi)), # estima o sigmasq no optim
        #                   nugget = nugget, p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(betahat, each = nrow(coord)), sigma = s22$varcov, log = T)
        
        ####################
        
        # Se tau NÃO é fixo (qualquer valor inicial) - fix.mean = T: betahat = 0 --> ele estima o v = tau²/sig² no optim
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, nugget = tausq*(ssres/(n*p)),
        #                   cov.pars = list(c(ssres/(n*p), phi)), p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(0, length(z)), sigma = s22$varcov, log = T) #  --> fix.mean = T
        # ## ou
        # s22 <- (ssres/(n*p))*kronecker(SigmaBq, v$varcov)
        # mvtnorm::dmvnorm(z, mean = rep(0, length(z)), sigma = s22, log = T) #  --> fix.mean = T

        # Se tau NÃO é fixo (qualquer valor inicial) - fix.mean = F: estima betahat.  Estima o v = tau²/sig² no optim
        # s22 <- CovSimpler(coords = coord, cov.model = cov.model, nugget = tausq*(ssres/(n*p)),
        #                   cov.pars = list(c(ssres/(n*p), phi)), p=2, SigmaB = SigmaB)
        # mvtnorm::dmvnorm(z, mean = rep(betahat, each = nrow(coord)), sigma = s22$varcov, log = T) #  --> fix.mean = T
        # ## ou
        # s22 <- (ssres/(n*p))*kronecker(SigmaBq, v$varcov)
        # mvtnorm::dmvnorm(z, mean = rep(betahat, each = nrow(coord)), sigma = s22, log = T) #  --> fix.mean = T
        
        
        if(negloglik > (.Machine$double.xmax^0.5) | negloglik == Inf | negloglik == -Inf)
          negloglik <- .Machine$double.xmax^0.5
    
      } else {  negloglik <- .Machine$double.xmax^0.5}
    } 

    if(temp.list$print.pars)
      cat(paste("log-likelihood = ", -negloglik, "\n"))
    
  }
  
  return(negloglik) 
}

vecdist <- function(x){as.vector(dist(x))}

