#' @title Calculates the log-likelihood value for the Simpler Model
#' @name loglikSimpler
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom mvtnorm dmvnorm
#' 
#' @export
#' 
loglikSimpler <- function(data, coords = NULL, dists.lowertri = NULL, 
                          cov.pars, cov.model, nugget = rep(0, length(cov.model)), 
                          est_mean = F, trend = rep("cte", length(cov.model)), SigmaB = NULL)
{
  ## função de uso externo!! 
  
  ## data --> dados (vetor, matriz ou data.frame). Se vetor p = 1, 
  ## se matriz ou data.frame, p = ncol(data)
  ## dist.matrix --> matriz de distâncias euclidianas - dist(coords)
  ## cov.pars --> lista com os parametros da matriz de Covariância - mesmo de CovSimpler.
  ## cov.model --> vetor de strings contendo o modelo da variável
  ## nugget = tau² --> deve ser um vetor de dimensão p, com p = número de respostas
  ## est_mean --> Se est_mean = T, o betahat será estimado pelos dados. 
  ## Se est_mean = F, o betahat é igual a zero e não é estimado.
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## trend deve ser uma lista com as especificações para a média de cada resposta. 
    ## o default é "cte" para todas as respostas.
  
  
  if (is.null(coords) && is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) && !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if(is.null(cov.pars))
    stop("cov.pars must be provided")
  
  if (!is.null(coords))
    n <- nrow(coords)
  if (!is.null(dists.lowertri))
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  
  if((is.vector(data) && (length(data) != n)) ||
     ((is.matrix(data) || is.data.frame(data)) && (nrow(data) != n)) )
    stop("coords and data have incompatible sizes")
  
  if(is.vector(data) || (ncol(data) == 1)){p <- 1}
  if(is.matrix(data) || is.data.frame(data)){p <- ncol(data)}
  
  if ((p>1) && is.null(SigmaB)) 
    stop("for p>1, SigmaB must be a positive definite matrix of dimension pxp")
  if(!is.null(SigmaB) && !matrixcalc::is.positive.definite(SigmaB))
    stop("SigmaB must be a positive definite matrix")
  if (!is.null(SigmaB) && (nrow(SigmaB)!= ncol(SigmaB))) 
    stop("SigmaB must be a SQUARE positive definite matrix of dimension pxp")
  if ((p>1) && !is.null(SigmaB) && (nrow(SigmaB)!=p)) 
    stop("SigmaB must be a positive definite matrix of dimension pxp")
  
  if(missing(trend)){trend <- rep(list("cte"), p)}
  if(!is.list(trend) || length(trend) != p)
      stop("trend must be a list of dimension p")
  
     
  ncov.model <- length(cov.model)
  
  if (length(nugget) != ncov.model) 
    stop("the length of the nugget vector must be equal to the length of the cov.model vector")
  
  params <- param_check_break(cov.pars = cov.pars, cov.model = cov.model, p = p)
  y <- as.numeric(data) #empilha os dados para todas as respostas
  
  
  if(ncov.model == 1)
  {     ## Cenario univariado 
    if(p == 1)
    {
      V <- varcov.spatial2(coords = coords, dists.lowertri = dists.lowertri, 
                           cov.model = cov.model, kappa = unlist(params$kappap), nugget = nugget, 
                           cov.pars = c(params$sigmap, params$phip), det = TRUE, sqrt.inv = FALSE)
  
      if(est_mean)
      {
        xmat <- unclass(trend.spatial(trend = trend[[1]], 
                                      geodata = list(coords = coords, data = data)))
        #beta.size <- ncol(xmat)
        ivx <- solve(V$varcov, xmat) # = solve(varcov)%*% xmat, onde xmat = D no caso de tendencia constante
        xivx <- crossprod(ivx, xmat) # = crossprod(solve(varcov)%*% D, D)  --> t(solve(varcov)%*% D) %*% D
        betahat <- .solve.geoR(xivx, crossprod(ivx, data)) 
        betahat_est <- as.numeric(betahat)
        
        #D <- as.matrix(rep(1, length = length(y))) ## 
        # = t(D) %*% solve(varcov) %*% D
        # crossprod(x, y)  --> t(x) %*% y 
        # crossprod(ivx, data[[i]]) = crossprod(solve(varcov)%*% D, y) = t(D) %*% solve(varcov) %*% y
        # .solve.geoR --> função interna do geoR que faz o solve mas dá os tratamentos adequados aos dados.
        # solve(x,y) = as.vector(solve(x)%*%y)
        res <- y - xmat %*% betahat # = (y-D %*% betahat)
      } else
      {
        res <- y # (considera betahat = 0)
      }
      
      ssres <- drop(crossprod(res, solve(V$varcov, res))) # = t(y-D %*% betahat)%*% solve(varcov) %*% (y-D %*% betahat)
      
      negloglik <- (n/2) * (log(2 * pi)) + V$log.det.to.half + 0.5 * ssres # method.lik == "ML"
      # mvtnorm::dmvnorm(as.numeric(data), mean=rep(0, n*p), sigma=Sigma2$varcov, log=TRUE)
      
    } else
      {# else (p > 1) --> multivariado: mesmo modelo de covariância para todas as respostas
      # calcula a cholesky apenas uma vez
      # entra com uma lista de 1 vetor de parâmetros, que é o mesmo para todas as variáveis.
      results <- varcov.spatial2(coords = coords, dists.lowertri = dists.lowertri, 
                      cov.model = cov.model, kappa = unlist(params$kappap), 
                      nugget = nugget, cov.pars = c(params$sigmap, params$phip), 
                      det = TRUE, sqrt.inv = TRUE)
      
      ## Bdiag of Inverse of each Cholesky matrix
      Q <- bdiag(lapply(1:p, function(x){t(results$sqrt.inverse)})) # sqrt.inverse = solve(chol(results$varcov))
      # Q deve ser triangular inferior 
      kron_prod <- kronecker(SigmaB, diag(n))
      invSigma <- crossprod(Q, solve(kron_prod, Q)) # = tQ %*% inv(kron_prod) %*% Q

      if(est_mean)
      {
        ## Aqui trend deve ser uma lista de dimensão p
        ## data deve ser matrix ou data.frame, p = ncol(data)
        ## considerando que trend pode ser diferente para cada variável (faz sentido??)
        
        xmat <- lapply(1:p, function(x){unclass(trend.spatial(trend = trend[[x]],
                                                              geodata = list(coords = coords,
                                                                             data = data[,x])))})
        beta.size <- sapply(1:p, function(x){ncol(xmat[[x]])})
        D <- bdiag(xmat) # agrupa as matrizes de covariáveis em uma matriz bloco diagonal (bonat (2020))
        
        tD_invSigma <- as.matrix(crossprod(D, invSigma)) # # = tD %*% invSigma
        betahat <- .solve.geoR(tD_invSigma%*%D, tD_invSigma%*%y)
        # = (tD %*% invSigma %*% D)^-1 %*% (tD %*% invSigma %*% D)
        
        betahat_est <- split(as.numeric(betahat), rep(1:length(beta.size), beta.size))
        names(betahat_est) <- paste("response", 1:p, sep = "")
        
        res <- y - D %*% betahat # = (y-D %*% betahat)
      } else {res <- y}
     
      # solve((t(D)%*%solve(Sigma2$varcov)%*%D), (t(D)%*%solve(Sigma2$varcov)%*%data))
      ## Log sum of cholesky diagonal elements (soma a diagonal de uma matriz e multiplica por p)
      Sigma_log_sum_chol_diag <- p*(results$log.det.to.half)
      ssres <- drop(crossprod(res, invSigma%*%res)) # = t(y-D %*% betahat)%*% solve(varcov) %*% (y-D %*% betahat)
      
      negloglik <- as.numeric(drop((n*p/2)*log(2*pi) + Sigma_log_sum_chol_diag + 
                               (n/2)*log(det(SigmaB)) + 0.5*ssres))
      # mvtnorm::dmvnorm(as.numeric(data), mean=rep(0, N), sigma=Sigma2$varcov, log=TRUE)
      }
    } else
      {# else (ncov.model > 1) --> multivariado --> ncov.model == p
        if(ncov.model == p)
        {# Cenario multivariado (modelos distintos de covariância para as respostas)
         # calcula a cholesky das covariancias marginais para as p variáveis
            
          results <- lapply(1:ncov.model, function(x) 
          {
            varcov.spatial2(coords = coords, dists.lowertri = dists.lowertri, 
                            kappa = params$kappap[[x]], nugget = nugget[x],
                            cov.model = cov.model[x], det = TRUE, sqrt.inv = TRUE,
                            cov.pars = c(params$sigmap[x], params$phip[x]))
           })         
          
          ## Bdiag of Inverse of each Cholesky matrix
          Q <- bdiag(lapply(1:p, function(x){t(results[[x]]$sqrt.inverse)})) # sqrt.inverse = solve(chol(results$varcov))
          # Q deve ser triangular inferior 
          kron_prod <- kronecker(SigmaB, diag(n))
          invSigma <- crossprod(Q, solve(kron_prod, Q)) # = tQ %*% inv(kron_prod) %*% Q
          
          if(est_mean)
          {
            ## Aqui trend deve ser uma lista de dimensão p
            ## data deve ser matrix ou data.frame, p = ncol(data)
            
            xmat <- lapply(1:p, function(x){unclass(trend.spatial(trend = trend[[x]],
                                                                  geodata = list(coords = coords,
                                                                                 data = data[,x])))})
            beta.size <- sapply(1:p, function(x){ncol(xmat[[x]])})
            D <- bdiag(xmat) # agrupa as matrizes de covariáveis em uma matriz bloco diagonal (bonat (2020))
            tD_invSigma <- as.matrix(crossprod(D, invSigma)) # # = tD %*% invSigma
            
            betahat <- .solve.geoR(as.matrix(tD_invSigma%*%D), tD_invSigma%*%y)
            # = (tD %*% invSigma %*% D)^-1 %*% (tD %*% invSigma %*% D)
            
            betahat_est <- split(as.numeric(betahat), rep(1:length(beta.size), beta.size))
            names(betahat_est) <- paste("response", 1:p, sep = "")
         
            res <- y - D %*% betahat # = (y-D %*% betahat)
          } else {res <- y}
          
          # solve((t(D)%*%solve(Sigma2$varcov)%*%D), (t(D)%*%solve(Sigma2$varcov)%*%data))
          ## Log sum of cholesky diagonal elements (soma a diagonal de uma matriz e multiplica por p)
          
          Sigma_log_sum_chol_diag <- sum(sapply(1:p, function(x){results[[x]]$log.det.to.half}))
          ssres <- drop(crossprod(res, invSigma%*%res)) # = t(y-D %*% betahat)%*% solve(varcov) %*% (y-D %*% betahat)
          
          negloglik <- as.numeric(drop((n*p/2)*log(2*pi) + Sigma_log_sum_chol_diag + 
                                         (n/2)*log(det(SigmaB)) + 0.5*ssres))
          
        } else
          { # (ncov.model > 1) e (ncov.model != p)
            stop("error: wrong model specification. The size of cov.model must be equal to p")
          }
      }
  
  if(est_mean)
  {
    return(list(betahat = betahat_est, loglik = -negloglik))
  } else
  {
    return(list(loglik = -negloglik))
  }
}

