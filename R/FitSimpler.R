#' @title Calculates the log-likelihood value for the Simpler Model
#' @name FitSimpler
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom Matrix solve diag crossprod tcrossprod bdiag
#' 
#' @export
#' 
FitSimpler <- function(ini.cov.pars, data, coords, cov.model, trend,
                       fix.nugget = TRUE, fix.kappa = TRUE, fix.mean = TRUE, 
                       nugget = rep(0, length(cov.model)), SigmaB = NULL,
                       limits = pars.limits2(), print.pars = FALSE, ...)
  #, mean = 0, method, hessian, logpars)
{
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## Estima os parâmetros da função de verossimilhança Simpler
  
  ## ini.cov.pars --> lista de vetores iniciais, onde cada vetor é relacionado aos parametros 
  #  de cada modelo. 
  # Ex. cov.model("exp", "matern") --> ini.cov.pars = list(c(sig2, phi), c(sig2, phi, kappa))
  
  ## data --> dados (vetor, matriz ou data.frame). Se vetor p = 1, 
  ## se matriz ou data.frame, p = ncol(data)
  ## coords --> matriz de coordenadas --> igual em varcov.spatial
  
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## trend --> deve ser uma lista com as especificações para a média de cada resposta. 
  ## o default é "cte" para todas as respostas.
  # fix.nugget, fix.kappa --> TRUE ou FALSE. Se TRUE, não estima estes parâmetros (para todas as respostas)
  ## nugget --> vetor de valores iniciais para o parâmetro nugget. Deve ter dimensão p
  # fix.mean, se TRUE, não estima a média e fixa igual a zero.
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp

  ## limits --> valores dos limites dos parâmetros na otimização --> mesmo de likfit
  
  ldots <- list()
  temp.list <- list()
  temp.list$print.pars <- print.pars
    
  if(is.vector(data) || (ncol(data) == 1)){p <- 1}
  if(is.matrix(data) || is.data.frame(data)){p <- ncol(data)}
  ## p deve ser igual a ncol(as.matrix(data))
  
  if(missing(trend)){trend <- rep(list("cte"), p)}
  
  ncov.model <- length(cov.model) ## número de variáveis
  params <- param_check_break(cov.pars = cov.pars, cov.model = cov.model, p = p)
  
  # for (i in 1:ncov.model)
  # {
  #   if (any(cov.model[i] == c("power", "gneiting.matern", "gencauchy"))) 
  #   {
  #     stop(paste("parameter estimation for", cov.model[i], "is not yet implemented"))
  #   }
  #   
  # }


  ####################################
  
  fixed.pars <- list(cov.model = cov.model)
  sigmas <- unlist(params$sigmap) # vetor de sigmas
  phis <- unlist(params$phip) # vetor de phis

  if (fix.nugget)
    fixed.pars$nugget <- nugget
   # fixed.pars$nugget <- nugget[which(fix.nugget == T)]
  
  if (fix.kappa) 
    fixed.pars$kappa <- kappa
   # fixed.pars$kappa <- kappa[which(fix.kappa == T)]
  

  method.lik <- "ML"
  n <- nrow(coords)
  z <- as.vector(data)
  
  temp.list$n <- n
  temp.list$p <- p
  temp.list$z <- z # length = n*p
  temp.list$coords <- as.matrix(coords)
  temp.list$cov.model <- cov.model
  temp.list$ncov.model <- ncov.model
  temp.list$fix.mean <- fix.mean
  temp.list$method.lik <- method.lik
  temp.list$SigmaB <- SigmaB
  
if(ncov.model == 1)
{     ## Cenario univariado
  if(p == 1)
  {
    tausq <- nugget
    ini <- params$phip
    fixed.values <- list()

    lower.optim <- c(limits$phi["lower"])
    upper.optim <- c(limits$phi["upper"])

    if (fix.nugget) {
      # não estima o nugget
      fixed.values$tausq <- nugget
    } else {
      # estima o nugget
      ini <- c(ini, nugget/params$sigmap)
      lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
      upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
    }

    if (fix.kappa) {
      # não estima o kappa
      fixed.values$kappa <- kappa
    } else
      {
        # estima o kappa
        ini <- c(ini, kappa)
        lower.optim <- c(lower.optim, limits$kappa["lower"])
        upper.optim <- c(upper.optim, limits$kappa["upper"])
      }

    if (fix.nugget && nugget > 0)
      {
        # não estima o nugget, e estima o sigma
        ini <- c(ini, params$sigmap)
        lower.optim <- c(lower.optim, limits$sigmasq["lower"])
        upper.optim <- c(upper.optim, limits$sigmasq["upper"])
      }

    names(ini) <- NULL
    beta.size <- 0
    # print(ini)

    # if (length(ini) == 1) # se houver só um parâmetro a ser estimado
    #   justone <- TRUE else justone <- FALSE

    if(!fix.mean)
    {
      xmat <- unclass(trend.spatial(trend = trend[[1]], geodata = list(coords = coords,
                                                                       data = data)))
      temp.list$xmat <- xmat
      # if (nrow(xmat) != n)
      #   stop("trend matrix has dimension incompatible with the data")
      beta.size <- temp.list$beta.size <- dim(xmat)[2]
      # número de colunas de xmat é o número de betas a serem estimados
    }

    ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa)

    ## temp.list$loglik.cte --> termo constante da log-verossimilhança
    npars <- beta.size + 2*p + p*(sum(unlist(ip) == FALSE)) # p de sigma e phi para cada resposta


    ## Usando o método "ML"
    if (ip$f.tausq && all(nugget > 0)) # se não estima o nugget. Esse termo fica igual na loglik usual
      temp.list$loglik.cte <- ((temp.list$n)/2) * (-log(2 * pi))  else
      # se tau é fixo e for zero ou se o tau não for fixo, ele usa sigma = 1 como valor inicial
      # na função FitSimpler_aux que será usada no optim
      temp.list$loglik.cte <- ((temp.list$n)/2) * (-log(2 * pi) + log(temp.list$n) - 1)
      # aqui ou o tau será estimado, ou ele é fixado em zero, neste caso algumas coisas mudam na loglik
      # o termo (n/2)*(log(temp.list$n) - 1), vem do termo -(1/2)(n log sigmasq + n) (model based geo - p. 112 - verossimilhança concentrada),
      # pois sigmasq = ssres/np.
      # logo, desenvolvendo -(1/2)(np log sigmasq + np), obtemos:
      # -(1/2)(np log (ssres/n) + n) = -(1/2)(n log (ssres) - nlog(n) + n)
      # Cujo termo constante é: (1/2)(nlog(n) - n)

    likGRF.dists.vec <- lapply(split(as.data.frame(coords), 1), vecdist)
    range.dist <- range(likGRF.dists.vec)
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)

    if (length(ini) == 1)
      {
        if (upper.optim == Inf)
          upper.optim <- 50 * max.dist
        lik.minim <- do.call("optimize", c(list(FitSimpler_aux,
                                                lower = lower.optim, upper = upper.optim, fp = fixed.values,
                                                ip = ip, temp.list = temp.list), ldots))
        lik.minim <- list(par = lik.minim$minimum, value = lik.minim$objective,
                          convergence = 0, message = "function optimize used")
      } else
        {
          MET <- pmatch(names(ldots), names(formals(optim)))

          if(is.na(MET) || all(names(formals(optim))[MET] != "method"))
            ldots$method <- "L-BFGS-B"

          if(!is.null(names(ldots))){
            names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- "method"
          }

          if(!is.null(ldots$method) && ldots$method == "L-BFGS-B"){
            ldots$lower <- lower.optim
            ldots$upper <- upper.optim
          }

          lik.minim <- do.call("optim", c(list(par = ini, fn = FitSimpler_aux,
                                               fp = fixed.values, ip = ip, temp.list = temp.list), ldots))
      }


    ## Atribuindo os parâmetros estimados:

    par.est <- lik.minim$par

    if (any(par.est < 0))
      par.est <- round(par.est, digits = 12)

    phi <- par.est[1]

    if (is.R())  ## pq isso? ele pode internamente no optim usar outra linguagem para otimizar?
      loglik.max <- -lik.minim$value else loglik.max <- -lik.minim$objective

    if (ip$f.tausq & !ip$f.kappa) {
      kappa <- par.est[2]
    }

    if (!ip$f.tausq & ip$f.kappa) {
      tausq <- par.est[2]
    }

    if (!ip$f.tausq & !ip$f.kappa) {
      tausq <- par.est[2]
      kappa <- par.est[3]
    }

    if (fix.nugget && nugget > 0) {
      sigmasq <- par.est[length(par.est)]
      if (sigmasq > 1e-12) tausq <- nugget/sigmasq
      check.sigmasq <- TRUE # após estimar o sigma ele estima o tau fazendo nugget/sigma?
    } else check.sigmasq <- FALSE # o sigma não foi estimado no optim, será estimado pelos dados

    ## Estimando a matriz de covariância
    if ((phi < 1e-12))
    {
      V <- diag(x = (1 + tausq), n) # fixou sigmaq = 1
    } else {
      if (check.sigmasq) # sigma foi estimado no optim
      {
        if (sigmasq < 1e-12)
        {
          if (!fix.nugget) # estimou o nugget no optim
            V <- diag(x = (1 + tausq), n) else
              V <- diag(x = sqrt(tausq), n) # não entendi
        } else V <- varcov.spatial2(coords = coords,
                                    cov.model = cov.model, kappa = kappa, nugget = tausq,
                                    cov.pars = c(1, phi))$varcov
      } else V <- varcov.spatial2(coords = coords,
                                  cov.model = cov.model, kappa = kappa, nugget = tausq,
                                  cov.pars = c(1, phi))$varcov
    }


    ## Estimando os betas

    if(!fix.mean)
    {
      ivvx <- solve(V, xmat)
      xivx <- crossprod(ivvx, xmat) # t(D)%*%inv(Sig)%*%D
      xivy <- crossprod(ivvx, data) # t(D)%*%inv(Sig)%*%y

      betahat <- .solve.geoR(xivx, xivy)
      betahat.var <- .solve.geoR(xivx) # variância de beta
      res <- as.vector(temp.list$z - xmat %*% betahat)

    } else
    {
      res <- as.vector(temp.list$z)
    }

    ivvy <- solve(V, data)
    yivy <- crossprod(data, ivvy) # t(y) %*%inv(Sig)%*%y

    ## Estimando o sigmasq
    if (!fix.nugget | (nugget < 1e-12))
    { # nugget é estimado no optim (não fixo) ou nulo
      if(!fix.mean)
      {
        ssres <- as.vector(yivy - 2 * crossprod(betahat, xivy) +
                             crossprod(betahat, xivx) %*% betahat)
        # = crossprod((data - xmat%*%betahat), solve(V, (data - xmat%*%betahat)))
      } else {ssres <- as.vector(crossprod(data, ivvy))}

      if (method.lik == "ML")
        sigmasq <- ssres/n else sigmasq <- ssres/(n - beta.size)
    }

    if (fix.nugget)
    {
      if (nugget > 0)
        tausq <- nugget
    }  else tausq <- tausq * sigmasq

    if(!fix.mean)
      if (sigmasq > 1e-12)
        betahat.var <- sigmasq * betahat.var

    if (phi < 0.001 * min.dist)
    {
      tausq <- tausq + sigmasq
      sigmasq <- 0
    }

    if (sigmasq < 1e-12)
      phi <- 0


    n.model.pars <- beta.size + 4
    par.su <- data.frame(status = rep(-9, n.model.pars))
    ind.par.su <- c(rep(0, beta.size), ip$f.tausq, 0, 0, ip$f.kappa)
    par.su$status <- ifelse(ind.par.su, "fixed", "estimated")

    if(!fix.mean)
    {
      par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa), digits = 4)

      if (beta.size == 1)
        beta.name <- "beta" else beta.name <- paste("beta", 0:(beta.size - 1), sep = "")

        row.names(par.su) <- c(beta.name, "tausq", "sigmasq", "phi", "kappa")
        par.su <- par.su[c((1:(n.model.pars - 3)), n.model.pars -
                             1, n.model.pars - 2, n.model.pars), ]

        lik.results <- list(cov.model = cov.model, nugget = tausq,
                            cov.pars = c(sigmasq, phi), sigmasq = sigmasq, phi = phi,
                            kappa = kappa, beta = as.vector(betahat), beta.var = betahat.var,
                            tausq = tausq, practicalRange = practicalRange(cov.model = cov.model,
                                                                           phi = phi, kappa = kappa),
                            method.lik = method.lik,
                            trend = trend, loglik = loglik.max, npars = npars, AIC = -2 * (loglik.max - npars),
                            BIC = -2 * (loglik.max - 0.5 * log(n) * npars),
                            parameters.summary = par.su, info.minimisation.function = lik.minim,
                            max.dist = max.dist, trend.matrix = xmat)

        if (length(lik.results$beta.var) == 1)
          lik.results$beta.var <- as.vector(lik.results$beta.var)

        if (length(lik.results$beta) > 1)
        {
          if (inherits(trend, "formula") || (length(class(trend)) >
                                             0 && any(class(trend) == "trend.spatial")))
            beta.names <- c("intercept", paste("covar", 1:(ncol(xmat) - 1), sep = ""))
          else if (trend == "1st")
            beta.names <- c("intercept", "x", "y")
          else if (trend == "2nd")
            beta.names <- c("intercept", "x", "y", "x2", "xy",
                            "y2")
          names(lik.results$beta) <- beta.names
        }

    } else
    {
      par.su$values <- round(c(tausq, sigmasq, phi, kappa), digits = 4)
      row.names(par.su) <- c("tausq", "sigmasq", "phi", "kappa")
      par.su <- par.su[c((1:(n.model.pars - 3)), n.model.pars -
                           1, n.model.pars - 2, n.model.pars), ]

      lik.results <- list(cov.model = cov.model, nugget = tausq,
                          cov.pars = c(sigmasq, phi), sigmasq = sigmasq, phi = phi,
                          kappa = kappa, tausq = tausq,
                          practicalRange = practicalRange(cov.model = cov.model,
                                                          phi = phi, kappa = kappa),
                          method.lik = method.lik, loglik = loglik.max, npars = npars,
                          AIC = -2 * (loglik.max - npars),
                          BIC = -2 * (loglik.max - 0.5 * log(n) * npars),
                          parameters.summary = par.su, info.minimisation.function = lik.minim,
                          max.dist = max.dist)
    }

  } else
  {# else (p > 1) --> multivariado: mesmo modelo de covariância para todas as respostas
    # calcula a cholesky apenas uma vez
    # entra com uma lista de 1 vetor de parâmetros, que é o mesmo para todas as variáveis.
    # neste caso estima apenas valor para cada parâmetro da matriz de covariância, estima as correlações entre as variáveis
    # e os betas correspondentes para cada variável
    
    # valores iniciais
    tausq <- nugget
    ini <- params$phip
    ini_rho <- SigmaB[lower.tri(SigmaB)]
    nrhos <- length(ini_rho)
    temp.list$nrhos <- nrhos
    
    fixed.values <- list()
    lower.optim <- c(limits$phi["lower"])
    upper.optim <- c(limits$phi["upper"])
    
    if (fix.nugget) {
      # não estima o nugget
      fixed.values$tausq <- nugget
    } else {
      # estima o nugget
      ini <- c(ini, nugget/params$sigmap)
      lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
      upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
    }
    
    if (fix.kappa) {
      # não estima o kappa
      fixed.values$kappa <- kappa
    } else 
    {
      # estima o kappa
      ini <- c(ini, kappa)
      lower.optim <- c(lower.optim, limits$kappa["lower"])
      upper.optim <- c(upper.optim, limits$kappa["upper"])
    }
    
    if (fix.nugget && nugget > 0)
    {
      # não estima o nugget, e estima o sigma
      ini <- c(ini, params$sigmap)
      lower.optim <- c(lower.optim, limits$sigmasq["lower"])
      upper.optim <- c(upper.optim, limits$sigmasq["upper"])
    }
    
    # se tau é fixo e igual a zero, ele estima o sigma2 usando a expressão do estimador de máxima verossimilhança: ssres/N
    # se tau é fixo e igual a zero, ele estima o sigma2 usando a expressão do estimador de máxima verossimilhança: ssres/N
    
    
    # estima os rhos
    ini <- c(ini, ini_rho)
    lower.optim <- c(lower.optim, rep(limits$rho["lower"], nrhos))
    upper.optim <- c(upper.optim, rep(limits$rho["upper"], nrhos))
    
    names(ini) <- NULL
    beta.size <- 0
    # print(ini)
    
    # if (length(ini) == 1) # se houver só um parâmetro a ser estimado 
    #   justone <- TRUE else justone <- FALSE
    
    if(!fix.mean)
    {
      xmat <- lapply(1:p, function(x){unclass(trend.spatial(trend = trend[[x]],
                                                    geodata = list(coords = coords,
                                                                   data = data[,x])))})
      D <- bdiag(xmat) # agrupa as matrizes de covariáveis em uma matriz bloco diagonal (bonat (2020))
      temp.list$D <- D
      
      # if (nrow(xmat) != n) 
      #   stop("trend matrix has dimension incompatible with the data")
      beta.size <- sapply(1:p, function(x){ncol(xmat[[x]])})
    } 
    
    ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa)
    
    ## temp.list$loglik.cte --> termo constante da log-verossimilhança
    npars <- sum(beta.size) + 2 + (sum(unlist(ip) == FALSE)) + nrhos # 2 de sigma e phi 
    # nrhos retorna o número de rhos
    
    ## Usando o método "ML"
    if (ip$f.tausq && nugget > 0) # se não estima o nugget. Esse termo fica igual na loglik usual
      temp.list$loglik.cte <- (-(temp.list$n*p)/2) * (log(2 * pi))  else 
        # se tau é fixo e for zero ou se o tau não for fixo, ele usa sigma = 1 como valor inicial
        # na função FitSimpler_aux que será usada no optim
        temp.list$loglik.cte <- ((temp.list$n*p)/2) * (-log(2 * pi) + log(temp.list$n*p) - 1)
    # aqui ou o tau será estimado, ou ele é fixado em zero, neste caso algumas coisas mudam na loglik
    # o termo (np/2)*(log(temp.list$np) - 1), vem do termo -(1/2)(np log sigmasq + np) (model based geo - p. 112 - verossimilhança concentrada),
    # pois sigmasq = ssres/np.
    # logo, desenvolvendo -(1/2)(np log sigmasq + np), obtemos:
    # -(1/2)(np log (ssres/np) + np) = -(1/2)(np log (ssres) - nplog(np) + np)
    # Cujo termo constante é: (1/2)(nplog(np) - np)
    
    likGRF.dists.vec <- lapply(split(as.data.frame(coords), 1), vecdist)
    range.dist <- range(likGRF.dists.vec)
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
    
    if (length(ini) == 1)
    {
      if (upper.optim == Inf)
        upper.optim <- 50 * max.dist
      lik.minim <- do.call("optimize", c(list(FitSimpler_aux,
                                              lower = lower.optim, upper = upper.optim, fp = fixed.values,
                                              ip = ip, temp.list = temp.list), ldots))
      lik.minim <- list(par = lik.minim$minimum, value = lik.minim$objective,
                        convergence = 0, message = "function optimize used")
    } else
    {
      MET <- pmatch(names(ldots), names(formals(optim)))
      
      if(all(is.na(MET)) || all(names(formals(optim))[MET] != "method"))
        ldots$method <- "L-BFGS-B"
      
      if(!is.null(names(ldots))){
        names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- "method"
      }
      
      if(!is.null(ldots$method) && ldots$method == "L-BFGS-B"){
        ldots$lower <- lower.optim
        ldots$upper <- upper.optim
      }
      
      #  FitSimpler_aux(pars = ini, fp = fixed.values, ip = ip, temp.list = temp.list)
      # pars = ini
      # fp = fixed.values;ip = ip;temp.list = temp.list
      
      lik.minim <- do.call("optim", c(list(par = ini, fn = FitSimpler_aux, 
                                           fp = fixed.values, ip = ip, temp.list = temp.list), ldots))
    }
    
    
    ## Atribuindo os parâmetros estimados:
    
    par.est <- lik.minim$par
    
    if (any(par.est < 0)) 
      par.est <- round(par.est, digits = 12)
    
    phi <- par.est[1]
    
    if (is.R())  ## pq isso? ele pode internamente no optim usar outra linguagem para otimizar?
      loglik.max <- -lik.minim$value else loglik.max <- -lik.minim$objective
    
    if (ip$f.tausq & !ip$f.kappa) {
      kappa <- par.est[2]
    }
    
    if (!ip$f.tausq & ip$f.kappa) {
      tausq <- par.est[2]
    }
    
    if (!ip$f.tausq & !ip$f.kappa) {
      tausq <- par.est[2]
      kappa <- par.est[3]
    }
    
    if (fix.nugget && nugget > 0) {
      sigmasq <- par.est[length(par.est) - nrhos]
      if (sigmasq > 1e-12) tausq <- nugget/sigmasq
      check.sigmasq <- TRUE # o sigma foi estimado no optim
    } else check.sigmasq <- FALSE # o sigma não foi estimado no optim, será estimado pelos dados
    
    rhosq <- par.est[(length(par.est) - nrhos +1):length(par.est)]
    
    ## Estimando a matriz de covariância
    if (phi < 1e-12) 
    {
      V <- list(varcov = diag(x=(1+tausq), n),
                sqrt.varcov = sqrt(tausq+1)*diag(n),
                sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
                log.det.to.half = (n/2) * log(1+tausq)) # fixou sigmaq = 1
    } else {
      if (check.sigmasq) # sigma foi estimado no optim
      {
        if (sigmasq < 1e-12) 
        {
          if (!fix.nugget) # estima o v = tau²/sig² no optim
            V <- list(varcov = diag(x=(1+tausq), n),
                    sqrt.varcov = sqrt(tausq+1)*diag(n),
                    sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
                    log.det.to.half = (n/2) * log(1+tausq)) else
              # V <- diag(x = sqrt(tausq), n) # considerou sigma=0
                      V <- diag(x=sqrt(tausq), ni) # não entendi
              # V <- list(varcov = diag(x=(1+tausq), n),
              #           sqrt.varcov = sqrt(tausq+1)*diag(n),
              #           sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
              #           log.det.to.half = (n/2) * log(1+tausq))     
        } else V <- varcov.spatial2(coords = coords, 
                                    cov.model = cov.model, kappa = kappa, nugget = tausq, 
                                    cov.pars = c(1, phi), sqrt.inv = T)
      } else V <- varcov.spatial2(coords = coords, 
                                  cov.model = cov.model, kappa = kappa, nugget = tausq, 
                                  cov.pars = c(1, phi), sqrt.inv = T)
    }    
    
    ## Estimando a matriz SigmaB:
    
    m1 <- matrix(1, nc = p, nr = p)
    m1[lower.tri(m1)] <- rhosq
    SigmaBq <- t(m1)
    SigmaBq[lower.tri(SigmaBq)] <- rhosq
    
    ## Estimando os betas
    
    kroninv <- kronecker(solve(SigmaBq), crossprod(V$sqrt.inverse))
    
    if(!fix.mean)
    {
      xiv <- crossprod(D, kroninv) #t(D)%*% inv(Sig)
      xivx <- as.matrix(xiv %*% D) #t(D)%*% inv(Sig) %*% D
      xivy <- as.matrix(xiv %*% z) # t(D) %*%inv(Sig)%*%y
      betahat <- .solve.geoR(xivx, xivy)
      #betahat.var <- .solve.geoR(xivx) # variância de beta
      res <- as.vector(z - as.matrix(D %*% betahat))
    } else {res <- z}

    ## Estimando o sigmasq
    if (!fix.nugget || (nugget < 1e-12)) 
    { # nugget é estimado no optim (não fixo) ou nulo
      ssres <- drop(crossprod(res, kroninv%*%res)) # e^t %*% kronecker(inv(Sigmab), inv(V)) %*% e
      # = crossprod((data - xmat%*%betahat), solve(V, (data - xmat%*%betahat)))
      
      if (method.lik == "ML") 
        sigmasq <- ssres/(n*p) else sigmasq <- ssres/(n*p - sum(beta.size)) ## confirmar!!
    }
    
    if (fix.nugget) 
    {
      if (nugget > 0) 
        tausq <- nugget
    }  else tausq <- tausq * sigmasq
    
    # if(!fix.mean)
    #   if (sigmasq > 1e-12) 
    #     betahat.var <- sigmasq * betahat.var
    
    if (phi < 0.001 * min.dist) 
    {
      tausq <- tausq + sigmasq
      sigmasq <- 0
    }
    
    if (sigmasq < 1e-12) 
      phi <- 0
    
    ## Apresentação dos resultados
    n.model.pars <- sum(beta.size) + 4 + nrhos # betas + (sig² + phi + kappa + nugget) + rhos
    par.su <- data.frame(status = rep(-9, n.model.pars))
    ind.par.su <- c(rep(fix.mean, sum(beta.size)), ip$f.tausq, 0, 0, ip$f.kappa, rep(0, nrhos))
    par.su$status <- ifelse(ind.par.su, "fixed", "estimated")
    ncomb_rho <- paste(combn(p,2)[1,], combn(p,2)[2,], sep = "")
    rho.name <- paste("rho", ncomb_rho, sep = "")
    
    if(!fix.mean)
    {
      par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa, rhosq), digits = 4)
      betahat.est <- split(betahat, rep(1:length(beta.size), beta.size))
      beta.name <- lapply(1:p, function(x){paste("y", x, ".", "b", 0:(beta.size[x] - 1), sep = "")})
      #names(beta.name) <- paste("var", 1:p, sep = "")
      
      # beta.size = c(1,2,3)
      # betahatt = c(6,7,2,3,6,9)
      # bbn = split(betahatt, rep(1:length(beta.size), beta.size))
      # bb.name <- lapply(1:p, function(x){paste("b", 0:(beta.size[x] - 1), sep = "")})
      # names(bbn ) <- bb.name
      
      row.names(par.su) <- c(unlist(beta.name), "tausq", "sigmasq", "phi", "kappa", rho.name)
      # par.su <- par.su[c((1:(n.model.pars - 3)), n.model.pars - 
      #                      1, n.model.pars - 2, n.model.pars), ]
      
      lik.results <- list(cov.model = cov.model, nugget = tausq, 
                          cov.pars = c(sigmasq, phi), sigmasq = sigmasq, phi = phi, 
                          kappa = kappa, beta = betahat.est, rho = rhosq, #beta.var = betahat.var, 
                          tausq = tausq, practicalRange = practicalRange(cov.model = cov.model, 
                                                                         phi = phi, kappa = kappa), 
                          method.lik = method.lik, 
                          trend = trend, loglik = loglik.max, npars = npars, AIC = -2 * (loglik.max - npars), 
                          BIC = -2 * (loglik.max - 0.5 * log(n) * npars), 
                          parameters.summary = par.su, info.minimisation.function = lik.minim, 
                          max.dist = max.dist, trend.matrix = xmat)
        
       # if (length(lik.results$beta.var) == 1) 
       #   lik.results$beta.var <- as.vector(lik.results$beta.var)
      
       # Renomeando os betas
       # for (i in 1:length(trend)) 
       #  {
       #    if (length(lik.results$beta) > 1) 
       #    {
       #      if (inherits(trend[[i]], "formula") || (length(class(trend[[i]])) > 
       #                                         0 && any(class(trend[[i]]) == "trend.spatial"))) 
       #        beta.names[[i]] <- c("intercept", paste("covar", 1:(ncol(xmat) - 1), sep = ""))
       #      else if (trend[[i]] == "1st") 
       #        beta.names[[i]] <- c("intercept", "x", "y")
       #      else if (trend[[i]] == "2nd") 
       #        beta.names[[i]] <- c("intercept", "x", "y", "x2", "xy", 
       #                        "y2")
       #      names(lik.results$beta) <- beta.names[[i]]
       #    }
       #  }
        
    } else 
    {
      par.su$values <- round(c(tausq, sigmasq, phi, kappa, rhosq), digits = 4)
      row.names(par.su) <- c("tausq", "sigmasq", "phi", "kappa", rho.name)
      
      lik.results <- list(cov.model = cov.model, nugget = tausq, 
                          cov.pars = c(sigmasq, phi), sigmasq = sigmasq, phi = phi, 
                          kappa = kappa, tausq = tausq, rho = rhosq, 
                          practicalRange = practicalRange(cov.model = cov.model, 
                                                          phi = phi, kappa = kappa), 
                          method.lik = method.lik, loglik = loglik.max, npars = npars,
                          AIC = -2 * (loglik.max - npars), 
                          BIC = -2 * (loglik.max - 0.5 * log(n) * npars), 
                          parameters.summary = par.su, info.minimisation.function = lik.minim, 
                          max.dist = max.dist)
      }
      
    
  }
} else
{# else (ncov.model > 1) --> multivariado: diferentes modelos para as respostas
  # calcula a cholesky p vezes

  # valores iniciais
  kappas <- unlist(params$kappap) # vetor de kappas
  tausq <- nugget
  ini <- params$phip
  ini_rho <- SigmaB[lower.tri(SigmaB)]
  nrhos <- length(ini_rho)
  temp.list$nrhos <- nrhos
  
  fixed.values <- list()
  lower.optim <- rep(c(limits$phi["lower"]), p)
  upper.optim <- rep(c(limits$phi["upper"]), p)
  
  if (fix.nugget) {
    # não estima o nugget
    fixed.values$tausq <- nugget
  } else {
    # estima o nugget
    ini <- c(ini, nugget/sigmas)
    lower.optim <- c(lower.optim, rep(limits$tausq.rel["lower"], length = length(nugget)))
    upper.optim <- c(upper.optim, rep(limits$tausq.rel["upper"], length = length(nugget)))
  }
  
  if (fix.kappa) {
    # não estima o kappa
    fixed.values$kappa <- kappa
  } else 
  {
    # estima o kappa
    ini <- c(ini, kappas)
    lower.optim <- c(lower.optim, rep(limits$kappa["lower"], length = length(kappas)))
    upper.optim <- c(upper.optim, rep(limits$kappa["upper"], length = length(kappas)))
  }
  
  ini <- c(ini, sigmas)
  lower.optim <- c(lower.optim, rep(limits$sigmasq["lower"], length = length(sigmas)))
  upper.optim <- c(upper.optim, rep(limits$sigmasq["upper"], length = length(sigmas)))
  # se tau é fixo e igual a zero, ele estima o sigma2 usando a expressão do estimador de máxima verossimilhança: ssres/N

  
  # estima os rhos
  ini <- c(ini, ini_rho)
  lower.optim <- c(lower.optim, rep(limits$rho["lower"], nrhos))
  upper.optim <- c(upper.optim, rep(limits$rho["upper"], nrhos))
  
  names(ini) <- NULL
  beta.size <- 0
  # print(ini)
  
  # if (length(ini) == 1) # se houver só um parâmetro a ser estimado 
  #   justone <- TRUE else justone <- FALSE
  
  if(!fix.mean)
  {
    xmat <- lapply(1:p, function(x){unclass(trend.spatial(trend = trend[[x]],
                                                          geodata = list(coords = coords,
                                                                         data = data[,x])))})
    D <- bdiag(xmat) # agrupa as matrizes de covariáveis em uma matriz bloco diagonal (bonat (2020))
    temp.list$D <- D
    
    # if (nrow(xmat) != n) 
    #   stop("trend matrix has dimension incompatible with the data")
    beta.size <- sapply(1:p, function(x){ncol(xmat[[x]])})
  } 
  
  ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa)
  
  ## temp.list$loglik.cte --> termo constante da log-verossimilhança
  npars <- sum(beta.size) + 2*p + (sum(unlist(ip) == FALSE)) + nrhos # 2 de sigma e phi para cada variável p
  # npars --> número de parâmetros a serem estimados
  # nrhos --> retorna o número de rhos
  
  ## Usando o método "ML"
  temp.list$loglik.cte <- -((temp.list$n*p)/2)*log(2 * pi) 
 
  likGRF.dists.vec <- lapply(split(as.data.frame(coords), 1), vecdist)
  range.dist <- range(likGRF.dists.vec)
  max.dist <- max(range.dist)
  min.dist <- min(range.dist)
  
  MET <- pmatch(names(ldots), names(formals(optim)))
    
    if(all(is.na(MET)) || all(names(formals(optim))[MET] != "method"))
      ldots$method <- "L-BFGS-B"
    
    if(!is.null(names(ldots))){
      names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- "method"
    }
    
    if(!is.null(ldots$method) && ldots$method == "L-BFGS-B"){
      ldots$lower <- lower.optim
      ldots$upper <- upper.optim
    }
    
    # FitSimpler_aux(pars = ini, fp = fixed.values, ip = ip, temp.list = temp.list)
    pars = ini
    fp = fixed.values;ip = ip;temp.list = temp.list
    
    lik.minim <- do.call("optim", c(list(par = ini, fn = FitSimpler_aux, 
                                         fp = fixed.values, ip = ip, temp.list = temp.list), ldots))

  
  
  ## Atribuindo os parâmetros estimados:
  
  par.est <- lik.minim$par
  
  if (any(par.est < 0)) 
    par.est <- round(par.est, digits = 12)
  
  phi <- par.est[1]
  
  if (is.R())  ## pq isso? ele pode internamente no optim usar outra linguagem para otimizar?
    loglik.max <- -lik.minim$value else loglik.max <- -lik.minim$objective
  
  if (ip$f.tausq & !ip$f.kappa) {
    kappa <- par.est[2]
  }
  
  if (!ip$f.tausq & ip$f.kappa) {
    tausq <- par.est[2]
  }
  
  if (!ip$f.tausq & !ip$f.kappa) {
    tausq <- par.est[2]
    kappa <- par.est[3]
  }
  
  if (fix.nugget && nugget > 0) {
    sigmasq <- par.est[length(par.est) - nrhos]
    if (sigmasq > 1e-12) tausq <- nugget/sigmasq
    check.sigmasq <- TRUE # o sigma foi estimado no optim
  } else check.sigmasq <- FALSE # o sigma não foi estimado no optim, será estimado pelos dados
  
  rhosq <- par.est[(length(par.est) - nrhos +1):length(par.est)]
  
  ## Estimando a matriz de covariância
  if (phi < 1e-12) 
  {
    V <- list(varcov = diag(x=(1+tausq), n),
              sqrt.varcov = sqrt(tausq+1)*diag(n),
              sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
              log.det.to.half = (n/2) * log(1+tausq)) # fixou sigmaq = 1
  } else {
    if (check.sigmasq) # sigma foi estimado no optim
    {
      if (sigmasq < 1e-12) 
      {
        if (!fix.nugget) # estima o v = tau²/sig² no optim
          V <- list(varcov = diag(x=(1+tausq), n),
                    sqrt.varcov = sqrt(tausq+1)*diag(n),
                    sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
                    log.det.to.half = (n/2) * log(1+tausq)) else
                      # V <- diag(x = sqrt(tausq), n) # considerou sigma=0
                      V <- diag(x=sqrt(tausq), ni) # não entendi
                    # V <- list(varcov = diag(x=(1+tausq), n),
                    #           sqrt.varcov = sqrt(tausq+1)*diag(n),
                    #           sqrt.inverse = (1/sqrt(tausq+1))*diag(n),
                    #           log.det.to.half = (n/2) * log(1+tausq))     
      } else V <- varcov.spatial2(coords = coords, 
                                  cov.model = cov.model, kappa = kappa, nugget = tausq, 
                                  cov.pars = c(1, phi), sqrt.inv = T)
    } else V <- varcov.spatial2(coords = coords, 
                                cov.model = cov.model, kappa = kappa, nugget = tausq, 
                                cov.pars = c(1, phi), sqrt.inv = T)
  }    
  
  ## Estimando a matriz SigmaB:
  
  m1 <- matrix(1, nc = p, nr = p)
  m1[lower.tri(m1)] <- rhosq
  SigmaBq <- t(m1)
  SigmaBq[lower.tri(SigmaBq)] <- rhosq
  
  ## Estimando os betas
  
  kroninv <- kronecker(solve(SigmaBq), crossprod(V$sqrt.inverse))
  
  if(!fix.mean)
  {
    xiv <- crossprod(D, kroninv) #t(D)%*% inv(Sig)
    xivx <- as.matrix(xiv %*% D) #t(D)%*% inv(Sig) %*% D
    xivy <- as.matrix(xiv %*% z) # t(D) %*%inv(Sig)%*%y
    betahat <- .solve.geoR(xivx, xivy)
    #betahat.var <- .solve.geoR(xivx) # variância de beta
    res <- as.vector(z - as.matrix(D %*% betahat))
  } else {res <- z}
  
  ## Estimando o sigmasq
  if (!fix.nugget || (nugget < 1e-12)) 
  { # nugget é estimado no optim (não fixo) ou nulo
    ssres <- drop(crossprod(res, kroninv%*%res)) # e^t %*% kronecker(inv(Sigmab), inv(V)) %*% e
    # = crossprod((data - xmat%*%betahat), solve(V, (data - xmat%*%betahat)))
    
    if (method.lik == "ML") 
      sigmasq <- ssres/(n*p) else sigmasq <- ssres/(n*p - sum(beta.size)) ## confirmar!!
  }
  
  if (fix.nugget) 
  {
    if (nugget > 0) 
      tausq <- nugget
  }  else tausq <- tausq * sigmasq
  
  # if(!fix.mean)
  #   if (sigmasq > 1e-12) 
  #     betahat.var <- sigmasq * betahat.var
  
  if (phi < 0.001 * min.dist) 
  {
    tausq <- tausq + sigmasq
    sigmasq <- 0
  }
  
  if (sigmasq < 1e-12) 
    phi <- 0
  
  ## Apresentação dos resultados
  n.model.pars <- sum(beta.size) + 4 + nrhos
  par.su <- data.frame(status = rep(-9, n.model.pars))
  ind.par.su <- c(rep(fix.mean, sum(beta.size)), ip$f.tausq, 0, 0, ip$f.kappa, rep(0, nrhos))
  par.su$status <- ifelse(ind.par.su, "fixed", "estimated")
  ncomb_rho <- paste(combn(p,2)[1,], combn(p,2)[2,], sep = "")
  rho.name <- paste("rho", ncomb_rho, sep = "")
  
  if(!fix.mean)
  {
    par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa, rhosq), digits = 4)
    betahat.est <- split(betahat, rep(1:length(beta.size), beta.size))
    beta.name <- lapply(1:p, function(x){paste("y", x, ".", "b", 0:(beta.size[x] - 1), sep = "")})
    #names(beta.name) <- paste("var", 1:p, sep = "")
    
    # beta.size = c(1,2,3)
    # betahatt = c(6,7,2,3,6,9)
    # bbn = split(betahatt, rep(1:length(beta.size), beta.size))
    # bb.name <- lapply(1:p, function(x){paste("b", 0:(beta.size[x] - 1), sep = "")})
    # names(bbn ) <- bb.name
    
    row.names(par.su) <- c(unlist(beta.name), "tausq", "sigmasq", "phi", "kappa", rho.name)
    par.su <- par.su[c((1:(n.model.pars - 3)), n.model.pars - 
                         1, n.model.pars - 2, n.model.pars), ]
    
    lik.results <- list(cov.model = cov.model, nugget = tausq, 
                        cov.pars = c(sigmasq, phi), sigmasq = sigmasq, phi = phi, 
                        kappa = kappa, beta = betahat.est, rho = rhosq, #beta.var = betahat.var, 
                        tausq = tausq, practicalRange = practicalRange(cov.model = cov.model, 
                                                                       phi = phi, kappa = kappa), 
                        method.lik = method.lik, 
                        trend = trend, loglik = loglik.max, npars = npars, AIC = -2 * (loglik.max - npars), 
                        BIC = -2 * (loglik.max - 0.5 * log(n) * npars), 
                        parameters.summary = par.su, info.minimisation.function = lik.minim, 
                        max.dist = max.dist, trend.matrix = xmat)
    
    # if (length(lik.results$beta.var) == 1) 
    #   lik.results$beta.var <- as.vector(lik.results$beta.var)
    
    # Renomeando os betas
    # for (i in 1:length(trend)) 
    #  {
    #    if (length(lik.results$beta) > 1) 
    #    {
    #      if (inherits(trend[[i]], "formula") || (length(class(trend[[i]])) > 
    #                                         0 && any(class(trend[[i]]) == "trend.spatial"))) 
    #        beta.names[[i]] <- c("intercept", paste("covar", 1:(ncol(xmat) - 1), sep = ""))
    #      else if (trend[[i]] == "1st") 
    #        beta.names[[i]] <- c("intercept", "x", "y")
    #      else if (trend[[i]] == "2nd") 
    #        beta.names[[i]] <- c("intercept", "x", "y", "x2", "xy", 
    #                        "y2")
    #      names(lik.results$beta) <- beta.names[[i]]
    #    }
    #  }
    
  } else 
  {
    par.su$values <- round(c(tausq, sigmasq, phi, kappa, rhosq), digits = 4)
    row.names(par.su) <- c("tausq", "sigmasq", "phi", "kappa", rho.name)
    
    par.su <- par.su[c((1:(n.model.pars - 3)), n.model.pars - 
                         1, n.model.pars - 2, n.model.pars), ]
    
    lik.results <- list(cov.model = cov.model, nugget = tausq, 
                        cov.pars = c(sigmasq, phi), sigmasq = sigmasq, phi = phi, 
                        kappa = kappa, tausq = tausq, rho = rhosq, 
                        practicalRange = practicalRange(cov.model = cov.model, 
                                                        phi = phi, kappa = kappa), 
                        method.lik = method.lik, loglik = loglik.max, npars = npars,
                        AIC = -2 * (loglik.max - npars), 
                        BIC = -2 * (loglik.max - 0.5 * log(n) * npars), 
                        parameters.summary = par.su, info.minimisation.function = lik.minim, 
                        max.dist = max.dist)
  }
  
}
  
      
  # attr(lik.results, "geodata") <- name.geodata
  return(lik.results)
  
}
