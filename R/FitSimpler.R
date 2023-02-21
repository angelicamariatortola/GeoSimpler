#' @title Calculates the log-likelihood value for the Simpler Model
#' @name FitSimpler
#' 
#' @description Computes the value of the log-likelihood for the
#' Simpler Covariance Model associated to a p-dimensional
#' multivariate Gaussian random field.
#' 
#' @importFrom Matrix chol solve diag crossprod tcrossprod bdiag
FitSimpler <- function(ini.cov.pars, data, coords, cov.model, trend = rep("cte", length(cov.model)),
                       fix.nugget = FALSE, nugget = rep(0, length(cov.model)), fix.kappa = TRUE,  
                       est_mean = F, SigmaB = NULL, limits = pars.limits(), print.pars = FALSE, ...)
  #, mean = 0, method, hessian, logpars)
{
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## ini.cov.pars --> lista de vetores iniciais, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  ## data --> dados (vetor, matriz ou data.frame). Se vetor p = 1, se matriz ou data.frame, p = ncol(data)
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## trend deve ser uma lista com as especificações para a média de cada resposta. 
  ## o default é "cte" para todas as respostas.
  ## est_mean --> Se est_mean = T, o betahat será estimado pelos dados. 
  ## Se est_mean = F, o betahat é igual a zero e não é estimado.
  
  
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp

  ldots <- list()
  temp.list <- list()
  temp.list$print.pars <- print.pars
  
  if(is.vector(data) || (ncol(data) == 1)){p <- 1}
  if(is.matrix(data) || is.data.frame(data)){p <- ncol(data)}
  ## p deve ser igual a ncol(as.matrix(data))
  
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
  
  if (fix.nugget) 
    fixed.pars$nugget <- nugget
  if (fix.kappa) 
    fixed.pars$kappa <- unlist(params$kappap)
  
  temp.list$cov.model <- cov.model
  method.lik <- "ML"
  temp.list$method.lik <- method.lik
  
  coords <- as.matrix(coords)
  data <- as.vector(data)
  n <- nrow(coords)
  
  
  # if(ncov.model == 1)
  # {     ## Cenario univariado 
  #   if(p == 1)
  #   {

      tausq <- nugget
      ini <- params$phip
      
      lower.optim <- c(limits$phi["lower"])
      upper.optim <- c(limits$phi["upper"])
      
      fixed.values <- list()
      kappa <- unlist(params$kappap)
      
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
      } else {
        # estima o kappa
        ini <- c(ini, kappa)
        lower.optim <- c(lower.optim, limits$kappa["lower"])
        upper.optim <- c(upper.optim, limits$kappa["upper"])
      }
      
      if (fix.nugget && nugget > 0) {
        # não estima o nugget, e estima o sigma
        ini <- c(ini, params$sigmap)
        lower.optim <- c(lower.optim, limits$sigmasq["lower"])
        upper.optim <- c(upper.optim, limits$sigmasq["upper"])
      }
      
      names(ini) <- NULL
      
      # if (length(ini) == 1) # se houver só um parâmetro a ser estimado 
      #   justone <- TRUE else justone <- FALSE
      
      if(est_mean)
      {
        xmat <- unclass(trend.spatial(trend = trend[[1]], geodata = list(coords = coords, 
                                                                         data = data)))
        temp.list$xmat <- xmat
        # if (nrow(xmat) != n) 
        #   stop("trend matrix has dimension incompatible with the data")
        beta.size <- temp.list$beta.size <- dim(xmat)[2] 
        # número de colunas de xmat é o número de betas a serem estimados
      } else beta.size <- 0
      
      ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa)
      temp.list$coords <- coords
      temp.list$z <- as.vector(data)
      temp.list$n <- n
      temp.list$est_mean <- est_mean
      ## temp.list$loglik.cte --> termo constante da log-verossimilhança
      npars <- beta.size + 2 + sum(unlist(ip) == FALSE) # 2 de sigma e phi
      
      
      ## Usando o método "ML"
      if (ip$f.tausq && (nugget > 0)) # se não estima o nugget. Esse termo fica igual na loglik usual
        temp.list$loglik.cte <- (temp.list$n/2) * (-log(2 * pi))  else 
        # se tau é fixo e for zero ou se o tau não for fixo, ele usa sigma = 1 como valor inicial
        # na função FitSimpler_aux que será usada no optim
        temp.list$loglik.cte <- (temp.list$n/2) * (-log(2 * pi) + log(temp.list$n) - 1)
        # aqui ou o tau será estimado, ou ele é fixado em zero, neste caso algumas coisas mudam na loglik
        # não entendi de onde saiu o termo (n/2)*(log(temp.list$n) - 1). Aproximou a cholesky para a raiz quadrada da diagonal?
      
      likGRF.dists.vec <- lapply(split(as.data.frame(coords), 1), vecdist)
      
      range.dist <- range(likGRF.dists.vec)
      max.dist <- max(range.dist)
      min.dist <- min(range.dist)
      
      if (length(ini) == 1) {
        if (upper.optim == Inf) 
          upper.optim <- 50 * max.dist
        lik.minim <- do.call("optimize", c(list(FitSimpler_aux, 
                                                lower = lower.optim, upper = upper.optim, fp = fixed.values, 
                                                ip = ip, temp.list = temp.list), ldots))
        lik.minim <- list(par = lik.minim$minimum, value = lik.minim$objective, 
                          convergence = 0, message = "function optimize used")
      } else {
       
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
      
      
  #   }
  # }
  
      
      
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
    if (sigmasq > 1e-12) 
      tausq <- nugget/sigmasq
    check.sigmasq <- TRUE # o sigma foi estimado no optim
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
  
  if(est_mean)
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
      if(est_mean)
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

  if(est_mean)
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
  
  if(est_mean)
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
  
  # attr(lik.results, "geodata") <- name.geodata
  return(lik.results)
  
}
