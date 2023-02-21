#' @title Simulation of Gaussian random fields with a Simpler Covariance Function
#' @name grfSimpler
#' 
#' @description Simulation of Gaussian random fields with a Simpler Covariance Function
#' 
#' @importFrom matrixcalc is.positive.definite
#' 
grfSimpler <- function(n, coords = NULL, nx, ny,
                       xlims = c(0, 1), ylims = c(0, 1), nsim = 1, 
                       cov.model = "exp", cov.pars, nugget = rep(0, length(cov.model)),
                       p = 1,  mean = rep(0, p), SigmaB = NULL)
{
  #### Função para uso externo
  
  # nx optional. Number of points in the X direction.
  # ny optional. Number of points in the Y direction.
  ## nsim --> número de simulações
  
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  ## SigmaB --> (opcional) matriz de correlação de dimensão pxp, 
  #  necessária se p=length(cov.model) > 1
  ## considerando grid regular. O usuário entra com as coordenadas em forma de matriz ou não
  ## method = "cholesky"
  # mean --> vetor de médias, uma para cada variável
  
  rseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  ncov.model <- length(cov.model)
  params <- param_check_break(cov.pars = cov.pars, cov.model = cov.model, p = p)
  
  
  # if (cov.model == "stable") 
  #   cov.model <- "powered.exponential"
  # if (cov.model == "matern" && kappa == 0.5) 
  #   cov.model <- "exponential"
  
  results <- list(coords = NULL, data = NULL)

  ######## Coordenadas #############
  if ((!missing(nx) && nx == 1) || (!missing(ny) && ny == 1) || 
      diff(xlims) == 0 || diff(ylims) == 0) 
    {
      sim1d <- TRUE
      #if (messages.screen) 
      cat("simulations in 1D\n")
    } else sim1d <- FALSE
  
  if(!is.null(coords))
  { # se o usuário fornece as coordenadas
    results$coords <- as.matrix(coords)
    x1vals <- sort(unique(round(results$coords[, 1], digits = 12)))
    x2vals <- sort(unique(round(results$coords[, 2], digits = 12)))
    
    cat("grf: simulation on a set of locations provided by the user\n")
  } else
  { # se o usuário não fornece as coordenadas. Precisamos criar a sequencia de x e y
    if (missing(nx)) 
      {
        if (sim1d) {nx <- ifelse(diff(xlims) == 0, 1, n)}
        else nx <- n
      }
    if (missing(ny))
      {
        if (sim1d) {ny <- ifelse(diff(ylims) == 0, 1, n)}
        else ny <- n
      }
    
    # Gerando a sequencia de pontos x e y
    xpts <- seq(xlims[1], xlims[2], length = nx)
    ypts <- seq(ylims[1], ylims[2], length = ny)
    
    # # Amplitude entre as coordenadas x e entre as coordenadas y
    # xspacing <- ifelse(length(xpts) == 1, 0, diff(xpts[1:2]))
    # yspacing <- ifelse(length(ypts) == 1, 0, diff(ypts[1:2]))
    
    # Criando o grid de locações 
    results$coords <- as.matrix(expand.grid(x = xpts, 
                                            y = ypts))
    
    # # Verifica se o espaçamento entre os pontos x é o mesmo que entre os pontos y
    # equal.spacing <- ifelse(abs(xspacing - yspacing) < 1e-12, TRUE, FALSE)
    
    cat(paste("grf: generating grid ", nx, " * ", ny, " with ", (nx * ny), " points\n"))
  }
  
  ################ Simulação #############
  
  n <- nrow(results$coords)
  
  # # Número de pontos distintos nas coordenadas
  # if (length(unique(round(results$coords[, 1], digits = 12))) ==  1 ||
  #     length(unique(round(results$coords[, 2], digits = 12))) ==  1) 
  #   {
  #     sim1d <- TRUE
  #   } else {sim1d <- FALSE}
  
  Sigma_cov <- CovSimpler(coords = results$coords, nugget = nugget,
                    cov.model = cov.model, cov.pars = cov.pars, p = p, 
                    SigmaB = SigmaB) 
      
  z <- matrix(rnorm(n*p * nsim), nrow = n*p, ncol = nsim)
  
  cholSz <- crossprod(chol(Sigma_cov$varcov), z) ## matriz de dimensão np x nsim

  sim <- lapply(1:nsim, function(x){matrix(cholSz[,x], nc = p)}) # quebra cada simulação em p colunas (uma resposta em cada coluna)

  results$data <- lapply(1:nsim, function(y){sapply(1:p, function(x){sim[[y]][,x] + mean[x]})})
  names(results$data) <- paste("sim", 1:nsim, sep = "")
  
  results[c("cov.model", "nugget", "cov.pars", ".Random.seed")] <- 
    list(cov.model, nugget, cov.pars, rseed)
  
  return(results)
}
