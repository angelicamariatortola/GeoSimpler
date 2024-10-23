varcov.spatial2 <- function (coords = NULL, dists.lowertri = NULL, cov.model = "matern", 
                             nugget = 0, cov.pars = stop("no cov.pars argument"), 
                             det = FALSE, scaled = FALSE, sqrt.inv = FALSE, ...) 
{
  #### Função para uso interno
  ## calcula matrix de covariância univariada. Adaptação da varcov.spatial que 
  ## retorna além da varcov, a sua decomposição de Cholesky.
  # A saída dessa função já retorna:
  # varcov --> matrix de covariância univariada
  # sqrt.varcov --> decomposição de cholesky de varcov
  # log.det.to.half --> log do determinante da decomposição de cholesky de varcov log(det(chol(varcov)))
  # sqrt.inverse --> inversa da decomposição de cholesky de varcov
  
  # Argumentos:
  # cov.pars --> deve ser um vetor contendo c(sig2, phi, kappa (se houver))
  
  if (is.null(coords) & is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) & !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if (!is.null(coords)) 
    n <- nrow(coords)
  if (!is.null(dists.lowertri)) 
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  
  tausq <- nugget
  
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    if(length(cov.pars) > 2)
    {
      kappa <- cov.pars[3:length(cov.pars)]
    }
  } else {
    stop("cov.pars must be a vector")
  }
  
  if (!is.null(coords)) 
    dists.lowertri <- as.vector(dist(coords))
  if (round(1e+12 * min(dists.lowertri)) == 0) 
    warning("Two or more pairs of data at coincident (or very close) locations. \nThis may cause crashes in some matrices operations.\n")
  varcov <- matrix(0, n, n)
  
  if (scaled) {
    if (all(phi < 1e-12)) 
      varcov <- diag(x = (1 + (tausq/sum(sigmasq))), n)
    else {
      if (is.vector(cov.pars)) 
        cov.pars.sc <- c(1, phi)
      else cov.pars.sc <- cbind(1, phi)
      covvec <- geoR::cov.spatial(obj = dists.lowertri, cov.model = cov.model, 
                            kappa = kappa, cov.pars = cov.pars.sc)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      remove("covvec")
      if (sum(sigmasq) < 1e-16) 
        diag(varcov) <- 1
      else diag(varcov) <- 1 + (tausq/sum(sigmasq))
    }
  } else {
    if (all(sigmasq < 1e-10) | all(phi < 1e-10)) {
      varcov <- diag(x = (tausq + sum(sigmasq)), n)
    }
    else {
      covvec <- geoR::cov.spatial(obj = dists.lowertri, cov.model = cov.model, 
                            kappa = kappa, cov.pars = c(sigmasq, phi))
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      #remove("covvec")
      diag(varcov) <- tausq + sum(sigmasq)
    }
  }
  
  if (det | sqrt.inv) 
  {
    varcov.sqrt <- try(chol(varcov), silent = TRUE)
    if (inherits(varcov.sqrt, "try-error")) 
    {
      print(varcov.sqrt[1])
      stop()
    } else 
    {
      if (det) 
        cov.logdeth <- sum(log(diag(varcov.sqrt)))
      if (sqrt.inv) 
        inverse.sqrt <- solve(varcov.sqrt)
    }
  }
  
  if (det | sqrt.inv) 
  {
    result <- list(varcov = varcov, sqrt.varcov = varcov.sqrt) 
    
    if(det)
    {
      result$log.det.to.half <- cov.logdeth
    }
    if (sqrt.inv) 
    {
      result$sqrt.inverse <- inverse.sqrt
    }
    
  } else
  {
    result <- list(varcov = varcov)
  }
  
  result$crash.parms <- NULL
  return(result)
}