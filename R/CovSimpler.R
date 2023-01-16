#' @title Calculates Value of the Simpler Covariance Function
#' @name CovSimpler
#' 
#' @description Computes the Simpler covariance function associated to a p-dimensional
#' multivariate Gaussian random field.
#' Different correlation functions are available to model marginal-covariance behaviors.
#' When \eqn{p = 1} the function reduces to the marginal-covariance function (univariate case).
#' 
#' @param varcov 	the covariance matrix.
#' @param log.det.chol.varcov the logarithmic determinant of the Cholesky decomposition of the covariance matrix.
#' @param func.decomp algorithm used for the decomposition of the covariance matrix in determinant calculation.
#' 
#' @importFrom matrixcalc is.positive.definite
#' @importFrom geoR matern
#' @importFrom Matrix chol crossprod tcrossprod Diagonal bdiag t
#' 
#' @export
#' 
CovSimpler <- function(coords = NULL, dists.lowertri = NULL,
                       cov.model, cov.pars, p, nugget = 0, 
                       chol.varcov = FALSE, SigmaB = NULL)
{
  #### Função para uso externo
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## dist.matrix --> matrix de distâncias (dimensão nxn)
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros 
  #  de cada modelo -->
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## nugget = tau²
  ## ele precisa entrar com o p, pois não entra com os dados???
  
  ## Exemplo:
  # cov.model = c("matern","exp")
  # cov.pars = list(c(sig2, phi, kappa),   -> parametros referentes a "matern"
  #                 c(sig2, phi))          -> parametros referentes a "exp"
  
  n <- nrow(dist.matrix) ## número de locações amostrais
  ncov.model <- length(cov.model)
  Sigma_cov <- matrix(0, n*p, n*p)
  func.det <- match.arg(func.det)
  
  if(ncov.model == 1)
  {
    pars_uni <- unlist(cov.pars) # apenas var, phi, kappa, ...
    sigmap <- pars_uni[1]
    
    ## Cenario univariado  
    if(p == 1)
    {
      Sigma_cov <- varcov.spatial(cov.model = cov.model[1], 
                                  dists.lowertri = dist(s100$coords[1:5,]),
                                  cov.pars = c(1,2), 
                                  func.inv = "cholesky")
      
      cov_marg(dist.matrix = dist.matrix, cov.model = cov.model[1],
                            cov.pars.uni = pars_uni)
      diag(Sigma_cov) <- nugget + sigmap
    
     # return(Sigma_cov)
    } else
    {# else (p > 1) --> multivariado: mesmo modelo de covariância para todas as respostas
    # calcula a cholesky apenas uma vez
    # entra com uma lista de 1 vetor de parâmetros, que é o mesmo para todas as variáveis.
  
      ct_uni <- chol(cov_marg(dist.matrix = dist.matrix, 
                              cov.model = cov.model[1], 
                              cov.pars.uni = pars_uni))
      
      ## Atribuindo a mesma ct_uni para todas as respostas
      ct <- lapply(1:p, function(x) {ct_uni}) 
      cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
      Sigma_cov <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
      Sigma_cov <- as.matrix(Sigma_cov)
      
      diag(Sigma_cov) <- nugget + sigmap
  
      # return(Sigma_cov)
    }
  } else
  {# else (ncov.model > 1) --> multivariado --> ncov.model == p
    if(ncov.model == p)
    {# Cenario multivariado (modelos distintos de covariância para as respostas)
     # calcula a cholesky das covariancias marginais para as p variáveis

      ct <- lapply(1:p, function(x) {Matrix::chol(cov_marg(dist.matrix = dist.matrix, 
                                                           cov.model = cov.model[x], 
                                                           cov.pars.uni = cov.pars[[x]]))})

      cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
      Sigma_cov <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
      Sigma_cov <- as.matrix(Sigma_cov)
      
      sigmap <- sapply(1:p, function(x){ cov.pars[[x]][1]})
      Matrix::diag(Sigma_cov) <- rep(sigmap + nugget, each = n)
      # neste caso, nugget é um vetor de dimensão p
      
      # return(Sigma_cov)
    }else
     { # (ncov.model > 1) e (ncov.model != p)
       stop("error: wrong model specification. The size of cov.model must be equal to p")
     }
  }
  
  if(chol.varcov)
  { # chol.varcov = T
    varcov.sqrt <- try(chol(Sigma_cov), silent = TRUE) # usando a decomposição de cholesky
   
    if (inherits(varcov.sqrt, "try-error")) 
    { # se der erro no calculo da cholesky, para.
      print(varcov.sqrt[1])
      stop()
    } else
    {  # se NÃO der erro no calculo da cholesky
      cov.logdeth <- sum(log(diag(varcov.sqrt))) # = sum(log(diag(R))) = 0.5*log|Sigma_cov| --> usa na verossimilhança
    }
    result <- list(varcov = Sigma_cov, 
                   chol.varcov = varcov.sqrt,
                   log.det.chol.varcov = cov.logdeth)
  } else
  { # chol.varcov = F
    result <- list(varcov = Sigma_cov)
  }
  
  return(result)
}

