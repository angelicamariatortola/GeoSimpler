#' @title Calculates Value of the Simpler Covariance Function
#' @name CovSimpler
#' 
#' @description Computes the Simpler covariance function associated to a p-dimensional
#' multivariate Gaussian random field.
#' Different correlation functions are available to model marginal-covariance behaviors.
#' When \eqn{p = 1} the function reduces to the marginal-covariance function (univariate case).
#' 
#' @importFrom matrixcalc is.positive.definite
#' @importFrom geoR matern
#' @importFrom Matrix chol crossprod tcrossprod Diagonal bdiag t
#' 
#' @export
#' 
CovSimpler <- function(dist.matrix, cov.model, cov.pars, p, nugget = NULL, SigmaB = NULL)
{
  #### Função para uso externo
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## dist.matrix --> matrix de distâncias (dimensão nxn)
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## nugget = tau²
  
  ## Exemplo:
  # cov.model = c("matern","exp")
  # cov.pars = list(c(sig2, phi, kappa),   -> parametros referentes a "matern"
  #                 c(sig2, phi))          -> parametros referentes a "exp"
  
  n <- nrow(dist.matrix) ## número de locações amostrais
  ncov.model <- length(cov.model)
  
  if(ncov.model == 1)
  {
    ## Cenario univariado  
    if(p == 1)
    {
      pars_uni <- unlist(cov.pars)
     
      Sigma <- cov_marg(dist.matrix = dist.matrix, cov.model = cov.model[1],
                        cov.pars.uni = pars_uni)
      
      if(!is.null(nugget)) # nugget != NULL
      {# neste caso nugget deve ser um valor numérico
        sigma <- pars_uni[1]
        Matrix::diag(Sigma) <- sigma + nugget
      }

      return(Sigma)
    }else
    {# else (p > 1) --> multivariado: mesmo modelo de covariância para todas as respostas
     # calcula a cholesky apenas uma vez
      pars_uni <- unlist(cov.pars)
      ct_uni <- Matrix::chol(cov_marg(dist.matrix = dist.matrix, 
                                      cov.model = cov.model, 
                                      cov.pars.uni = pars_uni))
    
      ## Atribuindo a mesma ct_uni para todas as respostas
      ct <- lapply(1:p, function(x) {ct_uni}) 
      cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
      Sigma <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
      Sigma <- as.matrix(Sigma)
      
      if(!is.null(nugget)) # nugget != NULL
      {# neste caso nugget deve ser um valor numérico (size = 1)
        sigma <- pars_uni[1]
        Matrix::diag(Sigma) <- rep(sigma + nugget, each = n*p)
      }
      
      return(Sigma)
    }
  }else
  {# else (ncov.model > 1) --> multivariado --> ncov.model == p
    if(ncov.model == p)
    {# Cenario multivariado (modelos distintos de covariância para as respostas)
     # calcula a cholesky das covariancias marginais para as p variáveis

      ct <- lapply(1:p, function(x) {Matrix::chol(cov_marg(dist.matrix = dist.matrix, 
                                                           cov.model = cov.model[x], 
                                                           cov.pars.uni = cov.pars[[x]]))})

      cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
      Sigma <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
      Sigma <- as.matrix(Sigma)
      
      if(!is.null(nugget)) # nugget != NULL
      {# neste caso nugget deve ser um vetor (size = length(cov.model))
        sigma <- sapply(1:p, function(x){ cov.pars[[x]][1]})
        Matrix::diag(Sigma) <- rep(sigma + nugget, each = n)
      }
      
      return(Sigma)
    }else
     { # ncov.model != p
       return("error: wrong model specification")
     }
  }
}

