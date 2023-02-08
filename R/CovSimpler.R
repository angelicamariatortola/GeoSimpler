#' @title Calculates Value of the Simpler Covariance Function
#' @name CovSimpler
#' 
#' @description Computes the Simpler covariance function associated to a p-dimensional
#' multivariate Gaussian random field.
#' Different correlation functions are available to model marginal-covariance behaviors.
#' When \eqn{p = 1} the function reduces to the marginal-covariance function (univariate case).
#' 
#' @param varcov 	the covariance matrix.
#' 
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix crossprod tcrossprod Diagonal bdiag t
#' 
#' @export
#' 
CovSimpler <- function(coords = NULL, dists.lowertri = NULL,
                       cov.model = "exp", cov.pars, p = 1, 
                       nugget = rep(0, length(cov.model)), 
                       SigmaB = NULL) 
{
  #### Função para uso externo
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## Calcula a matriz de covariância uni e multivariada.
  
  ## coords ou dists.lowertri devem ser dados --> igual em varcov.spatial
  ## cov.model --> vetor de strings contendo os modelos para cada variável.
  ## modelos de correlação da função geoR::cov.spatial
  ## cov.pars --> lista de vetores contendo os parâmetros var, phi e kappa (se houver), 
  ## onde cada vetor é relacionado aos parametros de cada modelo 
  ## Exemplo:
  # cov.model = c("matern","exp")
  # cov.pars = list(c(sig2, phi, kappa),   -> parametros referentes a "matern"
  #                 c(sig2, phi))          -> parametros referentes a "exp"
  ## nugget = tau² --> deve ser um vetor de dimensão p, com p = número de respostas
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  ## ele precisa entrar com o p, pois não entra com os dados???
  
  if (is.null(coords) && is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) && !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if(is.null(cov.pars))
    stop("cov.pars must be provided")
  
  if ((p>1) && is.null(SigmaB)) 
    stop("for p>1, SigmaB must be a positive definite matrix of dimension pxp")
  if(!is.null(SigmaB) && !matrixcalc::is.positive.definite(SigmaB))
    stop("SigmaB must be a positive definite matrix")
  if (!is.null(SigmaB) && (nrow(SigmaB)!= ncol(SigmaB))) 
    stop("SigmaB must be a SQUARE positive definite matrix of dimension pxp")
  if ((p>1) && !is.null(SigmaB) && (nrow(SigmaB)!=p)) 
    stop("SigmaB must be a positive definite matrix of dimension pxp")
  
  
  if (!is.null(coords))
    n <- nrow(coords)
  if (!is.null(dists.lowertri))
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  
  ncov.model <- length(cov.model)
  Sigma_cov <- matrix(0, n*p, n*p)
  params <- param_check_break(cov.pars = cov.pars, cov.model = cov.model, p = p)
  
  if (length(nugget) != ncov.model) 
    stop("the length of the nugget vector must be equal to the length of the cov.model vector")
  
  if(ncov.model == 1)
  {     ## Cenario univariado 
    if(p == 1)
    {
      results <- varcov.spatial2(coords = coords, dists.lowertri = dists.lowertri,
                        kappa = unlist(params$kappap), nugget = nugget,
                        cov.model = cov.model, 
                        cov.pars = c(params$sigmap, params$phip))
        # usei varcov.spatial2, para retornar além da varcov a cholesky da varcov.
        # Opção que no varcov.spatial do geoR não consegui obter
      Sigma_cov <- results$varcov
    } else
    {# else (p > 1) --> multivariado: mesmo modelo de covariância para todas as respostas
    # calcula a cholesky apenas uma vez
    # entra com uma lista de 1 vetor de parâmetros, que é o mesmo para todas as variáveis.
  
      results <- varcov.spatial2(coords = coords, dists.lowertri = dists.lowertri,
                                 kappa = unlist(params$kappap), nugget = nugget,
                                 cov.model = cov.model, det = T, # det = T para calcular a cholesky
                                 cov.pars = c(params$sigmap, params$phip))
      ## results$sqrt.varcov --> decomposição de Cholesky de Sigma_cov
      ct <- lapply(1:p, function(x) {results$sqrt.varcov}) ## Atribuindo a mesma ct_uni para todas as respostas
      cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
      Sigma_cov <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
      Sigma_cov <- as.matrix(Sigma_cov)
    }
  } else
  {# else (ncov.model > 1) --> multivariado --> ncov.model == p
    if(ncov.model == p)
    {# Cenario multivariado (modelos distintos de covariância para as respostas)
     # calcula a cholesky das covariancias marginais para as p variáveis

      results <- lapply(1:ncov.model, function(x) 
      {
        varcov.spatial2(coords = coords, dists.lowertri = dists.lowertri, det = T,
                        kappa = params$kappap[[x]], nugget = nugget[x],
                        cov.model = cov.model[x], 
                        cov.pars = c(params$sigmap[x], params$phip[x]))
        # usei varcov.spatial2, para retornar além da varcov a cholesky da varcov.
        # Opção que no varcov.spatial do geoR não consegui obter
      })
      
      ct <- lapply(1:p, function(x) {results[[x]]$sqrt.varcov})
      cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
      Sigma_cov <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
      Sigma_cov <- as.matrix(Sigma_cov)
    }else
     { # (ncov.model > 1) e (ncov.model != p)
       stop("error: wrong model specification. The size of cov.model must be equal to p")
     }
  }
 
  return(list(varcov = Sigma_cov))
}

