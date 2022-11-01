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
CovSimpler <- function(dist.matrix, cov.model, cov.pars, SigmaB = NULL)
{
  #### Função para uso externo
  ### Considerando que todas as variáveis foram observadas nas mesmas locações amostrais
  ## dist.matrix --> matrix de distâncias (dimensão nxn)
  ## cov.model --> vetor de strings contendo os modelos para cada variável
  ## cov.pars --> lista de vetores, onde cada vetor é relacionado aos parametros 
  #  de cada modelo
  ## SigmaB --> (necessária se p>1) matriz de correlação de dimensão pxp
  
  ## Exemplo:
  # cov.model = c("matern","exp")
  # cov.pars = list(c(var, phi, kappa),   -> parametros referentes a "matern"
  #                 c(var, phi))   -> parametros referentes a "exp"
  
  ##-- Time
  # time_start <- Sys.time()
  
  ##-- Tests
  # Verificando se:
  ## todos os argumentos foram dados
  # if(missing(dist.matrix)) stop("You must provide the distance matrix", call. = F)
  # if(missing(cov.model)) stop("You must provide the vector cov.model", call. = F)
  # if(missing(cov.pars)) stop("You must provide the list cov.pars", call. = F)
  # 
  # ## todos os argumentos são do tipo correto
  # if(!is.matrix(dist.matrix))
  # {stop("dist.matrix must be a matrix", call. = F)} 
  # if(!is.vector(cov.model))
  # {stop("cov.model must be a vector", call. = F)} 
  # if(!is.list(cov.pars))
  # {stop("cov.pars must be a list", call. = F)} 
  # if(all(!is.null(SigmaB)) && !is.matrix(SigmaB))
  # {stop("SigmaB must be a correlation matrix", call. = F)} 
  # 
  # ## todos os argumentos assumem os valores corretos
  # if(any(!cov.model %in% c("exp", "matern", "gaussian", "cauchy"))) 
  # {stop("cov.model must be a vector whose elements must be between the options: 
  #        'exp', 'matern', 'gaussian', 'cauchy'", call. = F)}
  # if(any(!unlist(cov.pars)>0)) 
  # {stop("You must provide a list with positive elements for cov.pars", call. = F)}
  # if(all(!is.null(SigmaB)) && !matrixcalc::is.positive.definite(SigmaB)) 
  #   stop("You must provide a positive definite correlation matrix for SigmaB", call. = F)
  # 
  # ## todos os argumentos estão com as dimensões corretas
  # if((length(cov.model) > 1) && all(is.null(SigmaB)))
  # {stop("SigmaB is missing", call. = F)}
  # if(length(cov.pars) !=  length(cov.model))
  # {stop("the size of the cov.pars list must be equal 
  #       to the size of the cov.model vector", call. = F)} 
  # if(all(!is.null(SigmaB)) && (dim(SigmaB)[1] != length(cov.model)))
  # {stop("the number of rows and columns of SigmaB must be equal 
  #       to the length of cov.model", call. = F)} 
  
  
  n <- nrow(dist.matrix) ## número de locações amostrais
  p <- length(cov.model) ## número de variáveis
  
  
  ## todos os modelos tem o numero correto de parâmetros
  # lapply(1:p, function(x) 
  # {
  #   if((cov.model[x] == "exp") && (length(cov.pars[[x]]) != 2)){
  #     stop("the number of parameters for the exponential model 
  #         must be equal to 2: variance and range, respectively", call. = F)} 
  #   if((cov.model[x] == "gaussian") && (length(cov.pars[[x]]) != 2)){
  #     stop("the number of parameters for the gaussian model 
  #         must be equal to 2: variance and range, respectively", call. = F)}
  #   if((cov.model[x] == "matern") && (length(cov.pars[[x]]) != 3)){
  #     stop("the number of parameters for the matern model 
  #         must be equal to 3: variance, range and smoothness, respectively", call. = F)}
  #   if((cov.model[x] == "cauchy") && (length(cov.pars[[x]]) != 3)){
  #     stop("the number of parameters for the cauchy model 
  #         must be equal to 3: variance, range and smoothness, respectively", call. = F)}
  # })
  
  
  ## Marginal-covariance (univariate)
  if(ncov.model==1 && (p == 1))
  {
    pars_uni <- unlist(cov.pars)
    Sigma <- cov_marg(dist.matrix = dist.matrix, cov.model = cov.model[1],
                      cov.pars.uni = pars_uni)
    return(Sigma)
  }
  
  ## Cross-covariance (multivariate --> p>1) 
  ## Modelos marginais de correlação podem ser iguais ou não (a variancia pode ser diferente)
  
  ## calcula a cholesky das matrizes de covariancia marginais para as p variáveis
  ct <- lapply(1:p, function(x) {Matrix::chol(cov_marg(dist.matrix = dist.matrix, 
               cov.model = cov.model[x], 
               cov.pars.uni = cov.pars[[x]]))})
  
  cros1 <- Matrix::crossprod(Matrix::bdiag(ct), kronecker(SigmaB, Matrix::Diagonal(n)))
  Sigma <- Matrix::tcrossprod(cros1, Matrix::t(Matrix::bdiag(ct)))
  Sigma <- as.matrix(Sigma)
  
  return(Sigma)
}
