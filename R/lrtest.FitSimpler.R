#' @title Realiza o teste de razão de verossimilhanças para modelos encaixados
#' @name lrtest.FitSimpler
#' 
#' @description função que calcula....
#' 
lrtest.FitSimpler <- function(fullmodel, reducedmodel)
{
  ## fullmodel --> lista ou vetor (do modelo completo) com o valor da log verossimilhança 
  # e o número de parâmetros do modelo
  ## reducedmodel --> lista ou vetor (do modelo reduzido) com o valor da log verossimilhança 
  # e o número de parâmetros do modelo
  
  reducedmodel_vec <- as.numeric(reducedmodel)
  fullmodel_vec <- as.numeric(fullmodel)
  
  method <- "Ratio likelihood test"
  
  x1 <- -2*(reducedmodel_vec[1] - fullmodel_vec[1])
  df <- fullmodel_vec[2]-reducedmodel_vec[2]
  p_value <- pchisq(x1, df, lower.tail = F)

  results <- list(x1 = x1, df = df, p_value = p_value)

  cat('\n') 
  cat('---------------------------------------------------\n')
  cat('  Ratio likelihood test \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("statistic:", results$x1, "\n")
  cat("df:", results$df, "\n")
  cat("p-value:", results$p_value, "\n")
  
  # cat(paste("statistic:", results$x1), paste("df:", results$df), 
  #     paste("p-value:", results$p_value), sep = ", ")
  cat('\n')
  
}