
## Para carregar as funções do pacote:
setwd("/media/angelica/HD/Pesquisa/Projetos/Modelagem de dados espaciais multivariados/GeoSimpler package/GeoSimpler")

library(devtools)
document()
load_all()

## para selecionar vários itens que serão comitados para o git no RStudio
# clicar na área git
# Ctrl + A --> seleciona tudo
# enter --> marca todas as caixas

## 1 - Funcão do GeoSimpler  |  depende (função)        | depende (pacote)

## cov_marg                 | cov.spatial                   | geoR

## CovSimpler                | crossprod                | Matrix
##                           | tcrossprod               | Matrix
##                           | t                        | Matrix
##                           | bdiag                    | Matrix
##                           | Diagonal                 | Matrix
##                           | chol                     | Matrix
##                           | cov_marg                 | GeoSimpler
