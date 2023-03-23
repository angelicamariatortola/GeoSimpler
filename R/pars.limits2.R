pars.limits2 <- function (phi = c(lower = 0, upper = +Inf),
                          sigmasq = c(lower = 0, upper = +Inf), 
                          nugget.rel = c(lower = 0, upper = +Inf), 
                          kappa = c(lower = 0, upper = +Inf), 
                          kappa2 = c(lower = 0, upper = +Inf), 
                          lambda = c(lower = -3, upper = 3), 
                          psiR = c(lower = 1, upper = +Inf), 
                          psiA = c(lower = 0, upper = 2 * pi), tausq.rel = nugget.rel,
                          rho = c(lower = -1, upper = 1)) 
{
  if (length(phi) != 2) 
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi")
  if (length(sigmasq) != 2) 
    stop("sigmasq must be a 2 components vector with lower and upper limits for the parameter sigmasq")
  if (length(tausq.rel) != 2) 
    stop("tausq.rel must be a 2 components vector with lower and upper limits for the parameter tausq.rel")
  if (length(kappa) != 2) 
    stop("kappa must be a 2 components vector with lower and upper limits for the parameter kappa")
  if (length(kappa2) != 2) 
    stop("kappa must be a 2 components vector with lower and upper limits for the parameter kappa")
  if (length(lambda) != 2) 
    stop("lambda must be a 2 components vector with lower and upper limits for the parameter lambda")
  if (length(psiR) != 2) 
    stop("psiR must be a 2 components vector with lower and upper limits for the parameter psiR")
  if (length(psiA) != 2) 
    stop("psiA must be a 2 components vector with lower and upper limits for the parameter psiA")
  if (length(rho) != 2) 
    stop("rho must be a 2 components vector with lower and upper limits for the parameter rho")
  
  if (phi[1] >= phi[2]) 
    stop("parameter phi: lower limit greater or equal upper limit")
  if (sigmasq[1] >= sigmasq[2]) 
    stop("parameter sigmasq: lower limit greater or equal upper limit")
  if (tausq.rel[1] >= tausq.rel[2]) 
    stop("parameter tausq.rel: lower limit greater or equal upper limit")
  if (kappa[1] >= kappa[2]) 
    stop("parameter kappa: lower limit greater or equal upper limit")
  if (kappa2[1] >= kappa2[2]) 
    stop("parameter kappa: lower limit greater or equal upper limit")
  if (lambda[1] >= lambda[2]) 
    stop("parameter lambda: lower limit greater or equal upper limit")
  if (psiR[1] >= psiR[2]) 
    stop("parameter psiR: lower limit greater or equal upper limit")
  if (psiA[1] >= psiA[2]) 
    stop("parameter psiA: lower limit greater or equal upper limit")
  if (rho[1] >= rho[2]) 
    stop("parameter rho: lower limit greater or equal upper limit")
  
  names(phi) <- names(sigmasq) <- names(tausq.rel) <- names(kappa) <- 
    names(kappa2) <- names(lambda) <- names(psiR) <- names(psiA) <-  names(rho) <- c("lower", "upper")
  
  return(list(phi = phi, sigmasq = sigmasq, tausq.rel = tausq.rel, 
              kappa = kappa, kappa2 = kappa2, lambda = lambda, psiR = psiR, 
              psiA = psiA, rho = rho))
}
