% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovSimpler.R
\name{CovSimpler}
\alias{CovSimpler}
\title{Calculates Value of the Simpler Covariance Function}
\usage{
CovSimpler(
  coords = NULL,
  dists.lowertri = NULL,
  cov.model = "exp",
  cov.pars,
  p = 1,
  nugget = rep(0, length(cov.model)),
  SigmaB = NULL
)
}
\arguments{
\item{varcov}{the covariance matrix.}
}
\description{
Computes the Simpler covariance function associated to a p-dimensional
multivariate Gaussian random field.
Different correlation functions are available to model marginal-covariance behaviors.
When \eqn{p = 1} the function reduces to the marginal-covariance function (univariate case).
}
