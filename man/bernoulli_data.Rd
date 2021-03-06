% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bernoulli_data}
\alias{bernoulli_data}
\title{Synthetic Bernoulli data}
\format{
A list with 300 elements, where each element element is an L x 2
  matrix of observations, where: \describe{ \item{1st column}{locations of
  obseravtions} \item{2nd column}{methylation level, 0 - unmethylated, 1 -
  methylated} }
}
\usage{
bernoulli_data
}
\value{
Synthetic Bernoulli methylation data.
}
\description{
A synthetic dataset containinig 300 entries from Bernoulli
  probit regression observations (e.g. single cell methylation data).
}
\seealso{
\code{\link{binomial_data}}, \code{\link{beta_data}},
  \code{\link{gaussian_data}}
}
\keyword{datasets}
