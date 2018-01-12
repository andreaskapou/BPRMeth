#' @title Processed ENCODE methylation data
#'
#' @description Small subset of ENCODE methylation data already in pre-processed
#'   format, which are used as a case study for the vignette.
#'
#' @format A list object containing methylation regions and annotation data.
#'   This in general would be the output of the
#'   \code{\link{create_region_object}} function. It has the following two
#'   objects: \itemize{ \item{ \code{met}: A list containing the methylation
#'   regions, each element of the list is a different genomic region.} \item{
#'   \code{anno}: Corresponding annotation data.}}
#'
#' @return Encode methylation data
#'
#' @seealso \code{\link{create_region_object}}, \code{\link{read_met}},
#'   \code{\link{read_anno}}, \code{\link{read_expr}}
"encode_met"


#' @title Processed ENCODE expression data
#'
#' @description Small subset of ENCODE expression data already in pre-processed
#'   format, which are used as a case study for the vignette.
#'
#' @format Expression data in format as the output of \code{\link{read_expr}}
#'   fucntion, i.e. a \code{data.table} with two columns: "id", "expr".
#'
#' @return Encode expression data
#'
#' @seealso \code{\link{create_region_object}}, \code{\link{read_met}},
#'   \code{\link{read_anno}}, \code{\link{read_expr}}
"encode_expr"


#' @title Synthetic Bernoulli data
#'
#' @description A synthetic dataset containinig 300 entries from Bernoulli
#'   probit regression observations (e.g. single cell methylation data).
#'
#' @format A list with 300 elements, where each element element is an L x 2
#'   matrix of observations, where: \describe{ \item{1st column}{locations of
#'   obseravtions} \item{2nd column}{methylation level, 0 - unmethylated, 1 -
#'   methylated} }
#'
#' @return Synthetic Bernoulli methylation data.
#'
#' @seealso \code{\link{binomial_data}}, \code{\link{beta_data}},
#'   \code{\link{gaussian_data}}
"bernoulli_data"


#' @title Synthetic Binomial data
#'
#' @description A synthetic dataset containinig 300 entries from Binomial probit
#'   regression observations (e.g. bulk methylation data).
#'
#' @format A list with 300 elements, where each element element is an L x 3
#'   matrix of observations, where: \describe{ \item{1st column}{locations of
#'   obseravtions} \item{2nd column}{total trials} \item{3rd column}{number of
#'   successes} }
#'
#' @return Synthetic Binomial methylation data.
#'
#' @seealso \code{\link{bernoulli_data}}, \code{\link{beta_data}},
#'   \code{\link{gaussian_data}}
"binomial_data"


#' @title Synthetic Beta data
#'
#' @description A synthetic dataset containinig 300 entries from Beta probit
#'   regression observations (e.g. Beta values from array methylation data).
#'
#' @format A list with 300 elements, where each element element is an L x 2
#'   matrix of observations, where: \describe{ \item{1st column}{locations of
#'   obseravtions} \item{2nd column}{methylation level } }
#'
#' @return Synthetic Beta methylation data.
#'
#' @section Details: The beta regression model is based on alternative
#'   parameterization of the beta density in terms of the mean and dispersion
#'   parameter:\url{https://cran.r-project.org/web/packages/betareg/}.
#'
#' @seealso \code{\link{binomial_data}}, \code{\link{bernoulli_data}},
#'   \code{\link{gaussian_data}}
"beta_data"


#' @title Synthetic Gaussian data
#'
#' @description A synthetic dataset containinig 300 entries from Gaussian (i.e.
#'   linear) regression observations (e.g. M-values from array methylation
#'   data).
#'
#' @format A list with 300 elements, where each element element is an L x 2
#'   matrix of observations, where: \describe{ \item{1st column}{locations of
#'   obseravtions} \item{2nd column}{methylation level} }
#'
#' @return Synthetic Bernoulli methylation data.
#'
#' @seealso \code{\link{binomial_data}}, \code{\link{beta_data}},
#'   \code{\link{bernoulli_data}}
"gaussian_data"


#' @title Synthetic bulk methylation data
#'
#' @description A synthetic dataset containinig 600 entries.
#'
#' @format A list with 600 elements, where each element element is an L x 3
#'   matrix of observations, where: \describe{ \item{1st column}{locations of
#'   obseravtions} \item{2nd column}{total trials} \item{3rd column}{number of
#'   successes} }
#'
#' @return Synthetic bulk methylation data
#'
#' @seealso \code{\link{gex_data}}, \code{\link{bernoulli_data}},
#'   \code{\link{binomial_data}}, \code{\link{beta_data}},
#'   \code{\link{gaussian_data}}
"meth_data"


#' @title Synthetic expression data
#'
#' @description Corresponding gene expression data for the
#'   \code{\link{meth_data}}
#'
#' @format A vector of length 600
#'
#' @return Synthetic gene expression data
#'
#' @seealso \code{\link{meth_data}}, \code{\link{bernoulli_data}},
#'   \code{\link{binomial_data}}, \code{\link{beta_data}},
#'   \code{\link{gaussian_data}}
"gex_data"
