#' Synthetic data for BPRMeth package
#'
#' A synthetic dataset containinig 600 entries.
#'
#' @format A list with 600 elements, where each element element is an L x 3
#'  matrix of observations, where:
#' \describe{
#'   \item{1st column}{locations of obseravtions}
#'   \item{2nd column}{total trials at corresponding locations}
#'   \item{3rd column}{number of successes at corresponding locations}
#' }
#'
#' @return Synthetic methylation data
"meth_data"

#' Synthetic data for mpgex package
#'
#' Corresponding gene expression data for the 'meth_data'
#'
#' @format A vector of length 600
#'
#' @return Synthetic gene expression data
"gex_data"
