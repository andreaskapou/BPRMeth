#' Create a polynomial object
#'
#' \code{create_polynomial_object} creates objects of class "polynomial". These
#' objects can be used as input to complex functions in order to perform
#' computations depending on the class of  basis function.
#'
#' @param M The degree of the polynomial object.
#'
#' @return An object of class 'polynomial'
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_rbf_object}}
#'
#' @examples
#' (obj <- create_polynomial_object(M = 2))
#'
#' @export
create_polynomial_object <- function(M = 1){
  # Check that M is numberic and integer
  assertthat::assert_that(is.numeric(M))
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  obj <- structure(list(M = M), class = "polynomial")
  return(obj)
}


#' Create an RBF object
#'
#' \code{create_rbf_object} creates objects of class "rbf". These objects can be
#' used as input to complex functions in order to perform computations depending
#' on the basis function.
#'
#' @param M The number of the RBFs.
#' @param gamma Inverse width of radial basis function.
#' @param mus Optional centers of the RBF.
#' @param eq_spaced_mus Logical, if TRUE, equally spaced centers are created,
#'   otherwise centers are created using \code{\link[stats]{kmeans}} algorithm.
#' @param whole_region Logical, indicating if the centers will be evaluated
#'   equally spaced on the whole region, or between the min and max of the
#'   observation values.
#'
#' @return An object of class 'rbf'.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_polynomial_object}}
#'
#' @examples
#' (obj <- create_rbf_object(M = 2))
#'
#' @export
create_rbf_object <- function(M = 2, gamma = NULL, mus = NULL, eq_spaced_mus = TRUE,
                                                         whole_region = TRUE){
  assertthat::assert_that(is.numeric(M))
  assertthat::assert_that(is.logical(eq_spaced_mus))
  assertthat::assert_that(is.logical(whole_region))
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  if (! is.null(gamma)){
    assertthat::assert_that(is.numeric(gamma))
    assertthat::assert_that(gamma > 0)
  }else{
    gamma <- M ^ 2 / ( abs(1) + abs(-1) ) ^ 2
  }
  if (! is.null(mus)){
    assertthat::assert_that(is.vector(mus))
    assertthat::assert_that(M == length(mus))
  }else{
    if (eq_spaced_mus){
      mus <- vector(mode = "numeric", M)
      if (whole_region){
        for (i in 1:M){
          mus[i] <- i * ( (1 - (-1)) / (M + 1) ) + (-1)
        }
      }
    }
  }
  obj <- structure(list(M = M,
                        mus = mus,
                        gamma = gamma,
                        eq_spaced_mus = eq_spaced_mus,
                        whole_region = whole_region),
                   class = "rbf")
  return(obj)
}


#------------------------------------------------------------


# (Internal) Apply polynomial basis function.
#
# Apply the polynomial basis function of degree M to the input X.
#
# @param X The input data, either a scalar, vector or matrix.
# @param M Integer, denoting the degree of the polynomial basis that will be
#   applied to the input X. M should not be negative.
#
# @return Input X, after being transformed from the polynomial basis function.
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{create_polynomial_object}}, \code{\link{rbf_basis}}
#
# @examples
# out <- polynomial_basis(X = 3, M = 2)
# #
# out <- polynomial_basis(X = c(2, 3, 4), M = 2)
#
.polynomial_basis <- function(X, M = 1){
  return(X ^ M)
}


# (Internal) Apply the radial basis function
#
# Apply the RBF function to the input X.
#
# @param X Input data.
# @param mus Centers from where we should compute the distance of the data X.
# @param gamma Inverse width of radial basis function.
#
# @return Input X, after being transformed from the RBF.
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{create_rbf_object}}, \code{\link{polynomial_basis}}
#
# @examples
# out <- rbf_basis(X = c(1,2), mus = c(1,1))
# #
# out <- rbf_basis(X = c(1,2), mus = c(1,1), gamma = 0.1)
#
.rbf_basis <- function(X, mus, gamma = 1){
  return(exp( (-1) * gamma * sum((X - mus)^2)))
}
