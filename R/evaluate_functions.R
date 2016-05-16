#' @name eval_functions
#' @rdname eval_functions
#'
#' @title Evaluate basis functions
#'
#' @description Method for evaluating an M basis function model with observation
#'   data \code{obs} and coefficients \code{w}.
#'
#' @param x The basis function object.
#' @param obs Observation data.
#' @param w Vector of length M, containing the coefficients of an
#'   M\eqn{^{th}}-order basis function.
#' @param ... Optional additional parameters
#'
#' @return The evaluated function values.
#'
#'   NOTE that the \code{eval_probit_function} computes the probit transformed
#'   basis function values.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}
NULL


#' @rdname eval_functions
#'
#' @examples
#' # Evaluate the probit transformed basis function values
#' x <- create_rbf_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_probit_function(x, obs, w)
#'
#' # -------------------------
#'
#' @export
eval_probit_function <- function(x, ...){
  f <- eval_function(x, ...)
  return(pnorm(f))
}


#' @rdname eval_functions
#'
#' @examples
#' # Evaluate the RBF basis function values
#' x <- create_rbf_object(M=2, mus = c(2,2.5))
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(x, obs, w)
#'
#' # -------------------------
#'
#' # Evaluate the Polynomial basis function values
#' x <- create_polynomial_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(x, obs, w)
#'
#' @export
eval_function <- function(x, ...){
  UseMethod("eval_function")
}


# Default function for the generic function 'eval_function'
eval_function.default <- function(x, ...){
  stop(paste("Object of class '", class(x), "' is not implemented.", sep = ""))
}


#' @rdname eval_functions
#'
#' @export
eval_function.rbf <- function(x, obs, w, ...){
  assertthat::assert_that(methods::is(x, "rbf"))
  assertthat::assert_that(is.vector(w))
  assertthat::assert_that(!is.null(x$mus))

  f <- rep(w[1], length(obs))
  obs <- as.matrix(obs)
  if (x$M > 0){
    for (i in 1:x$M){
      f <- f + w[i + 1] * apply(X      = obs,
                                MARGIN = 1,
                                FUN    = .rbf_basis,
                                mus    = x$mus[i],
                                gamma  = x$gamma)
    }
  }
  return(f)
}


#' @rdname eval_functions
#'
#' @export
eval_function.polynomial <- function(x, obs, w, ...){
  assertthat::assert_that(methods::is(x, "polynomial"))
  assertthat::assert_that(is.vector(obs))
  assertthat::assert_that(is.vector(w))

  f <- rep(w[1], length(obs))
  if (x$M > 0){
    for (i in 1:x$M){
      f <- f + w[i + 1] * .polynomial_basis(obs, i)
    }
  }
  return(f)
}
