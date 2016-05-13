#' Generic function for evaluating probit basis functions
#'
#' Method for evaluating the probit transformation of a basis function.
#'
#' @param x The basis function object.
#' @param ... Additional parameters that will be passed to more specific
#'   functions.
#'
#' @return The probit transformed basis function values.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{eval_function}}, \code{\link{eval_function.rbf}},
#'   \code{\link{eval_function.polynomial}}
#'
#' @examples
#' x <- create_rbf_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_probit_function(x, obs, w)
#'
#' @export
eval_probit_function <- function(x, ...){
  f <- eval_function(x, ...)
  return(pnorm(f))
}


#' Generic function for evaluating basis functions
#'
#' Method for evaluating a basis function object of class x.
#'
#' @inheritParams eval_probit_function
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{eval_probit_function}}, \code{\link{eval_function.rbf}},
#'   \code{\link{eval_function.polynomial}}
#'
#' @examples
#' x <- create_rbf_object(M=2, mus = c(2,2.5))
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(x, obs, w)
#'
#' #----------------
#'
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


#' Evaluate polynomial function
#'
#' Method for evaluating the polynomial function of degree M for observation
#' data obs and coefficients w.
#'
#' @param x The basis function object.
#' @param obs Input / observation data.
#' @param w Vector of length M, containing the coefficients of an Mth-order
#'   basis function.
#' @param ... Optional additional parameters
#'
#' @return The polynomial function values.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' x <- create_polynomial_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(x, obs, w)
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


#' Evaluate rbf function
#'
#' Method for evaluating the rbf function with M basis for observation data
#' obs and coefficients w.
#'
#' @inheritParams eval_function.polynomial
#'
#' @return The rbf function values.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' x <- create_rbf_object(M=2, mus = c(2,2.5))
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(x, obs, w)
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
