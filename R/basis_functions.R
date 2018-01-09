#' @name create_basis
#' @rdname create_basis
#' @aliases basis create_basis_function
#'
#' @title Create basis objects
#'
#' @description These functions create different basis objects, which can be
#'   used as input to complex functions in order to perform computations
#'   depending on the class of the basis function.
#'
#' @param M The number of the basis functions. In case of Fourier basis, this
#'   number should be even, since we need to have pairs of sines and cosines and
#'   the constant term is added by default.
#' @param period The period, that is the basis functions are periodic on a
#'   specific interval. Best choice is the range of the points used for
#'   regression.
#' @param gamma Inverse width of radial basis function.
#' @param mus Optional centers of the RBF.
#' @param eq_spaced_mus Logical, if TRUE, equally spaced centers are created,
#'   otherwise centers are created using \code{\link[stats]{kmeans}} algorithm.
#' @param whole_region Logical, indicating if the centers will be evaluated
#'   equally spaced on the whole region, or between the min and max of the
#'   observation values.
#'
#' @return A basis object of class 'rbf', 'polynomial' or 'fourier'.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{eval_functions}}, \code{\link{design_matrix}}
NULL


#' @rdname create_basis
#'
#' @examples
#' (obj <- create_rbf_object(M = 2))
#'
#' #---------------------------------
#'
#' @export
create_rbf_object <- function(M = 2, gamma = NULL, mus = NULL,
                              eq_spaced_mus = TRUE, whole_region = TRUE){
    # Check that M is numberic and integer
    assertthat::assert_that(is.numeric(M))
    assertthat::assert_that(is.logical(eq_spaced_mus))
    assertthat::assert_that(is.logical(whole_region))
    assertthat::assert_that(M %% 1 == 0)
    assertthat::assert_that(M > -1)
    if (!is.null(gamma)) {
        assertthat::assert_that(is.numeric(gamma))
        assertthat::assert_that(gamma > 0)
    }else {gamma <- M ^ 2 / (abs(1) + abs(-1)) ^ 2 }
    if (!is.null(mus)) {
        assertthat::assert_that(is.vector(mus))
        assertthat::assert_that(M == length(mus))
    }else{
        if (eq_spaced_mus) {
            mus <- vector(mode = "numeric", M)
            if (whole_region) {
                for (i in 1:M) { mus[i] <- i * ((1 - (-1)) / (M + 1) ) + (-1) }
            }
        }
    }
    obj <- structure(list(M = M, mus = mus, gamma = gamma,
                          eq_spaced_mus = eq_spaced_mus,
                          whole_region = whole_region),
                     class = "rbf")
    return(obj)
}


#' @rdname create_basis
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



#' @rdname create_basis
#'
#' @examples
#' (obj <- create_fourier_object(M = 2, period = 1))
#'
#' @export
create_fourier_object <- function(M = 2, period = 2){
    # Check that M is positive numberic and integer
    assertthat::assert_that(is.numeric(M))
    assertthat::assert_that(is.numeric(period))
    assertthat::assert_that(M %% 1 == 0)
    assertthat::assert_that(M > -1)
    assertthat::assert_that(period >= 0)
    # Check if we have even number of basis functions
    if (M %% 2 == 1 ) {
        message("Number of basis must be even; increasing by 1.")
        M <- M + 1
    }
    obj <- structure(list(M = M, period = period), class = "fourier")
    return(obj)
}



#------------------------------------------------------------


#' @name eval_functions
#' @rdname eval_functions
#'
#' @title Evaluate basis functions
#'
#' @description Method for evaluating an M basis function model with observation
#'   data \code{obs} and coefficients \code{w}.
#'
#' @param obj The basis function object.
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
#' @seealso \code{\link{create_basis}}, \code{\link{design_matrix}}
NULL


#' @rdname eval_functions
#'
#' @examples
#' # Evaluate the probit transformed basis function values
#' obj <- create_rbf_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_probit_function(obj, obs, w)
#'
#' # -------------------------
#'
#' @export
eval_probit_function <- function(obj, ...){
    f <- eval_function(obj, ...)
    return(pnorm(f))
}


#' @rdname eval_functions
#'
#' @examples
#' # Evaluate the RBF basis function values
#' obj <- create_rbf_object(M=2, mus = c(2,2.5))
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(obj, obs, w)
#'
#' # -------------------------
#'
#' # Evaluate the Polynomial basis function values
#' obj <- create_polynomial_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(obj, obs, w)
#'
#' # -------------------------
#'
#' # Evaluate the Fourier basis function values
#' obj <- create_fourier_object(M=2)
#' obs <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_function(obj, obs, w)
#'
#' @export
eval_function <- function(obj, ...){
    UseMethod("eval_function")
}


# Default function for the generic function 'eval_function'
eval_function.default <- function(obj, ...){
    stop(paste0("Object class '", class(obj), "' is not implemented."))
}


#' @rdname eval_functions
#'
#' @export
eval_function.rbf <- function(obj, obs, w, ...){
    assertthat::assert_that(methods::is(obj, "rbf"))
    assertthat::assert_that(is.vector(w))
    assertthat::assert_that(!is.null(obj$mus))
    f <- rep(w[1], length(obs))
    obs <- as.matrix(obs)
    if (obj$M > 0) {
        for (i in 1:obj$M) {
            f <- f + w[i + 1] * apply(X = obs, MARGIN = 1, FUN = .rbf_basis,
                                      mus = obj$mus[i], gamma = obj$gamma)
        }
    }
    return(f)
}


#' @rdname eval_functions
#'
#' @export
eval_function.polynomial <- function(obj, obs, w, ...){
    assertthat::assert_that(methods::is(obj, "polynomial"))
    assertthat::assert_that(is.vector(obs))
    assertthat::assert_that(is.vector(w))
    f <- rep(w[1], length(obs))
    if (obj$M > 0) {
        for (i in 1:obj$M) { f <- f + w[i + 1] * .polynomial_basis(obs, i) }
    }
    return(f)
}


#' @rdname eval_functions
#'
#' @export
eval_function.fourier <- function(obj, obs, w, ...){
    assertthat::assert_that(methods::is(obj, "fourier"))
    assertthat::assert_that(is.vector(obs))
    assertthat::assert_that(is.vector(w))
    # Create design matrix object
    H <- design_matrix(obj = obj, obs = obs)$H
    # Compute the inner product in order to get the predictions/evaluations
    f <- as.vector(H %*% w)
    return(f)
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
.rbf_basis <- function(X, mus, gamma = 1){
    return(exp( (-1) * gamma * sum( (X - mus) ^ 2) ))
}


# (Internal) Apply the Fourier basis function
#
# Apply the Fourier basis function to the input X.
#
# @param X Input data.
# @param M Integer denoting the degree of the basis function.
# @param period The period of the signal.
#
# @return Input X, after being transformed from the Fourier basis.
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{create_rbf_object}}, \code{\link{polynomial_basis}}
#
.fourier_basis <- function(X, M = 3, period = 2){
    # Compute base frequency
    omega <- 2 * pi / period
    if (M %% 2 == 1) { return(cos((M - 1)/2 * omega * X)) }
    else {return(sin(M/2 * omega * X)) }
}
