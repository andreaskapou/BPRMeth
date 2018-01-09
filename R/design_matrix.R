#' @name design_matrix
#' @rdname design_matrix
#' @aliases designmatrix des_matrix des_mat
#'
#' @title Generic function for creating design matrices
#'
#' @description These functions call the appropriate methods depending on the
#'   class of the object \code{obj} to create RBF, polynomial or Fourier design
#'   matrices.
#'
#' @param obj A basis function object.
#' @param obs A vector of observations.
#' @param ... Additional parameters.
#'
#' @return A design matrix object
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{eval_functions}}
#'
#' @export
NULL


#' @rdname design_matrix
#'
#' @examples
#' obj <- create_polynomial_object(M=2)
#' obs <- c(0,.2,.5)
#' poly <- design_matrix(obj, obs)
#'
#' #----------------
#'
#' obj <- create_rbf_object(M=2)
#' obs <- c(0,.2,.5)
#' rbf <- design_matrix(obj, obs)
#'
#' #----------------
#'
#' obj <- create_fourier_object(M=2)
#' obs <- c(0,.2,.5)
#' fourier <- design_matrix(obj, obs)
#'
#' @export
design_matrix <- function(obj, ...){
    UseMethod("design_matrix")
}


#' @rdname design_matrix
#'
design_matrix.default <- function(obj, ...){
    stop(paste0("Object type '", class(obj), "' is not implemented."))
}


#' @rdname design_matrix
#'
#' @export
design_matrix.polynomial <- function(obj, obs, ...){
    assertthat::assert_that(methods::is(obj, "polynomial"))
    assertthat::assert_that(is.vector(obs))
    N <- length(obs)  # Length of the dataset
    H <- matrix(1, nrow = N, ncol = obj$M + 1)
    # Compute X^(j)
    if (obj$M > 0) {for (j in 1:obj$M) {H[, j + 1] <- .polynomial_basis(obs,j)}}
    return(list(H = H, basis = obj))
}


#' @rdname design_matrix
#'
#' @export
design_matrix.rbf <- function(obj, obs, ...){
    assertthat::assert_that(methods::is(obj, "rbf"))
    assertthat::assert_that(is.vector(obs))
    N   <- length(obs)  # Length of the dataset
    if (obj$M == 0) { H <- matrix(1, nrow = N, ncol = 1); obj$mus <- 0;
    }else{
        if (is.null(obj$mus)) {
            if (obj$eq_spaced_mus) {
                obj$mus <- vector(mode = "numeric", obj$M)
                if (!obj$whole_region) {
                    # TODO: Keep this as functionality?
                    for (i in 1:obj$M) {
                        obj$mus[i] <- i*(max(obs) - min(obs))/(obj$M + 1) +
                            min(obs)
                    }
                }
            }else {
                repeat {
                    # TODO: Keep this as functionality?
                    km <- stats::kmeans(obs, obj$M, iter.max = 30, nstart = 10)
                    if (min(km$size) > 0) { break } # Accept non-empty clusters
                }
                obj$mus <- km$centers  # RBF centers
            }
        }
        # Convert the 'obs' vector to an N x 1 dimensional matrix
        obs <- as.matrix(obs)
        H <- matrix(1, nrow = N, ncol = obj$M + 1)
        for (j in 1:obj$M) {
            H[, j + 1] <- apply(obs,1,.rbf_basis,mus = obj$mus[j],
                                gamma = obj$gamma)
        }
    }
    return(list(H = H, basis = obj))
}


#' @rdname design_matrix
#'
#' @export
design_matrix.fourier <- function(obj, obs, ...){
    assertthat::assert_that(methods::is(obj, "fourier"))
    assertthat::assert_that(is.vector(obs))
    # Similar implementation to the FDA package.
    # Compute base frequency
    omega <- 2 * pi / obj$period
    # Initialize design matrix
    H <- matrix(1 / sqrt(2), nrow = length(obs), ncol = obj$M + 1)
    if (obj$M > 1) {
        j <- seq(2, obj$M, 2)           # Get the basis
        k <- j / 2
        evals <- outer(omega * obs, k)  # Compute outer product
        H[, j] <- sin(evals)            # Make the sine evaluations
        H[, j + 1] <- cos(evals)        # Make the cosine evaluations
    }
    # Divide by this constant to get actual magnitude
    # TODO: Check this step
    H <- H / sqrt(obj$period / 2)
    return(list(H = H, basis = obj))
}
