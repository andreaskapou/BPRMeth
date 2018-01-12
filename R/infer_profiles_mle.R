#' @name infer_profiles_mle
#' @rdname infer_profiles_mle
#' @aliases infer_profile_mle inference_mle
#'
#' @title Infer methylation profiles using MLE
#'
#' @description General purpose functions for inferring latent profiles for
#'   different observation models using maximum likelihood estimation (MLE).
#'   Current observation models are: 'bernoulli', 'binomial', 'beta' or
#'   'gaussian'. For most models we cannot obtain an analytically tractable
#'   solution, hence an optimization procedure is used. The
#'   \code{\link[stats]{optim}} package is used for performing optimization.
#'
#' @param X The input data, either a \code{\link[base]{matrix}} or a
#'   \code{\link[base]{list}} of elements of length N, where each element is an
#'   \code{L X C} matrix, where L are the total number of observations. The
#'   first column contains the input observations x (i.e. CpG locations). If
#'   "binomial" model then C=3, and 2nd and 3rd columns contain total number of
#'   trials and number of successes respectively. If "bernoulli" or "gaussian"
#'   model, then C=2 containing the output y (e.g. methylation level). If "beta"
#'   model, then C=3, where 2nd column contains output y and 3rd column the
#'   dispersion parameter.
#' @param H Optional, design matrix of the input data X. If NULL, H will be
#'   computed inside the function.
#' @param model Observation model name as character string. It can be either
#'   'bernoulli', 'binomial', 'beta' or 'gaussian'.
#' @param basis A 'basis' object. E.g. see \code{\link{create_basis}}. If NULL,
#'   will an RBF object will be created.
#' @param lambda The complexity penalty coefficient for ridge regression.
#' @param w A vector of initial parameters (i.e. coefficients of the basis
#'   functions).
#' @param beta_dispersion Dispersion parameter, only used for Beta distribution
#'   and will be the same for all observations.
#' @param opt_method The optimization method to be used. See
#'   \code{\link[stats]{optim}} for possible methods. Default is "CG".
#' @param opt_itnmax Optional argument giving the maximum number of iterations
#'   for the corresponding method. See \code{\link[stats]{optim}} for details.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param ... Additional parameters.
#'
#' @return An object of class \code{infer_profiles_mle_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   optimized parameter values. Each row of the matrix corresponds to each
#'   element of the list X; if X is a matrix, then N = 1. The columns are of the
#'   same length as the parameter vector w (i.e. number of basis functions). }
#'   \item{ \code{basis}: The basis object. } \item{\code{nll_feat}: NLL fit
#'   feature.} \item{\code{rmse_feat}: RMSE fit feature.}
#'   \item{\code{coverage_feat}: CpG coverage feature.}}
#'
#' @section Details: The beta regression model is based on alternative
#'   parameterization of the beta density in terms of the mean and dispersion
#'   parameter: \url{https://cran.r-project.org/web/packages/betareg/} . For
#'   modelling details for Binomial/Bernoulli observation model check the paper
#'   for BPRMeth:
#'   \url{https://academic.oup.com/bioinformatics/article/32/17/i405/2450762} .
#'
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{infer_profiles_vb}},
#'   \code{\link{infer_profiles_gibbs}}, \code{\link{create_region_object}}
NULL


#' @rdname infer_profiles_mle
#'
#' @examples
#' # Example of optimizing parameters for synthetic data using 3 RBFs
#' basis <- create_rbf_object(M=3)
#' out <- infer_profiles_mle(X = binomial_data, model = "binomial",
#'    basis = basis, is_parallel = FALSE, opt_itnmax = 10)
#'
#' #-------------------------------------
#'
#' basis <- create_rbf_object(M=3)
#' out <- infer_profiles_mle(X = beta_data, model = "beta",
#'    basis = basis, is_parallel = FALSE, opt_itnmax = 10)
#'
#' #-------------------------------------
#'
#' basis <- create_rbf_object(M=3)
#' out <- infer_profiles_mle(X = gaussian_data[[1]], model = "gaussian",
#'    basis = basis, is_parallel = FALSE, opt_itnmax = 10)
#'
#' @export
infer_profiles_mle <- function(X, model = NULL, basis = NULL, H = NULL,
                               lambda = .5, w = NULL, beta_dispersion = 5,
                               opt_method = "CG", opt_itnmax = 100,
                               is_parallel = FALSE, no_cores = NULL, ...){
    if (is.null(model)) { stop("Observation model not defined!") }
    # Create RBF basis object by default
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }
    if (is.null(H)) { # Create design matrix
        if (is.list(X)) { H <- lapply(X,function(x)
            design_matrix(basis,x[,1])$H)
        }else {H <- design_matrix(basis, X[,1])$H }
    }

    if (identical(model, "bernoulli") || identical(model, "binomial") ) {
        obj <- .infer_profiles_mle_bpr(X, H, basis, lambda, w, opt_method,
                                       opt_itnmax, is_parallel, no_cores)
    }else if (identical(model, "beta") ) {
        obj <- .infer_profiles_mle_beta(X, H, basis, lambda, w, beta_dispersion,
                                        opt_method, opt_itnmax, is_parallel,
                                        no_cores)
    }else if (identical(model, "gaussian") ) {
        obj <- .infer_profiles_mle_gaussian(X, H, basis, lambda, is_parallel,
                                            no_cores)
    }else{stop(paste0(model, " observation model not implented.")) }
    # Append general class infer_profiles
    class(obj) <- append(class(obj), c("infer_profiles_mle", "infer_profiles"))
    # Return object
    return(obj)
}


##------------------------------------------------------

# S3 method
.infer_profiles_mle_bpr <- function(X, ...){
    UseMethod(".infer_profiles_mle_bpr")
}

# Default function for the generic function 'infer_profiles_mle_bpr'
.infer_profiles_mle_bpr.default <- function(X, ...){
    stop("Object X should be either matrix or list!")
}

# Infer profiles Binomial/Bernoulli MLE list
.infer_profiles_mle_bpr.list <- function(X, H, basis, lambda, w, opt_method,
                                         opt_itnmax, is_parallel, no_cores,...){
    assertthat::assert_that(is.list(X)) # Check that X is a list object
    assertthat::assert_that(is.matrix(X[[1]]))
    i <- 0                  # Initialize so the CMD check passes without NOTES
    N <- length(X)          # Extract number of observations
    assertthat::assert_that(N > 0)
    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel)

    if (is_parallel) { # If parallel mode is ON create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        # Parallel optimization for each element of X, i.e. for each region i.
        res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- .infer_profiles_mle_bpr.matrix(
              X = X[[i]], H = H[[i]], basis = basis, lambda = lambda, w = w,
              opt_method = opt_method, opt_itnmax = opt_itnmax) })
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of X, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
           ex = {out_opt <- .infer_profiles_mle_bpr.matrix(
               X = X[[i]], H = H[[i]], basis = basis, lambda = lambda, w = w,
               opt_method = opt_method, opt_itnmax = opt_itnmax) })
    }

    # Matrix for storing optimized coefficients
    W <- sapply(res, function(x) x$W)
    if (is.matrix(W)) { W <- t(W)
    }else {W <- as.matrix(W) }
    colnames(W) <- paste("w", seq(1, NCOL(W)), sep = "")
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))

    # TODO : Remove this part
    # Matrix for storing the centers of RBFs if object class is 'rbf'
    # Mus <- NULL
    # if (methods::is(basis, "rbf")){
    #     if (is.null(basis$mus)){
    #         Mus <- sapply(lapply(res, function(x) x$basis), function(y) y$mus)
    #         if (is.matrix(Mus)){ Mus <- t(Mus)
    #         }else{ Mus <- as.matrix(Mus) }
    #         colnames(Mus) <- paste("mu", seq(1, NCOL(Mus)), sep = "")
    #     }
    # }

    # Store the object
    obj <- list(W = W, basis = basis, nll_feat = nll_feat,
                rmse_feat = rmse_feat, coverage_feat = coverage_feat)
    if (NCOL(X[[1]]) == 3) { class(obj) <- "infer_profiles_mle_binomial"
    }else {class(obj) <- "infer_profiles_mle_bernoulli" }
    return(obj)
}

# Infer profiles Binomial/Bernoulli MLE matrix
.infer_profiles_mle_bpr.matrix <- function(X, H, basis, lambda, w, opt_method,
                                           opt_itnmax, ...){
    if (is.null(w)) { w <- rep(0.5, NCOL(H)) }
    # Call optim function to perform minimization of the NLL of BPR function
    W <- c(optim(par = w, fn = bpr_log_likelihood, gr = bpr_gradient,
                   method = opt_method, control = list(maxit = opt_itnmax),
                   X = X, H = H, lambda = lambda, is_nll = TRUE)$par)
    # NLL fit feature
    nll_feat <- bpr_log_likelihood(w = W, X = X, H = H, lambda = lambda,
                                   is_nll = TRUE)
    # RMSE fit feature
    f_pred <- c(pnorm(H %*% W))
    if (NCOL(X) == 3) { f_true <- X[, 3] / X[, 2] } else {f_true <- X[, 2] }
    rmse_feat <- sqrt( mean( (f_pred - f_true)^2) )
    # CpG coverage feature
    coverage_feat <- NROW(X)
    # Make weight vector as row vector
    W <- matrix(W, nrow = 1)
    # Store the object
    obj <- list(W = W, basis = basis, nll_feat = nll_feat,
                rmse_feat = rmse_feat, coverage_feat = coverage_feat)
    if (NCOL(X) == 3) { class(obj) <- "infer_profiles_mle_binomial"
    }else {class(obj) <- "infer_profiles_mle_bernoulli" }
    return(obj)
}


##------------------------------------------------------

# S3 method
.infer_profiles_mle_beta <- function(X, ...){
    UseMethod(".infer_profiles_mle_beta")
}

# Default function for the generic function 'infer_profiles_mle_beta'
.infer_profiles_mle_beta.default <- function(X, ...){
    stop("Object X should be either matrix or list!")
}

# Infer profiles Beta MLE list
.infer_profiles_mle_beta.list <- function(X,H,basis,lambda,w,beta_dispersion,
                                          opt_method, opt_itnmax, is_parallel,
                                          no_cores, ...){
    assertthat::assert_that(is.list(X)) # Check that x is a list object
    i <- 0                  # Initialize so the CMD check passes without NOTES
    N <- length(X)          # Extract number of observations
    assertthat::assert_that(N > 0)
    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel)


    if (is_parallel) { # If parallel mode is ON create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        # Parallel optimization for each element of X, i.e. for each region i.
        res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- .infer_profiles_mle_beta.matrix(
              X = X[[i]], H = H[[i]], basis = basis, lambda = lambda, w = w,
              beta_dispersion = beta_dispersion, opt_method = opt_method,
              opt_itnmax = opt_itnmax) })
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of X, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
           ex = {out_opt <- .infer_profiles_mle_beta.matrix(
               X = X[[i]], H = H[[i]], basis = basis, lambda = lambda, w = w,
               beta_dispersion = beta_dispersion, opt_method = opt_method,
               opt_itnmax = opt_itnmax) })
    }

    # Matrix for storing optimized coefficients
    W <- sapply(res, function(x) x$W)
    if (is.matrix(W)) { W <- t(W)
    }else {W <- as.matrix(W) }
    colnames(W) <- paste("w", seq(1, NCOL(W)), sep = "")
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))

    # Store the object
    obj <- structure(list(W = W, basis = basis, nll_feat = nll_feat,
                          rmse_feat = rmse_feat, coverage_feat = coverage_feat),
                     class = "infer_profiles_mle_beta")
    return(obj)
}

# Infer profiles Beta MLE matrix
.infer_profiles_mle_beta.matrix <- function(X,H,basis,lambda,w,beta_dispersion,
                                            opt_method, opt_itnmax,...){
    # Concatenate the dispersion parameter
    X <- cbind(X, rep(beta_dispersion, NROW(X)))
    if (is.null(w)) { w <- rep(0.5, NCOL(H)) }

    # Call optim function to perform minimization of the NLL of BPR function
    W <- c(optim(par = w, fn = betareg_log_likelihood, gr = betareg_gradient,
                   method = opt_method, control = list(maxit = opt_itnmax),
                   X = X, H = H, lambda = lambda, is_nll = TRUE)$par)
    # NLL fit feature
    nll_feat <- betareg_log_likelihood(w = W, X = X, H = H, lambda = lambda,
                                       is_nll = TRUE)
    # RMSE fit feature
    f_pred <- c(pnorm(H %*% W))
    f_true <- X[, 2]
    rmse_feat <- sqrt( mean( (f_pred - f_true)^2) )
    # CpG coverage feature
    coverage_feat <- NROW(X)
    # Make weight vector as row vector
    W <- matrix(W, nrow = 1)
    # Store the object
    obj <- structure(list(W = W, basis = basis, nll_feat = nll_feat,
                          rmse_feat = rmse_feat, coverage_feat = coverage_feat),
                     class = "infer_profiles_mle_beta")
    return(obj)
}


##------------------------------------------------------

# S3 method
.infer_profiles_mle_gaussian <- function(X, ...){
    UseMethod(".infer_profiles_mle_gaussian")
}

# Default function for the generic function 'infer_profiles_mle_gaussian'
.infer_profiles_mle_gaussian.default <- function(X, ...){
    stop("Object X should be either matrix or list!")
}


# Infer profiles Gaussian MLE list
.infer_profiles_mle_gaussian.list <- function(X,H,basis,lambda,is_parallel,
                                              no_cores,...){
    assertthat::assert_that(is.list(X)) # Check that X is a list object
    i <- 0                  # Initialize so the CMD check passes without NOTES
    N <- length(X)          # Extract number of observations
    assertthat::assert_that(N > 0)
    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel)


    if (is_parallel) { # If parallel mode is ON create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        # Parallel optimization for each element of X, i.e. for each region i.
        res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
           ex = {out_opt <- .infer_profiles_mle_gaussian.matrix(
                  X = X[[i]], H = H[[i]], basis = basis, lambda = lambda) })
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of X, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
           ex = {out_opt <- .infer_profiles_mle_gaussian.matrix(
                  X = X[[i]], H = H[[i]], basis = basis, lambda = lambda) })
    }

    # Matrix for storing optimized coefficients
    W <- sapply(res, function(x) x$W)
    if (is.matrix(W)) { W <- t(W)
    }else {W <- as.matrix(W) }
    colnames(W) <- paste("w", seq(1, NCOL(W)), sep = "")
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))

    # Store the object
    obj <- structure(list(W = W, basis = basis, nll_feat = nll_feat,
                          rmse_feat = rmse_feat, coverage_feat = coverage_feat),
                     class = "infer_profiles_mle_gaussian")
    return(obj)
}

# Infer profiles Binomial MLE matrix
.infer_profiles_mle_gaussian.matrix <- function(X, H, basis, lambda, ...){
    if (lambda == 0) {
        qx <- qr(H)                 # Compute QR decomposition of H
        W <- c(solve.qr(qx, X[,2])) # Compute (H'H)^(-1)H'y
    }else{
        I <- diag(1, NCOL(H))  # Identity matrix
        # TODO: Should we do this or not for the bias term??
        # I[1,1]  <- 1e-10  # Do not change the intercept coefficient
        qx <- qr(lambda * I + t(H) %*% H)    # Compute QR decomposition
        W <- c(solve.qr(qx, t(H) %*% X[,2])) # Comp (lambda*I+H'H)^(-1)H'y
    }

    # NLL fit feature
    nll_feat <- lr_log_likelihood(w = W, X = X, H = H, lambda = lambda,
                                  is_nll = TRUE)
    # RMSE fit feature
    f_pred <- c(H %*% W)
    f_true <- X[, 2]
    rmse_feat <- sqrt( mean( (f_pred - f_true)^2) )
    # CpG coverage feature
    coverage_feat <- NROW(X)
    # Make weight vector as row vector
    W <- matrix(W, nrow = 1)
    # Store the object
    obj <- structure(list(W = W, basis = basis, nll_feat = nll_feat,
                          rmse_feat = rmse_feat, coverage_feat = coverage_feat),
                     class = "infer_profiles_mle_gaussian")
    return(obj)
}
