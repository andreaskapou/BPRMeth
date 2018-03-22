#' @name infer_profiles_vb
#' @rdname infer_profiles_vb
#' @aliases infer_profile_vb inference_vb
#'
#' @title Infer methylation profiles using VB
#'
#' @description General purpose functions for inferring latent profiles for
#'   different observation models using Variational Bayes (VB). Current
#'   observation models are: 'bernoulli', 'binomial' or 'gaussian'.
#'
#' @param X The input data, either a \code{\link[base]{matrix}} or a
#'   \code{\link[base]{list}} of elements of length N, where each element is an
#'   \code{L X C} matrix, where L are the total number of observations. The
#'   first column contains the input observations x (i.e. CpG locations). If
#'   "binomial" model then C=3, and 2nd and 3rd columns contain total number of
#'   trials and number of successes respectively. If "bernoulli" or "gaussian"
#'   model, then C=2 containing the output y (e.g. methylation level).
#' @param H Optional, design matrix of the input data X. If NULL, H will be
#'   computed inside the function.
#' @param model Observation model name as character string. It can be either
#'   'bernoulli', 'binomial', 'beta' or 'gaussian'.
#' @param basis A 'basis' object. E.g. see \code{\link{create_basis}}. If NULL,
#'   will an RBF object will be created.
#' @param w A vector of initial parameters (i.e. coefficients of the basis
#'   functions). If NULL, it will be initialized inside the function.
#' @param gaussian_l Noise precision parameter, only used when having "gaussian"
#'   observation model.
#' @param alpha_0 Hyperparameter: shape parameter for Gamma distribution. A
#'   Gamma distribution is used as prior for the precision parameter tau.
#' @param beta_0 Hyperparameter: rate parameter for Gamma distribution. A Gamma
#'   distribution is used as prior for the precision parameter tau.
#' @param vb_max_iter Integer denoting the maximum number of VB iterations.
#' @param epsilon_conv Numeric denoting the convergence threshold for VB.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_verbose Logical, print results during VB iterations.
#' @param ... Additional parameters.
#'
#' @return An object of class \code{infer_profiles_vb_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   optimized parameter values. Each row of the matrix corresponds to each
#'   element of the list X; if X is a matrix, then N = 1. The columns are of the
#'   same length as the parameter vector w (i.e. number of basis functions). }
#'   \item{ \code{W_Sigma}: A list with covariance matrices for each element row
#'   in W.} \item{ \code{basis}: The basis object. } \item{\code{nll_feat}: NLL
#'   fit feature.} \item{\code{rmse_feat}: RMSE fit feature.}
#'   \item{\code{coverage_feat}: CpG coverage feature.} \item{\code{lb_feat}:
#'   Lower Bound feature.}}
#'
#' @section Details: The modelling and mathematical details for inferring
#'   profiles using mean-field variational inference are explained here:
#'   \url{http://rpubs.com/cakapourani/} . More specifically: \itemize{\item{For
#'   Binomial/Bernoulli observation model check:
#'   \url{http://rpubs.com/cakapourani/variational-bayes-bpr}} \item{For
#'   Gaussian observation model check:
#'   \url{http://rpubs.com/cakapourani/variational-bayes-lr}}}
#'
#' @importFrom matrixcalc matrix.trace
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{infer_profiles_mle}},
#'   \code{\link{predict_expr}}, \code{\link{create_region_object}}
NULL


#' @rdname infer_profiles_vb
#'
#' @examples
#' # Example of inferring parameters for synthetic data using 3 RBFs
#' basis <- create_rbf_object(M=3)
#' out <- infer_profiles_vb(X = binomial_data, model = "binomial",
#'    basis = basis, is_parallel = FALSE, vb_max_iter = 10)
#'
#' #-------------------------------------
#'
#' basis <- create_rbf_object(M=3)
#' out <- infer_profiles_vb(X = gaussian_data, model = "gaussian",
#'    basis = basis, is_parallel = FALSE, vb_max_iter = 10)
#'
#' @export
infer_profiles_vb <- function(X, model = NULL, basis = NULL, H = NULL, w = NULL,
                              gaussian_l = 50, alpha_0 = .5, beta_0 = .1,
                              vb_max_iter = 100, epsilon_conv = 1e-5,
                              is_parallel = FALSE, no_cores = NULL,
                              is_verbose = FALSE, ...){
    if (is.null(model)) { stop("Observation model not defined!") }
    # Create RBF basis object by default
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }
    if (is.null(H)) { # Create design matrix
        if (is.list(X)) {H <- lapply(X,function(x) design_matrix(basis,x[,1])$H)
        }else {H <- design_matrix(basis, X[,1])$H }
    }

    if (identical(model, "bernoulli") || identical(model, "binomial") ) {
        obj <- .infer_profiles_vb_bpr(X, H, basis, w, alpha_0, beta_0,
                                      vb_max_iter,epsilon_conv,is_parallel,
                                      no_cores, is_verbose)
    }else if (identical(model, "gaussian") ) {
        obj <- .infer_profiles_vb_gaussian(X, H, basis, gaussian_l, alpha_0,
                                           beta_0, vb_max_iter, epsilon_conv,
                                           is_parallel, no_cores, is_verbose)
    }else {stop(paste0(model, " observation model not implented.")) }
    # Append general class
    class(obj) <- append(class(obj), c("infer_profiles_vb", "infer_profiles"))
    # Return object
    return(obj)
}


##------------------------------------------------------

# S3 method
.infer_profiles_vb_bpr <- function(X, ...){
    UseMethod(".infer_profiles_vb_bpr")
}

# Default function for the generic function 'infer_profiles_vb_bpr'
.infer_profiles_vb_bpr.default <- function(X, ...){
    stop("Object X should be either matrix or list!")
}

# Infer profiles Binomial/Bernoulli MLE list
.infer_profiles_vb_bpr.list <- function(X, H, basis, w, alpha_0, beta_0,
                                        vb_max_iter, epsilon_conv, is_parallel,
                                        no_cores, is_verbose, ...){
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
          ex = {out_opt <- .infer_profiles_vb_bpr.matrix(
           X = X[[i]], H = H[[i]], basis = basis, w = w, alpha_0 = alpha_0,
           beta_0 = beta_0, vb_max_iter = vb_max_iter,
           epsilon_conv = epsilon_conv, is_verbose = is_verbose) })
        # Stop parallel execution
        parallel::stopCluster(cl);
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of X, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- .infer_profiles_vb_bpr.matrix(
           X = X[[i]], H = H[[i]], basis = basis, w = w, alpha_0 = alpha_0,
           beta_0 = beta_0, vb_max_iter = vb_max_iter,
           epsilon_conv = epsilon_conv, is_verbose = is_verbose) })
    }

    # Matrix for storing optimized coefficients
    W <- sapply(res, function(x) x$W)
    W_Sigma <- lapply(res, function(x) x$W_Sigma[[1]])
    if (is.matrix(W)) { W <- t(W)
    }else{W <- as.matrix(W) }
    colnames(W) <- paste("w", seq(1, NCOL(W)), sep = "")
    alpha <- as.matrix(sapply(res, function(x) x$alpha))
    beta <- as.matrix(sapply(res, function(x) x$beta))
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))
    lb_feat <- as.matrix(sapply(res, function(x) x$lb_feat))

    # Store the object
    obj <- list(W = W, W_Sigma = W_Sigma, alpha = alpha, beta = beta,
                basis = basis, nll_feat = nll_feat, rmse_feat = rmse_feat,
                coverage_feat = coverage_feat, lb_feat = lb_feat)
    if (NCOL(X[[1]]) == 3) { class(obj) <- "infer_profiles_vb_binomial"
    }else {class(obj) <- "infer_profiles_vb_bernoulli" }
    return(obj)
}

# Infer profiles Binomial/Bernoulli MLE matrix
.infer_profiles_vb_bpr.matrix <- function(X, H, basis, w, alpha_0, beta_0,
                                          vb_max_iter, epsilon_conv,
                                          is_verbose, ...){
    # TODO: Make this more memory efficient
    H_tmp <- H # Copy design matrix
    if (NCOL(X) == 3) { # If we have Binomial distributed data
        T_i <- as.vector(X[, 2]) # Total number of reads for each CpG
        m_i <- as.vector(X[, 3]) # Methylated reads for each CpG
        y <- vector(mode = "integer") # Output vector y of length (Nx1)
        for (i in 1:NROW(X)) {y <- c(y,rep(1, m_i[i]),rep(0, T_i[i] - m_i[i]))}
        # Create extended design matrix Xx of dimension (N x M)
        H <- as.matrix(H[rep(1:NROW(H), T_i), ])
    }else {y <- X[, 2] } # Output y when having Bernoulli data

    N <- NROW(H)                   # Total observations
    D <- NCOL(H)                   # Number of covariates
    LB <- c(-1e+40)                # Store the lower bounds
    E_z <- vector(mode = "numeric", length = N) # Compute expectation E[z]

    HH    <- crossprod(H, H)   # Compute H'H
    alpha <- alpha_0 + 0.5*D   # Compute 'shape' of Gamma(\tau|\alpha,\beta)
    beta  <- beta_0            # Initialize 'rate' of Gamma(\tau|\alpha,\beta)
    S <- solve(diag(D) + HH)   # Initialize covariance of q(w | m, S)
    if (is.null(w)) {m <- rep(0, D)}else {m <- w} # Initialize mean of q(w|m,S)

    # Iterate to find optimal parameters
    for (i in 2:vb_max_iter) {
        # Update mean of q(z)
        mu <- H %*% m
        # Ensure that \mu is not large enough, for numerical stability
        mu[mu > 6] <- 6; mu[mu < -6] <- -6;

        mu_1 <- mu[y == 1]  # Keep data where y == 1
        mu_0 <- mu[y == 0]  # Keep data where y == 0
        # Compute expectation E[z]
        E_z[y == 1] <- mu_1 + dnorm(-mu_1) / (1 - pnorm(-mu_1))
        E_z[y == 0] <- mu_0 - dnorm(-mu_0) / pnorm(-mu_0)
        S <- solve(diag(alpha/beta, D) + HH) # Posterior covariance of q(w)
        m <- S %*% (crossprod(H, E_z))             # Posterior mean of q(w)
        beta <- beta_0 + 0.5 * c(crossprod(m) + matrix.trace(S))
        # Check beta parameter for numerical issues
        if (beta > 10*alpha) { beta <- 10*alpha }

        # Compute lower bound
        lb_p_zw_qw <- -0.5*matrix.trace(HH %*% (tcrossprod(m,m) + S)) +
            0.5*crossprod(mu) + sum(y*log(1 - pnorm(-mu)) +
                                        (1 - y)*log(pnorm(-mu)))
        lb_pw <- -0.5*D*log(2*pi) + 0.5*D*(digamma(alpha) - log(beta)) -
            alpha/(2*beta)*(t(m) %*% (m) + matrix.trace(S))
        lb_qw <- -0.5*log(det(S)) - 0.5*D*(1 + log(2*pi))
        lb_pa <- alpha_0*log(beta_0) + (alpha_0 - 1)*(digamma(alpha) -
            log(beta)) - beta_0*(alpha/beta) - lgamma(alpha_0)
        lb_qa <- -lgamma(alpha) + (alpha - 1)*digamma(alpha) + log(beta) - alpha

        LB <- c(LB, lb_p_zw_qw + lb_pw + lb_pa - lb_qw - lb_qa)

        # Show VB difference
        if (is_verbose) { cat("It:",i,"\tLB:\t",LB[i],"\tDiff:\t",
                              LB[i] - LB[i - 1],"\n")}
        # Check if lower bound decreases
        if (LB[i] < LB[i - 1]) { warning("Lower bound decreases!\n") }
        # Check for convergence
        if (abs(LB[i] - LB[i - 1]) < epsilon_conv) { break }
        # Check if VB converged in the given maximum iterations
        if (is_verbose) {
            if (i == vb_max_iter) {warning("VB did not converge!\n")}
        }
    }

    # NLL fit feature
    nll_feat <- bpr_log_likelihood(w = m, X = X, H = H_tmp, lambda = .1,
                                   is_nll = TRUE)
    # RMSE fit feature
    f_pred <- c(pnorm(H_tmp %*% m))
    if (NCOL(X) == 3) { f_true <- X[, 3] / X[, 2] } else{f_true <- X[, 2] }
    rmse_feat <- sqrt( mean( (f_pred - f_true)^2) )
    # CpG coverage feature
    coverage_feat <- NROW(X)
    # Lower bound feature
    lb_feat <- LB[length(LB)]
    # Make weight vector as row vector
    m <- matrix(m, nrow = 1)
    # Store covariance matrix in list object
    S <- list(S)
    # Store the object
    obj <- list(W = m, W_Sigma = S, alpha = alpha, beta = beta, basis = basis,
                nll_feat = nll_feat, rmse_feat = rmse_feat,
                coverage_feat = coverage_feat, lb_feat = lb_feat)
    if (NCOL(X) == 3) { class(obj) <- "infer_profiles_vb_binomial"
    }else {class(obj) <- "infer_profiles_vb_bernoulli" }
    return(obj)
}


##------------------------------------------------------


# S3 method
.infer_profiles_vb_gaussian <- function(X, ...){
    UseMethod(".infer_profiles_vb_gaussian")
}

# Default function for the generic function 'infer_profiles_vb_gaussian'
.infer_profiles_vb_gaussian.default <- function(X, ...){
    stop("Object X should be either matrix or list!")
}


# Infer profiles Gaussian MLE list
.infer_profiles_vb_gaussian.list <- function(X, H, basis, gaussian_l,
                                             alpha_0, beta_0, vb_max_iter,
                                             epsilon_conv, is_parallel,
                                             no_cores, is_verbose, ...){
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
           ex = {out_opt <- .infer_profiles_vb_gaussian.matrix(
               X = X[[i]], H = H[[i]], basis = basis, gaussian_l = gaussian_l,
               alpha_0 = alpha_0, beta_0 = beta_0, vb_max_iter = vb_max_iter,
               epsilon_conv = epsilon_conv, is_verbose = is_verbose) })
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of X, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
           ex = {out_opt <- .infer_profiles_vb_gaussian.matrix(
               X = X[[i]], H = H[[i]], basis = basis, gaussian_l = gaussian_l,
               alpha_0 = alpha_0, beta_0 = beta_0, vb_max_iter = vb_max_iter,
               epsilon_conv = epsilon_conv, is_verbose = is_verbose) })
    }

    # Matrix for storing optimized coefficients
    W <- sapply(res, function(x) x$W)
    W_Sigma <- lapply(res, function(x) x$W_Sigma[[1]])
    if (is.matrix(W)) { W <- t(W)
    }else{W <- as.matrix(W) }
    colnames(W) <- paste("w", seq(1, NCOL(W)), sep = "")
    alpha <- as.matrix(sapply(res, function(x) x$alpha))
    beta <- as.matrix(sapply(res, function(x) x$beta))
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))
    lb_feat <- as.matrix(sapply(res, function(x) x$lb_feat))

    # Store the object
    obj <- structure(list(W = W, W_Sigma = W_Sigma, alpha = alpha,
                          beta = beta, basis = basis, gaussian_l = gaussian_l,
                          nll_feat = nll_feat, rmse_feat = rmse_feat,
                          coverage_feat = coverage_feat, lb_feat = lb_feat),
                     class = "infer_profiles_vb_gaussian")
    return(obj)
}

# Infer profiles Binomial MLE matrix
.infer_profiles_vb_gaussian.matrix <- function(X, H, basis, gaussian_l, alpha_0,
                                               beta_0, vb_max_iter,epsilon_conv,
                                               is_verbose, ...){
    D <- NCOL(H)              # Number of covariates
    I <- NROW(H)              # Number of observations
    HH <- crossprod(H, H)     # Compute H'H
    Hy <- crossprod(H, X[,2]) # Compute H'y
    yy <- c(crossprod(X[,2])) # Compute y'y
    LB <- c(-Inf)             # Store the lower bounds

    alpha <- alpha_0 + D / 2    # Compute \alpha parameter of Gamma
    E_a <- alpha_0 / beta_0     # Initialize expectation of precision parameter

    # Iterate to find optimal parameters
    for (i in 2:vb_max_iter) {
        S <- solve(diag(E_a, D) + gaussian_l * HH) # Covariance Gaussian factor
        m <- gaussian_l * S %*% Hy                 # Mean of Gaussian factor
        E_ww <- c(crossprod(m) + matrix.trace(S))  # Expectation of E[w'w]
        beta <- beta_0 + 0.5 * E_ww  # Update \beta parameter of Gamma factor
        # Check beta parameter for numerical issues
        if (beta > 3*alpha) { beta <- 3*alpha }
        E_a <- alpha / beta          # Compute expectation of E[a]

        # Compute lower bound
        lb_py <- 0.5 * I * log(gaussian_l / (2 * pi)) - 0.5 * gaussian_l * yy +
            gaussian_l * crossprod(m, Hy) - 0.5 * gaussian_l *
            matrix.trace(HH %*% (tcrossprod(m, m) + S))
        lb_pw <- -0.5 * D * log(2*pi) + 0.5 * D * (digamma(alpha) - log(beta)) -
            0.5 * E_a * E_ww
        lb_pa <- alpha_0*log(beta_0) + (alpha_0 - 1) * (digamma(alpha) -
            log(beta)) - beta_0*E_a - log(gamma(alpha_0))
        lb_qw <- -0.5 * log(det(S)) - 0.5 * D * (1 + log(2*pi))
        lb_qa <- -lgamma(alpha) + (alpha - 1)*digamma(alpha) + log(beta) - alpha

        LB <- c(LB, lb_py + lb_pw + lb_pa - lb_qw - lb_qa)

        # Show VB difference
        if (is_verbose) { cat("It:",i,"\tLB:\t",LB[i],"\tDiff:\t",
                              LB[i] - LB[i - 1],"\n")}
        # Check if lower bound decreases
        if (LB[i] < LB[i - 1]) { warning("Lower bound decreases!\n")}
        # Check for convergence
        if (abs(LB[i] - LB[i - 1]) < epsilon_conv) { break }
        # Check if VB converged in the given maximum iterations
        if (is_verbose) {
            if (i == vb_max_iter) {warning("VB did not converge!\n")}
        }
    }

    # NLL fit feature
    nll_feat <- lr_log_likelihood(w = m, X = X, H = H, lambda = .1,
                                  is_nll = TRUE)
    # RMSE fit feature
    f_pred <- c(H %*% m)
    f_true <- X[, 2]
    rmse_feat <- sqrt( mean( (f_pred - f_true)^2) )
    # CpG coverage feature
    coverage_feat <- NROW(X)
    # Lower bound feature
    lb_feat <- LB[length(LB)]
    # Make weight vector as row vector
    m <- matrix(m, nrow = 1)
    # Store covariance matrix in list object
    S <- list(S)
    # Store the object
    obj <- structure(list(W = m, W_Sigma = S, alpha = alpha, beta = beta,
                          basis = basis, gaussian_l = gaussian_l,
                          nll_feat = nll_feat, rmse_feat = rmse_feat,
                          coverage_feat = coverage_feat, lb_feat = lb_feat),
                     class = "infer_profiles_vb_gaussian")
    return(obj)
}
