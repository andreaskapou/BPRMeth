#' @name infer_profiles_gibbs
#' @rdname infer_profiles_gibbs
#' @aliases infer_profile_gibbs inference_gibbs
#'
#' @title Infer methylation profiles using Gibbs sampling
#'
#' @description General purpose functions for inferring latent profiles for
#'   different observation models using Gibbs sampling. Currently implemented
#'   observation models are: 'bernoulli' and 'binomial' and the auxiliary
#'   variable approach is used.
#'
#' @param X The input data, either a \code{\link[base]{matrix}} or a
#'   \code{\link[base]{list}} of elements of length N, where each element is an
#'   \code{L X C} matrix, where L are the total number of observations. The
#'   first column contains the input observations x (i.e. CpG locations). If
#'   "binomial" model then C=3, and 2nd and 3rd columns contain total number of
#'   trials and number of successes respectively. If "bernoulli" then C=2
#'   containing the output y (e.g. methylation level).
#' @param H Optional, design matrix of the input data X. If NULL, H will be
#'   computed inside the function.
#' @param model Observation model name as character string. It can be either
#'   'bernoulli' or 'binomial'.
#' @param basis A 'basis' object. E.g. see \code{\link{create_basis}}. If NULL,
#'   will an RBF object will be created.
#' @param w A vector of initial parameters (i.e. coefficients of the basis
#'   functions). If NULL, it will be initialized inside the function.
#' @param mu_0 The prior mean hyperparameter vector for w.
#' @param cov_0 The prior covariance hyperparameter matrix for w.
#' @param gibbs_nsim Total number of simulations for the Gibbs sampler.
#' @param gibbs_burn_in Burn in period of the Gibbs sampler.
#' @param store_gibbs_draws Logical indicating if we should keep the whole MCMC
#'   chain for further analysis.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param ... Additional parameters.
#'
#' @return An object of class \code{infer_profiles_gibbs_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   posterior mean of the parameters w. Each row of the matrix corresponds to
#'   each element of the list X; if X is a matrix, then N = 1. The columns are
#'   of the same length as the parameter vector w (i.e. number of basis
#'   functions). } \item{ \code{W_sd}: An Nx(M+1) matrix with the posterior
#'   standard deviation (sd) of the parameters W.} \item{ \code{basis}: The
#'   basis object. } \item{\code{nll_feat}: NLL fit feature.}
#'   \item{\code{rmse_feat}: RMSE fit feature.} \item{\code{coverage_feat}: CpG
#'   coverage feature.} \item{\code{W_draws}: Optional, draws of the Gibbs
#'   sampler. } }
#'
#' @section Details: The modelling and mathematical details for inferring
#'   profiles using Gibbs sampling are explained here:
#'   \url{http://rpubs.com/cakapourani/} . More specifically: \itemize{\item{For
#'   Binomial observation model check:
#'   \url{http://rpubs.com/cakapourani/bayesian-bpr-model}} \item{For Bernoulli
#'   observation model check:
#'   \url{http://rpubs.com/cakapourani/bayesian-bpr-model}}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{infer_profiles_mle}},
#'   \code{\link{infer_profiles_vb}}, \code{\link{create_region_object}}
NULL


#' @rdname infer_profiles_gibbs
#'
#' @examples
#' # Example of inferring parameters for synthetic data using 3 RBFs
#' basis <- create_rbf_object(M=3)
#' out <- infer_profiles_gibbs(X = binomial_data, model = "binomial",
#'    basis = basis, is_parallel = FALSE, gibbs_nsim = 10, gibbs_burn_in = 5)
#'
#' #-------------------------------------
#'
#' @export
infer_profiles_gibbs <- function(X, model = NULL, basis = NULL, H = NULL,
                                 w = NULL, mu_0 = NULL, cov_0 = NULL,
                                 gibbs_nsim = 500, gibbs_burn_in = 100,
                                 store_gibbs_draws = FALSE, is_parallel = FALSE,
                                 no_cores = NULL, ...){
    if (is.null(model)) { stop("Observation model not defined!") }
    # Create RBF basis object by default
    if (is.null(basis)) {
        warning("Basis object not defined. Using as default M = 3 RBFs.\n")
        basis <- create_rbf_object(M = 3)
    }
    if (is.null(H)) { # Create design matrix
        if (is.list(X)) {
            # Remove rownames
            X <- lapply(X, function(x) { rownames(x) <- NULL; return(x) })
            H <- lapply(X, function(x) design_matrix(basis,x[,1])$H)
        }else {
            rownames(X) <- NULL
            H <- design_matrix(basis, X[,1])$H
        }
    }
    if (is.null(mu_0)) { mu_0 <- rep(0, basis$M + 1)}
    if (is.null(cov_0)) { cov_0 <- diag(4, basis$M + 1)}

    if (identical(model, "bernoulli") || identical(model, "binomial") ) {
        obj <- .infer_profiles_gibbs_bpr(X, H, basis, w, mu_0, cov_0,gibbs_nsim,
                                         gibbs_burn_in, store_gibbs_draws,
                                         is_parallel, no_cores)
    }else {stop(paste0(model, " observation model not implented.")) }
    # Append general class
    class(obj) <- append(class(obj),
                         c("infer_profiles_gibbs", "infer_profiles"))
    # Return object
    return(obj)
}


##------------------------------------------------------

# S3 method
.infer_profiles_gibbs_bpr <- function(X, ...){
    UseMethod(".infer_profiles_gibbs_bpr")
}

# Default function for the generic function 'infer_profiles_gibbs_bpr'
.infer_profiles_gibbs_bpr.default <- function(X, ...){
    stop("Object X should be either matrix or list!")
}

# Infer profiles Binomial/Bernoulli Gibbs list
.infer_profiles_gibbs_bpr.list <- function(X, H, basis, w, mu_0, cov_0,
                                           gibbs_nsim, gibbs_burn_in,
                                           store_gibbs_draws, is_parallel,
                                           no_cores, ...){
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
          ex = {out_opt <- .infer_profiles_gibbs_bpr.matrix(
              X = X[[i]], H = H[[i]], basis = basis, w = w, mu_0 = mu_0,
              cov_0 = cov_0, gibbs_nsim = gibbs_nsim,
              gibbs_burn_in = gibbs_burn_in,
              store_gibbs_draws = store_gibbs_draws)})
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of X, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
           ex = {out_opt <- .infer_profiles_gibbs_bpr.matrix(
               X = X[[i]], H = H[[i]], basis = basis, w = w, mu_0 = mu_0,
               cov_0 = cov_0, gibbs_nsim = gibbs_nsim,
               gibbs_burn_in = gibbs_burn_in,
               store_gibbs_draws = store_gibbs_draws) })
    }
    # Matrix for storing optimized coefficients
    W <- sapply(res, function(x) x$W)
    W_sd <- sapply(res, function(x) x$W_sd)
    W_draws <- lapply(res, function(x) x$W_draws)
    if (is.matrix(W)) { W <- t(W); W_sd <- t(W_sd)
    }else {W <- as.matrix(W); W_sd <- as.matrix(W_sd) }
    colnames(W) <- paste("w", seq(1, NCOL(W)), sep = "")
    colnames(W_sd) <- paste("w", seq(1, NCOL(W_sd)), sep = "")
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))

    # Store the object
    obj <- list(W = W, W_sd = W_sd, basis = basis, nll_feat = nll_feat,
                rmse_feat = rmse_feat, coverage_feat = coverage_feat,
                W_draws = W_draws)
    if (NCOL(X[[1]]) == 3) { class(obj) <- "infer_profiles_gibbs_binomial"
    }else {class(obj) <- "infer_profiles_gibbs_bernoulli" }
    return(obj)
}

# Infer profiles Binomial/Bernoulli Gibbs matrix
.infer_profiles_gibbs_bpr.matrix <- function(X, H, basis, w, mu_0, cov_0,
                                             gibbs_nsim, gibbs_burn_in,
                                             store_gibbs_draws, ...){
    if (is.null(w)) { w <- rep(0.5, NCOL(H)) }
    H_gibbs <- H
    if (NCOL(X) == 3) { # If we have Binomial distributed data
        T_i <- as.vector(X[, 2]) # Total number of reads for each CpG
        m_i <- as.vector(X[, 3]) # Methylated reads for each CpG
        N <- sum(T_i)            # Sum of total trials for each observation i
        y <- vector(mode = "integer") # Create extended vector y of length (Nx1)
        for (i in 1:NROW(X)) {y <- c(y,rep(1, m_i[i]),rep(0, T_i[i] - m_i[i]))}
        # Create extended design matrix Xx of dimension (N x M)
        H_gibbs <- as.matrix(H[rep(1:NROW(H), T_i), ])
    }else {# If we have Bernoulli distributed data
        N <- length(X[ ,1])  # Sum of total trials for each observation i
        y <- X[, 2]          # Output y
    }
    # Gibbs sampling algorithm
    w_chain <- .gibbs_bpr(H = H_gibbs, y = y, N = N, w = w, mu_0 = mu_0,
                          cov_0 = cov_0, gibbs_nsim = gibbs_nsim)
    # If we learn a constant function
    if (basis$M == 0) {
        W <- mean(w_chain[-(1:gibbs_burn_in)])
        W_sd <- stats::sd(w_chain[-(1:gibbs_burn_in)])
        if (store_gibbs_draws) { W_draws <- w_chain[-(1:gibbs_burn_in)]
        }else {W_draws <- NULL }

    }else{
        W <- colMeans(w_chain[-(1:gibbs_burn_in), ])
        W_sd <- apply(w_chain[-(1:gibbs_burn_in), ], 2, stats::sd)
        if (store_gibbs_draws) { W_draws <- w_chain[-(1:gibbs_burn_in), ]
        }else {W_draws <- NULL }
    }
    # NLL fit feature
    nll_feat <- bpr_log_likelihood(w = W, X = X, H = H, lambda = .1,
                                   is_nll = TRUE)
    # RMSE fit feature
    f_pred <- c(pnorm(H %*% W))
    if (NCOL(X) == 3) { f_true <- X[, 3] / X[, 2] }else {f_true <- X[, 2] }
    rmse_feat <- sqrt( mean( (f_pred - f_true)^2) )
    # CpG coverage feature
    coverage_feat <- NROW(X)
    # Make weight vector as row vector
    W <- matrix(W, nrow = 1)
    W_sd <- matrix(W_sd, nrow = 1)
    # Store the object
    obj <- list(W = W, W_sd = W_sd, basis = basis, nll_feat = nll_feat,
                rmse_feat = rmse_feat, coverage_feat = coverage_feat,
                W_draws = W_draws)
    if (NCOL(X) == 3) { class(obj) <- "infer_profiles_gibbs_binomial"
    }else {class(obj) <- "infer_profiles_gibbs_bernoulli" }
    return(obj)
}


# Perform Gibbs sampling for BPR model
#
# \code{.gibbs_bpr} performs Gibbs sampling for BPR model using auxiliary
#  variable approach.
#
# @param H Design matrix
# @param y Output y
# @param N Total number of trials
# @param w Maximum likelihood estimate for initial value
# @param mu_0 Prior mean vector for parameter w
# @param cov_0 Prior covariance matrix for parameter w
# @param gibbs_nsim Number of Gibbs simulations
#
# @return The chain with the sampled w values from the posterior
#
.gibbs_bpr <- function(H, y, N, w, mu_0, cov_0, gibbs_nsim){
    N1  <- sum(y)  # Number of successes
    N0  <- N - N1  # Number of failures
    # Matrix storing samples of the \w parameter
    w_chain <- matrix(0, nrow = gibbs_nsim, ncol = length(w))
    w_chain[1, ] <- w
    prec_0 <- solve(cov_0) # Compute posterior variance of w
    V <- solve(prec_0 + crossprod(H, H))
    z <- rep(0, N)  # Initialize latent variable Z, from truncated normal

    # Check all the cases when you might have totally methylated or unmethylated
    # read, then we cannot create 0 samples.
    if (N0 == 0) {
        for (t in 2:gibbs_nsim) {
            mu_z <- H %*% w # Update Mean of z
            # Draw latent variable z from its full conditional: z|w,y,X
            z[y == 1] <- truncnorm::rtruncnorm(N1, mean = mu_z[y == 1],
                                               sd = 1, a = 0, b = Inf)
            Mu <- V %*% (prec_0 %*% mu_0 + crossprod(H, z)) # Posterior mean w
            w <- c(mvtnorm::rmvnorm(1, Mu, V)) # Draw w from conditional: w|z,X
            w_chain[t, ] <- w         # Store the w draws
        }
    }else if (N1 == 0) {
        for (t in 2:gibbs_nsim) {
            mu_z <- H %*% w # Update Mean of z
            # Draw latent variable z from its full conditional: z|w,y,X
            z[y == 0] <- truncnorm::rtruncnorm(N0, mean = mu_z[y == 0],
                                               sd = 1, a = -Inf, b = 0)
            Mu <- V %*% (prec_0 %*% mu_0 + crossprod(H, z)) # Posterior mean w
            w <- c(mvtnorm::rmvnorm(1, Mu, V)) # Draw w from conditional: w|z,X
            w_chain[t, ] <- w         # Store the w draws
        }
    }else{
        for (t in 2:gibbs_nsim) {
            mu_z <- H %*% w # Update Mean of z
            # Draw latent variable z from its full conditional: z|w,y,X
            z[y == 0] <- truncnorm::rtruncnorm(N0, mean = mu_z[y == 0],
                                               sd = 1, a = -Inf, b = 0)
            z[y == 1] <- truncnorm::rtruncnorm(N1, mean = mu_z[y == 1],
                                               sd = 1, a = 0, b = Inf)
            Mu <- V %*% (prec_0 %*% mu_0 + crossprod(H, z)) # Posterior mean w
            w <- c(mvtnorm::rmvnorm(1, Mu, V)) # Draw w from conditional: w|z,X
            w_chain[t, ] <- w         # Store the w draws
        }
    }
    return(w_chain)
}
