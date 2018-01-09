#' @name cluster_profiles_mle
#' @rdname cluster_profiles_mle
#' @aliases cluster_profile_mle cluster_mle
#'
#' @title Cluster methylation profiles using EM
#'
#' @description General purpose functions for clustering latent profiles for
#'   different observation models using maximum likelihood estimation (MLE) and
#'   the EM algorithm. Initially, it performs parameter checking, and
#'   initializes main parameters, such as mixing proportions, basis function
#'   coefficients, then the EM algorithm is applied and finally model selection
#'   metrics are calculated, such as BIC and AIC.
#'
#' @param X The input data, which has to be a \code{\link[base]{list}} of
#'   elements of length N, where each element is an \code{L X C} matrix, where L
#'   are the total number of observations. The first column contains the input
#'   observations x (i.e. CpG locations). If "binomial" model then C=3, and 2nd
#'   and 3rd columns contain total number of trials and number of successes
#'   respectively. If "bernoulli" or "gaussian" model, then C=2 containing the
#'   output y (e.g. methylation level). If "beta" model, then C=3, where 2nd
#'   column contains output y and 3rd column the dispersion parameter.
#' @param K Integer denoting the total number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A MxK matrix of the initial parameters, where each column consists
#'   of the basis function coefficients for each corresponding cluster. If NULL,
#'   will be assigned with default values.
#' @param gaussian_sigma Initial standard deviation of the noise term, only used
#'   when having "gaussian" observation model.
#' @param em_max_iter Integer denoting the maximum number of EM iterations.
#' @param epsilon_conv Numeric denoting the convergence threshold for EM.
#' @param init_opt_itnmax Optimization iterations for obtaining the initial EM
#'   parameter values.
#' @param is_verbose Logical, print results during EM iterations.
#' @inheritParams infer_profiles_mle
#'
#' @return An object of class \code{cluster_profiles_mle_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An (M+1) X K matrix with the
#'   optimized parameter values for each cluster. Each column of the matrix
#'   corresponds a different cluster k.} \item{ \code{pi_k}: Mixing proportions.
#'   } \item{ \code{r_nk}: An (N X K) responsibility matrix of each
#'   observations being explained by a specific cluster. }  \item{ \code{basis}:
#'   The basis object. } \item{\code{nll}: The negative log likelihood vector.}
#'   \item{\code{labels}: Cluster assignment labels.} \item{\code{bic}: Bayesian
#'   Information Criterion metric.} \item{\code{aic}: Akaike Information
#'   Criterion metric.} \item{\code{icl}: Integrated Complete Likelihood
#'   criterion metric.} \item{\code{gaussian_sigma}: Optimized standard
#'   deviation for gaussian observation model.} }
#'
#' @section Details: The beta regression model is based on alternative
#'   parameterization of the beta density in terms of the mean and dispersion
#'   parameter: \url{https://cran.r-project.org/web/packages/betareg/}. For
#'   modelling details for Binomial/Bernoulli observation model check the paper
#'   for BPRMeth:
#'   \url{https://academic.oup.com/bioinformatics/article/32/17/i405/2450762} .
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{cluster_profiles_vb}}
#'   \code{\link{infer_profiles_vb}}, \code{\link{infer_profiles_mle}},
#'   \code{\link{infer_profiles_gibbs}}, \code{\link{create_met_region}}
NULL


#' @rdname cluster_profiles_mle
#'
#' @examples
#' # Example of optimizing parameters for synthetic data using 3 RBFs
#'
#' basis <- create_rbf_object(M=3)
#' out <- cluster_profiles_mle(X = binomial_data, model = "binomial",
#'   basis=basis, em_max_iter = 5, opt_itnmax = 5, init_opt_itnmax=5,
#'   is_parallel = FALSE)
#'
#' #-------------------------------------
#'
#' basis <- create_rbf_object(M=3)
#' out <- cluster_profiles_mle(X = gaussian_data, model = "gaussian",
#'   basis=basis, em_max_iter = 5, opt_itnmax = 5, init_opt_itnmax=5,
#'   is_parallel = FALSE)
#'
#' @export
cluster_profiles_mle <- function(X, K = 3, model = NULL, basis = NULL, H = NULL,
                                 pi_k = NULL, lambda = .5, beta_dispersion = 5,
                                 gaussian_sigma = rep(.2, K), w = NULL,
                                 em_max_iter = 50, epsilon_conv = 1e-4,
                                 opt_method = "CG", opt_itnmax = 50,
                                 init_opt_itnmax = 30, is_parallel = FALSE,
                                 no_cores = NULL, is_verbose = FALSE, ...){
    assertthat::assert_that(is.list(X))  # Check that X is a list object
    if (is.null(model)) { stop("Observation model not defined!") }
    tmp <- .em_checks(X, H, K, pi_k, model, basis, lambda, beta_dispersion, w,
                      opt_method, init_opt_itnmax, is_parallel, no_cores)
    w <- tmp$w; basis <- tmp$basis; pi_k <- tmp$pi_k

    # TODO: Make this faster, EM_checks also computes the design matrix.
    if (is.null(H)) { # Create design matrix
        H <- lapply(X,function(x) design_matrix(basis,x[,1])$H)
    }

    if (identical(model, "bernoulli") || identical(model, "binomial") ) {
        obj <- .cluster_profiles_mle_bpr(X = X, H = H, K = K, pi_k = pi_k,
            basis = basis, lambda = lambda, w = w, em_max_iter = em_max_iter,
            epsilon_conv = epsilon_conv, opt_method = opt_method,
            opt_itnmax = opt_itnmax, is_parallel = is_parallel,
            no_cores = no_cores, is_verbose = is_verbose)
    }else if (identical(model, "beta") ) {
        obj <- .cluster_profiles_mle_beta(X = X, H = H, K = K, pi_k = pi_k,
            basis = basis, lambda = lambda, beta_dispersion = beta_dispersion,
            w = w, em_max_iter = em_max_iter, epsilon_conv = epsilon_conv,
            opt_method = opt_method, opt_itnmax = opt_itnmax,
            is_parallel = is_parallel, no_cores = no_cores,
            is_verbose = is_verbose)
    }else if (identical(model, "gaussian")) {
        obj <- .cluster_profiles_mle_gaussian(X = X, H = H, K = K, pi_k = pi_k,
            basis = basis, lambda = lambda, gaussian_sigma = gaussian_sigma,
            w = w, em_max_iter = em_max_iter, epsilon_conv = epsilon_conv,
            is_verbose = is_verbose)
    }else {stop(paste0(model, " observation model not implented.")) }

    # Add names to the estimated parameters for clarity
    names(obj$pi_k) <- paste0("cluster", 1:K)
    colnames(obj$W) <- paste0("cluster", 1:K)
    colnames(obj$r_nk) <- paste0("cluster", 1:K)
    # Get hard cluster assignments for each observation
    obj$labels <- apply(X = obj$r_nk, MARGIN = 1,
                        FUN = function(x) which(x == max(x, na.rm = TRUE)))
    # Total params
    total_params <- (K - 1) + K * NROW(w)
    # Bayesian Information Criterion
    obj$bic <- 2 * utils::tail(obj$nll, n = 1) + total_params * log(length(X))
    # Akaike Iformation Criterion
    obj$aic <- 2 * utils::tail(obj$nll, n = 1) + 2 * total_params
    # Integrated Complete Likelihood criterion
    entropy <- (-1) * sum(obj$r_nk * log(obj$r_nk), na.rm = TRUE)
    obj$icl <- obj$bic + entropy
    # Append general class
    class(obj) <- append(class(obj),
                         c("cluster_profiles_mle", "cluster_profiles"))
    # Return object
    return(obj)
}


##------------------------------------------------------

# Cluster profiles EM Binomial or Bernoulli
.cluster_profiles_mle_bpr <- function(X, H, K, pi_k, basis, lambda, w,
                                      em_max_iter, epsilon_conv, opt_method,
                                      opt_itnmax, is_parallel, no_cores,
                                      is_verbose){
    assertthat::assert_that(is.list(X)) # Check that X is a list object
    k <- 0                  # Initialize so the CMD check passes without NOTES
    N <- length(X)          # Extract number of observations
    assertthat::assert_that(N > 0)
    weighted_pdf <- matrix(0, nrow = N, ncol = K) # Store weighted PDFs
    nll <- c(1e+40) # Initialize and store NLL for each EM iteration

    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel,
                                max_cores = K)
    # If parallel mode is ON create cluster object
    if (is_parallel) {
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
    }

    # Run EM algorithm until convergence
    for (t in 1:em_max_iter) {
        #
        # E-Step -----------------------------------------------
        for (k in 1:K) {
            # For each element in X, evaluate the log likelihood
            weighted_pdf[, k] <- vapply(X = 1:N, FUN = function(y)
                bpr_log_likelihood(w = w[, k], X = X[[y]], H = H[[y]],
                                   lambda = lambda, is_nll = FALSE),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
            weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
        }
        # Calculate probs using the logSumExp trick for numerical stability
        Z <- apply(weighted_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        r_nk <- exp(weighted_pdf - Z)
        # Evaluate and store the NLL
        nll  <- c(nll, (-1) * sum(Z))
        #
        # M-Step -----------------------------------------------
        # Compute sum of posterior probabilities for each cluster
        N_k <- colSums(r_nk)
        # Update mixing proportions for each cluster
        pi_k <- N_k / N
        # Update basis function coefficient vector w for each cluster
        if (is_parallel) { # If parallel mode is ON
            # Parallel optimization for each cluster k
            w <- foreach::"%dopar%"(obj = foreach::foreach(k = 1:K,
                                                           .combine = cbind),
             ex = {out <- optim(par = w[, k], fn = sum_weighted_bpr_lik,
              gr = sum_weighted_bpr_grad, method = opt_method,
              control = list(maxit = opt_itnmax), X_list = X, H_list = H,
              r_nk = r_nk[, k], lambda = lambda, is_nll = TRUE)$par})
        }else{
            # Sequential optimization for each clustrer k
            w <- foreach::"%do%"(obj = foreach::foreach(k = 1:K,
                                                        .combine = cbind),
             ex = {out <- optim(par = w[, k], fn = sum_weighted_bpr_lik,
              gr = sum_weighted_bpr_grad, method = opt_method,
              control = list(maxit = opt_itnmax), X_list = X, H_list = H,
              r_nk = r_nk[, k], lambda = lambda, is_nll = TRUE)$par})
        }

        if (is_verbose) {
            cat("It:",t,"\tNLL:\t",nll[t + 1],"\tDiff:\t",
                nll[t] - nll[t + 1],"\n")
        }
        if (nll[t + 1] > nll[t]) {message("NLL increases!\n") }
        # Check for convergence
        if (nll[t] - nll[t + 1] < epsilon_conv) { break }
        if (K == 1) { w <- as.matrix(w) }
    }
    if (is_parallel) {
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }
    if (K == 1) { w <- as.matrix(w) }
    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter) { warning("EM did not converge!\n") }

    # Store the object
    obj <- list(W = w, pi_k = pi_k, r_nk = r_nk, basis = basis, nll = nll)
    if (NCOL(X[[1]]) == 3) { class(obj) <- "cluster_profiles_mle_binomial"
    }else {class(obj) <- "cluster_profiles_mle_bernoulli" }
    return(obj)
}


##------------------------------------------------------


# Cluster profiles EM Beta
.cluster_profiles_mle_beta <- function(X, H, K, pi_k, basis, lambda,
                                       beta_dispersion, w, em_max_iter,
                                       epsilon_conv, opt_method, opt_itnmax,
                                       is_parallel, no_cores, is_verbose){
    assertthat::assert_that(is.list(X)) # Check that X is a list object
    k <- 0                  # Initialize so the CMD check passes without NOTES
    N <- length(X)          # Extract number of observations
    assertthat::assert_that(N > 0)
    weighted_pdf <- matrix(0, nrow = N, ncol = K) # Store weighted PDFs
    nll <- c(1e+40) # Initialize and store NLL for each EM iteration

    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel,
                                max_cores = K)
    # If parallel mode is ON create cluster object
    if (is_parallel) {
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl) }
    # Add dispersion parameter in observation matrix X
    for (i in 1:N) {
        X[[i]] <- cbind(X[[i]], rep(beta_dispersion, NROW(X[[i]])))
    }

    # Run EM algorithm until convergence
    for (t in 1:em_max_iter) {
        #
        # E-Step -----------------------------------------------
        for (k in 1:K) {
            # For each element in X, evaluate the log likelihood
            weighted_pdf[, k] <- vapply(X = 1:N, FUN = function(y)
                betareg_log_likelihood(w = w[, k], X = X[[y]], H = H[[y]],
                                       lambda = lambda, is_nll = FALSE),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
            weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
        }
        # Calculate probs using the logSumExp trick for numerical stability
        Z <- apply(weighted_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        r_nk <- exp(weighted_pdf - Z)
        # Evaluate and store the NLL
        nll  <- c(nll, (-1) * sum(Z))
        #
        # M-Step -----------------------------------------------
        # Compute sum of posterior probabilities for each cluster
        N_k <- colSums(r_nk)
        # Update mixing proportions for each cluster
        pi_k <- N_k / N
        # Update basis function coefficient vector w for each cluster
        if (is_parallel) { # If parallel mode is ON
            # Parallel optimization for each cluster k
            w <- foreach::"%dopar%"(obj = foreach::foreach(k = 1:K,
                                                           .combine = cbind),
             ex = {out <- optim(par = w[, k], fn = sum_weighted_betareg_lik,
              gr = sum_weighted_betareg_grad, method = opt_method,
              control = list(maxit = opt_itnmax), X_list = X, H_list = H,
              r_nk = r_nk[, k], lambda = lambda, is_nll = TRUE)$par})
        }else{
            # Sequential optimization for each clustrer k
            w <- foreach::"%do%"(obj = foreach::foreach(k = 1:K,
                                                        .combine = cbind),
             ex = {out <- optim(par = w[, k], fn = sum_weighted_betareg_lik,
              gr = sum_weighted_betareg_grad, method = opt_method,
              control = list(maxit = opt_itnmax), X_list = X, H_list = H,
              r_nk = r_nk[, k], lambda = lambda, is_nll = TRUE)$par})
        }

        if (is_verbose) {
            cat("It:",t,"\tNLL:\t",nll[t + 1],"\tDiff:\t",
                nll[t] - nll[t + 1],"\n")
        }
        if (nll[t + 1] > nll[t]) {message("NLL increases!\n")}
        # Check for convergence
        if (nll[t] - nll[t + 1] < epsilon_conv) { break }
        if (K == 1) { w <- as.matrix(w) }
    }
    if (is_parallel) {
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }
    if (K == 1) { w <- as.matrix(w) }
    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter) { warning("EM did not converge!\n") }

    # Store the object
    obj <- structure(list(W = w, pi_k = pi_k, r_nk = r_nk, basis = basis,
                          nll = nll), class = "cluster_profiles_mle_beta")
    return(obj)
}


##------------------------------------------------------


# Cluster profiles EM Binomial or Bernoulli
.cluster_profiles_mle_gaussian <- function(X, H, K, pi_k, basis, lambda,
                                           gaussian_sigma, w, em_max_iter,
                                           epsilon_conv, is_verbose){
    assertthat::assert_that(is.list(X)) # Check that X is a list object
    k <- 0                  # Initialize so the CMD check passes without NOTES
    N <- length(X)          # Extract number of observations
    assertthat::assert_that(N > 0)
    weighted_pdf <- matrix(0, nrow = N, ncol = K) # Store weighted PDFs
    nll <- c(1e+40) # Initialize and store NLL for each EM iteration

    # Precompute total observations for efficiency
    L_n <- vector("numeric", length = N)
    for (n in 1:N) { L_n[n] <- NROW(X[[n]]) }

    # Run EM algorithm until convergence
    for (t in 1:em_max_iter) {
        #
        # E-Step -----------------------------------------------
        for (k in 1:K) {
            # For each element in x, evaluate the BPR log likelihood
            weighted_pdf[, k] <- vapply(X = 1:N, FUN = function(y)
                sum(dnorm(x = X[[y]][,2], mean = H[[y]] %*% w[, k],
                          sd = gaussian_sigma[k], log = TRUE)) -
                    lambda * t(w[,k]) %*% w[,k],
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
            weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
        }
        # Calculate probs using the logSumExp trick for numerical stability
        Z <- apply(weighted_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        r_nk <- exp(weighted_pdf - Z)
        # Evaluate and store the NLL
        nll  <- c(nll, (-1) * sum(Z))
        #
        # M-Step -----------------------------------------------
        # Compute sum of posterior probabilities for each cluster
        N_k <- colSums(r_nk)
        # Update mixing proportions for each cluster
        pi_k <- N_k / N
        # Iterate over each cluster
        for (k in 1:K) {
            # Update basis function coefficient vector w for each cluster
            tmp_HH <- matrix(0, ncol = basis$M + 1, nrow = basis$M + 1)
            tmp_Hy <- vector("numeric", length = basis$M + 1)
            for (n in 1:N) {
                tmp_HH <- tmp_HH + crossprod(H[[n]]) * r_nk[n, k]
                tmp_Hy <- tmp_Hy + crossprod(H[[n]], X[[n]][,2])*r_nk[n,k]
            }
            w[, k] <- solve(tmp_HH + lambda) %*% tmp_Hy
            # Update variance of regression model for each cluster
            tmp <- sum(vapply(X = 1:N, FUN = function(y)
                crossprod(X[[y]][,2] - H[[y]] %*% w[, k]) * r_nk[y, k],
                FUN.VALUE = numeric(1),  USE.NAMES = FALSE))
            gaussian_sigma[k] <- tmp / r_nk[, k] %*% L_n
        }

        if (is_verbose) {
            cat("It:",t,"\tNLL:\t",nll[t + 1],"\tDiff:\t",
                nll[t] - nll[t + 1],"\n")
        }
        if (nll[t + 1] > nll[t]) {message("NLL increases!\n")}
        # Check for convergence
        if (nll[t] - nll[t + 1] < epsilon_conv) { break }
        if (K == 1) { w <- as.matrix(w) }
    }
    if (K == 1) { w <- as.matrix(w) }
    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter) { warning("EM did not converge!\n") }

    # Store the object
    obj <- structure(list(W = w, pi_k = pi_k, r_nk = r_nk, basis = basis,
                          nll = nll, gaussian_sigma = gaussian_sigma),
                     class = "cluster_profiles_mle_gaussian")
    return(obj)
}
