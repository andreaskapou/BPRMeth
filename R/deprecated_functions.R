#' @title (DEPRECATED) Predict gene expression from methylation profiles
#'
#' @description (DEPRECATED) \code{bpr_predict_wrap} is a function that wraps
#'   all the necessary subroutines for performing prediction on gene expression
#'   levels. Initially, it optimizes the parameters of the basis functions so as
#'   to learn the methylation profiles. Then, uses the learned parameters /
#'   coefficients of the basis functions as input features for performing
#'   regression in order to predict the corresponding gene expression levels.
#'
#' @param formula An object of class \code{\link[stats]{formula}}, e.g. see
#'   \code{\link[stats]{lm}} function. If NULL, the simple linear regression
#'   model is used.
#' @param x The binomial distributed observations, which has to be a list of
#'   elements of length N, where each element is an L x 3 matrix of
#'   observations, where 1st column contains the locations. The 2nd and 3rd
#'   columns contain the total trials and number of successes at the
#'   corresponding locations, repsectively. See
#'   \code{\link{process_haib_caltech_wrap}} on a possible way to get this data
#'   structure.
#' @param y Corresponding gene expression data for each element of the list x.
#' @param model_name A string denoting the regression model. Currently,
#'   available models are: \code{"svm"}, \code{"randomForest"}, \code{"rlm"},
#'   \code{"mars"} and \code{"lm"}.
#' @param w Optional vector of initial parameter / coefficient values.
#' @param basis Optional basis function object, default is an 'rbf' object, see
#'   \code{\link{create_rbf_object}}.
#' @param fit_feature Return additional feature on how well the profile fits the
#'   methylation data. Either NULL for ignoring this feature or one of the
#'   following: 1) "RMSE" for returning the fit of the profile using the RMSE as
#'   measure of error or 2) "NLL" for returning the fit of the profile using the
#'   Negative Log Likelihood as measure of error.
#' @param cov_feature Logical, whether to return an additional feature for the
#'   CpG coverage across the promoter region.
#' @param train_ind Optional vector containing the indices for the train set.
#' @param train_perc Optional parameter for defining the percentage of the
#'   dataset to be used for training set, the remaining will be the test set.
#' @param is_summary Logical, print the summary statistics.
#' @inheritParams bpr_optimize
#'
#' @return A 'bpr_predict' object which, in addition to the input parameters,
#'   consists of the following variables: \itemize{ \item{ \code{W_opt}: An
#'   Nx(M+1) matrix with the optimized parameter values. Each row of the matrix
#'   corresponds to each element of the list x. The columns are of the same
#'   length as the parameter vector w (i.e. number of basis functions). } \item{
#'   \code{Mus}: An N x M matrix with the RBF centers if basis object is
#'   \code{\link{create_rbf_object}}, otherwise NULL.} \item{train}: The
#'   training data. \item{test}: The test data. \item \code{gex_model}: The
#'   fitted regression model. \item \code{train_pred} The predicted values for
#'   the training data. \item \code{test_pred} The predicted values for the test
#'   data. \item \code{train_errors}: The training error metrics. \item
#'   \code{test_errors}: The test error metrics.}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{bpr_optimize}}, \code{\link{create_basis}}
#'
#' @examples
#' obs <- meth_data
#' y   <- gex_data
#' basis <- create_rbf_object(M = 5)
#' out   <- bpr_predict_wrap(x = obs, y = y, basis = basis,
#'                           is_parallel = FALSE, opt_itnmax = 3)
#'
#' @export
bpr_predict_wrap <- function(formula = NULL, x, y, model_name = "svm", w = NULL,
                             basis = NULL, train_ind = NULL, train_perc = 0.7,
                             fit_feature = "RMSE", cov_feature = TRUE,
                             opt_method = "CG", opt_itnmax = 100,
                             is_parallel = TRUE, no_cores = NULL,
                             is_summary = TRUE){

    # Check that x is a list object
    assertthat::assert_that(is.list(x))

    # Learn methylation profiles for each gene promoter region
    message("Learning methylation profiles ...\n")
    out_opt <- bpr_optim(x           = x,
                         w           = w,
                         basis       = basis,
                         opt_method  = opt_method,
                         opt_itnmax  = opt_itnmax,
                         is_parallel = is_parallel,
                         no_cores    = no_cores)

    W <- out_opt$W_opt
    if (out_opt$basis$M != 0) {
        if (cov_feature) { W <- cbind(W, out_opt$coverage_feat)}
        if (identical(fit_feature, "RMSE")) { W <- cbind(W, out_opt$rmse_feat)}
        if (identical(fit_feature, "NLL")) { W <- cbind(W, out_opt$nll_feat)}
    }

    # Create training and test sets
    message("Partitioning to test and train data ...\n")
    dataset <- .partition_data(x          = W,
                               y          = y,
                               train_ind  = train_ind,
                               train_perc = train_perc)

    # Train regression model from methylation profiles
    message("Training linear regression model ...\n")
    train_model <- inner_train_model_expr(formula = formula,
        model_name = model_name, train = dataset$train, is_summary = is_summary)

    # Predict gene expression from methylation profiles
    message("Making predictions ...\n")
    predictions <- inner_predict_model_expr(model = train_model$model,
        test = dataset$test, is_summary = is_summary)
    message("Done!\n\n")

    # Create 'bpr_predict' object
    obj <- structure(list(formula      = formula,
                          W_opt        = out_opt$W_opt,
                          model_name   = model_name,
                          basis        = out_opt$basis,
                          train_ind    = dataset$train_ind,
                          train_perc   = train_perc,
                          fit_feature  = fit_feature,
                          cov_feature  = cov_feature,
                          opt_method   = opt_method,
                          opt_itnmax   = opt_itnmax,
                          Mus          = out_opt$Mus,
                          train        = dataset$train,
                          test         = dataset$test,
                          gex_model    = train_model$model,
                          train_pred   = train_model$train_pred,
                          test_pred    = predictions$test_pred,
                          train_errors = train_model$train_errors,
                          test_errors  = predictions$test_errors),
                     class = "bpr_predict")
    return(obj)
}



#' @name bpr_optimize
#' @rdname bpr_optimize
#' @aliases bpr_optimise
#'
#' @title (DEPRECATED) Optimize BPR negative log likelihood function
#'
#' @description (DEPRECATED) The function bpr_optimize minimizes the negative
#'   log likelihood of the BPR function. Since it cannot be evaluated
#'   analytically, an optimization procedure is used. The
#'   \code{\link[stats]{optim}} packages is used for performing optimization.
#'
#' @param x The input object, either a \code{\link[base]{matrix}} or a
#'   \code{\link[base]{list}}.
#' @param ... Additional parameters.
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}.
#' @param lambda The complexity penalty coefficient for ridge regression.
#' @param opt_method The optimization method to be used. See
#'   \code{\link[stats]{optim}} for possible methods. Default is "CG".
#' @param opt_itnmax Optional argument giving the maximum number of iterations
#'   for the corresponding method. See \code{\link[stats]{optim}} for details.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 2.
#'
#' @return Depending on the input object \code{x}: \itemize{\item{If \code{x} is
#'   a \code{\link[base]{list}}:}  An object containing the following elements:
#'   \itemize{ \item{ \code{W_opt}: An Nx(M+1) matrix with the optimized
#'   parameter values. Each row of the matrix corresponds to each element of the
#'   list x. The columns are of the same length as the parameter vector w (i.e.
#'   number of basis functions). } \item{ \code{Mus}: An N x M matrix with the
#'   RBF centers if basis object is \code{\link{create_rbf_object}}, otherwise
#'   NULL.} \item{ \code{basis}: The basis object. } \item{ \code{w}: The
#'   initial values of the parameters w. } } \item{If \code{x} is a
#'   \code{\link[base]{matrix}}:} An object containing the following elements:
#'   \itemize{ \item{ \code{w_opt}: Optimized values for the coefficient vector
#'   w. The length of the result is the same as the length of the vector w. }
#'   \item{ \code{basis}: The basis object. } } \item{If calling
#'   \code{bpr_optim_fast} just the optimal weight matrix W_opt.} }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{eval_functions}}
NULL


#' @rdname bpr_optimize
#'
#' @examples
#' # Example of optimizing parameters for synthetic data using default values
#' data <- meth_data
#' out_opt <- bpr_optim(x = data, is_parallel = FALSE, opt_itnmax = 3)
#'
#' #-------------------------------------
#'
#' @export
bpr_optim <- function(x, ...){
    UseMethod("bpr_optim")
}


# Default function for the generic function 'bpr_optim'
bpr_optim.default <- function(x, ...){
    stop("Object x should be either matrix or list!")
}


#' @rdname bpr_optimize
#'
#' @examples
#' # Example of optimizing parameters for synthetic data using 3 RBFs
#' ex_data <- meth_data
#' basis <- create_rbf_object(M=3)
#' out_opt <- bpr_optim(x = ex_data, is_parallel = FALSE, basis = basis,
#'                      opt_itnmax = 3)
#'
#' #-------------------------------------
#'
#' @export
bpr_optim.list <- function(x, w = NULL, basis = NULL, lambda = 1/2,
                           opt_method = "CG", opt_itnmax = 100,
                           is_parallel = TRUE, no_cores = NULL, ...){
    assertthat::assert_that(is.list(x)) # Check that x is a list object
    i <- 0  # Initialize so the CMD check on R passes without NOTES
    N <- length(x)          # Extract number of observations
    assertthat::assert_that(N > 0)
    # Perform checks for initial parameter values
    out <- .do_checks(w = w, basis = basis); w <- out$w; basis <- out$basis
    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel)

    if (is_parallel) { # If parallel mode is ON
        # Create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
        # Parallel optimization for each element of x, i.e. for each region i.
        res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- bpr_optim.matrix(x = x[[i]], w = w, basis = basis,
            lambda = lambda, opt_method = opt_method, opt_itnmax = opt_itnmax)})
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of x, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- bpr_optim.matrix(x = x[[i]], w = w, basis = basis,
            lambda = lambda, opt_method = opt_method, opt_itnmax = opt_itnmax)})
    }

    # Matrix for storing optimized coefficients
    W_opt <- sapply(res, function(x) x$W_opt)
    if (is.matrix(W_opt)) { W_opt <- t(W_opt)
    }else {W_opt <- as.matrix(W_opt) }
    colnames(W_opt) <- paste("w", seq(1, NCOL(W_opt)), sep = "")
    nll_feat <- as.matrix(sapply(res, function(x) x$nll_feat))
    rmse_feat <- as.matrix(sapply(res, function(x) x$rmse_feat))
    coverage_feat <- as.matrix(sapply(res, function(x) x$coverage_feat))

    # Matrix for storing the centers of RBFs if object class is 'rbf'
    Mus <- NULL
    if (methods::is(basis, "rbf")) {
        if (is.null(basis$mus)) {
            Mus <- sapply(lapply(res, function(x) x$basis), function(y) y$mus)
            if (is.matrix(Mus)) { Mus <- t(Mus)
            }else {Mus <- as.matrix(Mus) }
            colnames(Mus) <- paste("mu", seq(1, NCOL(Mus)), sep = "")
        }
    }
    return(list(W_opt = W_opt, Mus = Mus, basis = basis, nll_feat = nll_feat,
                rmse_feat = rmse_feat, coverage_feat = coverage_feat, w = w))
}


#' @rdname bpr_optimize
#'
#' @examples
#' # Example of of specific promoter region using 2 RBFs
#' basis <- create_rbf_object(M=2)
#' w <- c(0.1, 0.1, 0.1)
#' data <- meth_data[[1]]
#' out_opt <- bpr_optim(x = data, w = w, basis = basis, fit_feature = "NLL",
#'                      opt_itnmax = 3)
#'
#' @importFrom stats optim
#'
#' @export
bpr_optim.matrix <- function(x, w = NULL, basis = NULL, lambda = 1/2,
                             opt_method = "CG", opt_itnmax = 100, ...){
    obs <- as.vector(x[, 1]) # Vector for storing CpG locations relative to TSS
    # Perform checks for initial parameter values
    out <- .do_checks(w = w, basis = basis); w <- out$w; basis <- out$basis
    # Create design matrix H
    des_mat <- design_matrix(obj = basis, obs = obs)
    H <- des_mat$H
    basis <- des_mat$basis

    # Call optim function to perform minimization of the NLL of BPR function
    w_opt <- optim(par = w, fn = bpr_log_likelihood, gr = bpr_gradient,
                   method = opt_method, control = list(maxit = opt_itnmax),
                   X = x, H = H, lambda = lambda, is_nll = TRUE)$par

    # NLL fit
    nll_feat <- bpr_log_likelihood(w = w_opt,X = x,H = H,lambda = lambda,
                                   is_nll = TRUE)
    # RMSE fit
    f_pred <- as.vector(pnorm(H %*% w_opt))
    if (NCOL(x) == 3) { f_true <- x[, 3] / x[, 2]
    }else{f_true <- x[, 2] }
    rmse_feat <- sqrt( mean( (f_pred - f_true) ^ 2) )
    # CpG coverage feature
    coverage_feat <- length(obs)

    return(list(W_opt = w_opt, basis = basis, nll_feat = nll_feat,
                rmse_feat = rmse_feat, coverage_feat = coverage_feat))
}


# @rdname bpr_optimize
#
# @examples
# # Example of optimizing parameters for synthetic data using 3 RBFs
# ex_data <- meth_data
# basis <- create_rbf_object(M=3)
# H <- list()
# for (i in 1:length(ex_data)){
#   H[[i]] <- design_matrix(obj = basis, obs = as.vector(ex_data[[i]][, 1]))$H
# }
# w <- rep(0.5, basis$M + 1)
# out_opt <- bpr_optim_fast(x = ex_data, H = H, w = w, is_parallel = FALSE,
#                           opt_itnmax = 3)
#
# @export
.bpr_optim_fast <- function(x, H, w = NULL, lambda = 1/6, opt_method = "CG",
                           opt_itnmax = 100, is_parallel = TRUE,
                           no_cores = NULL, ...){

    assertthat::assert_that(is.list(x)) # Check that x is a list object
    N <- length(x)                      # Extract number of observations
    if (is.null(w)) { w <- rep(0.5, NCOL(H[[1]])) }
    i <- 0  # Initialize so the CMD check on R passes without NOTES
    no_cores <- .parallel_cores(no_cores = no_cores, is_parallel = is_parallel)
    # If parallel mode is ON
    if (is_parallel) {
        # Create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)

        # Parallel optimization for each element of x, i.e. for each region i.
        res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- optim(par = w, fn = bpr_log_likelihood,
           gr = bpr_gradient, method = opt_method,
           control = list(maxit = opt_itnmax), X = x[[i]], H = H[[i]],
           lambda  = lambda, is_nll = TRUE)$par })
        # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }else{
        # Sequential optimization for each element of x, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
          ex = {out_opt <- optim(par = w, fn = bpr_log_likelihood,
           gr = bpr_gradient, method = opt_method,
           control = list(maxit = opt_itnmax), X = x[[i]], H = H[[i]],
           lambda = lambda, is_nll = TRUE)$par })
    }

    # Matrix for storing optimized coefficients
    W_opt <- sapply(res, function(x) x)
    if (is.matrix(W_opt)) { W_opt <- t(W_opt) }
    else {W_opt <- as.matrix(W_opt) }
    colnames(W_opt) <- paste("w", seq(1, NCOL(W_opt)), sep = "")
    return(list(W_opt = W_opt))
}



##-------------------------------------------------------------------



#' @title (DEPRECATED) Cluster methylation profiles
#'
#' @description (DEPRECATED) \code{bpr_cluster_wrap} is a wrapper function that
#'   clusters methylation profiles using the EM algorithm. Initially, it
#'   performs parameter checking, and initializes main parameters, such as
#'   mixing proportions, basis function coefficients, then the EM algorithm is
#'   applied and finally model selection metrics are calculated, such as BIC and
#'   AIC.
#'
#' @param x The binomial distributed observations, which has to be a list of
#'   elements of length N, where each element is an L x 3 matrix of
#'   observations, where 1st column contains the locations. The 2nd and 3rd
#'   columns contain the total trials and number of successes at the
#'   corresponding locations, repsectively. See
#'   \code{\link{process_haib_caltech_wrap}} on a possible way to get this data
#'   structure.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Mixing proportions.
#' @param w A MxK matrix, where each column consists of the basis function
#'   coefficients for each corresponding cluster.
#' @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}.
#' @param em_max_iter Integer denoting the maximum number of EM iterations.
#' @param epsilon_conv Numeric denoting the convergence parameter for EM.
#' @param init_opt_itnmax Optimization iterations for obtaining the initial EM
#'   parameter values.
#' @param is_verbose Logical, print results during EM iterations.
#' @inheritParams bpr_optimize
#'
#' @return A 'bpr_cluster' object which, in addition to the input parameters,
#'   consists of the following variables: \itemize{ \item{\code{pi_k}: Fitted
#'   mixing proportions.} \item{\code{w}: A MxK matrix with the fitted
#'   coefficients of the basis functions for each cluster k.} \item{\code{NLL}:
#'   The Negative Log Likelihood after the EM algorithm has finished.}
#'   \item{\code{r_nk}: Posterior probabilities of each promoter region
#'   belonging to each cluster.} \item{\code{labels}: Hard clustering
#'   assignments of each observation/promoter region.} \item{\code{BIC}:
#'   Bayesian Information Criterion metric.} \item{\code{AIC}: Akaike
#'   Information Criterion metric.} \item{\code{ICL}: Integrated Complete
#'   Likelihood criterion metric.} }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' ex_data <- meth_data
#' data_clust <- bpr_cluster_wrap(x = ex_data, em_max_iter = 3, opt_itnmax = 3,
#'                                init_opt_itnmax = 5, is_parallel = FALSE)
#'
#' @export
bpr_cluster_wrap <- function(x, K = 3, pi_k = NULL, w = NULL, basis = NULL,
                             em_max_iter = 100, epsilon_conv = 1e-04,
                             lambda = 1/2, opt_method = "CG", opt_itnmax = 100,
                             init_opt_itnmax = 100, is_parallel = TRUE,
                             no_cores = NULL, is_verbose = FALSE){
    assertthat::assert_that(is.list(x)) # Check that x is a list object
    N <- length(x)  # Extract number of observations
    assertthat::assert_that(N > 0)

    # Perform checks for initial parameter values
    out <- .do_EM_checks(x = x, K = K, pi_k = pi_k, w = w, basis = basis,
                         lambda = lambda, opt_method = opt_method,
                         init_opt_itnmax = init_opt_itnmax,
                         is_parallel = is_parallel, no_cores = no_cores)
    w <- out$w; basis <- out$basis; pi_k <- out$pi_k

    # Apply EM algorithm to cluster similar methylation profiles
    message("Clustering methylation profiles via EM ...\n")
    bpr_cluster <- .bpr_EM(x = x, K = K, pi_k = pi_k, w = w, basis = basis,
                           em_max_iter = em_max_iter,
                           epsilon_conv = epsilon_conv,
                           lambda = lambda, opt_method = opt_method,
                           opt_itnmax = opt_itnmax, is_parallel = is_parallel,
                           no_cores = no_cores, is_verbose = is_verbose)
    message("Finished clustering!\n\n")

    # Add names to the estimated parameters for clarity
    names(bpr_cluster$pi_k) <- paste0("clust", 1:K)
    colnames(bpr_cluster$w) <- paste0("clust", 1:K)

    # Get hard cluster assignments for each observation
    bpr_cluster$labels <- apply(X = bpr_cluster$r_nk, MARGIN = 1,
                        FUN = function(x) which(x == max(x, na.rm = TRUE)))

    # Perform model selection
    total_params <- (K - 1) + K * NROW(w)
    # Bayesian Information Criterion
    bpr_cluster$BIC <- 2 * utils::tail(bpr_cluster$NLL, n = 1) +
        total_params * log(N)
    # Akaike Iformation Criterion
    bpr_cluster$AIC <- 2 * utils::tail(bpr_cluster$NLL, n = 1) + 2*total_params
    # Integrated Complete Likelihood criterion
    entropy <- (-1) * sum(bpr_cluster$r_nk * log(bpr_cluster$r_nk),
                          na.rm = TRUE)
    bpr_cluster$ICL <- bpr_cluster$BIC + entropy
    # Store initial max optimization iterations
    bpr_cluster$init_opt_itnmax <- init_opt_itnmax
    class(bpr_cluster) <- "bpr_cluster"
    return(bpr_cluster)
}


# EM algorithm for BPR mixture model
#
# \code{.bpr_EM} implements the EM algorithm for performing clustering on DNA
#  methylation profiles, where the observation model is the Binomial
#  distributed Probit Regression function.
#
# @param x A list of elements of length N, where each element is an L x 3
# matrix of observations, where 1st column contains the locations. The 2nd
# and 3rd columns contain the total trials and number of successes at the
# corresponding locations, repsectively.
# @param K Integer denoting the number of clusters K.
# @param pi_k Vector of length K, denoting the mixing proportions.
# @param w A MxK matrix, where each column contains the basis function
# coefficients for the corresponding cluster.
# @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}.
# @param em_max_iter Integer denoting the maximum number of EM iterations.
# @param epsilon_conv Numeric denoting the convergence parameter for EM.
# @param opt_method The optimization method to be used. See
#  \code{\link[stats]{optim}} for possible methods. Default is "CG".
# @param opt_itnmax Optional argument giving the maximum number of iterations
#  for the corresponding method. See \code{\link[stats]{optim}} for details.
# @param is_parallel Logical, indicating if code should be run in parallel.
# @param no_cores Number of cores to be used, default is max_no_cores - 2.
# @param is_verbose Logical, print results during EM iterations
#
# @importFrom stats optim
.bpr_EM <- function(x, K = 2, pi_k = NULL, w = NULL, basis = NULL,
                    em_max_iter = 100, epsilon_conv = 1e-05, lambda = 1/2,
                    opt_method = "CG", opt_itnmax = 100, is_parallel = TRUE,
                    no_cores = NULL, is_verbose = FALSE){
    N <- length(x)  # Extract number of observations
    k <- 0
    weighted_pdf <- matrix(0, nrow = N, ncol = K) # Store weighted PDFs
    NLL <- c(1e+40) # Initialize and store NLL for each EM iteration
    # Create design matrix for each observation
    des_mat <- lapply(X = x,FUN = function(y) design_matrix(obj = basis,
                                                            obs = y[,1])$H)

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
            # For each element in x, evaluate the BPR log likelihood
            weighted_pdf[, k] <- vapply(X = 1:N, FUN = function(y)
                bpr_log_likelihood(w = w[, k], X = x[[y]], H = des_mat[[y]],
                                   lambda = lambda, is_nll = FALSE),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
            weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
        }

        # Calculate probs using the logSumExp trick for numerical stability
        Z <- apply(weighted_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        r_nk <- exp(weighted_pdf - Z)
        # Evaluate and store the NLL
        NLL  <- c(NLL, (-1) * sum(Z))

        #
        # M-Step -----------------------------------------------
        # Compute sum of posterior probabilities for each cluster
        N_k <- colSums(r_nk)
        # Update mixing proportions for each cluster
        pi_k <- N_k / N

        # Update basis function coefficient vector w for each cluster
        # If parallel mode is ON
        if (is_parallel) {
            # Parallel optimization for each cluster k
            w <- foreach::"%dopar%"(obj = foreach::foreach(k = 1:K,
                                                           .combine = cbind),
            ex = {out <- optim(par = w[, k], fn = sum_weighted_bpr_lik,
            gr = sum_weighted_bpr_grad, method = opt_method,
            control = list(maxit = opt_itnmax), X_list = x, H_list = des_mat,
            r_nk = r_nk[, k], lambda = lambda, is_nll = TRUE)$par})
        }else{
            # Sequential optimization for each clustrer k
            w <- foreach::"%do%"(obj = foreach::foreach(k = 1:K,
                                                        .combine = cbind),
             ex = {out <- optim(par = w[, k], fn = sum_weighted_bpr_lik,
             gr = sum_weighted_bpr_grad, method = opt_method,
             control = list(maxit = opt_itnmax), X_list = x, H_list = des_mat,
             r_nk = r_nk[, k], lambda = lambda, is_nll = TRUE)$par})
        }

        if (is_verbose) {
            cat("It:\t",t,"\tNLL:\t",NLL[t + 1],"\tDiff:\t",
                NLL[t] - NLL[t + 1],"\n")
        }
        if (NLL[t + 1] > NLL[t]) { message("NLL increases!\n") }
        # Check for convergence
        if (NLL[t] - NLL[t + 1] < epsilon_conv) { break }
        if (K == 1) { w <- as.matrix(w) }
    }
    if (is_parallel) {
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }
    if (K == 1) { w <- as.matrix(w) }
    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter) { warning("EM did not converge!\n") }

    obj <- structure(list(w = w, pi_k = pi_k, r_nk = r_nk, NLL = NLL,
                          basis = basis, K = K, N = N, lambda = lambda,
                          em_max_iter = em_max_iter, opt_method = opt_method,
                          opt_itnmax = opt_itnmax),
                     class = "bpr_EM")
    return(obj)
}


# @title Cluster methylation profiles fast method
#
# @description Efficient implementation for clustering methylation profiles
#   using the EM algorithm, given a design matrix H.
#
# @param x Observations, stored in a list object.
# @param H The design matrix H
# @param K Integer denoting the number of clusters K.
# @param pi_k Vector of length K, denoting the mixing proportions.
# @param w A MxK matrix, where each column consists of the basis function
#    coefficients for each corresponding cluster.
# @param em_max_iter Integer denoting the maximum number of EM iterations.
# @param epsilon_conv Numeric denoting the convergence parameter for EM.
# @param is_verbose Logical, print results during EM iterations.
# @inheritParams bpr_optimize
.bpr_EM_fast <- function(x, H, K = 2, pi_k = rep(1/K, K), w = NULL,
                        em_max_iter = 100, epsilon_conv = 1e-05, lambda = 1/2,
                        opt_method = "CG", opt_itnmax = 100, is_verbose=FALSE){

    N <- length(x)  # Extract number of observations
    weighted_pdf <- matrix(0, nrow = N, ncol = K) # Store weighted PDFs
    NLL <- c(1e+40) # Initialize and store NLL for each EM iteration
    # Run EM algorithm until convergence
    for (t in 1:em_max_iter) {
        # E-Step -----------------------------------------------
        for (k in 1:K) { # For each element in x, evaluate the BPR likelihood
            weighted_pdf[, k] <- vapply(X = 1:N, FUN = function(y)
                bpr_log_likelihood(w = w[, k], X = x[[y]], H = H[[y]],
                                   lambda = lambda, is_nll = FALSE),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
            weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
        }
        # Calculate probs using the logSumExp trick for numerical stability
        Z <- apply(weighted_pdf, 1, .log_sum_exp)
        r_nk <- exp(weighted_pdf - Z) # Get responsibilities
        NLL  <- c(NLL, (-1) * sum(Z))      # Evaluate and store the NLL

        # M-Step -----------------------------------------------
        #
        N_k <- colSums(r_nk)  # Sum of responsibilities for each cluster
        pi_k <- N_k / N            # Update mixing proportions for each cluster

        # Update basis function coefficient vector w for each cluster k
        for (k in 1:K) {
            w[, k] <- optim(par = w[, k], fn = sum_weighted_bpr_lik,
                            gr = sum_weighted_bpr_grad, method = opt_method,
                            control = list(maxit = opt_itnmax), X_list = x,
                            H_list = H, r_nk = r_nk[, k],
                            lambda = lambda, is_nll = TRUE)$par
        }
        if (is_verbose) {
            cat("It:\t",t,"\tNLL:\t",NLL[t + 1],"\tDiff:\t",
                NLL[t] - NLL[t + 1],"\n")
        }
        if (NLL[t + 1] > NLL[t]) {
            message("Negative Log Likelihood increases!\n")
        }
        if (K == 1) { w <- as.matrix(w) }
        # Check for convergence
        if (NLL[t] - NLL[t + 1] < epsilon_conv) { break }
    }
    # Store the object
    obj <- structure(list(w = w,pi_k = pi_k, r_nk = r_nk,
                          NLL = NLL[length(NLL)]), class = "bpr_EM_fast")
    return(obj)
}


# Internal function to make all the appropriate type checks.
.do_EM_checks <- function(x, K = 2, pi_k,  w, basis, lambda = 1/2,
                          opt_method = "CG", init_opt_itnmax = 100,
                          is_parallel = TRUE, no_cores = NULL){
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }
    if (is.null(w)) {
        w <- rep(0.5, basis$M + 1)

        # Optimize the BPR function for each element in x
        out_opt <- bpr_optim(x = x, w = w, basis = basis, fit_feature = NULL,
                             cov_feature = FALSE, lambda = lambda,
                             method = opt_method, itnmax = init_opt_itnmax,
                             is_parallel = is_parallel, no_cores = no_cores)
        W_opt <- out_opt$W_opt # Keep only the optimized coefficients
        cl <- stats::kmeans(W_opt, K, nstart = 25) # Use K-means
        C_n <- cl$cluster   # Get the mixture components
        w <- t(cl$centers)  # Mean for each cluster

        if (is.null(pi_k)) { # Mixing proportions
            N <- length(x)
            pi_k <- as.vector(table(C_n) / N)
        }
    }
    if (is.null(pi_k)) { pi_k <- rep(1 / K, K) }
    if (NROW(w) != (basis$M + 1) ) {
        stop("Coefficient vector should be M+1, M: number of basis functions!")
    }
    return(list(w = w, basis = basis, pi_k = pi_k))
}



#' @title (DEPRECATED) Plot the fit of methylation profiles across a region
#'
#' @description (DEPRECATED) \code{plot_fitted_profiles} is a simple function
#'   for plotting the methylation data across a give region, together with the
#'   fit of the methylation profiles.
#'
#' @param region Promoter region number
#' @param X Methylation data observations
#' @param fit_prof Fitted profile
#' @param fit_mean Fitted mean function
#' @param title Title of the plot
#' @param ... Additional parameters
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{old_plot_cluster_prof}}
#'
#' @examples
#' # Fit methylation profiles using 3 RBFs
#' obs <- meth_data
#' y   <- gex_data
#' basis <- create_rbf_object(M = 3)
#' out   <- bpr_predict_wrap(x = obs, y = y, basis = basis,
#'                           is_parallel = FALSE, opt_itnmax = 5)
#'
#' # Create the plot
#' old_plot_fitted_profiles(region = 16, X = meth_data, fit_prof = out)
#'
#' @export
old_plot_fitted_profiles <- function(region, X, fit_prof, fit_mean = NULL,
                                     title = "Gene promoter", ...){

    graphics::par(cex = 1.05, mai = c(1.37,1.37,.7,.3) )
    x <- X[[region]][,1]
    y <- X[[region]][,3]/X[[region]][,2]
    xs <- seq(from = -1, to = 1, by = 0.01)
    graphics::plot(x, y, col = "blue2", pch = 21, ylim = c(0,1),
                   xlim = c(-1,1), lwd = 0.8, xlab = NA, ylab = NA,
                   cex.axis = 1.1, xaxt = "n")
    graphics::mtext(side = 1, "genomic region", line = 3, cex = 1.2)
    graphics::mtext(side = 2, "methylation level", line = 3, cex = 1.2)
    graphics::axis(side = 1, at = c(-1, 0, 1), labels = c("-7kb","TSS","+7kb"))
    graphics::title(main = title, line = 1, cex.main = 1.4)
    if (!is.null(fit_mean)) {
        graphics::lines(x = xs, y = eval_probit_function(fit_mean$basis, xs,
                    fit_mean$W_opt[region, ]), col = 'coral', lwd = 2, lty = 2)
    }
    graphics::lines(x = xs, y = eval_probit_function(fit_prof$basis, xs,
                    fit_prof$W_opt[region, ]), col = 'red2', lwd = 2)
}

#' @title (DEPRECATED) Plot of clustered methylation profiles
#'
#' @description (DEPRECATED) \code{plot_cluster_prof} creates a plot of cluster
#'   methylation profiles, where each colour denotes a different cluster.
#'
#' @param bpr_cluster_obj The output of the \code{bpr_cluster_wrap} function.
#' @param main_lab The title of the plot
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{old_plot_fitted_profiles}},
#'   \code{\link{old_boxplot_cluster_gex}}
#'
#' @examples
#' # Cluster methylation profiles using 4 RBFs
#' obs <- meth_data
#' basis <- create_rbf_object(M = 4)
#' res   <- bpr_cluster_wrap(x = obs, K = 3, em_max_iter = 2, opt_itnmax = 3,
#'                           init_opt_itnmax = 2, is_parallel = FALSE)
#'
#' # Create the plot
#' old_plot_cluster_prof(bpr_cluster_obj = res)
#'
#' @export
old_plot_cluster_prof <- function(bpr_cluster_obj,
                                  main_lab = "Clustered methylation profiles"){
    graphics::par(mar = c(4.2, 4.1, 3.1, 2), xpd = TRUE)
    cols <- c("darkolivegreen4", "cornflowerblue",
              "coral", "firebrick","#E69F00")
    xs <- seq(-1,1,len = 2000) # create some values
    graphics::plot(x = xs,
                   y = eval_probit_function(bpr_cluster_obj$basis, xs,
                                            bpr_cluster_obj$w[, 1]),
                   xlim = c(-1, 1), ylim = c(0, 1),
                   type = "l", col = cols[1], lwd = 4,
                   xlab = "promoter region",
                   ylab = "methylation level",
                   main = main_lab)
    K <- 5
    if (bpr_cluster_obj$K < 5) { K <- bpr_cluster_obj$K }
    for (k in 2:K) {
        graphics::lines(x = xs, y = eval_probit_function(bpr_cluster_obj$basis,
                           xs, bpr_cluster_obj$w[, k]), col = cols[k], lwd = 4)
    }
}

#' @title (DEPRECATED) Boxplot of clustered expression levels
#'
#' @description (DEPRECATED) \code{boxplot_cluster_gex} creates a boxplot of
#'   clustered gene expression levels which depend on the clustered methylation
#'   profiles. Each colour denotes a different cluster.
#'
#' @param bpr_cluster_obj The output of the \code{bpr_cluster_wrap} function.
#' @param gex The vector of gene expression data for each promoter region.
#' @param main_lab The title of the plot
#'
#' @return The figure to be plotted in the device.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{old_plot_cluster_prof}},
#'   \code{\link{old_plot_fitted_profiles}}
#'
#' @examples
#' # Cluster methylation profiles using 4 RBFs
#' obs <- meth_data
#' basis <- create_rbf_object(M = 4)
#' res   <- bpr_cluster_wrap(x = obs, K = 3, em_max_iter = 2, opt_itnmax = 3,
#'                           init_opt_itnmax = 2, is_parallel = FALSE)
#'
#' # Create the plot
#' old_boxplot_cluster_gex(bpr_cluster_obj = res, gex = gex_data)
#'
#' @export
old_boxplot_cluster_gex <- function(bpr_cluster_obj, gex,
                                    main_lab = "Gene expression levels"){
    graphics::par(mar = c(4.2, 4.1, 3.1, 5.5), xpd = TRUE)

    cols <- c("darkolivegreen4","cornflowerblue","coral","firebrick","#E69F00")
    gex_list <- list()
    for (k in 1:bpr_cluster_obj$K) {
        gex_list[[k]] <- gex[which(bpr_cluster_obj$labels == k)]
    }
    graphics::boxplot(gex_list, col = cols[1:bpr_cluster_obj$K], notch = TRUE,
                      xlab = "Cluster K", ylab = "expression level",
                      main = main_lab)
    graphics::legend("right", inset = c(-0.18,0),
                     legend = seq(1,bpr_cluster_obj$K),
                     lty = 1, lwd = 3, col = cols[1:bpr_cluster_obj$K],
                     title = "Cluster")
}

