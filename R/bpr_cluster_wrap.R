#' Cluster methylation profiles
#'
#' \code{bpr_cluster_wrap} is a wrapper function that clusters methylation
#' profiles using the EM algorithm. Initially, it performs parameter checking,
#' and initializes main parameters, such as mixing proportions, basis function
#' coefficients, then the EM algorithm is applied and finally model selection
#' metrics are calculated, such as BIC and AIC.
#'
#' @param x The binomial distributed observations, which has to be a list of
#'   elements of length N, where each element is an L x 3 matrix of
#'   observations, where 1st column contains the locations. The 2nd and 3rd
#'   columns contain the total trials and number of successes at the
#'   corresponding locations, repsectively. See
#'   \code{\link{process_haib_caltech_wrap}} on a possible way to get this data
#'   structure.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
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
#'   \item{\code{post_prob}: Posterior probabilities of each promoter region
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
#' data_clust <- bpr_cluster_wrap(x = ex_data, em_max_iter = 3,
#'                    is_parallel = FALSE, opt_itnmax = 10)
#'
#' @export
bpr_cluster_wrap <- function(x, K = 3, pi_k = NULL, w = NULL, basis = NULL,
                             em_max_iter = 100, epsilon_conv = 1e-04,
                             opt_method = "CG", opt_itnmax = 100,
                             init_opt_itnmax = 100, is_parallel = TRUE,
                             no_cores = NULL, is_verbose = FALSE){
  # Check that x is a list object
  assertthat::assert_that(is.list(x))

  # Extract number of observations
  N <- length(x)
  assertthat::assert_that(N > 0)

  # Perform checks for initial parameter values
  out <- .do_EM_checks(x = x,
                       K = K,
                       pi_k = pi_k,
                       w = w,
                       basis = basis,
                       opt_method = opt_method,
                       init_opt_itnmax = init_opt_itnmax,
                       is_parallel = is_parallel,
                       no_cores    = no_cores)
  w     <- out$w
  basis <- out$basis
  pi_k  <- out$pi_k

  # Apply EM algorithm to cluster similar methylation profiles
  message("Clustering methylation profiles via EM ...\n")
  bpr_cluster <- .bpr_EM(x = x,
                         K = K,
                         pi_k = pi_k,
                         w = w,
                         basis = basis,
                         em_max_iter = em_max_iter,
                         epsilon_conv = epsilon_conv,
                         opt_method = opt_method,
                         opt_itnmax = opt_itnmax,
                         is_parallel = is_parallel,
                         no_cores    = no_cores,
                         is_verbose = is_verbose)
  message("Finished clustering!\n\n")

  # Add names to the estimated parameters for clarity
  names(bpr_cluster$pi_k) <- paste0("clust", 1:K)
  colnames(bpr_cluster$w) <- paste0("clust", 1:K)

  # Get hard cluster assignments for each observation
  bpr_cluster$labels <- apply(X      = bpr_cluster$post_prob,
                              MARGIN = 1,
                              FUN    = function(x)
                                which(x == max(x, na.rm = TRUE)))

  # Perform model selection
  total_params <- (K - 1) + K * NROW(w)

  # Bayesian Information Criterion
  bpr_cluster$BIC <- 2 * utils::tail(bpr_cluster$NLL, n = 1) +
                          total_params * log(N)

  # Akaike Iformation Criterion
  bpr_cluster$AIC <- 2 * utils::tail(bpr_cluster$NLL, n = 1) + 2 * total_params

  # Integrated Complete Likelihood criterion
  entropy <- (-1) * sum(bpr_cluster$post_prob * log(bpr_cluster$post_prob),
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
#
.bpr_EM <- function(x, K = 2, pi_k = NULL, w = NULL, basis = NULL,
                    em_max_iter = 100, epsilon_conv = 1e-05, opt_method = "CG",
                    opt_itnmax = 100, is_parallel = TRUE, no_cores = NULL,
                    is_verbose = FALSE){

  # Extract number of observations
  N <- length(x)

  # Store weighted PDFs
  weighted_pdf <- matrix(0, nrow = N, ncol = K)

  # Initialize and store NLL for each EM iteration
  NLL <- c(1e+40)

  # If parallel mode is ON
  if (is_parallel){
    # If number of cores is not given
    if (is.null(no_cores)){
      no_cores <- parallel::detectCores() - 2
    }else{
      if (no_cores >= parallel::detectCores()){
        no_cores <- parallel::detectCores() - 1
      }
    }
    if (is.na(no_cores)){
      no_cores <- 2
    }
    if (no_cores > K){
      no_cores <- K
    }
    # Create cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
  }

  if (is_parallel){
    # Create design matrix for each observation
    des_mat <- parallel::mclapply(X   = x,
                                  FUN = function(y)
                                    .design_matrix(x = basis, obs = y[, 1]),
                                  mc.cores = no_cores)
  }else{
    # Create design matrix for each observation
    des_mat <- lapply(X   = x,
                      FUN = function(y)
                        .design_matrix(x = basis, obs = y[, 1]))
  }

  # Run EM algorithm until convergence
  for (t in 1:em_max_iter){

    #
    # E-Step -----------------------------------------------
    #
    # Compute weighted pdfs for each cluster
    for (k in 1:K){
      # For each element in x, evaluate the BPR log likelihood
      weighted_pdf[, k] <- vapply(X   = 1:N,
                                  FUN = function(y)
                                    .bpr_likelihood(w = w[, k],
                                                    H = des_mat[[y]]$H,
                                                    data = x[[y]][, 2:3],
                                                    is_NLL = FALSE),
                                  FUN.VALUE = numeric(1),
                                  USE.NAMES = FALSE)
      weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z <- apply(weighted_pdf, 1, .log_sum_exp)
    # Get actual posterior probabilities, i.e. responsibilities
    post_prob <- exp(weighted_pdf - Z)
    # Evaluate and store the NLL
    NLL  <- c(NLL, (-1) * sum(Z))

    #
    # M-Step -----------------------------------------------
    #
    # Compute sum of posterior probabilities for each cluster
    N_k <- colSums(post_prob)
    # Update mixing proportions for each cluster
    pi_k <- N_k / N

    # Update basis function coefficient vector w for each cluster
    # If parallel mode is ON
    if (is_parallel){
      # Parallel optimization for each cluster k
      w <- foreach::"%dopar%"(obj = foreach::foreach(k = 1:K,
                                                     .combine = cbind),
                        ex  = {
                          out <- optim(par       = w[, k],
                                       fn        = .sum_weighted_bpr_lik,
                                       gr        = .sum_weighted_bpr_grad,
                                       method    = opt_method,
                                       control   = list(maxit = opt_itnmax),
                                       x         = x,
                                       des_mat   = des_mat,
                                       post_prob = post_prob[, k],
                                       is_NLL    = TRUE)$par
                        })
    }else{
      # Sequential optimization for each clustrer k
      w <- foreach::"%do%"(obj = foreach::foreach(k = 1:K,
                                                  .combine = cbind),
                       ex  = {
                         out <- optim(par       = w[, k],
                                      fn        = .sum_weighted_bpr_lik,
                                      gr        = .sum_weighted_bpr_grad,
                                      method    = opt_method,
                                      control   = list(maxit = opt_itnmax),
                                      x         = x,
                                      des_mat   = des_mat,
                                      post_prob = post_prob[, k],
                                      is_NLL    = TRUE)$par
                       })
    }

    if (is_verbose){
      cat("It:\t", t, "\tNLL:\t", NLL[t + 1],
          "\tNLL_diff:\t", NLL[t] - NLL[t + 1], "\n")
    }
    if (NLL[t + 1] > NLL[t]){
      stop("Negative Log Likelihood increases - Stopping EM!\n")
    }
    # Check for convergence
    if (NLL[t] - NLL[t + 1] < epsilon_conv){
      break
    }
  }
  if (is_parallel){
    # Stop parallel execution
    parallel::stopCluster(cl)
  }
  # Check if EM converged in the given maximum iterations
  if (t == em_max_iter){
    warning("EM did not converge with the given maximum iterations!\n")
  }

  obj <- structure(list(K = K,
                        N = N,
                        w = w,
                        pi_k = pi_k,
                        em_max_iter = em_max_iter,
                        opt_method = opt_method,
                        opt_itnmax = opt_itnmax,
                        NLL = NLL,
                        basis = basis,
                        post_prob = post_prob),
                   class = "bpr_EM")
  return(obj)
}



# Internal function to make all the appropriate type checks.
.do_EM_checks <- function(x, K = 2, pi_k = NULL,  w = NULL, basis = NULL,
                          opt_method = "CG", init_opt_itnmax = 100,
                          is_parallel = TRUE, no_cores = NULL){
  if (is.null(basis)){
    basis <- create_rbf_object(M = 3)
  }
  if (is.null(w)){
    w <- rep(0.5, basis$M + 1)

    # Optimize the BPR function for each element in x
    out_opt <- bpr_optim(x           = x,
                         w           = w,
                         basis       = basis,
                         fit_feature = NULL,
                         cpg_dens_feat = FALSE,
                         method      = opt_method,
                         itnmax      = init_opt_itnmax,
                         is_parallel = is_parallel,
                         no_cores    = no_cores)

    # Keep only the optimized coefficients
    W_opt <- out_opt$W_opt

    # Use Kmeans with random starts
    cl <- stats::kmeans(W_opt, K, nstart = 25)
    # Get the mixture components
    C_n <- cl$cluster
    # Mean for each cluster
    w <- t(cl$centers)

    # Mixing proportions
    if (is.null(pi_k)){
      N <- length(x)
      pi_k <- as.vector(table(C_n) / N)
    }
  }
  if (is.null(pi_k)){
    pi_k <- rep(1 / K, K)
  }
  if (NROW(w) != (basis$M + 1) ){
    stop("Coefficients vector should be M+1, M: number of basis functions!")
  }
  return(list(w = w, basis = basis, pi_k = pi_k))
}
