#' EM algorithm for BPR mixture model
#'
#' \code{bpr_EM} implements the EM algorithm for performing clustering on DNA
#'  methylation profiles, where the observation model is the Binomial
#'  distributed Probit Regression function, \code{\link{bpr_likelihood}}.
#'
#' @param x A list of elements of length N, where each element is an L x 3
#' matrix of observations, where 1st column contains the locations. The 2nd
#' and 3rd columns contain the total trials and number of successes at the
#' corresponding locations, repsectively.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A MxK matrix, where each column contains the basis function
#' coefficients for the corresponding cluster.
#' @param basis A 'basis' object. E.g. see \code{\link{polynomial.object}}
#' @param em_max_iter Integer denoting the maximum number of EM iterations.
#' @param epsilon_conv Numeric denoting the convergence parameter for EM.
#' @param opt_method The optimization method to be used. See
#'  \code{\link[stats]{optim}} for possible methods. Default is 'CG'.
#' @param opt_itnmax Optional argument giving the maximum number of iterations
#'  for the corresponding method. See \code{\link[stats]{optim}} for details.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_verbose Logical, print results during EM iterations
#'
#' @importFrom stats optim
#'
#' @export
bpr_EM <- function(x, K = 2, pi_k = NULL, w = NULL, basis = NULL,
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
                                    design_matrix(x = basis, obs = y[ ,1]),
                                  mc.cores = no_cores)
  }else{
    # Create design matrix for each observation
    des_mat <- lapply(X   = x,
                      FUN = function(y)
                        design_matrix(x = basis, obs = y[ ,1]))
  }

  # Run EM algorithm until convergence
  for (t in 1:em_max_iter){

    #
    # E-Step -----------------------------------------------
    #
    # Compute weighted pdfs for each cluster
    for (k in 1:K){
      # For each element in x, evaluate the BPR log likelihood
      weighted_pdf[ ,k] <- vapply(X   = 1:N,
                                  FUN = function(y)
                                    bpr_likelihood(w = w[ ,k],
                                                   H = des_mat[[y]]$H,
                                                   data = x[[y]][ ,2:3],
                                                   is_NLL = FALSE),
                                  FUN.VALUE = numeric(1),
                                  USE.NAMES = FALSE)
      weighted_pdf[ ,k] <- log(pi_k[k]) + weighted_pdf[ ,k]
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z <- apply(weighted_pdf, 1, log_sum_exp)

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
                        out <- optim(par       = w[ ,k],
                                     fn        = sum_weighted_bpr_lik,
                                     gr        = sum_weighted_bpr_grad,
                                     method    = opt_method,
                                     control   = list(maxit = opt_itnmax),
                                     x         = x,
                                     des_mat   = des_mat,
                                     post_prob = post_prob[ ,k],
                                     is_NLL    = TRUE)$par
                              })
    }else{
      # Sequential optimization for each clustrer k
      w <- foreach::"%do%"(obj = foreach::foreach(k = 1:K,
                                                  .combine = cbind),
                           ex  = {
                     out <- optim(par       = w[ ,k],
                                  fn        = sum_weighted_bpr_lik,
                                  gr        = sum_weighted_bpr_grad,
                                  method    = opt_method,
                                  control   = list(maxit = opt_itnmax),
                                  x         = x,
                                  des_mat   = des_mat,
                                  post_prob = post_prob[ ,k],
                                  is_NLL    = TRUE)$par
                           })
    }

    if (is_verbose){
      cat("It:\t",t, "\tNLL:\t", NLL[t + 1],
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
