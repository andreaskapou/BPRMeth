#' Cluster similar methylation profiles
#'
#' \code{mpgex_cluster} is a wrapper function that clusters similar methylation
#' profiles using the EM algorithm. Initially, it performs parameter checking,
#' and initializes main parameters, such as mixing proportions, basis function
#' coefficients, then the EM algorithm is applied and finally model selection
#' metrics are calculated, such as BIC and AIC.
#'
#' @param init_opt_itnmax Optimization iterations for obtaining the initial EM
#'  parameter values.
#' @inheritParams bpr_EM
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' ex_data <- meth_data
#' mpgex_clust <- bpr_cluster_wrap(x = ex_data, em_max_iter = 3, is_parallel = FALSE, opt_itnmax = 10)
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
  mpgex <- bpr_EM(x = x,
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
  message("Done!\n\n")

  # Add names to the estimated parameters for clarity
  names(mpgex$pi_k) <- paste0("clust", 1:K)
  colnames(mpgex$w) <- paste0("clust", 1:K)

  # Get hard cluster assignments for each observation
  mpgex$labels <- apply(X      = mpgex$post_prob,
                        MARGIN = 1,
                        FUN    = function(x)
                          which(x == max(x, na.rm = TRUE)))

  # Perform model selection
  total_params <- (K - 1) + K * NROW(w)

  # Bayesian Information Criterion
  mpgex$BIC <- 2 * utils::tail(mpgex$NLL, n = 1) + total_params * log(N)

  # Akaike Iformation Criterion
  mpgex$AIC <- 2 * utils::tail(mpgex$NLL, n = 1) + 2 * total_params

  # Integrated Complete Likelihood criterion
  entropy <- (-1) * sum(mpgex$post_prob * log(mpgex$post_prob),
                        na.rm = TRUE)
  mpgex$ICL <- mpgex$BIC + entropy

  # Store initial max optimization iterations
  mpgex$init_opt_itnmax <- init_opt_itnmax

  class(mpgex) <- "mpgex_cluster"

  return(mpgex)
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
