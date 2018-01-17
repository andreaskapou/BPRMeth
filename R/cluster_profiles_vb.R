#' @name cluster_profiles_vb
#' @rdname cluster_profiles_vb
#' @aliases cluster_profile_vb cluster_vb
#'
#' @title Cluster methylation profiles using VB
#'
#' @description General purpose functions for clustering latent profiles for
#'   different observation models using Variational Bayes (VB) EM-like
#'   algorithm.
#'
#' @param X The input data, which has to be a \code{\link[base]{list}} of
#'   elements of length N, where each element is an \code{L X C} matrix, where L
#'   are the total number of observations. The first column contains the input
#'   observations x (i.e. CpG locations). If "binomial" model then C=3, and 2nd
#'   and 3rd columns contain total number of trials and number of successes
#'   respectively. If "bernoulli" or "gaussian" model, then C=2 containing the
#'   output y (e.g. methylation level).
#' @param K Integer denoting the total number of clusters K.
#' @param delta_0 Parameter vector of the Dirichlet prior on the mixing
#'   proportions pi.
#' @param w Optional, an (M+1)xK matrix of the initial parameters, where each
#'   column consists of the basis function coefficients for each corresponding
#'   cluster k. If NULL, will be assigned with default values.
#' @inheritParams infer_profiles_vb
#'
#' @return An object of class \code{cluster_profiles_vb_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An (M+1) X K matrix with the
#'   optimized parameter values for each cluster, M are the number of basis
#'   functions. Each column of the matrix corresponds a different cluster k.}
#'   \item{ \code{W_Sigma}: A list with the covariance matrices of the posterior
#'   parmateter W for each cluster k.} \item{ \code{r_nk}: An (N X K)
#'   responsibility matrix of each observations being explained by a specific
#'   cluster. } \item{ \code{delta}: Optimized Dirichlet paramter for the mixing
#'   proportions. } \item{ \code{alpha}: Optimized shape parameter of Gamma
#'   distribution. } \item{ \code{beta}: Optimized rate paramter of the Gamma
#'   distribution } \item{ \code{basis}: The basis object. } \item{\code{lb}:
#'   The lower bound vector.} \item{\code{labels}: Cluster assignment labels.}
#'   \item{ \code{pi_k}: Expected value of mixing proportions.} }
#'
#' @section Details: The modelling and mathematical details for clustering
#'   profiles using mean-field variational inference are explained here:
#'   \url{http://rpubs.com/cakapourani/} . More specifically: \itemize{\item{For
#'   Binomial/Bernoulli observation model check:
#'   \url{http://rpubs.com/cakapourani/vb-mixture-bpr}} \item{For Gaussian
#'   observation model check: \url{http://rpubs.com/cakapourani/vb-mixture-lr}}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{cluster_profiles_mle}}
#'   \code{\link{infer_profiles_vb}}, \code{\link{infer_profiles_mle}},
#'   \code{\link{infer_profiles_gibbs}}, \code{\link{create_region_object}}
NULL


#' @rdname cluster_profiles_vb
#'
#' @examples
#' # Example of optimizing parameters for synthetic data using 3 RBFs
#' basis <- create_rbf_object(M=3)
#' out <- cluster_profiles_vb(X = binomial_data, model = "binomial",
#'   basis=basis, vb_max_iter = 10)
#'
#' #-------------------------------------
#'
#' basis <- create_rbf_object(M=3)
#' out <- cluster_profiles_vb(X = gaussian_data, model = "gaussian",
#'   basis=basis, vb_max_iter = 10)
#'
#' @export
cluster_profiles_vb <- function(X, K = 3, model = NULL, basis = NULL, H = NULL,
                                delta_0 = rep(1/K, K), w = NULL,
                                gaussian_l = .5, alpha_0 = 1e-1, beta_0 = 1e-1,
                                vb_max_iter = 100, epsilon_conv = 1e-4,
                                is_verbose = FALSE, ...){
    assertthat::assert_that(is.list(X))  # Check that X is a list object
    if (is.null(model)) { stop("Observation model not defined!") }
    # Create RBF basis object by default
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }
    if (is.null(H)) { # Create design matrix
        H <- lapply(X = X, FUN = function(x) design_matrix(basis, x[, 1])$H)
    }
    if (is.null(w)) {
        # Optimal weight vector for each region using MLE
        w <- infer_profiles_mle(X = X, model = model, basis = basis, H = H,
                                lambda = .5, opt_itnmax = 25)$W
        # Use Kmeans with random starts
        cl <- stats::kmeans(w, K, nstart = 25)
        # Mean for each cluster
        w <- matrix(t(cl$centers), ncol = K, nrow = basis$M + 1)
    }
    # If we have Binomial distributed data
    if (identical(model, "binomial") || NCOL(X[[1]]) == 3) {
        T_n = M_n = y <- list()
        for (n in 1:length(X)) {
            T_n[[n]] <- c(X[[n]][, 2]) # Total number of reads for each CpG
            M_n[[n]] <- c(X[[n]][, 3]) # Methylated reads for each CpG
            y[[n]] <- vector(mode = "integer") # Output vector y
            for (i in 1:NROW(X[[n]])) {
                y[[n]] <- c(y[[n]], rep(1, M_n[[n]][i]),
                            rep(0, T_n[[n]][i] - M_n[[n]][i]))
            }
            # Create extended design matrix from binary observations
            H[[n]] <- as.matrix(H[[n]][rep(1:nrow(H[[n]]), T_n[[n]]), ])
        }
    }else{
        y <- lapply(X = X, FUN = function(x) x[,2])  # Extract responses y_{n}
    }

    if (identical(model, "bernoulli") || identical(model, "binomial") ) {
        obj <- .cluster_profiles_vb_bpr(H = H, y = y, K = K, model = model,
            basis = basis, w = w, delta_0 = delta_0, alpha_0 = alpha_0,
            beta_0 = beta_0, vb_max_iter = vb_max_iter,
            epsilon_conv = epsilon_conv, is_verbose = is_verbose)
    }else if (identical(model, "gaussian")) {
        obj <- .cluster_profiles_vb_gaussian(H = H, y = y, K = K, basis = basis,
            w = w, delta_0 = delta_0, gaussian_l = gaussian_l,
            alpha_0 = alpha_0, beta_0 = beta_0, vb_max_iter = vb_max_iter,
            epsilon_conv = epsilon_conv, is_verbose = is_verbose)
    }else {stop(paste0(model, " observation model not implented.")) }

    # Add names to the estimated parameters for clarity
    names(obj$delta) <- paste0("cluster", 1:K)
    names(obj$pi_k) <- paste0("cluster", 1:K)
    colnames(obj$W) <- paste0("cluster", 1:K)
    colnames(obj$r_nk) <- paste0("cluster", 1:K)
    # Get hard cluster assignments for each observation
    obj$labels <- apply(X = obj$r_nk, MARGIN = 1,
                        FUN = function(x) which(x == max(x, na.rm = TRUE)))
    # Subtract \ln(K!) from Lower bound to get a better estimate
    obj$lb <- obj$lb - log(factorial(K))
    # Append general class
    class(obj) <- append(class(obj),
                         c("cluster_profiles_vb", "cluster_profiles"))
    # Return object
    return(obj)
}


##------------------------------------------------------

# Cluster profiles EM Binomial or Bernoulli
.cluster_profiles_vb_bpr <- function(H, y, K, model, basis, w, delta_0, alpha_0,
                                     beta_0, vb_max_iter, epsilon_conv,
                                     is_verbose){
    assertthat::assert_that(is.list(H)) # Check that H is a list object
    D <- basis$M + 1        # Number of covariates
    N <- length(H)          # Extract number of observations
    assertthat::assert_that(N > 0)
    LB <- c(-Inf)           # Store the lower bounds
    r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = N, ncol = K)
    E_ww <- vector("numeric", length = K)

    # Compute H_{n}'H_{n}
    HH <- lapply(X = H, FUN = function(h) crossprod(h))
    y_1 <- lapply(1:N, function(n) which(y[[n]] == 1))
    y_0 <- lapply(1:N, function(n) which(y[[n]] == 0))
    E_z = E_zz <- y

    # Compute shape parameter of Gamma
    alpha_k <- rep(alpha_0 + D/2, K)
    # TODO: Remove this
    # Use Kmeans with random starts
    # cl  <- stats::kmeans(w, K, nstart = 25)
    # Create a hot 1-K encoding
    # C_n <- as.matrix(Matrix::sparseMatrix(1:N, cl$cluster, x = 1))
    # Mean for each cluster
    m_k <- w
    # Covariance of each cluster
    S_k <- lapply(X = 1:K, function(k) solve(diag(2, D)))
    # Scale of precision matrix
    beta_k   <- rep(beta_0, K)
    # Dirichlet parameter
    delta_k  <- delta_0
    # Expectation of log Dirichlet
    e_log_pi <- digamma(delta_k) - digamma(sum(delta_k))
    mk_Sk    <- lapply(X = 1:K, function(k) tcrossprod(m_k[, k]) + S_k[[k]])
    # Update \mu, E[z] and E[z^2]
    mu <- lapply(1:N, function(n) c(H[[n]] %*% rowMeans(m_k)))
    # mu <- lapply(1:N, function(n) c(H[[n]] %*% rowSums(sapply(X = 1:K,
    #                    function(k) C_n[n,k]*m_k[,k]))))
    # Ensure that \mu is not 0 or 1
    for (my_iter in 1:length(mu)) {
        mu[[my_iter]][which(mu[[my_iter]] > (1 - 1e-15))] <- 1 - 1e-15
        mu[[my_iter]][which(mu[[my_iter]] < 1e-15)] <- 1e-15
    }
    E_z <- lapply(1:N, function(n) .update_Ez(E_z = E_z[[n]], mu = mu[[n]],
                                              y_0 = y_0[[n]], y_1 = y_1[[n]]))
    E_zz <- lapply(1:N, function(n) 1 + mu[[n]] * E_z[[n]])

    # Show progress bar
    if (is_verbose) { pb <- utils::txtProgressBar(min = 2, max = vb_max_iter,
                                                  style = 3) }
    # Run VB algorithm until convergence
    for (i in 2:vb_max_iter) {
        ##-------------------------------
        # Variational E-Step
        ##-------------------------------
        for (k in 1:K) {
            log_rho_nk[,k] <- e_log_pi[k] + sapply(1:N,
              function(n) -0.5*crossprod(E_zz[[n]]) + m_k[,k] %*% t(H[[n]]) %*%
                  E_z[[n]] - 0.5*matrix.trace(HH[[n]] %*% mk_Sk[[k]]))
        }
        # Calculate responsibilities using logSumExp for numerical stability
        log_r_nk <- log_rho_nk - apply(log_rho_nk, 1, .log_sum_exp)
        r_nk     <- exp(log_r_nk)

        ##-------------------------------
        # Variational M-Step
        ##-------------------------------
        # Update Dirichlet parameter
        delta_k <- delta_0 + colSums(r_nk)
        # TODO: Compute expected value of mixing proportions
        pi_k <- (delta_0 + colSums(r_nk)) / (K * delta_0 + N)
        for (k in 1:K) {
            # Update covariance for Gaussian
            w_HH <- lapply(X = 1:N, function(n) HH[[n]]*r_nk[n,k])
            S_k[[k]] <- solve(diag(alpha_k[k]/beta_k[k], D) + .add_func(w_HH))
            # Update mean for Gaussian
            w_Hz <- lapply(X = 1:N,function(n) t(H[[n]]) %*% E_z[[n]]*r_nk[n,k])
            m_k[,k] <- S_k[[k]] %*% .add_func(w_Hz)
            # Update \beta_k parameter for Gamma
            E_ww[k] <- crossprod(m_k[,k]) + matrix.trace(S_k[[k]])
            beta_k[k]  <- beta_0 + 0.5*E_ww[k]
        }
        # Update \mu, E[z] and E[z^2]
        if (D == 1) {
            mu <- lapply(X = 1:N, FUN = function(n) c(H[[n]] %*%
                    sum(sapply(X = 1:K, function(k) r_nk[n,k]*m_k[,k]))))
        }else {
            mu <- lapply(X = 1:N, FUN = function(n) c(H[[n]] %*%
                    rowSums(sapply(X = 1:K, function(k) r_nk[n,k]*m_k[,k]))))
        }
        # Ensure that \mu is not 0 or 1
        for (my_iter in 1:length(mu)) {
            mu[[my_iter]][which(mu[[my_iter]] > (1 - 1e-15))] <- 1 - 1e-15
            mu[[my_iter]][which(mu[[my_iter]] < 1e-15)] <- 1e-15
        }
        E_z <- lapply(X = 1:N, FUN = function(n) .update_Ez(E_z = E_z[[n]],
                 mu = mu[[n]], y_0 = y_0[[n]], y_1 = y_1[[n]]))
        E_zz <- lapply(X = 1:N, FUN = function(n) 1 + mu[[n]] * E_z[[n]])
        # Update expectations over \ln\pi
        e_log_pi  <- digamma(delta_k) - digamma(sum(delta_k))
        # Compute expectation of E[a]
        E_alpha <- alpha_k / beta_k
        # TODO: Perform model selection using MLE of mixing proportions
        # e_log_pi <- colMeans(r_nk)

        ##-------------------------------
        # Variational lower bound
        ##-------------------------------
        mk_Sk    <- lapply(X = 1:K, function(k) tcrossprod(m_k[, k]) + S_k[[k]])
        lb_pz_qz <- sum(sapply(1:N, function(n) 0.5*crossprod(mu[[n]]) +
                       sum(y[[n]] * log(1 - pnorm(-mu[[n]])) + (1 - y[[n]]) *
                       log(pnorm(-mu[[n]]))) - 0.5*sum(sapply(1:K, function(k)
                       r_nk[n,k] * matrix.trace(HH[[n]] %*% mk_Sk[[k]]) )) ))
        lb_p_w   <- sum(-0.5*D*log(2*pi) + 0.5*D*(digamma(alpha_k) -
                       log(beta_k)) - 0.5*E_alpha*E_ww)
        lb_p_c   <- sum(r_nk %*% e_log_pi)
        lb_p_pi  <- sum((delta_0 - 1)*e_log_pi) + lgamma(sum(delta_0)) -
                       sum(lgamma(delta_0))
        lb_p_tau <- sum(alpha_0*log(beta_0) + (alpha_0 - 1)*(digamma(alpha_k) -
                       log(beta_k)) - beta_0*E_alpha - lgamma(alpha_0))
        lb_q_c   <- sum(r_nk*log_r_nk)
        lb_q_pi  <- sum((delta_k - 1)*e_log_pi) + lgamma(sum(delta_k)) -
            sum(lgamma(delta_k))
        lb_q_w   <- sum(-0.5*log(sapply(X = 1:K,function(k) det(S_k[[k]]))) -
                            0.5*D*(1 + log(2*pi)))
        lb_q_tau <- sum(-lgamma(alpha_k) + (alpha_k - 1)*digamma(alpha_k) +
                            log(beta_k) - alpha_k)
        # Sum all parts to compute lower bound
        LB <- c(LB, lb_pz_qz + lb_p_c + lb_p_pi + lb_p_w + lb_p_tau - lb_q_c -
                   lb_q_pi - lb_q_w - lb_q_tau)

        # Show VB difference
        if (is_verbose) {
            cat("It:\t",i,"\tLB:\t",LB[i],"\tDiff:\t",LB[i] - LB[i - 1],"\n")
            cat("Z: ",lb_pz_qz,"\tC: ",lb_p_c - lb_q_c,
                "\tW: ", lb_p_w - lb_q_w,"\tPi: ", lb_p_pi - lb_q_pi,
                "\tTau: ",lb_p_tau - lb_q_tau,"\n")
        }
        # Check if lower bound decreases
        if (LB[i] < LB[i - 1]) { warning("Warning: Lower bound decreases!\n") }
        # Check for convergence
        if (abs(LB[i] - LB[i - 1]) < epsilon_conv) { break }
        # Check if VB converged in the given maximum iterations
        if (i == vb_max_iter) {warning("VB did not converge!\n")}
        if (is_verbose) { utils::setTxtProgressBar(pb, i) }
    }
    if (is_verbose) { close(pb) }

    # Store the object
    obj <- list(W = m_k, W_Sigma = S_k, r_nk = r_nk, delta = delta_k,
                alpha = alpha_k, beta = beta_k, basis = basis, pi_k = pi_k,
                lb = LB)
    if (identical(model, "binomial")) {
        class(obj) <- "cluster_profiles_vb_binomial"
    }else {class(obj) <- "cluster_profiles_vb_bernoulli" }
    return(obj)
}


##------------------------------------------------------


# Cluster profiles EM Binomial or Bernoulli
.cluster_profiles_vb_gaussian <- function(H, y, K, basis, w, delta_0,
                                          gaussian_l, alpha_0, beta_0,
                                          vb_max_iter, epsilon_conv,
                                          is_verbose){
    assertthat::assert_that(is.list(H)) # Check that H is a list object
    D <- basis$M + 1        # Number of covariates
    N <- length(H)          # Extract number of observations
    assertthat::assert_that(N > 0)
    LB <- c(-Inf)           # Store the lower bounds
    r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = N, ncol = K)
    E_ww <- vector("numeric", length = K)

    # Compute H_{n}'H_{n}
    HH <- lapply(X = H, FUN = function(h) crossprod(h))
    # Compute y_{n}'y_{n}
    yy <- unlist(lapply(X = y, FUN = function(y) c(crossprod(y))))
    # Compute H_{n}'y_{n}
    Hy <- lapply(X = 1:N, FUN = function(i) crossprod(H[[i]], y[[i]]))
    # Extract total observations
    len_y <- unlist(lapply(X = y, FUN = function(y) length(y)))

    # Compute shape parameter of Gamma
    alpha_k <- rep(alpha_0 + D/2, K)
    # Mean for each cluster
    m_k <- w
    # Covariance of each cluster
    S_k <- lapply(X = 1:K, function(k) solve(diag(2, D)))
    # Scale of precision matrix
    beta_k   <- rep(beta_0, K)
    # Dirichlet parameter
    delta_k  <- delta_0
    # Expectation of log Dirichlet
    e_log_pi <- digamma(delta_k) - digamma(sum(delta_k))
    mk_Sk    <- lapply(X = 1:K, function(k) tcrossprod(m_k[, k]) + S_k[[k]])

    # Show progress bar
    if (is_verbose) { pb <- utils::txtProgressBar(min = 2, max = vb_max_iter,
                                                  style = 3) }
    # Run VB algorithm until convergence
    for (i in 2:vb_max_iter) {
        ##-------------------------------
        # Variational E-Step
        ##-------------------------------
        for (k in 1:K) {
            log_rho_nk[,k] <- e_log_pi[k] + gaussian_l*sapply(1:N, function(n)
                m_k[,k] %*% Hy[[n]] - 0.5*matrix.trace(HH[[n]] %*% mk_Sk[[k]]))
        }

        # Calculate responsibilities using logSumExp for numerical stability
        log_r_nk <- log_rho_nk - apply(log_rho_nk, 1, .log_sum_exp)
        r_nk     <- exp(log_r_nk)

        ##-------------------------------
        # Variational M-Step
        ##-------------------------------
        # Update Dirichlet parameter
        delta_k <- delta_0 + colSums(r_nk)
        # TODO: Compute expected value of mixing proportions
        pi_k <- (delta_0 + colSums(r_nk)) / (K * delta_0 + N)
        for (k in 1:K) {
            # Update covariance for Gaussian
            w_HH <- lapply(X = 1:N, function(n) HH[[n]]*r_nk[n,k])
            S_k[[k]] <- solve(diag(alpha_k[k]/beta_k[k], D) +
                                  gaussian_l * .add_func(w_HH))
            # Update mean for Gaussian
            w_Hy <- lapply(X = 1:N, function(n) Hy[[n]]*r_nk[n,k])
            m_k[,k] <- gaussian_l * S_k[[k]] %*% .add_func(w_Hy)

            # Update \beta_k parameter for Gamma
            E_ww[k] <- crossprod(m_k[,k]) + matrix.trace(S_k[[k]])
            beta_k[k]  <- beta_0 + 0.5*E_ww[k]
        }
        # Update expectations over \ln\pi
        e_log_pi  <- digamma(delta_k) - digamma(sum(delta_k))
        # Compute expectation of E[a]
        E_alpha <- alpha_k / beta_k

        ##-------------------------------
        # Variational lower bound
        ##-------------------------------
        mk_Sk    <- lapply(X = 1:K, function(k) tcrossprod(m_k[, k]) + S_k[[k]])
        lb_p_y   <- -0.5*sum(len_y)*log(2*pi*(1/gaussian_l)) - 0.5*gaussian_l *
            sum(yy) + sum(sapply(1:K, function(k) gaussian_l*(sum(sapply(1:N,
            function(n) r_nk[n,k]*(m_k[,k] %*% Hy[[n]] -
            0.5*matrix.trace(HH[[n]] %*% mk_Sk[[k]])))))))
        lb_p_w   <- sum(-0.5*D*log(2*pi) + 0.5*D*(digamma(alpha_k) -
            log(beta_k)) - 0.5*E_alpha*E_ww)
        lb_p_c   <- sum(r_nk %*% e_log_pi)
        lb_p_pi  <- sum((delta_0 - 1)*e_log_pi) + lgamma(sum(delta_0)) -
            sum(lgamma(delta_0))
        lb_p_tau <- sum(alpha_0*log(beta_0) + (alpha_0 - 1)*(digamma(alpha_k) -
            log(beta_k)) - beta_0*E_alpha - lgamma(alpha_0))
        lb_q_c   <- sum(r_nk*log_r_nk)
        lb_q_pi  <- sum((delta_k - 1)*e_log_pi) + lgamma(sum(delta_k)) -
            sum(lgamma(delta_k))
        lb_q_w   <- sum(-0.5*log(sapply(X = 1:K,function(k) det(S_k[[k]]))) -
            0.5*D*(1 + log(2*pi)))
        lb_q_tau <- sum(-lgamma(alpha_k) + (alpha_k - 1)*digamma(alpha_k) +
            log(beta_k) - alpha_k)
        # Sum all parts to compute lower bound
        LB <- c(LB, lb_p_y + lb_p_c + lb_p_pi + lb_p_w + lb_p_tau - lb_q_c -
            lb_q_pi - lb_q_w - lb_q_tau)

        # Show VB difference
        if (is_verbose) {
            cat("It:\t",i,"\tLB:\t",LB[i],"\tDiff:\t",LB[i] - LB[i - 1],"\n")
            cat("Y: ",lb_p_y,"\tC: ",lb_p_c - lb_q_c,
                "\tW: ", lb_p_w - lb_q_w,"\tPi: ", lb_p_pi - lb_q_pi,
                "\tTau: ",lb_p_tau - lb_q_tau,"\n")
        }
        # Check if lower bound decreases
        if (LB[i] < LB[i - 1]) { warning("Warning: Lower bound decreases!\n") }
        # Check for convergence
        if (abs(LB[i] - LB[i - 1]) < epsilon_conv) { break }
        # Check if VB converged in the given maximum iterations
        if (i == vb_max_iter) {warning("VB did not converge!\n")}
        if (is_verbose) { utils::setTxtProgressBar(pb, i) }
    }
    if (is_verbose) { close(pb) }

    # Store the object
    obj <- structure(list(W = m_k, W_Sigma = S_k, r_nk = r_nk, delta = delta_k,
                          alpha = alpha_k, beta = beta_k, basis = basis,
                          gaussian_l = gaussian_l, pi_k = pi_k, lb = LB),
                     class = "cluster_profiles_vb_gaussian")
    return(obj)
}


# Compute E[z]
.update_Ez <- function(E_z, mu, y_1, y_0){
    E_z[y_1] <- mu[y_1] + dnorm(-mu[y_1]) / (1 - pnorm(-mu[y_1]))
    E_z[y_0] <- mu[y_0] - dnorm(-mu[y_0]) / pnorm(-mu[y_0])
    return(E_z)
}
