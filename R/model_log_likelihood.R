#' @name model_log_likelihood
#' @rdname model_log_likelihood
#' @aliases obs_model_likelihood obs_log_likelihood model_likelihood
#'
#' @title Compute observation model log-likelihood
#'
#' @description These functions evaluate the model log-likelihood and gradient
#'   for different observation models. Available models are "bpr" (i.e.
#'   "bernoulli" or "binomial"), "beta" and "lr" (i.e. "gaussian"). There are
#'   also functions to compute the sum and weighted sum of the observation model
#'   likelihoods, e.g. required for the EM algorithm. These functions are
#'   written in C++ for efficiency.
#'
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param X An \code{L X C} matrix, where L are the total number of
#'   observations. The first column contains the input observations x (i.e. CpG
#'   locations). If "binomial" model then C=3, and 2nd and 3rd columns contain
#'   total number of trials and number of successes respectively. If "bernoulli"
#'   or "gaussian" model, then C=2 containing the output y (e.g. methylation
#'   level). if "beta" model, then C=3, where 2nd column contains output y and
#'   3rd column the dispersion parameter. Each row corresponds to each row of
#'   the design matrix \code{H}.
#' @param H The \code{L x M} matrix design matrix, where L is the number of
#'   observations and M the number of basis functions.
#' @param X_list A list of elements of length N, where each element is an
#'   \code{L x K} matrix of observations X.
#' @param H_list A list of elements of length N, where each element contains the
#'   \code{L x M} design matrices \code{H}.
#' @param r_nk A vector of length N containing the posterior probabilities (i.e.
#'   responsibilities) for each element of \code{X_list}.
#' @param lambda The complexity penalty coefficient for penalized regression.
#' @param is_nll Logical, indicating if the Negative Log Likelihood should be
#'   returned.
#'
#' @return Returns the log likelihood or gradient of the observation model.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{eval_functions}}, \code{\link{infer_profiles_mle}}
NULL

#' @rdname model_log_likelihood
#'
#' @export
lr_log_likelihood <- function(w, X, H, lambda = .5, is_nll = FALSE){
    fit <- H %*% w                       # Predictions
    residuals <- c(X[,2]) - fit          # Residuals
    df <- NROW(H) - NCOL(H)              # Degrees of freedom
    sigma <- sqrt(sum(residuals^2) / df) # Standard deviation of residuals
    # Compute the log likelihood
    res <- sum(dnorm(x = X[,1], mean = fit, sd = sigma, log = TRUE)) -
        lambda * t(w) %*% w
    # If we required the Negative Log Likelihood
    if (is_nll) {res <- (-1) * res }
    return(res)
}

# (Internal) BPR log likelihood function
#
# \code{bpr_likelihood} evaluates the Binomial distributed Probit regression log
# likelihood function for a given set of coefficients, observations and a design
# matrix.
#
# @section Mathematical formula: The Binomial distributed Probit Regression log
#   likelihood function is computed by the following formula: \deqn{log p(y | f,
#   w) = \sum_{l=1}^{L} log Binom(m_{l} | t_{l}, \Phi(w^{t}h(x_{l})))} where
#   h(x_{l}) are the basis functions.
#
# @param w A vector of parameters (i.e. coefficients of the basis functions)
# @param H The \code{L x M} matrix design matrix, where L is the number of
#   observations and M the number of basis functions.
# @param data An \code{L x 2} matrix containing in the 1st column the total
#   number of trials and in the 2nd the number of successes. Each row
#   corresponds to each row of the design matrix.
# @param lambda The complexity penalty coefficient for penalized regression.
# @param is_nll Logical, indicating if the Negative Log Likelihood should be
#   returned.
# @return the log likelihood
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
.old_bpr_likelihood <- function(w, H, data, lambda = 1/2, is_nll = FALSE){
    # If the data is matrix, we have Binomial observations
    if (is.matrix(data)) { total <- data[, 1]; succ  <- data[, 2]
    # ... otherwise we have Bernoulli observations.
    }else {succ <- data; total <- rep(1, length(succ)) }
    # Predictions of the target variables
    # Compute the cdf of N(0,1) distribution (i.e. probit function)
    Phi <- pnorm(H %*% w)
    # In extreme cases where probit is 0 or 1, subtract a tiny number
    # so we can evaluate the log(0) when computing the Binomial
    Phi[which(Phi > (1 - 1e-15))] <- 1 - 1e-15
    Phi[which(Phi < 1e-15)] <- 1e-15
    # Compute the log likelihood
    res <- sum(dbinom(x = succ, size = total, prob = Phi, log = TRUE)) -
        lambda * t(w) %*% w
    #M <- length(w)
    #if (M > 1){ res <- res - lambda * t(w[2:M]) %*% w[2:M] }
    if (is_nll) { res <- (-1)*res } # If we require the Negative Log Likelihood
    return(res)
}


# (Internal) Gradient of the BPR log likelihood function
#
# \code{bpr_gradient} computes the gradient w.r.t the coefficients w of the
# Binomial distributed Probit regression log likelihood function.
#
# @inheritParams bpr_likelihood
# @return the gradient vector of the log likelihood w.r.t the vector of
#   coefficients w
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
.old_bpr_gradient <- function(w, H, data, lambda = 1/2, is_nll = FALSE){
    # If the data is matrix, we have Binomial observations
    if (is.matrix(data)) { total <- data[, 1]; succ  <- data[, 2]
    # ... otherwise we have Bernoulli observations.
    }else {succ <- data; total <- rep(1, length(succ)) }
    g <- as.vector(H %*% w)  # Predictions of the target variables
    Phi <- pnorm(g) # Compute the cdf of N(0,1) (i.e. probit function)
    # In extreme cases where probit is 0 or 1, subtract a tiny number
    # so we can evaluate the log(0) when computing the Binomial
    Phi[which(Phi > (1 - 1e-15))] <- 1 - 1e-15
    Phi[which(Phi < 1e-15)] <- 1e-15
    N <- dnorm(g)  # Compute the density of a N(0,1) distribution
    N[which(N < 1e-15)] <- 1e-15

    # Compute the gradient vector w.r.t the coefficients w
    gr <- (N*(succ*(1/Phi) - (total - succ)*(1/(1 - Phi)))) %*% H - 2*lambda*w
    #M <- length(w)
    #if (M > 1){ gr[2:M] <- gr[2:M] - 2 * lambda * w[2:M] }
    if (is_nll) { gr <- (-1) * gr } # If we require the Negative Log Likelihood
    return(gr)
}


# (Internal) Sum of weighted BPR log likelihoods
#
# \code{sum_weighted_bpr_lik} computes the sum of the BPR log likelihoods for
# each elements of x, and then weights them by the corresponding posterior
# probabilities. This function is mainly used for the M-Step of the EM algorithm
# \code{\link{bpr_EM}}.
#
# @param w A vector of parameters (i.e. coefficients of the basis functions)
# @param x A list of elements of length N, where each element is an L x 3 matrix
#   of observations, where 1st column contains the locations. The 2nd and 3rd
#   columns contain the total trials and number of successes at the
#   corresponding locations, repsectively.
# @param des_mat A list of length N, where each element contains the \code{L x
#   M} design matrices, where L is the number of observations and M the number
#   of basis functions.
# @param r_nk A vector of length N containing the posterior probabilities
#   fore each element of list x, respectively.
# @param lambda The complexity penalty coefficient for penalized regression.
# @param is_nll Logical, indicating if the Negative Log Likelihood should be
#   returned.
#
# @return The weighted sum of BPR log likelihoods
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
.old_sum_weighted_bpr_lik <- function(w,x,des_mat,r_nk,lambda,is_nll=TRUE){
    N <- length(x)
    # For each element in x, evaluate the BPR log likelihood
    res <- vapply(X = 1:N, FUN = function(y) .old_bpr_likelihood(w = w,
      H = des_mat[[y]], data = x[[y]][,2:3], lambda = lambda, is_nll = is_nll),
        FUN.VALUE = numeric(1), USE.NAMES = FALSE)
    # Return the dot product of the result and the posterior probabilities
    return(r_nk %*% res)
}


# (Internal) Sum of weighted gradients of the BPR log likelihood
#
# \code{sum_weighted_bpr_grad} computes the sum of the gradients of BPR log
# likelihood for each elements of x, and then weights them by the corresponding
# posterior probabilities. This function is mainly used for the M-Step of the EM
# algorithm \code{\link{bpr_EM}}.
#
# @inheritParams sum_weighted_bpr_lik
# @return A vector with weighted sum of the gradients of BPR log likelihood.
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
.old_sum_weighted_bpr_grad <- function(w,x,des_mat,r_nk,lambda,is_nll=TRUE){
    N <- length(x)
    # For each element in x, evaluate the gradient of the BPR log likelihood
    res <- vapply(X = 1:N, FUN = function(y) .old_bpr_gradient(w = w,
        H = des_mat[[y]], data = x[[y]][,2:3], lambda = lambda,
          is_nll = is_nll),
          FUN.VALUE = numeric(length(w)), USE.NAMES = FALSE)
    # Return the dot product of the result and the posterior probabilities
    return(r_nk %*% t(res))
}
