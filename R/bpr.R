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
# @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#   returned.
#
# @return the log likelihood
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{bpr_gradient}}, \code{\link{design_matrix}}
#
# @importFrom stats pnorm dbinom
#
.bpr_likelihood <- function(w, H, data, is_NLL = FALSE){
    total <- data[, 1]
    succ  <- data[, 2]

    # Predictions of the target variables
    # Compute the cdf of N(0,1) distribution (i.e. probit function)
    Phi <- pnorm(H %*% w)

    # In extreme cases where probit is 0 or 1, subtract a tiny number
    # so we can evaluate the log(0) when computing the Binomial
    Phi[which(Phi > (1 - 1e-289))] <- 1 - 1e-289
    Phi[which(Phi < 1e-289)] <- 1e-289

    # Compute the log likelihood
    res <- sum(dbinom(x = succ, size = total, prob = Phi, log = TRUE)) -
        1 / 2 * t(w) %*% w

    # If we required the Negative Log Likelihood
    if (is_NLL){
        res <- (-1) * res
    }
    return(res)
}


# (Internal) Gradient of the BPR log likelihood function
#
# \code{bpr_gradient} computes the gradient w.r.t the coefficients w of the
# Binomial distributed Probit regression log likelihood function.
#
# @inheritParams bpr_likelihood
#
# @return the gradient vector of the log likelihood w.r.t the vector of
#   coefficients w
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{bpr_likelihood}}, \code{\link{design_matrix}}
#
# @importFrom stats pnorm dnorm
#
.bpr_gradient <- function(w, H, data, is_NLL = FALSE){
    total <- data[, 1]
    succ  <- data[, 2]

    # Predictions of the target variables
    g <- as.vector(H %*% w)
    # Compute the cdf of N(0,1) distribution (i.e. probit function)
    Phi <- pnorm(g)

    # In extreme cases where probit is 0 or 1, subtract a tiny number
    # so we can evaluate the log(0) when computing the Binomial
    Phi[which(Phi > (1 - 1e-289))] <- 1 - 1e-289
    Phi[which(Phi < 1e-289)] <- 1e-289

    # Compute the density of a N(0,1) distribution
    N <- dnorm(g)
    N[which(N < 1e-289)] <- 1e-289

    # Compute the gradient vector w.r.t the coefficients w
    gr <- (N * (succ * (1 / Phi) - (total - succ) * (1 / (1 - Phi)))) %*% H - w

    # If we required the Negative Log Likelihood
    if (is_NLL){
        gr <- (-1) * gr
    }
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
# @param post_prob A vector of length N containing the posterior probabilities
#   fore each element of list x, respectively.
# @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#   returned.
#
# @return The weighted sum of BPR log likelihoods
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{bpr_likelihood}}, \code{\link{bpr_EM}}
#
.sum_weighted_bpr_lik <- function(w, x, des_mat, post_prob, is_NLL = TRUE){
    N <- length(x)

    # TODO: Create tests
    # For each element in x, evaluate the BPR log likelihood
    res <- vapply(X   = 1:N,
                  FUN = function(y) .bpr_likelihood(w = w,
                                                    H = des_mat[[y]]$H,
                                                    data = x[[y]][, 2:3],
                                                    is_NLL = is_NLL),
                  FUN.VALUE = numeric(1),
                  USE.NAMES = FALSE)

    # Return the dot product of the result and the posterior probabilities
    return(post_prob %*% res)
}


# (Internal) Sum of weighted gradients of the BPR log likelihood
#
# \code{sum_weighted_bpr_grad} computes the sum of the gradients of BPR log
# likelihood for each elements of x, and then weights them by the corresponding
# posterior probabilities. This function is mainly used for the M-Step of the EM
# algorithm \code{\link{bpr_EM}}.
#
# @inheritParams sum_weighted_bpr_lik
#
# @return A vector with weighted sum of the gradients of BPR log likelihood.
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{bpr_gradient}}, \code{\link{bpr_EM}}
#
.sum_weighted_bpr_grad <- function(w, x, des_mat, post_prob, is_NLL = TRUE){
    N <- length(x)

    # TODO: Create tests
    # For each element in x, evaluate the gradient of the BPR log likelihood
    res <- vapply(X   = 1:N,
                  FUN = function(y) .bpr_gradient(w = w,
                                                  H = des_mat[[y]]$H,
                                                  data = x[[y]][, 2:3],
                                                  is_NLL = is_NLL),
                  FUN.VALUE = numeric(length(w)),
                  USE.NAMES = FALSE)

    # Return the dot product of the result and the posterior probabilities
    return(post_prob %*% t(res))
}
