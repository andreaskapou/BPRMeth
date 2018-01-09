#' Generating Bernoulli probit regression data for a given region
generate_bernoulli_data <- function(basis, w, max_cov = 30, xmin = -500,
                                    xmax = 500){
    # L is the number of observations found in the ith region
    L <- rbinom(n = 1, size = max_cov, prob = .8)
    x <- matrix(0, nrow = L, ncol = 2)
    # Randomly sample locations for the observations
    obs <- sort(sample(xmin:xmax, L))
    # Scale to (-1,1)
    x[, 1] <- BPRMeth:::.minmax_scaling(data = obs, xmin = xmin, xmax = xmax,
                                        fmin = -1, fmax = 1)
    H      <- BPRMeth::design_matrix(basis, x[, 1])$H
    p_suc  <- pnorm(H %*% w) + rnorm(NROW(H), mean = 0, sd = 0.05)
    p_suc[which(p_suc > (1 - 1e-10))] <- 1 - 1e-10
    p_suc[which(p_suc < 1e-10)] <- 1e-10
    x[, 2] <- rbinom(NROW(H), 1, p_suc)
    return(x)
}


#' Generating Bimomial probit regression data for a given region
generate_binomial_data <- function(basis, w, max_cov = 30, xmin = -500,
                                   xmax = 500){
    # L is the number of observations found in the ith region
    L <- rbinom(n = 1, size = max_cov, prob = .8)
    x <- matrix(0, nrow = L, ncol = 3)
    # Randomly sample locations for the observations
    obs <- sort(sample(xmin:xmax, L))
    # Scale to (-1,1)
    x[, 1] <- BPRMeth:::.minmax_scaling(data = obs, xmin = xmin, xmax = xmax,
                                        fmin = -1, fmax = 1)
    H      <- BPRMeth::design_matrix(basis, x[, 1])$H
    p_suc  <- pnorm(H %*% w) + rnorm(NROW(H), mean = 0, sd = 0.05)
    p_suc[which(p_suc > (1 - 1e-10))] <- 1 - 1e-10
    p_suc[which(p_suc < 1e-10)] <- 1e-10
    x[, 2] <- rbinom(NROW(H), 30, 0.8)
    x[, 3] <- rbinom(NROW(H), x[,2], p_suc)
    return(x)
}


#' Generating Beta probit regression data for a given region
generate_beta_data <- function(basis, w, beta_dispersion = 5, max_cov = 30,
                               xmin = -500, xmax = 500){
    # L is the number of observations found in the ith region
    L <- rbinom(n = 1, size = max_cov, prob = .8)
    x <- matrix(0, nrow = L, ncol = 2)
    # Randomly sample locations for the observations
    obs <- sort(sample(xmin:xmax, L))
    # Scale to (-1,1)
    x[, 1] <- BPRMeth:::.minmax_scaling(data = obs, xmin = xmin, xmax = xmax,
                                        fmin = -1, fmax = 1)
    H      <- BPRMeth::design_matrix(basis, x[, 1])$H
    p_suc  <- pnorm(H %*% w) + rnorm(NROW(H), mean = 0, sd = 0.05)
    p_suc[which(p_suc > (1 - 1e-10))] <- 1 - 1e-10
    p_suc[which(p_suc < 1e-10)] <- 1e-10

    shape1 <- p_suc * beta_dispersion
    shape2 <- beta_dispersion * (1 - p_suc)
    x[, 2] <- rbeta(NROW(H), shape1 = shape1, shape2 = shape2)
    return(x)
}


#' Generating Gaussian regression data for a given region
generate_gaussian_data <- function(basis, w, max_cov = 30, xmin = -500,
                                   xmax = 500){
    # L is the number of observations found in the ith region
    L <- rbinom(n = 1, size = max_cov, prob = .8)
    x <- matrix(0, nrow = L, ncol = 2)
    # Randomly sample locations for the observations
    obs <- sort(sample(xmin:xmax, L))
    # Scale to (-1,1)
    x[, 1] <- BPRMeth:::.minmax_scaling(data = obs, xmin = xmin, xmax = xmax,
                                        fmin = -1, fmax = 1)
    H <- BPRMeth::design_matrix(basis, x[, 1])$H
    # Compute mean of data
    mu <- H %*% w
    # Randomly generate data from this function
    x[, 2] <- rnorm(NROW(H), mu, 0.5)
    return(x)
}

#' Generating probit regression data
generate_data <- function(N=300, type="bernoulli", K=3, pi_k=c(0.45, 0.35, 0.2),
                          w=NULL, basis=BPRMeth::create_rbf_object(M = 3),
                          beta_dispersion = 5, max_cov=50, xmin=-500, xmax=500){
    set.seed(123)
    #require(BPRMeth)
    x <- vector("list", length = N)
    if (is.null(w)) {
        w <- list(w1 = c(-1, -1, 0.9, 3),
                  w2 = c(0.1, -2.4, 3, -2),
                  w3 = c(0.4, 0.7, .7, -2.8))
    }
    # Choose cluster to generate the data
    C_true <- t(rmultinom(N, 1, pi_k))
    C_true <- unlist(apply(C_true, 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
    for (i in 1:N) {
        if (identical(type, "bernoulli")) {
            x[[i]] <- generate_bernoulli_data(basis = basis, w = w[[C_true[i]]],
                                              max_cov = max_cov, xmin = xmin,
                                              xmax = xmax)
        }else if (identical(type, "binomial")) {
            x[[i]] <- generate_binomial_data(basis = basis, w = w[[C_true[i]]],
                                             max_cov = max_cov, xmin = xmin,
                                             xmax = xmax)
        }else if (identical(type, "beta")) {
            x[[i]] <- generate_beta_data(basis = basis, w = w[[C_true[i]]],
                                         beta_dispersion = beta_dispersion,
                                         max_cov = max_cov, xmin = xmin,
                                         xmax = xmax)
        }else if (identical(type, "gaussian")) {
            x[[i]] <- generate_gaussian_data(basis = basis, w = w[[C_true[i]]],
                                             max_cov = max_cov, xmin = xmin,
                                             xmax = xmax)
        }
    }
    return(x)
}


set.seed(1)
bernoulli_data <- generate_data(N = 300, type = "bernoulli")
binomial_data  <- generate_data(N = 300, type = "binomial")
beta_data      <- generate_data(N = 300, type = "beta")
gaussian_data  <- generate_data(N = 300, type = "gaussian",
                                w = list(w1 = c(-1.1, -2, 2.7, -2),
                                         w2 = c(2, 1, -2, 1),
                                         w3 = c(0.4, -0.7, 1.7, 2.8)))
devtools::use_data(bernoulli_data, binomial_data, beta_data, gaussian_data,
                   overwrite = TRUE)
