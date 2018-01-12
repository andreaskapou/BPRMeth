# Compute predictive distribution of inferred profiles using VB
.predictive_infer_profile <- function(model, x, region = 1){
    # Predictive mean
    H <- design_matrix(model$basis, x)$H
    W_pred <- c(H %*% model$W[region, ])
    if (methods::is(model, "infer_profiles_vb_binomial") ||
        methods::is(model, "infer_profiles_vb_bernoulli")) {
        # Predictive variance
        W_sd_pred <- sqrt(1 + diag(H %*% model$W_Sigma[[region]] %*% t(H)))
        W_pred <- pnorm(W_pred / W_sd_pred)
    }else if (methods::is(model, "infer_profiles_vb_gaussian")) {
        W_sd_pred <- sqrt(1 / model$gaussian_l +
                              diag(H %*% model$W_Sigma[[region]] %*% t(H)))
    }else {stop("Not predictive distribution for this object") }
    return(list(W_pred = W_pred, W_sd_pred = W_sd_pred))
}


# Compute predictive distribution of clustered profiles using VB
.predictive_cluster_profile <- function(model, x){
    K <- NCOL(model$W) # NUmber of clusters
    H <- design_matrix(model$basis, x)$H # Design matrix
    W_pred = W_sd_pred <- matrix(0, ncol = K, nrow = NROW(H))
    for (k in 1:K) {
        # Predictive mean
        W_pred[,k] <- c(H %*% model$W[,k])
        if (methods::is(model, "cluster_profiles_vb_binomial") ||
            methods::is(model, "cluster_profiles_vb_bernoulli")) {
            # Predictive variance
            W_sd_pred[,k] <- sqrt(1 + diag(H %*% model$W_Sigma[[k]] %*% t(H)))
            W_pred[,k] <- pnorm(W_pred[,k] / W_sd_pred[,k])
        }else if (methods::is(model, "cluster_profiles_vb_gaussian")) {
            W_sd_pred <- sqrt(1 / model$gaussian_l +
                                  diag(H %*% model$W_Sigma[[k]] %*% t(H)))
        }else {stop("Not predictive distribution for this object") }
    }
    return(list(W_pred = W_pred, W_sd_pred = W_sd_pred))
}

# Compute the min-max scaling
#
# \code{.minmax_scaling} normalizes a given vector using the the min-max
# scaling method. More formally:
# \deqn{scaled = \frac{data -x_{min}}{x_{max} - x_{min}} \times (f_{max} -
#  f_{min}) + f_{min}}
#
# @param data Vector with numeric data to be scaled.
# @param xmin Optional minimum value, otherwise \code{min(data)} will be used.
# @param xmax Optional maximum value, otherwise \code{max(data)} will be used.
# @param fmin Optional minimum range value, default is -1.
# @param fmax Optional maximum range value, default is 1.
#
# @return The scaled data in the given range, default is between (-1, 1). If
#  xmin = xmax the input vector \code{data} is returned.
#
.minmax_scaling <- function(data, xmin = NULL, xmax = NULL,
                            fmin = -1, fmax = 1){
    if (is.null(xmin)) { xmin <- min(data) }
    if (is.null(xmax)) { xmax <- max(data) }
    if ( (xmin - xmax) == 0) { return(data) }
    minmax <- (data - xmin) / (xmax - xmin)
    minmax_scaled <- minmax * (fmax - fmin) + fmin
    return(minmax_scaled)
}


# Centre CpG locations
#
# \code{centre_loc} centres CpG locations relative to middle point
#   of feature of interest.
#
# @param region CpG locations
# @param centre Genomic region central point
# @param strand_direction Strand direction
#
# @return Centred location data
#
.do_centre_loc <- function(region, centre, strand_direction){
    assertthat::assert_that(is.character(strand_direction))
    tmp <- region - centre
    # If '-' strand, swap CpG locations
    if (identical(strand_direction, "-")) { tmp <- (-tmp) }
    return(tmp)
}

# Partition data in train and test set
#
# \code{.partition_data} partition data randomly in train and test sets.
#
# @param x Input / Independent variables
# @param y Dependent variables
# @param train_ind Index of traininig data, if NULL a random one is generated.
# @param train_perc Percentage of training data when partitioning.
#
# @return A list containing the train, test data and index of training data.
#
.partition_data <- function(x, y, train_ind = NULL, train_perc = 0.7){
    # Convert both x and y to matrices
    x <- data.matrix(x); colnames(x) <- NULL
    y <- data.matrix(y); colnames(y) <- NULL
    if (is.null(train_ind)) {
        pivot <- NROW(x) * train_perc
        train_ind <- sample(NROW(x), round(pivot))
    }
    train <- data.frame(x = x[train_ind, , drop = FALSE],
                        y = y[train_ind, , drop = FALSE])
    test <- data.frame(x = x[-train_ind, , drop = FALSE],
                       y = y[-train_ind, , drop = FALSE])
    return(list(train = train, test = test, train_ind = train_ind))
}


# @title Number of parallel cores
#
# @description Function for creating the number of parallel cores that will be
#   used during EM.
# @param no_cores Number of cores given as input
# @param is_parallel Logical, did we require parallel computations
# @param M Total number of sources
#
.parallel_cores <- function(no_cores=NULL, is_parallel=FALSE, max_cores = NULL){
    if (is_parallel) { # If parallel mode is ON
        # If number of cores is not given
        if (is.null(no_cores)) { no_cores <- parallel::detectCores() - 1
        } else{if (no_cores >= parallel::detectCores()) {
            no_cores <- parallel::detectCores() - 1 }
        }
        if (is.na(no_cores)) { no_cores <- 2 }
        if (!is.null(max_cores)) {
            if (no_cores > max_cores) { no_cores <- max_cores }
        }
    }
    return(no_cores)
}


# Calculate error metrics
#
# \code{.calculate_errors} calculates error metrics so as to assess the
#  performance of a model, e.g. linear regression.
# @param x Actual values.
# @param y Predicted values.
# @param summary Logical, indicating if the erorr metrics should be printed.
#
# @return A list containing the following components:
# \itemize{
#  \item \code{mae} mean absolute error.
#  \item \code{mse}  mean squared error (the variance).
#  \item \code{rmse} root mean squared error (std dev).
#  \item \code{mape} mean absolute percentage error.
#  \item \code{rstd} relative standard deviation.
# }
#
.calculate_errors <- function(x, y, summary = FALSE){
    # TODO Compute actual errors using the right DoF!!
    R <- list()
    if (!is.numeric(x) || !is.numeric(y))
        stop("Arguments 'x' and 'y' must be numeric vectors.")
    if (length(x) != length(y))
        stop("Vectors 'x' and ' y' must have the same length.")
    error  <- y - x
    # mean absolute error
    R$mae  <- mean(abs(error))
    mae_f  <- formatC(R$mae, digits = 4, format = "f")
    # mean squared error (the variance?!)
    R$mse  <- mean(error ^ 2)
    mse_f  <- formatC(R$mse, digits = 4, format = "f")
    # root mean squared error (std. dev.)
    R$rmse <- sqrt(R$mse)
    rmse_f <- formatC(R$rmse, digits = 4, format = "f")
    # mean absolute percentage error
    R$mape <- mean(abs(error / x))
    # relative standard deviation
    R$rstd <- R$rmse / mean(x)
    # Compute r-squared
    R$rsq  <- 1 - (sum(error ^ 2) / sum( (x - mean(x)) ^ 2) )
    rsq_f  <- formatC(R$rsq, digits = 4, format = "f")
    # Pearson Correlation Coefficient
    R$pcc  <- stats::cor(x, y)
    pcc_f  <- formatC(R$pcc, digits = 4, format = "f")
    if (summary) {
        cat("-- Error Terms ----\n")
        cat(" MAE:  ", mae_f, "\n")
        cat(" MSE:  ", mse_f, "\n")
        cat(" RMSE: ", rmse_f, "\n")
        cat(" R-sq: ", rsq_f, "\n")
        cat(" PCC:  ", pcc_f, "\n")
        cat("-------------------\n\n")
    }
    if (summary) { invisible(R)
    }else {return(R) }
}


# Internal function to make all the appropriate type checks.
.do_checks <- function(w = NULL, basis = NULL){
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }
    if (is.null(w)) { w <- rep(0.5, basis$M + 1) }
    if (length(w) != (basis$M + 1) ) {
        stop("Coefficient vector does not match basis functions!")
    }
    return(list(w = w, basis = basis))
}

# Internal function to make all the appropriate type checks for EM.
.em_checks <- function(X, H, K=3, pi_k=NULL, model, basis=NULL, lambda=.5,
                       beta_dispersion=5, w=NULL, opt_method="CG",
                       init_opt_itnmax=100, is_parallel=FALSE, no_cores=NULL){
    if (is.null(w)) { # If we haven't initialized coefficient matrix
        out <- infer_profiles_mle(X = X, model = model, basis = basis, H = H,
          lambda = lambda, w = w, beta_dispersion = beta_dispersion,
          opt_method = opt_method, init_opt_itnmax = init_opt_itnmax,
          is_parallel = is_parallel, no_cores = no_cores)
        basis <- out$basis # Extract basis object
        W <- out$W # Extract the optimized coefficients
        cl <- stats::kmeans(W, K, nstart = 25) # Use K-means
        C_n <- cl$cluster   # Get the mixture components
        w <- t(cl$centers)  # Mean for each cluster
        # Mixing proportions
        if (is.null(pi_k)) { pi_k <- as.vector(table(C_n) / length(X)) }
    }else{
        if (is.null(basis)) { basis <-  create_rbf_object(M = NROW(w)) }
        if (NROW(w) != (basis$M + 1) ) {
            stop("Coefficients w should match number of basis functions!")
        }
        if (is.null(pi_k)) { pi_k <- rep(1 / K, K) }
    }
    return(list(w = w, basis = basis, pi_k = pi_k))
}

# Compute stable Log-Sum-Exp
#
# \code{.log_sum_exp} computes the log sum exp trick for avoiding numeric
# underflow and have numeric stability in computations of small numbers.
#
# @param x A vector of observations
#
# @return The logs-sum-exp value
#
# @references
#  \url{https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/}
#
.log_sum_exp <- function(x) {
    # Computes log(sum(exp(x))
    offset <- max(x)
    s <- log(sum(exp(x - offset))) + offset
    i <- which(!is.finite(s))
    if (length(i) > 0) { s[i] <- offset }
    return(s)
}

## A general-purpose adder:
.add_func <- function(x) Reduce(f = "+", x = x)

# Extract FPKM from string
#
# \code{.extract_fpkm} Extracts FPKM value from a string
#
# @param x a string containing FPKM information
#
# @return The FPKM numeric value
#
.extract_fpkm <- function(x){
    # TODO test when no FPKM is available
    fpkm <- gsub(".* FPKM ([^;]+);.*", "\\1", x)
    return(as.numeric(fpkm))
}


# Extract gene name from string
#
# \code{.extract_gene_name} Extracts gene name from a string
#
# @param x a string containing gene name information
#
# @return The gene name as a string
#
.extract_gene_name <- function(x){
    # TODO test when no gene name is available
    gene_name <- gsub(".* gene_name ([^;]+);.*", "\\1", x)
    return(gene_name)
}


# Discard selected chromosomes
#
# \code{.discard_chr} Discards selected chromosomes
#
# @param x The HTS data stored in a data.table object
# @param chr_discarded A vector with chromosome names to be discarded.
#
# @return The reduced HTS data.
#
.discard_chr <- function(x, chr_discarded = NULL){
    assertthat::assert_that(methods::is(x, "data.table"))
    if (!is.null(chr_discarded)) {
        message("Removing selected chromosomes ...")
        for (i in 1:length(chr_discarded)) {
            x <- x[x$chr != chr_discarded[i]]
        }
    }
    return(x)
}


# Discard BS-Seq noisy reads
#
# \code{.discard_bs_noise_reads} discards low coverage and (really) high reads
#  from BS-Seq experiments. These reads can be thought as noise of the
#  experiment.
#
# @param bs_data A GRanges object containing the BS-Seq data.
# @param min_bs_cov The minimum number of reads mapping to each CpG site.
# @param max_bs_cov The maximum number of reads mapping to each CpG site.
#
# @return The reduced GRanges object without noisy observations
#
.discard_bs_noise_reads <- function(bs_data, min_bs_cov = 2, max_bs_cov = 1000){
    message("Discarding noisy reads ...")
    bs_data <- subset(bs_data, bs_data$total >= min_bs_cov)
    bs_data <- subset(bs_data, bs_data$total <= max_bs_cov)
    return(bs_data)
}

