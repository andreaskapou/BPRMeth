#' Predict gene expression from methylation profiles
#'
#' \code{mpgex_regr} is a function that wraps all the necessary subroutines
#' for performing predictions on gene expressions. Initially, it optimizes the
#' parameters of the basis functions so as to learn the methylation profiles.
#' Then uses the learned parameters / coefficients of the basis functions as
#' input features for performing linear regression in order to predict/regress
#' the corresponding gene expression data.
#'
#' @param formula An object of class \code{\link[stats]{formula}} needed when
#'  calling the \code{\link[stats]{lm}} function for performing linear
#'  regression. If NULL, the simple linear regression method is used.
#' @param x The binomial distributed observations, which has to be a list
#'  where each element is an L x 3 dimensional matrix.
#' @param y Corresponding gene expression data for each element of the list x
#' @param model_name A charcter denoting the regression model.
#' @param w Optional vector of initial parameter / coefficient values.
#' @param basis Optional basis function object, default is
#'  \code{\link{polynomial.object}}
#' @param train_ind Optional vector containing the indices for the
#'  train set.
#' @param train_perc Optional parameter for defining the percentage of the
#'  dataset to be used for training set, the remaining will be the test set.
#' @param fit_feature Additional feature on how well the profile fits the
#'  methylation data.
#' @param cpg_dens_feat Additional feature for the CpG density across the
#'  promoter region.
#' @param opt_method Parameter for defining the method to be used in the
#'  optimization procedure, see \code{\link[stats]{optim}}.
#' @param opt_itnmax Optional parameter for defining the max number of
#'  iterations of the optimization procedure, see \code{\link[stats]{optim}}.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return An mpgex object consisting of the following elements:
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_likelihood}},
#'  \code{\link{polynomial.object}}, \code{\link{rbf.object}},
#'  \code{\link{design_matrix}}
#'
#' @examples
#' obs <- meth_data
#' y   <- gex_data
#' basis <- rbf.object(M = 5)
#' out   <- mpgex_regr(x = obs, y = y, basis = basis, is_parallel = FALSE,
#'                     opt_itnmax = 50)
#'
#' @export
mpgex_regr <- function(formula = NULL, x, y, model_name = "svm", w = NULL,
                       basis = NULL, train_ind = NULL, train_perc = 0.7,
                       fit_feature = NULL, cpg_dens_feat = FALSE,
                       opt_method = "CG", opt_itnmax = 500,
                       is_parallel = TRUE, no_cores = NULL, is_summary = TRUE){

  # Check that x is a list object
  assertthat::assert_that(is.list(x))

  # Learn methylation profiles for each gene promoter region
  message("Learning methylation profiles ...\n")
  out_opt <- bpr_optim(x           = x,
                       w           = w,
                       basis       = basis,
                       fit_feature = fit_feature,
                       cpg_dens_feat = cpg_dens_feat,
                       opt_method  = opt_method,
                       opt_itnmax  = opt_itnmax,
                       is_parallel = is_parallel,
                       no_cores    = no_cores)

  # Create training and test sets
  message("Partitioning to test and train data ...\n")
  dataset <- partition_data(x          = out_opt$W_opt,
                            y          = y,
                            train_ind  = train_ind,
                            train_perc = train_perc)

  # Train regression model from methylation profiles
  message("Training linear regression model ...\n")
  train_model <- train_model_gex(formula    = formula,
                                 model_name = model_name,
                                 train      = dataset$train,
                                 is_summary = is_summary)

  # Predict gene expression from methylation profiles
  message("Making predictions ...\n")
  predictions <- predict_model_gex(model      = train_model$gex_model,
                                   test       = dataset$test,
                                   is_summary = is_summary)
  message("Done!\n\n")

  # Create 'mpgex_regr' object
  obj <- structure(list(formula      = formula,
                        model_name   = model_name,
                        opt_method   = opt_method,
                        opt_itnmax   = opt_itnmax,
                        train_ind    = dataset$train_ind,
                        gex_model    = train_model$gex_model,
                        train_pred   = train_model$train_pred,
                        test_pred    = predictions$test_pred,
                        train_errors = train_model$train_errors,
                        test_errors  = predictions$test_errors,
                        fit_feature  = fit_feature,
                        cpg_dens_feat = cpg_dens_feat,
                        opt_method   = opt_method,
                        opt_itnmax   = opt_itnmax,
                        train        = dataset$train,
                        test         = dataset$test,
                        basis        = out_opt$basis,
                        W_opt        = out_opt$W_opt,
                        Mus          = out_opt$Mus),
                   class = "mpgex_regr")
  return(obj)
}
