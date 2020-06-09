#' @title Predict gene expression from methylation profiles
#'
#' @description \code{bpr_predict_wrap} is a function that wraps all the
#'   necessary subroutines for performing prediction on gene expression levels.
#'   Initially, it optimizes the parameters of the basis functions so as to
#'   learn the methylation profiles. Then, uses the learned parameters /
#'   coefficients of the basis functions as input features for performing
#'   regression in order to predict the corresponding gene expression levels.
#'
#' @param formula An object of class \code{\link[stats]{formula}}, e.g. see
#'   \code{\link[stats]{lm}} function. If NULL, the simple linear regression
#'   model is used.
#' @param prof_obj Inferred profiles object. This in general will be the output
#'   of 'infer_profiles_'(inference_meth) function.
#' @param expr Gene expression data with two columns in \code{data.frame} or
#'   \code{data.table} format. 1st column will have gene IDs and should have
#'   column name "id", 2nd column will have expression levels.
#' @param anno Annotation data as a \code{GRanges} object.
#' @param model_name A string denoting the regression model. Currently,
#'   available models are: \code{"svm"}, \code{"randomForest"}, \code{"rlm"},
#'   \code{"mars"}, \code{"gp"}, and \code{"lm"}.
#' @param fit_feature Use additional feature of how well the profile fits the
#'   methylation data. Either NULL for ignoring this feature or one of the
#'   following: 1) "RMSE" or 2) "NLL" which will be used as input features for
#'   predicting expression.
#' @param cov_feature Logical, whether to use coverage as input feature for
#'   predictions.
#' @param train_ind Optional vector containing the indices for the train set.
#' @param train_perc Optional parameter for defining the percentage of the
#'   dataset to be used for training set, the remaining will be the test set.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A 'predict_expr' object which consists of the following variables:
#'   \itemize{ \item{train}: The training data. \item{test}: The test data.
#'   \item \code{model}: The fitted regression model. \item \code{train_pred}
#'   The predicted values for the training data. \item \code{test_pred} The
#'   predicted values for the test data. \item \code{train_errors}: The training
#'   error metrics. \item \code{test_errors}: The test error metrics.}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{infer_profiles_mle}}, \code{\link{infer_profiles_vb}},
#'   \code{\link{infer_profiles_gibbs}}, \code{\link{create_basis}},
#'
#' @examples
#' # Fit methylation profiles using 5 RBFs
#' basis <- create_rbf_object(M = 5)
#' prof <- infer_profiles_vb(X = encode_met$met, model = "binomial",
#'     basis = basis, is_parallel = FALSE, vb_max_iter = 5)
#' # Predict expression
#' pred_obj <- predict_expr(prof_obj = prof, expr = encode_expr,
#'   anno = encode_met$anno, model_name = "lm", is_summary = FALSE)
#'
#' @export
predict_expr <- function(formula = NULL, prof_obj, expr, anno,model_name = "lm",
                         train_ind = NULL, train_perc = 0.7,
                         fit_feature = "RMSE", cov_feature = TRUE,
                         is_summary = TRUE){
    assertthat::assert_that(methods::is(prof_obj, "infer_profiles"))
    W <- prof_obj$W
    if (prof_obj$basis$M != 0) {
        if (cov_feature) { W <- cbind(W, prof_obj$coverage_feat)}
        if (identical(fit_feature, "RMSE")) { W <- cbind(W, prof_obj$rmse_feat)}
        if (identical(fit_feature, "NLL")) { W <- cbind(W, prof_obj$nll_feat)}
    }
    assertthat::assert_that(methods::is(anno, "GRanges"))
    if (NROW(W) != NROW(anno)) {
        stop("anno object does not match length of prof_obj object.")
    }
    # Append feature ID so we can merge with expression data later
    W <- data.table::data.table(cbind(anno$id, W)); colnames(W)[1] <- "id"
    # Make sure that the first column name of expression matrix is called "id"
    if (!identical(colnames(expr)[1], "id")) {
        stop("First column of expression matrix must have name 'id'.")
    }
    # Merge met profiles with expression data by id
    W <- merge(W, expr, by = "id")
    # Convert to a data.frame
    W <- data.table::setDF(W)
    # Create training and test sets
    dataset <- .partition_data(x = W[, 2:(ncol(W) - 1), drop = FALSE],
                               y = W[, ncol(W), drop = FALSE],
                               train_ind  = train_ind, train_perc = train_perc)
    # Train regression model from methylation profiles
    train <- inner_train_model_expr(formula = formula, model_name = model_name,
                                    train = dataset$train,
                                    is_summary = is_summary)
    # Predict gene expression from methylation profiles
    pred <- inner_predict_model_expr(model = train$model, test = dataset$test,
                                     is_summary = is_summary)
    # Create 'bpr_predict' object
    obj <- structure(list(model = train$model, train_pred = train$train_pred,
                          test_pred = pred$test_pred,
                          train_errors = train$train_errors,
                          test_errors  = pred$test_errors,
                          train = dataset$train, test = dataset$test),
                     class = "predict_expr")
    return(obj)
}


#' @title (INNER) Train expression model from methylation profiles
#'
#' @description This function trains a regression model for
#'   predicting gene expression levels by taking as input the higher order
#'   methylation features extracted from specific genomic regions.
#'
#' @param formula An object of class \code{\link[stats]{formula}}, e.g. see
#'   \code{\link[stats]{lm}} function. If NULL, the simple linear model is used.
#' @param model_name A string denoting the regression model. Currently,
#'   available models are: \code{"svm"}, \code{"randomForest"}, \code{"rlm"},
#'   \code{"mars"}, \code{"gp"}, and \code{"lm"}.
#' @param train The training data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements: \itemize{ \item
#'   \code{formula}: The formula that was used. \item \code{gex_model}: The
#'   fitted model. \item \code{train_pred} The predicted values for the training
#'   data. \item \code{train_errors}: The training error metrics. }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{predict_expr}}
#'
#' @examples
#' # Create synthetic data
#' train_data <- data.frame(x = rnorm(20), y=rnorm(20, 1, 4))
#' res <- inner_train_model_expr(formula = y~., train = train_data)
#'
#' # Using a different model
#' res <- inner_train_model_expr(model_name = "randomForest",
#'                               train = train_data)
#'
#' @importFrom stats formula predict
#'
#' @export
inner_train_model_expr <- function(formula = NULL, model_name = "lm", train,
                       is_summary = TRUE){
    if (is.null(formula)) { formula <- y ~ . }
    if (model_name == "mars") {
        model <- earth::earth(formula = formula, data = train)
    }else if (model_name == "randomForest") {
        model <- randomForest::randomForest(formula = formula, data = train)
    }else if (model_name == "rlm") {
        model <- MASS::rlm(formula = formula, data = train)
    }else if (model_name == "svm") {
        model <- e1071::svm(formula = formula, data = train, kernel = "radial",
                            cross = 10)
    }else if (model_name == "gp") {
        model <- kernlab::gausspr(x = formula, data = train)
        train_pred <- c(kernlab::predict(object = model,
                newdata = train[, 1:(NCOL(train) - 1), drop = FALSE],
                type = "response"))
    }else {model <- stats::lm(formula = formula, data = train) }
    # Make predictions
    if (!identical(model_name, "gp")) {
        train_pred <- c(predict(object = model, type = "response"))
    }
    if (length(train_pred) != length(train$y)) {
        warning("The regression model returned NAs"); train_errors <- NULL
    }else{
        # Calculate model errors
        if (is_summary) message("-- Train Errors --")
        train_errors <- .calculate_errors(x = train$y, y = train_pred,
                                          summary = is_summary)
    }
    return(list(model = model, train_pred = train_pred,
                train_errors = train_errors))
}



#' @title (INNER) Predict expression
#'
#' @description This functions makes predictions of gene expression levels using
#'   a model trained on methylation features extracted from genomic regions.
#'
#' @param model The fitted regression model, i.e. the output of
#'   \code{\link{inner_train_model_expr}}.
#' @param test The testing data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements: \itemize{ \item
#'   \code{test_pred}: The predicted values for the test data. \item
#'   \code{test_errors}: The test error metrics. }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{predict_expr}}
#'
#' @examples
#' # Create synthetic data
#' train_data <- data.frame(x = rnorm(20), y=rnorm(20, 1, 4))
#' test_data <- data.frame(x = rnorm(40), y=rnorm(20, 1, 3))
#'
#' # Train the model
#' train_model <- inner_train_model_expr(formula = y~., model_name="svm",
#'                                       train = train_data)
#'
#' # Make predictions
#' res <- inner_predict_model_expr(model = train_model$model, test = test_data)
#'
#' @importFrom stats lm predict
#'
#' @export
inner_predict_model_expr <- function(model, test, is_summary = TRUE){
    # Convert to a data.frame
    test <- as.data.frame(test)
    # Make predictions
    if (is(model, "gausspr")) {
        test_pred <- c(kernlab::predict(object = model,
         newdata = test[, 1:(NCOL(test) - 1), drop = FALSE], type = "response"))
    }else{
        test_pred <- c(predict(object = model,
         newdata = test[, 1:(NCOL(test) - 1), drop = FALSE], type = "response"))
    }
    # Calculate model errors
    if (is_summary) message("-- Test Errors --")
    test_errors <- .calculate_errors(x = test$y, y = test_pred,
                                     summary = is_summary)
    return(list(test_pred = test_pred, test_errors = test_errors))
}
