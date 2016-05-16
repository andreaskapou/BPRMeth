#' Train gene expression model from methylation profiles
#'
#' \code{train_model_gex} trains a regression model for predicting gene
#' expression levels by taking as input the higher order methylation features
#' extracted from specific genomic regions.
#'
#' @param formula An object of class \code{\link[stats]{formula}}, e.g. see
#'   \code{\link[stats]{lm}} function. If NULL, the simple linear regression
#'   model is used.
#' @param model_name A string denoting the regression model. Currently,
#'   available models are: \code{"svm"}, \code{"randomForest"}, \code{"rlm"},
#'   \code{"mars"} and \code{"lm"}.
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
#' @seealso \code{\link{predict_model_gex}}
#'
#' @examples
#' # Create synthetic data
#' train_data <- data.frame(x = rnorm(20), y=rnorm(20, 1, 4))
#' res <- train_model_gex(formula = y~., train = train_data)
#'
#' # Using a different model
#' res <- train_model_gex(model_name = "randomForest", train = train_data)
#'
#' @importFrom stats formula predict
#'
#' @export
train_model_gex <- function(formula = NULL, model_name = "svm", train,
                            is_summary = TRUE){
  if (is.null(formula)){
    formula <- y ~ .
  }
  if (model_name == "mars"){
    model <- earth::earth(formula = formula,
                          data = train)
  }else if (model_name == "randomForest"){
    model <- randomForest::randomForest(formula = formula,
                                        data = train)
  }else if (model_name == "rlm"){
    model <- MASS::rlm(formula = formula,
                       data = train)
  }else if (model_name == "svm"){
    model <- e1071::svm(formula = formula,
                        data = train,
                        kernel = "radial",
                        cross = 10)
  }else{
    model <- stats::lm(formula = formula,
                       data = train)
  }

  # Make predictions
  train_pred <- predict(object = model,
                        type   = "response")

  if (length(train_pred) != length(train$y)){
    warning("The regression model returned NAs")
    train_errors <- NULL
  }else{
    # Calculate model errors
    if (is_summary) message("-- Train Errors --")
    train_errors <- .calculate_errors(x = train$y,
                                      y = train_pred,
                                      summary = is_summary)
  }

  out <- list(formula      = formula,
              gex_model    = model,
              train_pred   = train_pred,
              train_errors = train_errors)
  return(out)
}



#' Predict gene expression model from methylation profiles
#'
#' \code{predict_model_gex} makes predictions of gene expression levels using a
#' model trained on higher order methylation features extracted from specific
#' genomic regions.
#'
#' @param model The fitted regression model, i.e. the output of
#'   \code{\link{train_model_gex}}.
#' @param test The testing data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements: \itemize{ \item
#'   \code{test_pred}: The predicted values for the test data. \item
#'   \code{test_errors}: The test error metrics. }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{train_model_gex}}
#'
#' @examples
#' # Create synthetic data
#' train_data <- data.frame(x = rnorm(20), y=rnorm(20, 1, 4))
#' test_data <- data.frame(x = rnorm(20), y=rnorm(20, 1, 3))
#'
#' # Train the model
#' train_model <- train_model_gex(formula = y~., train = train_data)
#'
#' # Make predictions
#' res <- predict_model_gex(model = train_model$gex_model, test = test_data)
#'
#' @importFrom stats lm predict
#'
#' @export
predict_model_gex <- function(model, test, is_summary = TRUE){
  # Convert to a data.frame
  test <- as.data.frame(test)

  # Make predictions
  test_pred <- predict(object  = model,
                       newdata = test[, 1:(NCOL(test) - 1), drop = FALSE],
                       type    = "response")

  # Calculate model errors
  if (is_summary) message("-- Test Errors --")
  test_errors <- .calculate_errors(x = test$y,
                                   y = test_pred,
                                   summary = is_summary)

  out <- list(test_pred   = test_pred,
              test_errors = test_errors)
  return(out)
}
