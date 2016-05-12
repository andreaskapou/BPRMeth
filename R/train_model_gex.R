#' Train model for gene expression data from methylation profiles
#'
#' \code{train_model_gex} trains a regression model for gene expression data
#'  from methylation profiles.
#'
#' @param formula An object of class \code{\link[stats]{formula}} needed when
#'  calling the the model, e.g. \code{\link[stats]{lm}} function for performing
#'  linear regression. If NULL, the simple linear regression model is used.
#' @param model_name A charcter denoting the regression model
#' @param train The training data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements:
#' \itemize{
#'  \item \code{formula} the formula that was used.
#'  \item \code{gex_model} the fitted model.
#'  \item \code{train_pred} the training predicted values.
#'  \item \code{train_errors} the training error metrics.
#' }
#'
#' @seealso \code{\link{predict_model_gex}}
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
                        kernel="radial",
                        cross=10)
  }else{
    model <- stats::lm(formula = formula,
                       data = train)
  }

  # Make predictions
  train_pred <- predict(object = model,
                        type   = "response")

  if (length(train_pred) != length(train$y)){
    warning("The regression model returned NAs")
    train_errors = NULL
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
