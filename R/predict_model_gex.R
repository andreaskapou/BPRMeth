#' Predict gene expression data from methylation profiles
#'
#' \code{predict_model_gex} predicts gene expression data from methylation
#'  profiles, using a regression model. The same method can be used to predict
#'  differential expression from differential methylation profiles.
#'
#' @param model The regression model learned during training
#' @param test The test data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements:
#' \itemize{
#'  \item \code{test_pred} the test predicted values.
#'  \item \code{test_errors} the test error metrics.
#' }
#'
#' @seealso \code{\link{train_model_gex}}
#'
#' @importFrom stats lm predict
#'
#' @export
predict_model_gex <- function(model, test, is_summary = TRUE){
  # Convert to a data.frame
  test <- as.data.frame(test)

  # Make predictions
  test_pred <- predict(object  = model,
                       newdata = test[ , 1:(NCOL(test) - 1), drop = FALSE],
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
