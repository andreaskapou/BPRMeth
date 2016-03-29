#' Fitting linear models using Basis Functions
#'
#' \code{blm} is used to fit linear models using basis functions such as
#' Radial Basis Functions (RBFs) and Polynomial Basis Functions.
#'
#' @param x The observations.
#' @param y The response.
#' @param basis Basis function object e.g. \code{\link{rbf.object}}.
#' @param lambda Optional parameter for performing penalized least squares.
#' @param return.all Optional logical, indicating if all the metrics should be
#'  computed (mainly for efficiency).
#'
#' @return An object of class "blm" is a list containing the following
#'  components:
#' \itemize{
#'  \item \code{coefficients} a named vector of coefficients.
#'  \item \code{residuals} the residuals, that is response minus fitted values.
#'  \item \code{fitted.values} the fitted mean values.
#'  \item \code{df.residuals} the residual degrees of freedom.
#'  \item \code{sigma} the standard deviation of the residuals.
#'  \item \code{vcov} The covariance matrix.
#'  \item \code{basis} The basis object used.
#'  \item \code{lambda} The regularization parameter.
#'  \item \code{call} the matched call.
#' }
#'
#' @seealso \code{\link{predict.blm}}, \code{\link{polynomial.object}},
#'  \code{\link{rbf.object}}, \code{\link{summary.blm}}
#'
#' @export
blm <- function(x, y, basis, lambda = 0, return.all = TRUE){
  est <- list(basis = basis, lambda = lambda)
  # Create the design matrix
  des_mat <- design_matrix(x = basis, obs = x)
  H <- des_mat$H
  if (lambda == 0){
    # Compute QR decomposition of H
    qx <- qr(H)
    #  Compute (H'H)^(-1)H'y
    est$coefficients <- solve.qr(qx, y)
  }else{
    I <- diag(1, NCOL(H))  # Identity matrix
    I[1,1]  <- 1e-10  # Do not change the intercept coefficient
    # Compute QR decomposition
    qx <- qr(lambda * I + t(H) %*% H)
    # Compute (lambda * I + H'H)^(-1)H'y
    est$coefficients <- solve.qr(qx, t(H) %*% y)
  }

  # Standard deviation of residuals
  est$sigma <- sqrt(sum(est$residuals ^ 2) / est$df.residuals)
  if (return.all){
    # Fitted values
    est$fitted.values <- as.vector(H %*% est$coefficients)
    # Residuals
    est$residuals <- y - est$fitted.values
    # Degrees of freedom
    est$df.residuals <- NROW(H) - NCOL(H)
    # Compute sigma^2 * (x'x)^(-1)
    est$vcov <- est$sigma ^ 2 * chol2inv(qx$qr)
    colnames(est$vcov) <- rownames(est$vcov) <- colnames(H)
  }
  est$call <- match.call()
  class(est) <- "blm"

  return(est)
}


#' Make predictions using the Basis Linear Model
#'
#' S3 method for class blm, that makes predictions for \code{newdata} using
#' the Basis Linear Model. This funciton is similar to \code{\link[stats]{lm}},
#' however currently it does not provide the \code{\link[stats]{formula}}
#' functionality.
#'
#' @param object Object of class \code{\link{blm}}.
#' @param newdata An optional data frame in which to look for variables with
#'  which to predict. If omitted, the fitted values are used.
#' @param ... Optional additional parameters.
#'
#' @return \code{predict.blm} produces a vector of predictions.
#'
#' @export
predict.blm <- function(object, newdata = NULL, ...){
  if(is.null(newdata)){
    y <- stats::fitted(object)
  }
  else{
    # Create the design matrix
    H <- design_matrix(x = object$basis, obs = newdata)$H
    y <- as.vector(H %*% stats::coef(object))
  }
  return(y)
}


#' Print the output of the Basis Linear Model
#'
#' S3 method for class blm, that prints the output of the Basis Linear Model
#' in a similar way to \code{\link[stats]{lm}}.
#'
#' @param x Object of class \code{\link{blm}}.
#' @param ... Optional additional parameters.
#'
#' @export
print.blm <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}


#' Summary output of the Basis Linear Model
#'
#' S3 method for class blm, that computes the summary of the Basis Linear
#' Model in a similar way to \code{\link[stats]{lm}}.
#'
#' @param object Object of class \code{\link{blm}}.
#' @param ... Optional additional parameters.
#'
#' @return A summary.blm object
#'
#' @export
summary.blm <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  tval <- stats::coef(object) / se

  TAB <- cbind(Estimate = stats::coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2 * stats::pt(-abs(tval), df = object$df))
  colnames(TAB) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  res <- structure(list(call = object$call,
                        coefficients = TAB),
                   class = "summary.blm")
  return(res)
}


#' Print summary output of the Basis Linear Model
#'
#' S3 method for class summary.blm, that prints the summary of the Basis Linear
#' Model in a similar way to \code{\link[stats]{lm}}.
#'
#' @param x Object of class \code{\link{summary.blm}}.
#' @param ... Optional additional parameters.
#'
#' @export
print.summary.blm <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
}
