#' Influence Functions from Generalized Linear Models
#'
#' \code{IF.glm} returns a matrix of influence functions from either a \code{glm} model object (run with \code{x=TRUE}),
#'  or a \code{glm.fit} object (in which case the user must parse the design matrix, \code{x}).
#'  The resulting matrix has one row for each observation, and one column for each GLM model parameter.
#'
#' @param object A \code{glm} model object (run with \code{x=TRUE}), or a \code{glm.fit} object
#' @param x The model matrix. Not required for a \code{glm} model object (run with \code{x=TRUE}).
#' @return Matrix of influence function values.
#' @export

IF.glm <- function(object,x=object$x){ 
  NROW(object$y)*(x*object$prior.weights*(object$y-object$fitted.values))%*%chol2inv(object$qr$qr)
}