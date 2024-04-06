#' Treatment Effects for Causal Inference
#'
#' \code{teffects} is used to fit treatment effect estimands of a binary
#' exposure, on an outcome, given a set of variables which adjust for
#' confounding. Estimands included are: the Average Treatment Effect
#' (\code{"ATE"}), the Average Treatment Effect in the Treated
#' (\code{"ATT"}), and the Average Treatment Effect in the Controls
#' (\code{"ATC"}).
#'
#' Several methods are implemented to adjust for confounding.
#' Methods based on regression adjustment
#' (\code{"\link[=fit_ra]{RA}","\link[=fit_aipw]{AIPW}"}) require
#' that a glm outcome model is specified (See \code{\link[stats]{glm}}).
#' Methods base on the propensity score
#' (\code{"\link[=fit_ipw]{IPW}","\link[=fit_aipw]{AIPW}","\link[=fit_psmatch]{PSMatch}"})
#' require that a binomial glm model is sepcified for the exposure
#' (propensity score model).
#' When models are not required, dummy models must speficied, to parse
#' the appropriate exposure/ outcome variable.
#' E.g. when using the method \code{"RA"}, then the exposure model might be
#' \code{A~1}, where \code{A} is the binary exposure variable of interest.
#' The right hand side of these dummy models is ignored.
#'
#' @param exposure.formula an object of class "\code{\link[stats]{formula}}"
#' (or one that can be coerced to that class) where the left hand side of the
#' formula contains the binary exposure variable of interest.
#' @param outcome.formula an object of class "\code{\link[stats]{formula}}"
#' (or one that can be coerced to that class) where the left hand side of the
#' formula contains the outcome variable of interest.
#' The right hand side should not contain the exposure.
#' @param outcome.family family function for the outcome model, can be can be a
#' character string naming a family function, a family function or the result of
#' a call to a family function. (See \link[stats]{family} for details of family
#' functions.)
#' For methods \code{"IPW"} and \code{"PSMatching"} no outcome model is
#' required.
#' @param treatment.effect The treatment effect estimand. Can be either
#' \code{"ATE","ATT","ATC"} or \code{"All"} to fit all three.
#' The method, \code{"PSMatch"}, can only fit the ATE estimand.
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link{as.data.frame}}  to a data frame)
#' containing the variables in the model. If not found in data, the
#' variables are taken from environment(formula), typically the environment from
#' which \code{teffects} is called.
#' @param weights an optional vector of ‘observation weights’ to be used in
#' the fitting process. Should be \code{NULL} or a numeric vector.
#' @param method The mediation fitting method to be used. Can be either
#' \code{"AIPW","IPW","RA"} or \code{"PSMatch"}
#' @param ... arguments to be passed to the matching function when using the
#' \code{"PSMatch"} method.
#' See \link[Matching]{Match} for details.
#'
#' @return An object of class \code{teffect} with effect estimates,
#' estimated standard errors, Wald based test statistics and estimated
#' potential outcome means.
#' @examples
#' # Example on Generated data
#' N <- 50
#' X <- rnorm(N) # confounder
#' A <- rbinom(N, 1, plogis(X)) # treatment variable
#' Y <- X + 0.5 * A # continuous outcome
#' Z <- rbinom(N, 1, plogis(Y)) # binary outcome
#' df <- data.frame(X = X, A = A, Y = Y, Z = Z)
#'
#' # Fit AIPW by default
#' teffect(A ~ X, Y ~ X, data = df)
#' teffect(A ~ X, Z ~ X, data = df, outcome.family = "binomial")
#'
#' # Fit IPW
#' teffect(A ~ X, Y ~ 1, data = df, method = "IPW")
#' teffect(A ~ X, Z ~ 1, data = df, method = "IPW")
#'
#' # Fit RA
#' teffect(A ~ 1, Y ~ X, data = df, method = "RA")
#' teffect(A ~ 1, Z ~ X, data = df, outcome.family = "binomial", method = "RA")
#'
#' # Fit PSMatch
#' teffect(A ~ X, Y ~ 1, data = df, method = "PSMatch")
#' teffect(A ~ X, Z ~ 1, data = df, method = "PSMatch")
#' @export
teffect <- function(exposure.formula, outcome.formula,
                    outcome.family = "gaussian",
                    treatment.effect = "ATE",
                    data, weights, method = "AIPW", ...) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c("exposure.formula", "outcome.formula", "data", "weights"), names(mf), 0L
  )

  # Get exposure model frame if required
  xf <- mf[c(1, m[c(1, 3, 4)])]
  xf[[1L]] <- quote(stats::model.frame)
  xf$drop.unused.levels <- TRUE
  names(xf)[2] <- "formula"
  xf <- eval(xf, parent.frame())
  X <- stats::model.response(xf, "any")
  weights <- as.vector(stats::model.weights(xf))
  xt <- attr(xf, "terms") # allow model.frame to have updated it
  Zx <- stats::model.matrix(xt, xf)

  yf <- mf[c(1, m[c(2, 3, 4)])]
  yf[[1L]] <- quote(stats::model.frame)
  yf$drop.unused.levels <- TRUE
  names(yf)[2] <- "formula"
  yf <- eval(yf, parent.frame())
  Y <- stats::model.response(yf, "any")
  yt <- attr(yf, "terms") # allow model.frame to have updated it
  Zy <- stats::model.matrix(yt, yf)

  ## Check exposure and outcome families
  if (method %in% c("IPW", "AIPW", "PSMatch")) {
    xfam <- "binomial"
  } else {
    xfam <- "none"
  }
  if (method %in% c("RA", "AIPW")) {
    ## family
    if (is.character(outcome.family)) {
      outcome.family <- get(
        outcome.family,
        mode = "function", envir = parent.frame()
      )
    }
    if (is.function(outcome.family)) Ofam <- outcome.family()
    ofam <- Ofam$family
    if (is.null(ofam)) {
      print(family)
      stop("'family' not recognized")
    }
  } else {
    ofam <- "none"
  }

  ## Check Weights
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(weights) && any(weights < 0)) {
    stop("negative weights not allowed")
  }

  ## exposure family
  if (method == "AIPW") {
    a <- teffectsR::fit_aipw(
      X, Y, Zx, Zy, Ofam,
      treatment.effect = treatment.effect, weights = weights
    )
  } else if (method == "IPW") {
    if (treatment.effect != "ATE") {
      warning("Method: 'IPW' only compatible with ATE. Estimating ATE. ")
    }
    a <- teffectsR::fit_ipw(X, Y, Zx, weights = weights)
  } else if (method == "RA") {
    a <- teffectsR::fit_ra(
      X, Y, Zy, Ofam,
      treatment.effect = treatment.effect, weights = weights
    )
  } else if (method == "PSMatch") {
    a <- teffectsR::fit_psmatch(X, Y, Zx,
      treatment.effect = treatment.effect,
      weights = weights, ...
    )
  } else {
    stop("Method not recognized")
  }


  a$Method <- method
  a$call <- mf
  a$exposure.family <- xfam
  a$outcome.family <- ofam
  class(a) <- "teffect"
  a
}


#' @export
print.teffect <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call),
      sep = "\n",
      collapse = "\n"
    ), "\n\n",
    sep = ""
  )
  cat(
    "Fitting Method: \t", x$Method,
    "\nExposure GLM family:\t", x$exposure.family,
    "\nOutcome  GLM family:\t", x$outcome.family
  )

  cat("\n\nTreatment Effects:\n")

  df <- data.frame(
    Estimate = x$coef, Std.Error = x$std.err,
    Wald.value = x$Wald,
    Wald.pval = stats::pchisq(x$Wald, df = 1, lower.tail = FALSE)
  )

  stats::printCoefmat(df,
    digits = 6, signif.stars = TRUE, na.print = "NA",
    tst.ind = 3, P.values = TRUE, has.Pvalue = TRUE
  )
}
