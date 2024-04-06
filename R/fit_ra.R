#' Fit Regression Adjusment
#'
#' Fitting function called by \link[teffectsR]{teffect} when the method is set
#' to \code{'RA'}.
#' Regression adjusment fits a \link[stats]{glm} model for the outcome to each
#' treatment subgroup.
#' Potential outcomes are derived using these outcome models.
#'
#' @param X numeric vector of (0,1) specifiying the treatment variable.
#' @param Y numeric vector sepcifying the outcome variable
#' @param Zy design matrix for the outcome model to be fitted to each treatment
#' subgroup.
#' @param Ofam family function for the outcome model.
#' See \link[stats]{family} for details of family functions.
#' @param treatment.effect the treament effect estimand of interest.
#' Can be either \code{"ATE","ATT","ATC"} or \code{"All"} to fit all three.
#' @param weights an optional numeric vector of ‘observation weights’ to be used
#' in the fitting process.
#'
#' @return List of fit parameters, which is used to derive an object of class
#' \code{teffects} when called by \link[teffectsR]{teffect}.
#' @examples
#' # generate some data
#' N <- 50
#' X <- rnorm(N) # confounder
#' A <- rbinom(N, 1, plogis(X)) # treatment variable
#' Y <- X + 0.5 * A # continuous outcome
#' Z <- rbinom(N, 1, plogis(Y)) # binary outcome
#' df <- data.frame(X = X, A = A, Y = Y, Z = Z)
#'
#' teffect(A ~ 1, Y ~ X, data = df, method = "RA")
#' teffect(A ~ 1, Z ~ X, data = df, outcome.family = "binomial", method = "RA")
#' @importFrom stats gaussian family
#' @export
fit_ra <- function(X, Y, Zy, Ofam = gaussian(), treatment.effect = "ATE", weights = rep(1, N)) {
  ## Using math formulas similar to https://boris.unibe.ch/130362/1/jann-2019-influencefunctions.pdf

  link <- Ofam$linkinv # invlink function for outcome model
  dlink <- Ofam$mu.eta # derivative of mean wrt to linear predictor
  N <- NROW(X)
  if (is.null(weights)) {
    weights <- rep.int(1, N)
  }

  Y.mod1 <- stats::glm.fit(Zy, Y, family = Ofam, weights = weights * X)
  Y.mod0 <- stats::glm.fit(Zy, Y, family = Ofam, weights = weights * (1 - X))

  IFb1 <- teffectsR::IF.glm(Y.mod1, Zy) # Influence function of outcome model
  IFb0 <- teffectsR::IF.glm(Y.mod0, Zy) # Influence function of outcome model

  eta1 <- Y.mod1$linear.predictor
  eta0 <- Y.mod0$linear.predictor

  Po1 <- link(eta1)
  Po0 <- link(eta0)

  P.wts <- list(ATE = rep.int(1, N), ATT = X, ATC = 1 - X)

  if (!treatment.effect %in% c("ATE", "ATT", "ATC", "All")) {
    stop("'treatment.effect' not recognized")
  }

  if (treatment.effect != "All") P.wts <- P.wts[treatment.effect]

  run <- lapply(P.wts, function(wt) {
    wt <- as.vector(weights) * wt
    N.wt <- sum(wt)

    Po1.mean <- sum(Po1 * wt) / N.wt
    Po0.mean <- sum(Po0 * wt) / N.wt

    xbar1 <- (dlink(eta1) * wt) %*% Zy / N.wt
    xbar0 <- (dlink(eta0) * wt) %*% Zy / N.wt

    IF1 <- (Po1 - Po1.mean) * wt * N / N.wt + IFb1 %*% t(xbar1)
    IF0 <- (Po0 - Po0.mean) * wt * N / N.wt + IFb0 %*% t(xbar0)

    TE <- Po1.mean - Po0.mean
    Var <- sum((IF1 - IF0)^2) / N^2

    a <- list()
    a$coefs <- TE
    a$std.err <- sqrt(Var)
    a$Wald <- TE^2 / Var
    a$PO.means <- c(Po1 = Po1.mean, Po0 = Po0.mean)
    a$PO.std.err <- c(Po1 = sqrt(sum(IF1^2)) / N, Po0 = sqrt(sum(IF0^2)) / N)
    a
  })

  run <- simplify2array(run)
  a <- list(
    coef = as.double(run["coefs", ]),
    std.err = as.double(run["std.err", ]),
    Wald = as.double(run["Wald", ]),
    Po.means = run["PO.means", ],
    Po.std.err = run["PO.std.err", ]
  )
  names(a$coef) <- names(a$std.err) <- names(a$Wald) <- names(P.wts)
  a
}
