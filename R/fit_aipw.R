#' Fit Augmented Inverse Probability Weighting
#'
#' Fitting function called by \link[teffectsR]{teffect} when the method is set
#' to \code{'AIPW'}.
#' Augmented Inverse Probability Weighting combines the IPW method and the RA
#' method.
#' A (binomial \link[stats]{glm}) exposure model is fitted which provides
#' observation weights.
#' Also a (\link[stats]{glm}) outcome model is fitted to each treatment
#' subgroup, and used to predict potential outcomes.
#'
#' @param X numeric vector of (0,1) specifiying the treatment variable.
#' @param Y numeric vector sepcifying the outcome variable
#' @param Zx design matrix for the exposure (propensity score) model.
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
#' teffect(A ~ X, Y ~ X, data = df, method = "AIPW")
#' teffect(A ~ X, Z ~ X, data = df, outcome.family = "binomial", method = "AIPW")
#' 
#' @importFrom stats gaussian quasibinomial family
#' @export
fit_aipw <- function(X, Y, Zx, Zy, Ofam = gaussian(), treatment.effect = "ATE", weights = rep(1, N)) {
  link <- Ofam$linkinv # invlink function for outcome model
  dlink <- Ofam$mu.eta # derivative of mean wrt to linear predictor
  N <- NROW(X)
  if (is.null(weights)) {
    weights <- rep.int(1, N)
  }

  X.mod <- stats::glm.fit(Zx, X, family = quasibinomial(), weights = weights)
  Y.mod1 <- stats::glm.fit(Zy, Y, family = Ofam, weights = weights * X)
  Y.mod0 <- stats::glm.fit(Zy, Y, family = Ofam, weights = weights * (1 - X))

  IFa <- teffectsR::IF.glm(X.mod, Zx) # Influence function of ps model
  IFb1 <- teffectsR::IF.glm(Y.mod1, Zy) # Influence function of outcome model
  IFb0 <- teffectsR::IF.glm(Y.mod0, Zy) # Influence function of outcome model

  ps <- X.mod$fitted.values
  eta1 <- Y.mod1$linear.predictor
  eta0 <- Y.mod0$linear.predictor

  mu1 <- link(eta1)
  mu0 <- link(eta0)

  X_ps <- X / ps
  cX_ps <- (1 - X) / (1 - ps)
  dX_ps <- -X_ps * (1 - ps) # derivative of weights
  dcX_ps <- cX_ps * ps

  # Get potential outcomes
  Po1 <- X_ps * Y + (1 - X_ps) * mu1
  Po0 <- cX_ps * Y + (1 - cX_ps) * mu0

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

    xbar1a <- ((Y - mu1) * dX_ps * wt) %*% Zx / N.wt
    xbar0a <- ((Y - mu0) * dcX_ps * wt) %*% Zx / N.wt
    xbar1b <- ((1 - X_ps) * dlink(eta1) * wt) %*% Zy / N.wt
    xbar0b <- ((1 - cX_ps) * dlink(eta0) * wt) %*% Zy / N.wt

    IF1 <- (Po1 - Po1.mean) * wt * N / N.wt + IFa %*% t(xbar1a) + IFb1 %*% t(xbar1b)
    IF0 <- (Po0 - Po0.mean) * wt * N / N.wt + IFa %*% t(xbar0a) + IFb0 %*% t(xbar0b)
    Var <- sum((IF1 - IF0)^2) / N^2
    TE <- Po1.mean - Po0.mean

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
