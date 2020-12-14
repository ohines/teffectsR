#' Fit Inverse Probability Weighting
#'
#' Fitting function called by \link[teffectsR]{teffect} when the method is set to \code{'IPW'}.
#' Inverse Probability Weighting fits a (binomial \link[stats]{glm}) model for the exposure (propensity score model).
#' Predictions from this model are used to weight the observations and derive potential outcome means.
#' This implementation normalizes the propesnity score weights.
#' This is sometimes called the "sample bounded" IPW estimator.
#' 
#' @param X numeric vector of (0,1) specifiying the treatment variable.
#' @param Y numeric vector sepcifying the outcome variable
#' @param Zx design matrix for the exposure (propensity score) model.
#' @param weights an optional numeric vector of ‘observation weights’ to be used in the fitting process.
#' 
#' @return List of fit parameters, which is used to derive an object of class \code{teffects} when called by \link[teffectsR]{teffects}.
#' @examples
#' #generate some data
#' N = 50
#' X = rnorm(N)              #confounder
#' A = rbinom(N,1,plogis(X)) #treatment variable
#' Y = X+0.5*A               #continuous outcome
#' Z = rbinom(N,1,plogis(Y)) #binary outcome
#' df = data.frame(X=X,A=A,Y=Y,Z=Z)
#'
#' teffect(A~X,Y~1,data=df,method="IPW")
#' teffect(A~X,Z~1,data=df,method="IPW")
#'
#' @export
fit_ipw <- function(X,Y,Zx,weights=rep(1,N)){
  N = NROW(X)
  if (is.null(weights)){
    weights <- rep.int(1, N)
  }
  wt = as.vector(weights)
  N.wt = sum(wt)
  
  X.mod = glm.fit(Zx,X,family=quasibinomial(),weights=wt)
  IFa  = IF.glm(X.mod,Zx) #Influence function of ps model
  
  ps   = X.mod$fitted.values
  X_ps = X/ps
  cX_ps = (1-X)/(1-ps)
  dX_ps = -X_ps*(1-ps) #derivative of weights
  dcX_ps = cX_ps*ps
  N.wt1 = sum(wt*X_ps)
  N.wt0 = sum(wt*cX_ps)
  
  Po1.mean = sum(Y*X_ps*wt)/N.wt1
  Po0.mean = sum(Y*cX_ps*wt)/N.wt0
  
  xbar1a = ((Y-Po1.mean)*dX_ps*wt)%*%Zx/N.wt1
  xbar0a = ((Y-Po0.mean)*dcX_ps*wt)%*%Zx/N.wt0
  
  IF1 = (Y - Po1.mean)*wt*X_ps*N/N.wt1 + IFa%*%t(xbar1a) 
  IF0 = (Y - Po0.mean)*wt*cX_ps*N/N.wt0 + IFa%*%t(xbar0a)  
  
  TE = Po1.mean-Po0.mean
  Var = sum((IF1-IF0)^2)/N^2
  
  a = list()
  a$coefs <- c(ATE=TE)
  a$std.err <- c(ATE=sqrt(Var))
  a$Wald <- c(ATE= TE^2/Var)
  a$PO.means <- c(Po1 = Po1.mean, Po0 = Po0.mean)
  a$PO.std.err<- c(Po1 = sqrt(sum(IF1^2))/N, Po0 = sqrt(sum(IF0^2))/N )
  return(a)
}
