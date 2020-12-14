#' Fit Propensity Score Matching
#'
#' Fitting function called by \link[teffectsR]{teffects} when the method is set to \code{'PSMatch'}.
#' Propensity Score Matching fits a (binomial \link[stats]{glm}) model for the exposure (propensity score model).
#' Observations from each treatment group are matched based on predictions from this model.
#' See \code{\link[Matching]{Match}} for details on Matching.
#'
#' @param X numeric vector of (0,1) specifiying the treatment variable.
#' @param Y numeric vector sepcifying the outcome variable
#' @param Zx design matrix for the exposure (propensity score) model.
#' @param Zy design matrix for the outcome model to be fitted to each treatment subgroup.
#' @param Ofam family function for the outcome model. See \link[stats]{family} for details of family functions.
#' @param treatment.effect the treament effect estimand of interest. 
#' Can be either \code{"ATE","ATT"} or \code{"ATC"}.
#' @param weights an optional numeric vector of ‘observation weights’ to be used in the fitting process.
#' @param ... arguments to be passed to the matching function. See \link[Matching]{Match} for details.
#'
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
#' teffect(A~X,Y~1,data=df,method="PSMatch")
#' teffect(A~X,Z~1,data=df,method="PSMatch")
#' @export
fit_psmatch<- function(X,Y,Zx,treatment.effect="ATE",weights=rep(1,N),...){
  
  if(treatment.effect == "All"){
    P.wts = list(ATE="ATE",ATT="ATT",ATC="ATC")
  }else if(treatment.effect %in% c("ATE","ATT","ATC")){
    P.wts =  list(TE=treatment.effect)
    names(P.wts) <- treatment.effect
  }else{
    stop(gettextf("treatment.effect must be 'ATE', 'ATT', 'ATC' or 'All'."))
  }
  
  if (is.null(weights)){
    weights <- rep.int(1, N)
  }
  
  ps = glm.fit(Zx,X,family=quasibinomial(),weights=weights)$fitted.values

  run <- lapply(P.wts, function(TE){
    Mod_match = Matching::Match(Y=Y,Tr=X,X=ps,
                                estimand = TE,
                                weights = weights,
                                ...)
    
    if(Mod_match$ndrops.matches>0){
      warning(gettextf("Observations dropped while estimating '%s' due to caliper = '%i'",TE,
                       Mod_match$ndrops.matches))
    }
    a = list()
    a$coefs <- c(TE=Mod_match$est)
    a$std.err <- c(Mod_match$se)
    a$Wald <- c(TE=(Mod_match$est/Mod_match$se)^2)
    names(a$coefs) <- names(a$std.err) <- names(a$Wald) <- treatment.effect
    
    a$PO.means <- c(Po1 = NA, Po0 = NA)
    a$PO.std.err<- c(Po1 = NA, Po0 = NA )
    (a)
  })
    
  run <- simplify2array(run)
  a <- list(coef = as.double(run["coefs",]),
            std.err = as.double(run["std.err",]),
            Wald = as.double(run["Wald",]),
            Po.means = run["PO.means",],
            Po.std.err = run["PO.std.err",])
  names(a$coef) <- names(a$std.err) <- names(a$Wald) <- names(P.wts)
  return(a)
  
}
