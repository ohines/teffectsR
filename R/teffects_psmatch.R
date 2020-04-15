#' Treatment Effects - Propensity Score Matching
#'
#' Uses the \code{\link[Matching]{Match}} function to perform matching.
#'
#' @param PS.formula \code{\link[stats]{formula}} object for the Propensity Score logistic regression model.
#' must have the binary treatment variable as a response variable.
#' @param O.variable the name of the outcome variable as a string.
#' @param data  data frame to be used for fitting
#' @param te The treatment effect to be fitted. Can be \code{"ate"} for the average treatment effect or
#' \code{"atet"} for the average treatment effect in the treated population
#' @param caliper Optional scalar denoting the caliper which is passed to the \code{\link[Matching]{Match}} function.
#' Observations which are outside of the caliper are dropped, with a warning.
#'
#'
#' @return An object of class \code{teffects}
#' @examples
#' #generate some data
#' N = 1000
#' X = rnorm(N)              #confounder
#' A = rbinom(N,1,plogis(X)) #treatment variable
#' Y = X+0.5*A               #continuous outcome
#' Z = rbinom(N,1,plogis(Y)) #binary outcome
#' df = data.frame(X=X,A=A,Y=Y,Z=Z)
#'
#' teffects_psmatch(A~X,"Y",df)
#' teffects_psmatch(A~X,"Z",df)
#' @export
teffects_psmatch<- function(PS.formula,O.variable,data,te="ate",caliper=NULL){
  cl = match.call()
  PSmodel = glm(PS.formula,data=data,family="binomial",y=T,x=T)
  ps = predict(PSmodel,type='response')

  if(te=="ate"){
    te_pass <- "ATE"
  }else if(te=="atet"){
    te_pass <- "ATT"
  }else{
    stop(gettextf("te must be 'ate' or 'atet' not: '%s'",te))
  }

  #Use the Matching package to do all the hard work
  Mod_match = Matching::Match(Y=data[,O.variable]
                              ,Tr=PSmodel$y,X=ps,estimand = te_pass,caliper=caliper)

  if(Mod_match$ndrops.matches>0){
    warning(gettextf("Number of observations dropped due to caliper = '%i'",
                     Mod_match$ndrops.matches))
  }

  a <- structure(list(),class="teffect")
  a$Rfamily = 'matching'
  a$PSfamily = 'binomial'
  a$method = 'PS Match'
  a$call = cl
  a$te=te
  a$TE = Mod_match$est
  a$POmeans = list(Y1 =NA, Y0 = NA,
                   Y1_sd = NA, Y0_sd = NA)
  a$sd = Mod_match$se
  return(a)

}
