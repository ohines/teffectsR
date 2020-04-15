#' Treatment Effects - Inverse Probability Weighting
#'
#' Short Description
#'
#' @param PS.formula \code{\link[stats]{formula}} object for the Propensity Score logistic regression model.
#' must have the binary treatment variable as a response variable.
#' @param O.variable the name of the outcome variable as a string.
#' @param data  data frame to be used for fitting
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
#' teffects_ipw(A~X,"Y",df)
#' teffects_ipw(A~X,"Z",df)
#'
#' @export
teffects_ipw <- function(PS.formula,O.variable,data){
  te="ate"
  cl = match.call()
  PSmodel = glm(PS.formula,data=data,family="binomial",y=T,x=T)
  ps = predict(PSmodel,type='response')
  D = PSmodel$y #Treatment variable
  X = PSmodel$x
  Y = df[,O.variable]
  N = length(D)
  IFb = N*(X*(D-ps))%*%chol2inv(PSmodel$qr$qr) #Influence function of ps model

  w1 = D/ps
  w0 = (1-D)/(1-ps)
  w1 = N*w1/sum(w1)
  w0 = N*w0/sum(w0)
  Po1 = mean(Y*w1) #potential oucomes
  Po0 = mean(Y*w0)

  dw1 = -w1*(1-ps) #derivative of weights
  dw0 = w0*ps

  xbar1 = ((Y-Po1)*dw1)%*%X/N
  xbar0 = ((Y-Po0)*dw0)%*%X/N
  IF1 = IFb%*%t(xbar1) + (Y-Po1)*w1
  IF0 = IFb%*%t(xbar0) + (Y-Po0)*w0

  a <- structure(list(),class="teffect")
  a$Rfamily = 'weighted mean'
  a$PSfamily = 'binomial'
  a$method = 'IPW'
  a$call = cl
  a$te=te
  a$TE = Po1-Po0
  a$POmeans = list(Y1 =Po1, Y0 = Po0,
                   Y1_sd = sqrt(mean(IF1^2)/N), Y0_sd = sqrt(mean(IF0^2)/N))
  a$sd = sqrt(mean((IF1 - IF0)^2)/N)
  return(a)
}
