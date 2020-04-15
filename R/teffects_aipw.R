#' Treatment Effects - Augmented Inverse Probability Weighting
#'
#' Short Description
#'
#' @param PS.formula \code{\link[stats]{formula}} object for the Propensity Score logistic regression model.
#' must have the binary treatment variable as a response variable.
#' @param O.formula \code{\link[stats]{formula}} object for the Outcome model to be fitted to the treated and untreated strata separately.
#' @param O.family link function for the Outcome model, can be either \code{"gaussian"} or \code{"binomial"}
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
#' teffects_aipw(A~X,Y~X,"gaussian",df)
#' teffects_aipw(A~X,Z~X,"binomial",df)
#' @export
teffects_aipw <- function(PS.formula,O.formula,O.family,data){
cl = match.call()

PSmodel = glm(PS.formula,data=data,family='binomial',y=T,x=T)
ps = predict(PSmodel,type="response")
D = PSmodel$y
N = length(D)

Ymod1 = glm(O.formula,data=data[D==1,],family=O.family)
Ymod0 = glm(O.formula,data=data[D==0,],family=O.family)

Y = model.frame(O.formula,data=data)[,1]
X = model.matrix(O.formula,data=data)

eta1 = predict(Ymod1,data) #the linear predictor
eta0 = predict(Ymod0,data)

if(O.family=='binomial'){
  link<- function(eta) as.vector(plogis(eta))
  dlink<- function(eta) link(eta)*(1-link(eta))
}else if(O.family=='gaussian'){
  link<- function(eta) as.vector(eta)
  dlink<- function(eta) rep.int(1,length(eta))
}

mu1 = link(eta1)
mu0 = link(eta0)

w1 = D/ps
w0 = (1-D)/(1-ps)
dw1 = -w1*(1-ps) #derivative of weights
dw0 = w0*ps

#Get potential outcomes
Po1 = w1*Y + (1-w1)*mu1
Po0 = w0*Y + (1-w0)*mu0
Po1.mean = mean(Po1)
Po0.mean = mean(Po0)

xbar1a = ((Y-mu1)*dw1)%*%PSmodel$x/N
xbar0a = ((Y-mu0)*dw0)%*%PSmodel$x/N
xbar1b = ((1-w1)*dlink(eta1))%*%X/N
xbar0b = ((1-w0)*dlink(eta0))%*%X/N

IFa = N*(PSmodel$x*(D-ps))%*%chol2inv(PSmodel$qr$qr) #Influence function of ps model
IFb1 = sum(D)*(X*D*(Y-mu1))%*%chol2inv(Ymod1$qr$qr) #Influence function of outcome model
IFb0 = sum(1-D)*(X*(1-D)*(Y-mu0))%*%chol2inv(Ymod0$qr$qr) #Influence function of outcome model

IF1 = IFa%*%t(xbar1a) + IFb1%*%t(xbar1b) + Po1 - Po1.mean
IF0 = IFa%*%t(xbar0a) + IFb0%*%t(xbar0b) + Po0 - Po0.mean

a <- structure(list(),class="teffect")
a$Rfamily = O.family
a$PSfamily = 'binomial'
a$method = 'AIPW'
a$call = cl
a$te="ate"
a$TE = Po1.mean-Po0.mean
a$POmeans = list(Y1 =Po1.mean, Y0 = Po0.mean,
                 Y1_sd = sqrt(mean(IF1^2)/N), Y0_sd = sqrt(mean(IF0^2)/N))
a$sd = sqrt(mean((IF1 - IF0)^2)/N)
return(a)

}
