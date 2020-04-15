#' Treatment Effects - Regression Adjusment
#'
#' Short Description
#'
#' @param O.formula \code{\link[stats]{formula}} object for the Outcome model to be fitted to the treated and untreated strata separately.
#' @param T.variable the name of the treatment variable as a string.
#' @param data  data frame to be used for fitting
#' @param te The treatment effect to be fitted. Can be \code{"ate"} for the average treatment effect or
#' \code{"atet"} for the average treatment effect in the treated population or
#' \code{"atc"} for the average treatment effect in the control population.
#' @param O.family link function for the Outcome model, can be either \code{"gaussian"} or \code{"binomial"}
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
#' teffects_ra(Y~X,"A",O.family="gaussian",data=df)
#' teffects_ra(Z~X,"A",O.family="binomial",data=df)
#' @export
teffects_ra<- function(O.formula,T.variable,data,te="ate",O.family="gaussian"){
cl = match.call()
df0 = data[data[,T.variable]==0,]
df1 = data[data[,T.variable]==1,]
model0 <-glm(O.formula,data=df0,family=O.family,x=T)
model1 <-glm(O.formula,data=df1,family=O.family,x=T)

if(O.family=='binomial'){
  link<- function(eta) as.vector(plogis(eta))
  dlink<- function(eta) link(eta)*(1-link(eta))
}else if(O.family=='gaussian'){
  link<- function(eta) as.vector(eta)
  dlink<- function(eta) rep.int(1,length(eta))
}

## Using formulas similar to https://boris.unibe.ch/130362/1/jann-2019-influencefunctions.pdf
N = nrow(data)
X = rbind(model0$x,model1$x)
Y = c(model0$y,model1$y)
D = c(rep.int(0,nrow(df0)),rep.int(1,nrow(df1))) #Treatment variable
b0 = model0$coefficients
b1 = model1$coefficients

IF <- function(X,Y,b,R,P,link,dlink,w=rep.int(1,length(Y))){
  N = length(Y)
  eta = X%*%b
  Wp = sum(w*P) #weights of variables used for prediction
  Wr = sum(w*R) #weights of regressors
  PO =link(eta) #potential outcome
  PO.mean = sum(P*w*PO)/Wp #Potential outcome mean
  Q = Wr*solve(crossprod(X,X*w*R*dlink(eta)))
  U =  X*(Y-link(eta))

  Xp = (P*w*dlink(eta))%*%X/Wp

  IF = (U%*%Q%*%t(Xp))*R*w/Wr + P*w*(PO-PO.mean)/Wp
  return(list(PO.mean = PO.mean,
              IF = N*IF))
}

if(te=='atet'){
  P <- D
}else if(te=='atc'){
  P <- (1-D)
}else if(te=='ate'){
  P <- rep.int(1,N)
}else{
  stop(gettextf("te must be 'ate', 'atc', or 'atet' not: '%s'",te))}

IF1 = IF(X,Y,b1,D,P ,link,dlink)
IF0 = IF(X,Y,b0,1-D,P,link,dlink)

a <- structure(list(),class="teffect")
a$Rfamily = O.family
a$PSfamily = 'none'
a$method = 'Regression Adjustment'
a$call = cl
a$te=te
a$TE = IF1$PO.mean - IF0$PO.mean
a$POmeans = list(Y1 =IF1$PO.mean, Y0 = IF0$PO.mean,
                 Y1_sd = sqrt(mean(IF1$IF^2)/N), Y0_sd = sqrt(mean(IF0$IF^2)/N))
a$sd = sqrt(mean((IF1$IF - IF0$IF)^2)/N)
return(a)
}
