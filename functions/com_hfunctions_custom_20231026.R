#### These are functions modeling the decay in detection rate with radial distance in a line-transect survey
#### x = perpendicular distance
#### y = longitudinal forward distance
#### All exctracted from LT2D package

list.h = c("h.RE","h.yTRE","h.yTEE","h.HBDF","h.HBDFg0","h.IP3","h.IP","h.EP3","h.IP4","h.EP4","h.SS","h.okamura","h.const")

h.RE=function(y,x,b) { ## Radial exponential (half normal), noted h1 in Borchers et al. 2016
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(theta[1]*(y^2+x^2)^(-(theta[2]+2)/2))
}
attributes(h.RE)=list(fName="h.RE")

h.yTRE=function(y,x,b) { ## Radial exponential with longitudinal (y) translation for h(0)<1
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@references Borchers DL and Cox MJ, (2016). Distance sampling detection functions: 2D or not 2D? Biometrics
   if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b)
  return(theta[1]*((y+theta[3])^2+x^2)^(-(theta[2]+2)/2))
  # #theta=c(exp(b[1:2]),b[3])
  # theta1=theta[1]^(1/theta[2])
  # return(((x/theta1)^2+((y+theta[3])/theta1)^2)^(-theta[2]/2))
}
attributes(h.yTRE)=list(fName="h.yTRE")

h.yTEE=function(y,x,b) { ## Radial exponential with ellipse shape and h(0)<0
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@references Borchers DL and Cox MJ, (2016). Distance sampling detection functions: 2D or not 2D? Biometrics
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 4.")
  }
  theta=exp(b)
  thetax=theta[1]^(1/theta[2])
  thetay=theta[4]^(1/theta[2])
  return(((x/thetax)^2+((y+theta[3])/thetay)^2)^(-theta[2]/2))
}
attributes(h.yTEE)=list(fName="h.yTEE")


h.HBDF=function(y,x,b) {  ## Hayes and Buckland's (1983) function modelling that animal are more easily sighted if the observer is moving towards them, noted h2 in Borchers et al. 2016
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(theta[1]*y*(y^2+x^2)^(-(theta[2]+3)/2))
}
attributes(h.HBDF)=list(fName="h.HBDF")

h.HBDFg0=function(y,x,b) { ## h.HBDF with h(0)<1
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[1:2])
  dF=function(y,theta) theta[1]*y*(y^2+x^2)^(-(theta[2]+3)/2)
  g0=plogis(b[3])
  return(dF(y,theta)*g0)
}
attributes(h.HBDFg0)=list(fName="h.HBDFg0")


h.IP3=function(y,x,b) { ## Three-parameter inverse power hazard detection function  with h(0)<1
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-modulated Poisson process models for animal availability" Biometrics (in press).
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x/theta[2])^2+(y/theta[2])^2))^(theta[3]+1)
  return(p)
}
attributes(h.IP3)=list(fName="h.IP3")


h.IP=function(y,x,b) {  ## Inverse power hazard function
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x)^2+(y)^2))^(theta[2]+1)
  return(p)
}
attributes(h.IP)=list(fName="h.IP")


h.EP3=function(y,x,b) {  ## Three-parameter exponential power hazard detection function with h(0)<1
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-modulated Poisson process models for animal availability" Biometrics (in press).
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[2:3])
  dF=function(y,x,theta) exp(-(x^theta[2]+y^theta[2])/(theta[1]^theta[2]))
  g0=plogis(b[1])
  return(dF(y,x,theta)*g0)
}
attributes(h.EP3)=list(fName="h.EP3")


h.IP4=function(y,x,b) {  ## Four-parameter inverse power hazard detection function with h(0)<1
#'Has form h(y,x)=theta[1]*(1/(1+(x/theta[2])^2+(y/theta[4])^2))^(theta[3]+1).
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-modulated Poisson process models for animal availability" Biometrics (in press).
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x/theta[2])^2+(y/theta[4])^2))^(theta[3]+1)
  return(p)
}
attributes(h.IP4)=list(fName="h.IP4")


h.EP4=function(y,x,b) { ## Four-parameter exponential power hazard detection function  with h(0)<1
#'Has form h(y,x)=theta[1]*exp(-(x^theta[3]+y^theta[3])/(theta[2]^theta[3])).
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-modulated Poisson process models for animal availability" Biometrics (in press).
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[2:4])
  dF=function(y,x,theta) exp(-((x/theta[1])^theta[2]+(y/theta[3])^theta[2]))
  g0=plogis(b[1])
  return(dF(y,x,theta)*g0)
}
attributes(h.EP4)=list(fName="h.EP4")


h.SS=function(y,x,b=c(0,0)) {
#'@references Skaug, Hans J., and Tore Schweder. "Hazard models for line transect surveys with independent observers." Biometrics 55.1 (1999): 29-36.
  fName='h.exp2'
  mu=exp(b[1])
  sigma=exp(b[2])
  gama=2
  hr=mu*exp(-(x^gama + y^gama)/sigma^gama)
  return(hr)
}
attributes(h.SS)=list(fName="h.SS")


h.okamura=function(y,x,b=c(0,0)) {
#'@references Okamura, Hiroshi et al. (2003) Abundance Estimation of Diving Animals by the Double-Platoform Line Transect Method, Biometrics 59(3):512-520.
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  sigma.x=exp(b[1])
  sigma.y=exp(b[2])
  return(exp(-(x/sigma.x + y/sigma.y)))
}
attributes(h.okamura)=list(fName="h.okamura")


h.const=function(y,x,b=1) { ## Constant rate
  return(rep(b[1],length(y)))
}
attributes(h.const)=list(fName="h.const")

