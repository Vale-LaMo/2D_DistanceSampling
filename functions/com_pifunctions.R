#### These are functions modeling the change in animal density with perpendicular distance to a line-transect
#### x = perpendicular distance
#### w = perpendicular truncation distance
#### All excepted pi.sigmo are extracted from the LT2D package

list.pi = c("pi.HN","pi.CHN","pi.TN","pi.const","pi.sigmo", "pi.sigmoI")

pi.HN=function(x,logphi,w) { ## Half-normal
  hnF=function(x,logphi) exp(-x^2/(2*exp(logphi[1])^2))
  return(hnF(x,logphi)/integrate(hnF,0,w,logphi)$value)
}
attributes(pi.HN)=list(fName="pi.HN")


pi.CHN=function(x,logphi,w){ ## Complementary Half-normal
  chnF=function(x,logphi) 1-exp(-(x-logphi[1])^2/(2*exp(logphi[2])^2))
  return(chnF(x,logphi)/integrate(chnF,0,w,logphi)$value)
}
attributes(pi.CHN)=list(fName="pi.CHN")


pi.TN=function(x,logphi,w){ ## Truncated normal
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  mu=logphi[1]
  sigma=exp(logphi[2])
  f=dnorm(x,mean=mu,sd=sigma)
  denom=(pnorm(w,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))
  if(denom>0) f=f/denom else f=0
  return(f)
}
attributes(pi.TN)=list(fName="pi.TN")


pi.const=function(x,logphi=NULL,w){ ## Uniform
  return(rep(1/w,length(x)))
}
attributes(pi.const)=list(fName="pi.const")


pi.sigmo = function(x,logphi,w) {  ## sigmoid density function that sums to 1 over [0,w]
  a = plogis(logphi[1])*3*w
  b = plogis(logphi[2])*w 
  return(1/(a*log(exp((w-b)/a)+1)-a*log(exp(-b/a)+1))*1/(exp(-(x-b)/a)+1))
}      
attributes(pi.sigmo)=list(fName="pi.sigmo")

pi.sigmoI = function(x,logphi,w) {  ## sigmoid density function that sums to 1 over [0,w] and is >0 for distance zero
  a = plogis(logphi[1])*3*w
  b = plogis(logphi[2])*w 
  c = plogis(logphi[3])
  return(1/(c*w + a*log(exp((w-b)/a)+1)-a*log(exp(-b/a)+1))*(c + 1/(exp(-(x-b)/a)+1)))
}      
attributes(pi.sigmoI)=list(fName="pi.sigmoI")

