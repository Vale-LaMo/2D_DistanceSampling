## Utilities and likelihood for 2D line transect distance analysis
## Mostly extracted from LT2D package with very limited typo and readability edits
## Designed to be used within the following protocol: 1) choose relevant combinations of functions h and pi from the available list provided in com_hfunctions and com_pifunctions 2) fit to data using fityx ; 3) select best model using AIC or other ; 4) plot along the perpendicular axis using plotfit.x ; 5) optional not implemented: compute SE using a parametric bootstrap within MNorm(fit$par,fit$vcov)

## "Waiting function"
#   Calculates the pdf of the 'waiting distance' \eqn{f(y,x)=h(y,x)*\exp(-\int_y^{ystart} h(t,x) dt)}.
fyx=function(y,x,b,hr,ystart,nint=100,eps=1e-12)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  n=length(x)
  f=intval=rep(NA,n)
  hr=match.fun(hr)
  for(i in 1:n) {
    y0=max(y[i],eps)
    dy=(ystart-y0)/nint/2                           # for crude integration
    yy=seq(y0,ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    h=hr(yy,rep(x[i],nint),b)
    int=sum(hr(yy,rep(x[i],nint),b)*dy*2)  # crude integration
    intval[i]=exp(-int)
    #    int=integrate(f=hr,lower=max(y[i],ylo),upper=ystart,x=x[i],b=b)
    #    intval[i]=exp(-int$value)
  }
  hrval=hr(y,x,b)
  bads=which(hrval>=.Machine$double.xmax) # identify infinite hazards
  if(length(bads)>0) { # infinite hazard so p(detect)=0
    f[bads]=.Machine$double.xmax
    f[-bads]=hr(y[-bads],x[-bads],b)*intval[-bads]
  }else{
    f=hr(y,x,b)*intval
  }
return(f)
}

## "Survival function"
#  Calculates probability to "survive undetected" until position y, for a given perpendicular distance x and given forward distance range ymax.
#  y is a vector
#  hfun is a function from "com_hfunctions.R" (got rid of original option to give a char string)
#  Numerical integration (if using h.HBDF there is an analytical solution worked out but not used here)
Sy=function(x,y,ymax,b,hfun) {
  n=length(x)
  if(length(y)!=n) stop("Lengths of x and y must be the same.")
  pS=rep(NA,n)
  h=match.fun(hfun)
  for(i in 1:n)
    pS[i]=exp(-integrate(h,y[i],ymax,x=x[i],b=b,subdivisions = 1000L)$value)
  return(pS)
}

## Perpendicular detection function
#  Calculates the perpendicular detection function, \eqn{p(x)}, for a given hazard.
px=function(x,b,hr,ystart,nint=100) 1-Sy(x,rep(0.0001,length(x)),ystart,b,hr)

## Product of p and pi
p.pi.x=function(x,b,hr,ystart,pi.x,logphi,w) px(x,b,hr,ystart)*pi.x(x,logphi,w)

## Negative log-likelihood  ####################################################
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param pars c(b,logphi); hazard rate and density log-parameters in a vector (see details).  
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal /!\ Must to ensure the hazard function has decayed to almost zero by ystart.
#'@param pi.x perpendicular distance density distribution
#'@param w perpendicular truncation distance.
#'@param length.b length of the hazard rate parameter vector

negloglik.yx=function(pars,y,x,hr,ystart,pi.x,w,length.b=2,debug=FALSE)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  if(debug) print(pars)
  n=length(y)
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  # unpack parameters /!\ pars must be passed as c(b,logphi), typically a length-4 vector
  b=pars[1:length.b]
  if(attributes(pi.x)$fName=="pi.const") logphi=NULL else logphi=pars[(1+length.b):length(pars)]
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(fyx(y,x,b,hr,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom
  ##llik[is.nan(llik)]=-9e37
  if(debug) print(llik)
  return(-llik)
}

fityx=function(y,x,b,hr,ystart,pi.x,logphi,w,control=list(),hessian=FALSE,corrFlag=0.7,debug=FALSE,...)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  if(debug) print(pars)
  n=length(y)
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  if(attributes(pi.x)$fName=="pi.const") pars=b else pars=c(b,logphi)
  length.b=length(b)
  
  fit=optim(par=pars,fn=negloglik.yx,y=y,x=x,hr=hr,ystart=ystart,pi.x=pi.x,w=w,
            length.b=length.b,
            hessian=hessian,debug=debug,control=control,...)
  fit$error=FALSE
  if(fit$convergence!=0){
    warning('Convergence issue (code = ', fit$convergence,') . Check optim() help.')
    fit$error=TRUE
  }
  fit$hr=hr
  fit$pi.x=pi.x
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par[1:length.b]  # ***
  if(length.b!=length(pars)){
    fit$logphi=fit$par[(1+length.b):length(pars)]
    parnames=c(rep("b",length.b), rep("logphi",length(pars)-length.b))
  }else{
      fit$logphi=NA
      parnames=c(rep("b",length.b))
  }    # ***
  fit$AIC=2*fit$value+2*length(fit$par)
  fit$dat=data.frame(x=x,y=y) # attach data to fitted object
  
  ## Standard errors
  if(hessian){
    mNames=paste('b',1:length.b,sep='')
    if(!all(is.na(fit$logphi))) mNames=c(mNames,paste('logphi',1:length(fit$logphi),sep=''))
    fit$vcov=solve(fit$hessian)
    if(any(diag(fit$vcov)<=0)) fit$vcov=MASS::ginv(fit$hessian + diag(length(fit$par))*1e-6)
    if(any(diag(fit$vcov)<=0)) {
      warning('Failed to invert hessian.  Model covergance problem in fityx?')
      fit$error=TRUE
      fit$CVpar=rep(NA,length(fit$par))
    } else fit$SE=sqrt(diag(fit$vcov))
    fit$corr=cov2cor(fit$vcov)
    row.names(fit$corr)=parnames
    colnames(fit$corr)=parnames
    corr=fit$corr
    corr[upper.tri(corr,diag=TRUE)]=NA
    corrIND=which(abs(corr)>corrFlag,arr.ind=T)
    if(nrow(corrIND)){
      warning('absolute correlation exceeds ',corrFlag,' in parameter estimates: ',
              paste(paste(mNames[corrIND[,1]],mNames[corrIND[,2]],sep=' to '),collapse='; '))
      fit$error=TRUE}
  }
  return(fit)
}

## Plot results along the x axis
plotfit.x=function(x,est,nclass=10,nint=100,
                   plot=TRUE,title="",ymx=NULL,
                   ...)
{
  Nhat.yx=bias=NULL
  b=est$b; hr=match.fun(est$hr); ystart=est$ystart; pi.x=match.fun(est$pi.x)
  logphi=est$logphi; w=est$w
  x=x[x>=0 & x<=w]
  f.x=p.x.std=adbnTRUE=0
  gridx=seq(1e-10,w,length=100)

  p.xpifit=p.pi.x(gridx,b,hr,ystart,pi.x,logphi,w)
  mufit=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr
                  ,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value 
  f.xfit=p.xpifit/mufit
  p.xfit=px(gridx,b,hr,ystart,nint=nint)
  ptot=integrate(f=px,lower=0,upper=w,b=b,hr=hr,ystart=ystart)$value
  p.xfit.std=p.xfit/ptot
  adbn=pi.x(gridx,logphi,w)
  
  if(plot){
    breaks=seq(0,w,length=(nclass+1))
    hx=hist(x,breaks=breaks,plot=FALSE) # get hist bar heights
    ymax=max(f.xfit,p.xfit.std,adbn,f.x,p.x.std,adbnTRUE,hx$density,ymx) 
    ## Plot empirical distribution of sightings along x axis
    hx=hist(x,breaks=breaks,freq=FALSE,ylim=c(0,ymax),
            main=title,xlab="perpendicular distance (x)",ylab="pdf", ...)
    ## Overlay predicted distribution of sightings along x axis, f(x)
    lines(gridx,f.xfit,col="grey",lwd=2)
    # overlay fitted detection function p(x), scaled to have area=1
    lines(gridx,p.xfit.std,lty=2,col="black",lwd=2)
    # overlay fitted animal density function pi(x), scaled to have area=1
    lines(gridx,adbn,col="black",lwd=2)
    legend("topright",legend=c("f(x)","p(x)",expression(pi(x))),
                col=c("grey","black","black"),lwd=c(2,2,2),lty=c(1,2,1))
  }
  invisible(list(gridx=gridx,p.xpifit=p.xpifit,mufit=mufit,
                 f.xfit=f.xfit,p.xfit=p.xfit,ptot=ptot,p.xfit.std=p.xfit.std,adbn=adbn
                 ))
}


## Parametric bootstrap
boot <- function(fit, Nboot=1000) {
pHat = N = NULL
length.b = length(fit$b)
hr=match.fun(fit$hr); ystart=fit$ystart; 
pi.x=match.fun(fit$pi.x); w=fit$w
for (k in 1:Nboot) {
  theta = rmvnorm(1, fit$par, fit$vcov)
  b = theta[1:length.b]
  logphi = theta[(length.b+1):length(theta)]
  phat = integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,
                ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value
  pHat = c(pHat, phat)
  N = c(N, length(x)/phat)
}
phat.est = integrate(f=p.pi.x,lower=0,upper=w,b=fit$b,hr=hr,
                ystart=ystart,pi.x=pi.x,logphi=fit$logphi,w=w)$value
tabVU = data.frame(rbind( c(phat.est, quantile(pHat, c(0.025,0.985))),
             c(length(x)/phat.est, quantile(N, c(0.025,0.985)))
             ))
colnames(tabVU) = c("Est","CI.low","CI.upp")
rownames(tabVU) = c("p-hat","N")
round(tabVU, 3)
return(tabVU)
}


## Effective strip width p-hat   (do parametric bootstrap to compute SE)
phat=function(fit=NULL,w=NULL,hr=NULL,b=NULL,ystart=NULL,pi.x=NULL,logphi=NULL)
{
  int=NULL
  if(!is.null(fit)){
    upper=fit$w; b=fit$b; hr=match.fun(fit$hr)
    ystart=fit$ystart;pi.x=match.fun(fit$pi.x);logphi=fit$logphi;w=fit$w
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,
                ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value
  }
  return(int)
}

tryCatch.W.E <- function(expr)
{
    W <- NULL
    w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
				     warning = w.handler),
	 warning = W)
}

