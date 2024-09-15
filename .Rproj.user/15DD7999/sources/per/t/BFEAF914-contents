plotfit.x.red <- 
function(x,est,nclass=10,nint=100,
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
    lines(gridx,f.xfit,col="red",lwd=4)
    # overlay fitted detection function p(x), scaled to have area=1
    lines(gridx,p.xfit.std,lty=2,col="black",lwd=2)
    # overlay fitted animal density function pi(x), scaled to have area=1
    lines(gridx,adbn,col="black",lwd=2)
    legend("topright",legend=c("f(x)","p(x)",expression(pi(x))),
           col=c("red","black","black"),lwd=c(2,2,2),lty=c(1,2,1))
  }
  invisible(list(gridx=gridx,p.xpifit=p.xpifit,mufit=mufit,
                 f.xfit=f.xfit,p.xfit=p.xfit,ptot=ptot,p.xfit.std=p.xfit.std,adbn=adbn
  ))
}