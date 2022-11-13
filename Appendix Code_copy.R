#### Analysis of distance sampling data: live observations and dung piles
#### By simulating a plausible set of forward distances

library(mvtnorm)
library(tidyverse)
# if(!"devtools" %in% rownames(installed.packages())) {install.packages("devtools")}
# devtools::install_github('david-borchers/LT2D')
library(LT2D)

#### Load 2D distance functions
source("functions/com_hfunctions.R")
source("functions/com_pifunctions.R")
source("functions/com_likelihoodutilities.R")

#### Load dataset
data <- read.csv("data/monticola_FINAL_20122013.csv", sep=",")
data %>% 
  rename(x = distance) -> data ## rename perpendicular distance to x
data$transect_start_time = as.POSIXct(as.character(data$transect_start_time),format = "%d/%m/%Y %H:%M:%S")
data$transect_end_time = as.POSIXct(as.character(data$transect_end_time),format = "%d/%m/%Y %H:%M:%S")
data$obs_time = as.POSIXct(as.character(data$obs_time),format = "%d/%m/%Y %H:%M:%S")

#### Transect lengths 
transect.names = c("A","I","J","K","L")
transect.lengths = c(8, 6, 5.8, 6, 5.4)   # A I J K L
names(transect.lengths) =  transect.names
data = data[data$transect_id %in% transect.names,]
dim(data)
table(data$transect_id, as.POSIXlt(data$obs_time)$year)

#### Generate forward distance data
####    (scenario letters correspond to Fig. 3 in the main text)
scenario = "c"
if (scenario == "a")
  data$y = pmax(0.01, data$x + rnorm(length(data$x), -1+(data$obs_type=="Crotte"), 1))
# nel primo scenario la distanza yi è all'incirca dello stesso ordine di grandezza di xi, con qualche variabilità gaussiana
if (scenario == "b") {
  mx=max(data$x)
  data$y = pmax(0.01, (mx/2)*exp(-data$x/(mx/2)) + rnorm(length(data$x), -1+(data$obs_type=="Crotte"), 1))
# correlazione negativa tra yi e xi, più sei vicino più si allontanano, più grande yi più piccolo xi
  # la riga sopra è un modo per creare una correlazione negativa basata sull'esponenziale
}
if (scenario == "c") {
  # nel terzo scenario yi e xi sono indipendenti
  b = c(0.2,0.2)
  x<-runif(length(data$x), px(13-0.01,rev(b),h.RE,13,nint=100),1)
  f <- function(x, u) px(x,rev(b),h.RE,13,nint=100) - u
  my.uniroot <- function(x) uniroot(f, interval=c(0, 13), tol = 0.0001, u = x)$root
  data$y <- vapply(x, my.uniroot, numeric(1))
}


#### Divide into year- and cue-specific datasets
# crea 4 dataset dividendo in base al tipo di dato e all'anno
montiVU2012 = data[as.POSIXlt(data$obs_time)$year==112 & data$obs_type=="Vu" & data$obs_period=="Nocturne",]
montiCR2012 = data[as.POSIXlt(data$obs_time)$year==112 & data$obs_type=="Crotte" & data$obs_period=="Diurne",]
montiVU2013 = data[as.POSIXlt(data$obs_time)$year==113 & data$obs_type=="Vu" & data$obs_period=="Nocturne",]
montiCR2013 = data[as.POSIXlt(data$obs_time)$year==113 & data$obs_type=="Crotte" & data$obs_period=="Diurne",]
unlist(lapply(list(montiVU2012,montiCR2012,montiVU2013,montiCR2013), function(x) dim(x)[1]))


#### Fit 2D distance sampling model using multiple initial value to avoid local minima in the deviance
 y = c(montiVU2012$y , montiVU2013$y)
 x = c(montiVU2012$x , montiVU2013$x)
 hr = h.RE # h.yTRE non è compatibile con pi.sigmoI
 # funzionano h.RE, h.IP, h.SS, h.okamura
 pi.x = pi.sigmo
 # funzionano con h.RE: pi.sigmo, pi.CHN, pi.TN
 ystart = ceiling(max(y))
 w = ceiling(max(x))
 length.b = 2
 debug=FALSE

 FIT=list(); dev=NULL
 for (m in 1:10) {
   pars = rnorm(4, c(0.25,0.25,-4,-1), 3)
   # pars = rnorm(4, c(0.25,0.25,4,1), 3)
   tmp0 <- tryCatch.W.E (
     fityx(y,x,pars[1:length.b],
           hr,ystart,pi.x,pars[(length.b+1):length(pars)],w,control=list(),hessian=TRUE,corrFlag=0.7,debug=FALSE)
   )
   fit = NA
   if(! "error" %in% class(tmp0$value)) {
     fit <- tmp0$value
     fit$vcov <-  matrix(Matrix::nearPD(fit$vcov)$mat,4,4)
   }
   FIT[[m]] = fit
   if(is.na(fit[1])) dev=c(dev, 1e12) else dev = c(dev, fit$val)
 }
 fitVU = FIT[[which.min(dev)]]
 tabVU = matrix(NA,2,3)
 if(is.na(fitVU[1])) tabVU = matrix(NA,2,3) else {
 tmp1 <- tryCatch.W.E (boot(fitVU))
 if(! "error" %in% class(tmp1$value))  tabVU=tmp1$value
 }
 ### Population density
 tabVU[2,]/sum(transect.lengths)/(fitVU$w/1000) /2 ## divided by two to get the average over 2 years
 
 
 # vedere https://github.com/david-borchers/LT2D/blob/master/inst/FitsForPaper.r
 # linea nera sigmoide = distribuzione reale animali
 # linea grigia = detection osservata
 # linea tratteggiata = detection corretta tenendo conto della risposta comportamentale
 plotfit.x(x[x<=w],fitVU,nclass=20);rug(x[x<=w])
 fName = "h1"
 GoFx(fitVU,plot=TRUE)$pvals
 plotfit.y(y[x<=w],x,fitVU,nclass=20);rug(x=y[x<=w])
 plotfit.smoothfy(fitVU,nclass=32);rug(x=y[x<=w])
 GoFy(fitVU,plot=TRUE)$pvals # non funziona
 #EHSW:
 phatInterval(fitVU)
 phatInterval(fitVU)*w
 # p(0):
 p0.n=1-Sy(0,0,ystart,fitVU$b,h1);p0.n
 plotfit.smoothfy(fitVU,xmax=0.004)
 
 
 
 y = c(montiCR2012$y , montiCR2013$y) 
 x = c(montiCR2012$x , montiCR2013$x)
 hr = h.RE
 pi.x = pi.sigmo
 ystart = ceiling(max(y))
 w = ceiling(max(x))
 length.b = 2
 debug=FALSE

 FIT=list(); dev=NULL
 for (m in 1:10) {
   pars = rnorm(4, 0, 3)
   tmp0 <- tryCatch.W.E (
     fityx(y,x,pars[1:length.b],hr,ystart,pi.x,pars[(length.b+1):length(pars)],w,control=list(),hessian=TRUE,corrFlag=0.7,debug=FALSE)
   )
   fit = NA
   if(! "error" %in% class(tmp0$value)) {
     fit <- tmp0$value
     fit$vcov <-  matrix(Matrix::nearPD(fit$vcov)$mat,4,4)
   }
   FIT[[m]] = fit
   if(is.na(fit[1])) dev=c(dev, 1e12) else dev = c(dev, fit$val)
 }
 fitCR = FIT[[which.min(dev)]]
 tabCR = matrix(NA,2,3)
 if(is.na(fitCR[1])) tabCR = matrix(NA,2,3) else {
 tmp1 <- tryCatch.W.E (boot(fitCR))
 if(! "error" %in% class(tmp1$value))  tabCR=tmp1$value
 }
 
### Conversion Ratio dung/live
tabCR[2,]*(fitVU$w/fitCR$w)/tabVU[2,]
### Population density
tabVU[2,]/sum(transect.lengths)/(fitVU$w/1000) /2 ## divided by two to get the average over 2 years
### Fit details (parameters values, fitted detection and avoidance functions) are in the objects fitVU and fitCR for live and dung data respectively
tabVU[2,]/sum(transect.lengths)/(fitVU$w/1000) /2 -> res1


res1
names(fitVU)
class(fitVU$hr)



