
#### Load packages
library(tidyverse)
library(readxl)
library(mvtnorm)
# if(!"devtools" %in% rownames(installed.packages())) 
#   {install.packages("devtools")}
# devtools::install_github('david-borchers/LT2D')
library(LT2D)
library(Distance)

#### Load 2D distance functions
source("functions/com_hfunctions.R")
source("functions/com_pifunctions.R")
source("functions/com_likelihoodutilities.R")
source("functions/GoFy_mod.R") # custom GoFy function, modified by VLM 2022-11-11
source("functions/plotfit.x.red.R") # custom function, modified by VLM 2023-08-31, to have a red line instead of a grey one

params$trunc_perp_dist_perc <- 5
params$h.function <- "h.RE"
params$pi.function <- "pi.sigmo"
params$starting_values <-  c(0.25,0.25,-4,-1)
params$sd <-  6
params$n_models <-  200
params$n_hpars <- 2
params$n_pipars <- 2



#### Load dataset
data <- read.csv("data/SAMPLE.csv")
colnames(data) <- c("perp_dist", "forw_dist", "samplesize", "id")
data %>% 
  mutate(detected = 1) -> data
head(data)

#### Dealing with NA and non-detections
data_clean <- 
  data %>% 
  filter(detected != 0, # we only include actual observations in the dataset used to fit the detection function
         perp_dist != "NA", # we remove lines with NA distances
         forw_dist != "NA")
data_clean$forw_dist <- abs(data_clean$forw_dist) # we make sure all distances are positive (see Discussion for details)

plot(data_clean$perp_dist, data_clean$forw_dist,
     xlim=c(0,max(data_clean$perp_dist)),
     ylim=c(0,max(data_clean$forw_dist)),
     xlab = "Perpedincular distance",
     ylab = "Forward distance")

par(mfrow = c(1,2))
hist(data_clean$perp_dist, main = "", xlab = "Perpendicular distance (m)")
boxplot(data_clean$perp_dist, ylab = "Perpendicular distance (m)")



no_data <- round(params$trunc_perp_dist_perc*length(data_clean$perp_dist)/100,0) # no. data to be deleted
threshold <- sort(data_clean$perp_dist, decreasing = TRUE)[no_data+1] # threshold
data_trunc <- 
  data_clean %>% 
  filter(perp_dist <= threshold)

par(mfrow = c(1,2))
hist(data_trunc$forw_dist, main = "", xlab = "Forward distance (m)")
boxplot(data_trunc$forw_dist, ylab = "Forward distance (m)")

# ystart = max(data_trunc$forw_dist) # change this to the desired truncation distance if necessary, e.g.
ystart = params$trunc_forw_dist_m
data_trunc <- 
  data_trunc %>% 
  filter(forw_dist <= ystart)

#### Model fitting
y = data_trunc$forw_dist
x = data_trunc$perp_dist
hr = params$h.function # h.yTRE not compatible with pi.sigmoI
# these functions work: h.RE, h.IP, h.SS, h.okamura
pi.x = params$pi.function # perpendicular distance function used
# functions tested and working with h.RE: pi.sigmo, pi.CHN, pi.TN
ystart = ceiling(max(y))
w = ceiling(max(x))
length.b = params$n_hpars # pars for h function
length.logphi = params$n_pipars # pars for pi function
length.pars = length.b + length.logphi
debug=FALSE

FIT=list(); AICvalues=NULL
for (m in 1:params$n_models) {
  set.seed(m)
  pars = rnorm(length.pars, # tot no. pars 
               params$starting_values, params$sd) 
  set.seed(m)
  tmp0 <- tryCatch.W.E (
    fityx(y,x,pars[1:length.b],
          hr,ystart,pi.x,
          pars[(length.b+1):length(pars)],w,
          control=list(),
          hessian=TRUE,corrFlag=0.7,debug=FALSE)
  )
  fit = NA
  if(! "error" %in% class(tmp0$value)) {
    fit <- tmp0$value
    fit$vcov <-  matrix(Matrix::nearPD(fit$vcov)$mat,length.pars,length.pars)
  }
  FIT[[m]] = fit
  # if(is.na(fit[1])) dev=c(dev, 1e12) else dev = c(dev, fit$AIC)
  
  # with the funciton used in this analyses, we add the constraint that the pi.x pars should be negative
  # to maintain the sigmoid shape
  if(params$pi.function == "pi.sigmo" & params$h.function == "h.RE") {
    if(!is.na(fit[1])) {
      if(any(is.nan(fit$corr)) | any(fit$par[3:4] > 0)) {
        AICvalues=c(AICvalues, 1e12)
      } else {
        AICvalues=c(AICvalues, fit$AIC)
      }
    } else {
      AICvalues=c(AICvalues, 1e12)
    }
  } else {
    if(!is.na(fit[1])) {
      if(any(is.nan(fit$corr))) {
        # if(all(fit$b > 0)) {
        AICvalues=c(AICvalues, 1e12)
      } else {
        AICvalues=c(AICvalues, fit$AIC)
      }
    } else {
      AICvalues=c(AICvalues, 1e12)
    }
  }
}

data.frame(m = 1:params$n_models, modAIC = AICvalues) -> df.AIC
df.AIC %>% 
  arrange(modAIC) %>% 
  filter(modAIC <= min(df.AIC$modAIC) + 2) -> tab.AIC
tab.AIC

CV.phat.values <- vector("numeric", length(tab.AIC$m))
for (i in 1:length(tab.AIC$m)) {
  fName = params$h.function
  CV.phat.values[i] <- phatModels(list(FIT[[tab.AIC$m[i]]]))$CV.phat
  # LT2D::phatModels(modList = list(FIT[tab.AIC$m[i]]))$CV.phat
}
(tab.AIC %>%
    mutate(CV.phat = CV.phat.values) %>%
    filter(CV.phat == min(CV.phat)) -> best_mod)
fitVU = FIT[[best_mod$m]]

plotfit.x.red(x[x<=w],fitVU,nclass=20,nint=100);rug(x[x<=w])

fName = params$h.function
# GoF for perpendicular distances
GoFx(fitVU,plot=TRUE)$pvals

# GoF for forward distances
fName = params$pi.function
GoFy_mod(fitVU,plot=TRUE)$pvals
# plotfit.smoothfy(fitVU,nclass=32);rug(x=y[x<=w])
# plotfit.y(y[x<=w & y<=ystart],x,fitVU,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fitVU,xmax=199)

(LT2D::phatModels(modList = list(fitVU), # same as fitVU
                  n=length(na.omit(data_trunc$cluster_size))) -> stats_df_groups)
