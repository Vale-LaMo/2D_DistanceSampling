
####---- Load packages ----
library(tidyverse)
library(readxl)
library(mvtnorm)
# if(!"devtools" %in% rownames(installed.packages())) 
#   {install.packages("devtools")}
# devtools::install_github('david-borchers/LT2D')
library(LT2D)
library(Distance)
library(hrbrthemes)
library(ggpubr)

####---- Load 2D distance functions ----
source("functions/com_hfunctions.R")
source("functions/com_pifunctions.R")
source("functions/com_likelihoodutilities.R")
source("functions/GoFy_mod.R") # custom GoFy function, modified by VLM 2022-11-11
source("functions/plotfit.x.red.R") # custom function, modified by VLM 2023-08-31, to have a red line instead of a grey one

####---- Set pars ----
params <- vector("list",8)
params$trunc_perp_dist_perc <- 5
params$h.function <- "h.RE"
params$pi.function <- "pi.sigmo"
params$starting_values <-  c(0.25,0.25,-4,-1) # used for duiker
params$sd <-  6 # used for duiker
params$n_models <-  400
params$n_hpars <- 2
params$n_pipars <- 2


####---- Load and clean the dataset ----
# usati questi due file:
# data <- read.csv("data/duikersimula_9_7_2024.csv") # ok
# data <- read.csv("data/simula_20240730/impalasimula_bw_0_9.csv") # ok usata questa

# altri file che potrebbero tornare utili:
# data <- read.csv("data/impalasimula_9_7_2024.csv")
# data <- read.csv("data/simula_20240730/impalasimula_bw_1_2.csv")
# data <- read.csv("data/simula_20240730/impalasimula_bw_1_4.csv")
# data <- read.csv("data/simula_20240730/impalasimula_bw_1.csv")

colnames(data) <- c("perp_dist", "forw_dist", "samplesize", "replicate", "id")
data %>% 
  mutate(detected = 1,
         id_survey = paste(samplesize,replicate,sep="-")) -> data
head(data)


#### Dealing with NA and non-detections (not necessary for sims but kept for consistency with the original code)
data_clean_tot <- 
  data %>% 
  filter(detected != 0, # we only include actual observations in the dataset used to fit the detection function
         perp_dist != "NA", # we remove lines with NA distances
         forw_dist != "NA")
data_clean_tot$forw_dist <- abs(data_clean_tot$forw_dist) # we make sure all distances are positive (see Discussion for details)

####---- Simulations: For loop ----

stats_df_groups <- vector("list",length(unique(data_clean_tot$id_survey)))
fitVU <- vector("list",length(unique(data_clean_tot$id_survey)))
best_mod <- vector("list",length(unique(data_clean_tot$id_survey)))
cvs <- vector("list",length(unique(data_clean_tot$id_survey)))
# for (i in 1:3) {
for (i in 169:length(unique(data_clean_tot$id_survey))) {
  print(i)
  skip_to_next <- FALSE
  data_clean <- filter(data_clean_tot, id_survey == unique(id_survey)[i])
  
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
  
  ystart = max(data_trunc$forw_dist) # change this to the desired truncation distance if necessary, e.g.
  # ystart = params$trunc_forw_dist_m
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
  for (j in 1:length(tab.AIC$m)) {
    fName = params$h.function
    CV.phat.values[j] <- phatModels(list(FIT[[tab.AIC$m[j]]]))$CV.phat
    # LT2D::phatModels(modList = list(FIT[tab.AIC$m[i]]))$CV.phat
  }
  (tab.AIC %>%
      mutate(CV.phat = CV.phat.values) %>%
      filter(CV.phat == min(CV.phat)) -> best_mod[[i]])
  fitVU[[i]] = FIT[[best_mod[[i]]$m]]
  cvs[[i]] <- best_mod[[i]] %>% 
    mutate(id_survey = data_clean$id_survey[1],
           sample_size = data_clean$samplesize[1])
  
  # tryCatch(print(b), error = function(e) { skip_to_next <<- TRUE})
  
  # provo in extremis a togliere questo
  tryCatch(plotfit.x.red(x[x<=w],fitVU[[i]],nclass=20,nint=100),
           rug(x[x<=w]),
           error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }
  
  # if(skip_to_next) { fName = params$h.function
  # # GoF for perpendicular distances
  # GoFx(fitVU[[i]],plot=TRUE)$pvals
  # 
  # # # GoF for forward distances
  # # fName = params$pi.function
  # # GoFy_mod(fitVU[[i]],plot=TRUE)$pvals
  # # # plotfit.smoothfy(fitVU,nclass=32);rug(x=y[x<=w])
  # # # plotfit.y(y[x<=w & y<=ystart],x,fitVU,nclass=20);rug(x=y[x<=w])
  # # plotfit.smoothfy(fitVU[[i]],xmax=w)
  # 
  # (LT2D::phatModels(modList = list(fitVU[[i]]), # same as fitVU
  #                   n=length(na.omit(data_trunc$cluster_size))) %>%
  #     mutate(id_survey = data_clean$id_survey[1],
  #            sample_size = data_clean$samplesize[1]) -> stats_df_groups[[i]]) }
  
  # rug(x[x<=w])
  # plotfit.x.red(x[x<=w],fitVU[[i]],nclass=20,nint=100);rug(x[x<=w])

  fName = params$h.function
  # GoF for perpendicular distances
  GoFx(fitVU[[i]],plot=TRUE)$pvals

  # # GoF for forward distances
  # fName = params$pi.function
  # GoFy_mod(fitVU[[i]],plot=TRUE)$pvals
  # # plotfit.smoothfy(fitVU,nclass=32);rug(x=y[x<=w])
  # # plotfit.y(y[x<=w & y<=ystart],x,fitVU,nclass=20);rug(x=y[x<=w])
  # plotfit.smoothfy(fitVU[[i]],xmax=w)

  (LT2D::phatModels(modList = list(fitVU[[i]]), # same as fitVU
                    n=length(na.omit(data_trunc$cluster_size))) %>%
      mutate(id_survey = data_clean$id_survey[1],
             sample_size = data_clean$samplesize[1]) -> stats_df_groups[[i]])
  
}


# salvataggio dei risultati

# save(fitVU, file="output/simula/duikersimula_9_7_2024_fitVU.RData", compress = F)
# save(best_mod, file="output/simula/duikersimula_9_7_2024_bestmod.RData", compress = F)
# save(cvs, file="output/simula/duikersimula_9_7_2024_cvs.RData", compress = F)
# save(cvs_simula, file="output/simula/duikersimula_9_7_2024_cvs_simula.RData", compress = F)
# save(stats_df_groups, file="output/simula/duikersimula_9_7_2024_statsdfgroups.RData", compress = F)

# save(fitVU, file="output/simula/impalasimula_31_7_2024_fitVU.RData", compress = F)
# save(best_mod, file="output/simula/impalasimula_31_7_2024_bestmod.RData", compress = F)
# save(cvs, file="output/simula/impalasimula_31_7_2024_cvs.RData", compress = F)
# save(cvs_simula, file="output/simula/impalasimula_31_7_2024_cvs_simula.RData", compress = F)
# save(stats_df_groups, file="output/simula/impalasimula_31_7_2024_statsdfgroups.RData", compress = F)

####---- Results:  ----
load("output/simula/duikersimula_9_7_2024_cvs.RData") # dati di base dei cv
# do.call("rbind",stats_df_groups) %>% 
#   filter(CV.phat <= 1) -> cvs_simula # Remove non-sense CVs (> 1)
#   # mutate(id_survey = unique(data_clean_tot$id_survey),
#   #        sample_size = as.numeric(str_extract(id_survey, "[0-9]+")))-> cvs_simula
# do.call("rbind",stats_df_groups) -> cvs_simula
do.call("rbind",cvs) -> cvs_simula


library(Rmisc)
data_summary <- function(data, varname, groupnames){
  # require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      lcl = CI(x[[col]])[[3]],
      ucl = CI(x[[col]])[[1]])
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

library(plyr)
cvs_simula %>% 
  filter(CV.phat < 0.9) -> cvs_simula_clean # for duiker e impala seconda versione
# cvs_simula %>% 
#   filter(CV.phat < 0.7) -> cvs_simula_clean # for impala prima versione
cvs_simula_grouped <- data_summary(cvs_simula_clean, varname="CV.phat", 
                                   groupnames=c("sample_size"))
# Convert dose to a factor variable
cvs_simula_grouped$sample_size = as.factor(cvs_simula_grouped$sample_size)


# # sysfonts::font_add_google("Roboto Condensed")
# # import_roboto_condensed()
# # update_geom_font_defaults(font_rc_light)
# import_plex_sans()
# extrafont::font_import()


species.name <- "Impala"
species.name <- "Duiker"


## boxplot
ggplot(cvs_simula, aes(x = as.factor(sample_size), y = CV.phat)) +
  geom_boxplot() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  scale_x_discrete(expand=c(-1.1,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.05, 1)) +
  xlab("Sample size") +
  ylab("CV") +
  # theme(axis.text.y=element_blank()) +
  # labs(
  #   # title="IBM Plex Sans Test",
  #   subtitle=species.name,
  #   # caption="Source: hrbrthemes & IBM"
  # ) +
  theme_ipsum_ps(axis_title_size = 10, axis = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## plot with sd error bars
ggplot(cvs_simula_grouped, aes(x = sample_size, y = CV.phat)) +
  geom_point() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_errorbar(aes(ymin=CV.phat-sd, ymax=CV.phat+sd), width=.3,
                position=position_dodge(0)) +
  scale_x_discrete(expand=c(-1.05,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.05, 0.5)) +
  xlab("Sample size") +
  ylab("CV") +
  # theme(axis.text.y=element_blank()) +
  # labs(
  #   # title="IBM Plex Sans Test",
  #   subtitle=paste(species.name,": Plot with SD error bars",sep=""),
  #   # caption="Source: hrbrthemes & IBM"
  # ) +
  theme_ipsum_ps(axis_title_size = 10, axis = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
## plot with CI error bars
ggplot(cvs_simula_grouped, aes(x = sample_size, y = CV.phat)) +
  geom_point() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.3,
                position=position_dodge(0)) +
  scale_x_discrete(expand=c(-1.05,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.05, 0.5)) +
  xlab("Sample size") +
  ylab("CV") +
  # theme(axis.text.y=element_blank()) +
  # labs(
  #   # title="IBM Plex Sans Test",
  #   subtitle=paste(species.name,": Plot with 95% CI error bars",sep="")
  #   # caption="Source: hrbrthemes & IBM"
  # ) +
  theme_ipsum_ps(axis_title_size = 10, axis = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

table(cvs_simula_clean$sample_size)

# flush_ticks(p)

  smoothingSpline = smooth.spline(cvs_simula$sample_size, cvs_simula$CV.phat, spar=0.35)
  plot(cvs_simula$sample_size, cvs_simula$CV.phat)
  lines(smoothingSpline)

# cvs_simula %>% 
#   filter(sample_size != 25) %>%
#   filter(sample_size != 75) %>% 
#   filter(sample_size != 175) %>% 
#   filter(sample_size <= 200) -> cvs_simula2
# qplot(cvs_simula2$sample_size, cvs_simula2$CV.phat, geom='smooth', span =0.5) +
#   xlab("Sample size") +
#   ylab("CV") +
#   # scale_x_discrete(expand=c(0,0), breaks = c(50,100,150,200)) +
#   # scale_y_continuous(expand=c(0,0), limits=c(0, 0.2)) +
#   theme_ipsum_ps(axis_title_size = 15, grid = "XY")



# save(fitVU, file="output/simula/impala_fitVU.RData", compress = F)
# save(best_mod, file="output/simula/impala_bestmod.RData", compress = F)
# save(cvs, file="output/simula/impala_cvs.RData", compress = F)
# save(stats_df_groups, file="output/simula/impala_statsdfgroups.RData", compress = F)

# Error in integrate(h, y[i], ymax, x = x[i], b = b, subdivisions = 1000L) :
#   the integral is probably divergent
 