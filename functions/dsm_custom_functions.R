#' @author David Lawrence Miller
check.cols <- 
  
  function(ddf.obj, segment.data, observation.data, segment.area){
  
  ## check that the columns are there
  checks <- list(segment.data = c("Effort", "Sample.Label"),
                 observation.data = c("object", "Sample.Label", "size",
                                      "distance"))
  
  # don't need to check distance if we don't have a detection function
  if(is.null(ddf.obj)){
    checks$observation.data <- checks$observation.data[checks$observation.data!=
                                                         "distance"]
  }
  
  if(!is.null(segment.area)){
    checks$segment.data <- checks$segment.data[checks$segment.data != "Effort"]
  }
  
  for(i in seq_len(length(checks))){
    check.res <- checks[[i]] %in% names(get(names(checks)[[i]]))
    if(any(!check.res)){
      
      stop(paste0("Column(s) \"",
                  paste(checks[[i]][!check.res],collapse="\", \""),
                  "\" not found in ", names(checks)[[i]],
                  ".\n  Check ?\"dsm-data\"."))
    }
  }
  
  ## check that Sample.Label is unique
  if(length(segment.data$Sample.Label)!=length(
    unique(segment.data$Sample.Label))){
    warning("'Sample.Labels are non-unique in segment data!")
  }
  
  invisible()
  }




#' @importFrom stats aggregate
make.data_vlm <- 
  
  function(response, ddfobject, segdata, obsdata, group,
                          convert.units, availability, segment.area,
                          family){
  
  # probably want to do something smart here...
  seglength.name <- 'Effort'
  segnum.name <- 'Sample.Label'
  distance.name <- 'distance'
  cluster.name <- 'size'
  
  # avoid irritating "tibble" issues
  segdata <- data.frame(segdata)
  obsdata <- data.frame(obsdata)
  
  # matching doesn't work later if we don't ensure character labels
  obsdata$object <- as.character(obsdata$object)
  obsdata$Sample.Label <- as.character(obsdata$Sample.Label)
  segdata$Sample.Label <- as.character(segdata$Sample.Label)
  
  # Estimating group abundance/density
  if(group){
    obsdata[, cluster.name][obsdata[, cluster.name] > 0] <- 1
  }
  
  # for single ddfs, make a list of 1
  # if(any(c("ddf", "dsmodel")  %in% class(ddfobject))){
  ddfobject <- list(ddfobject)
  segdata$ddfobj <- 1
  obsdata$ddfobj <- 1
  # }else{
  #   if(!any("ddfobj" %in% names(segdata))){
  #     stop("If there are multiple detection functions there must be a column named \"ddfobj\" in the segment and observation data.frame, see ?\"dsm-data\"")
  #   }
  # }
  
  # iterate over the list of ddfs
  full_obsdata <- c()
  segdata$segment.area <- NA
  segdata$width <- NA
  
  # deal with availability
  if(response == "density.est"){
    if(!(length(availability) %in% c(1, nrow(obsdata)))){
      stop("Length of 'availability' must be 1 or 'obsdata' rows")
    }
    obsdata$availability <- availability
  }else if(response == "count"){
    if(!(length(availability) %in% c(1, nrow(segdata)))){
      stop("Length of 'availability' must be 1 or 'segdata' rows")
    }
    segdata$availability <- availability
  }else if(response == "abundance.est"){
    if(!(length(availability) %in% c(1, nrow(obsdata)))){
      stop("Length of 'availability' must be 1 or 'obsdata' rows")
    }
    obsdata$availability <- availability
  }
  
  for(i in seq_along(ddfobject)){
    
    this_ddf <- ddfobject[[i]]
    
    # make this a mrds ddf object if we had a Distance one
    if(all(class(this_ddf) == "dsmodel")){
      this_ddf <- this_ddf$ddf
      ddfobject[[i]] <- this_ddf
    }
    
    # grab the probabilities of detection
    # fitted.p <- fitted(this_ddf)
    fitted.p <- p.xfit.std
    
    # remove observations which were not in the detection function
    if(!("ddfobj" %in% names(obsdata))){
      stop("No ddfobj column in observation data")
    }
    # this_obsdata <- obsdata[obsdata[["ddfobj"]]==i, ]
    this_obsdata <- filter(obsdata, distance <= w)
    
    # Check that observations are between left and right truncation
    # warning only -- observations are excluded below
    # No truncation check for strip transects
    if(any(this_obsdata[, distance.name] > this_ddf$w)){
      # if(any(this_obsdata[, distance.name] > this_ddf$meta.data$width)){
      warning(paste("Some observations are outside of detection function", i,
                    "truncation!"))
    }
    
    # reorder the fitted ps, match the ordering in obsdata
    # put this in a column of this_obsdata
    fitted.p.df$object <- as.character(fitted.p.df$object)
    left_join(this_obsdata, fitted.p.df[,c("object","p")], by="object") -> this_obsdata
    # this_obsdata$p <- fitted.p[match(this_obsdata$object, names(fitted.p))]
    
    # calculate the "width" of the transect, make sure we get it right
    # if we are doing left truncation
    width <- this_ddf$w
    if(!is.null(this_ddf$meta.data$left)){
      width <- width - this_ddf$meta.data$left
    }
    segdata$width[segdata$ddfobj==i] <- width
    
    # what if there are no matches? Perhaps this is due to the object
    # numbers being wrong? (HINT: yes.)
    if(nrow(this_obsdata) == 0){
      stop(paste("No observations in detection function", i,
                 "matched those in observation table. Check the \"object\" column."))
    }
    
    # # make sure that the right columns are in the obsdata
    # if(this_ddf$method %in% c("io", "trial")){
    #   if(!("observer" %in% names(this_obsdata))){
    #     stop("obsdata must have a column named observer")
    #   }else{
    #     if((this_ddf$method == "trial") && !all(this_obsdata$observer==1)){
    #       stop("Only observer 1 data is needed for obsdata with trial mode")
    #     }
    #   }
    # }
    
    # # depending on the detection function (and its data format)
    # # we need a different subset of obsdata
    # # ds = all detections
    # # io = only unique obsns
    # # trial = only obs 1
    # if(this_ddf$method == "io"){
    #   if(any(duplicated(this_obsdata$object))){
    #     stop(paste0("Some object IDs are duplicated in observation data for detection function model ", i, " only unique IDs are required for observation data for io models"))
    #   }
    # }else if(this_ddf$method == "trial"){
    #   this_obsdata <- this_obsdata[this_obsdata[["observer"]] == 1, ]
    # }
    
    # bind this to the full data
    full_obsdata <- rbind(full_obsdata, this_obsdata)
    
    # set the segment area for this detection function in the segdata
    # if(this_ddf$meta.data$point){
    #   # here "Effort" is number of visits
    #   segdata$segment.area[segdata$ddfobj==i] <- pi *
    #     segdata$width[segdata$ddfobj==i]^2 *
    #     segdata[, seglength.name][segdata$ddfobj==i]
    # }else{
    # line transects
    segdata$segment.area[segdata$ddfobj==i] <- 2 *
      segdata[, seglength.name][segdata$ddfobj==i] *
      segdata$width[segdata$ddfobj==i]
  }
  
  
  # set the full observation data
  obsdata <- full_obsdata
  
  # set the segment area in the data
  if(!is.null(segment.area)){
    segdata$segment.area <- segment.area
  }
  
  
  ## Aggregate response values of the sightings over segments
  if(response == "density.est"){
    responsedata <- aggregate(obsdata[, cluster.name]/
                                (obsdata$p * obsdata$availability),
                              list(obsdata[, segnum.name]), sum)
    off.set <- "none"
  }else if(response == "count"){
    responsedata <- aggregate(obsdata[, cluster.name],
                              list(obsdata[, segnum.name]), sum)
    off.set <- "eff.area"
  }else if(response == "abundance.est"){
    responsedata <- aggregate(obsdata[, cluster.name]/
                                (obsdata$p * obsdata$availability),
                              list(obsdata[, segnum.name]), sum)
    off.set <- "area"
  }
  
  
  ## warn if any observations were not allocated
  responsecheck <- aggregate(obsdata[, cluster.name],
                             list(obsdata[, segnum.name]), sum)
  if(sum(obsdata[, cluster.name]) != sum(responsecheck[, 2])){
    message(paste0("Some observations were not allocated to segments!\n",
                   "Check that Sample.Labels match"))
  }
  
  # name the response data columns
  names(responsedata) <- c(segnum.name, response)
  
  # if the Sample.Labels don't match at all then we need to stop, nothing
  # can work as all the response values will be zero
  if(!any(segdata[, segnum.name] %in% responsedata[, segnum.name])){
    stop("No matches between segment and observation data.frame Sample.Labels!")
  }
  
  # Next merge the response variable with the segment records and any
  # response variable that is NA should be assigned 0 because these
  # occur due to 0 sightings
  dat <- merge(segdata, responsedata, by=segnum.name, all.x=TRUE)
  dat[, response][is.na(dat[, response])] <- 0
  
  # for the offsets with effective area, need to make sure that
  # the ps match the segments
  if(off.set == "eff.area"){
    dat$p <- NA
    for(i in seq_along(ddfobject)){
      this_ddf <- ddfobject[[i]]
      
      # get all the covariates in this model
      # df_vars <- all_df_vars(this_ddf)
      df_vars <- 0
      
      if("fake_ddf" %in% class(this_ddf)){
        # strip transect
        dat[dat$ddfobj == i, ]$p <- 1
      }else if(length(df_vars) == 0){
        # if there are no covariates, and all the fitted ps are the same
        # then just duplicate that value enough times for the segments
        dat[dat$ddfobj == i, ]$p <- rep(fitted(this_ddf)[1],
                                        nrow(dat[dat$ddfobj == i, ]))
      }else if(this_ddf$method %in% c("ds", "io", "trial")){
        # get all the covariates in this model
        # df_vars <- all_df_vars(this_ddf)
        df_vars <- NULL
        
        # check these vars are in the segment table
        if(!all(df_vars %in% colnames(dat))){
          stop(paste0("Detection function covariates are not in the segment data",
                      "\n  Missing: ", df_vars[!(df_vars %in% colnames(dat))]))
        }
        
        # make a data.frame to predict for
        nd <- dat[dat$ddfobj == i, ]#[, df_vars, drop=FALSE]
        nd$distance <- 0
        
        this_ddf$method <- "2DDS"
        if(this_ddf$method == "io"){
          nd <- rbind(nd, nd)
          nd$observer <- c(rep(1, nrow(nd)/2), rep(2, nrow(nd)/2))
          dat[dat$ddfobj == i, ]$p <- predict(this_ddf, newdata=nd)$fitted
        }else if(this_ddf$method == "trial"){
          nd$observer <- 1
          dat[dat$ddfobj == i, ]$p <- predict(this_ddf, newdata=nd)$fitted
        }else{
          dat[dat$ddfobj == i, ]$p <- predict(this_ddf, newdata=nd)$fitted
        }
      }else{
        stop("Only \"ds\", \"io\" and \"trial\" models are supported!")
      }
    } # end loop over detection functions
  }
  
  # check that none of the Effort values are zero
  if(any(dat[, seglength.name]==0)){
    stop(paste0("Effort values for segments: ",
                paste(which(dat[, seglength.name]==0), collapse=", "),
                " are 0."))
  }
  
  # correct segment area units
  dat$segment.area <- dat$segment.area*convert.units
  
  # calculate the offset
  #   area we just calculate the area
  #   effective area multiply by p (and availability)
  #   when density is response, offset should be 1 (and is ignored anyway)
  dat$off.set <- switch(off.set,
                        eff.area = dat$segment.area*dat$p*dat$availability,
                        area     = dat$segment.area,
                        none     = 1)
  
  # calculate the density (count/area)
  if(response == "abundance.est"){
    dat[, response] <- dat[, response]/(dat$segment.area)
  }
  
  
  # Set offset as log (or whatever link is) of area or effective area
  # dat$off.set <- family$linkfun(dat$off.set)
  dat$off.set <- log(dat$off.set)
  
  return(dat)
}


dsm_vale <- 
  
  function (formula, ddf.obj, segment.data, observation.data, engine = "gam", dat, 
                      convert.units = 1, family = quasipoisson(link = "log"), group = FALSE, 
                      control = list(keepData = TRUE), availability = 1, segment.area = NULL, 
                      weights = NULL, method = "REML", ...) 
{
  stopifnot(engine %in% c("gam", "bam", "glm", "gamm"))
  if (is.null(ddf.obj)) {
    stop("NULL detection functions no longer supported, see ?dummy_ddf")
  }
  if (all(class(ddf.obj) != "list")) {
    if (all(class(ddf.obj) == "dsmodel")) {
      ddf.obj <- ddf.obj$ddf
    }
  }
  else {
    if (length(ddf.obj) == 1) {
      ddf.obj <- ddf.obj[[1]]
    }
    for (i in seq_along(ddf.obj)) {
      if (all(class(ddf.obj[[i]]) == "dsmodel")) {
        ddf.obj[[i]] <- ddf.obj[[i]]$ddf
      }
    }
  }
  response <- as.character(formula)[2]
  if (response %in% c("presence", "D", "density", "Dhat", "N", 
                      "Nhat", "n")) {
    stop(paste("Response", response, "is deprecated, see ?dsm for details."))
  }
  possible.responses <- c("density.est", "count", "abundance.est")
  if (!(response %in% possible.responses)) {
    stop(paste("Model must be one of:", paste(possible.responses, 
                                              collapse = ", ")))
  }
  # check.cols(ddf.obj, segment.data, observation.data, segment.area)
  # dat <- make.data(response, ddf.obj, segment.data, observation.data, 
  # group, convert.units, availability, segment.area, family)
  if (!(response %in% c("density.est"))) {
    formula <- as.formula(paste(c(as.character(formula)[c(2, 
                                                          1, 3)], "+ offset(off.set)"), collapse = ""))
  }
  else {
    if (is.null(weights)) {
      weights <- dat$segment.area
    }
    else if (length(weights) == 1) {
      weights <- rep(weights, nrow(dat))
    }
  }
  if (engine == "gamm" && dsm_env$old_mgcv) {
    message("You are using mgcv version < 1.7-24, please update to at least 1.7-24 to avoid fitting problems.")
  }
  if (engine %in% c("glm", "gamm")) {
    control$keepData <- NULL
  }
  args <- list(formula = formula, family = family, data = dat, 
               weights = weights, control = control, method = method, 
               ...)
  if (engine == "glm") {
    args$method <- NULL
  }
  fit <- withCallingHandlers(do.call(engine, args), warning = "matrixnotposdef.handler")
  if ("knots" %in% names(match.call())) {
    fit$knots <- get(as.character(match.call()$knots))
  }
  if (engine == "gamm") {
    fit$gam$ddf <- ddf.obj
    fit$gam$data <- dat
    fit$gam$gamm.call.list <- list(formula = formula, family = family, 
                                   data = dat, control = control)
  }
  else {
    fit$ddf <- ddf.obj
  }
  class(fit) <- c("dsm", class(fit))
  return(fit)
  }


summary.dsm.var_vale <- 
  
  function(object, detfunc, alpha=0.05, boxplot.coef=1.5,
                                 bootstrap.subregions=NULL, fName = "h.RE"){
  
  # storage
  sinfo <- list()
  # save the alpha value for cis
  sinfo$alpha <- alpha
  
  ### analytical variance estimation (varprop and gam results)
  sinfo$varprop <- object$var.prop
  sinfo$saved <- object
  sinfo$bootstrap <- object$bootstrap
  
  # re run the variance calculation, putting everything together
  pd <- c()
  off <- c()
  for(i in seq_len(length(object$pred.data))){
    pd <- rbind(pd, object$pred.data[[i]])
    off <- rbind(off, object$off.set[[i]])
  }
  object$pred.data <- pd
  object$off.set <- as.vector(off)
  
  if(object$var.prop){
    var.prop <- dsm_var_prop(object$dsm.obj, object$pred.data,
                             object$off.set, object$seglen.varname,
                             object$type.pred)
  } else {
    var.prop <- dsm_var_gam(object$dsm.obj, object$pred.data,object$off.set,
                            object$seglen.varname, object$type.pred)
  }
  
  sinfo$se <- sqrt(var.prop$pred.var)
  
  # grab the predicted values
  if(length(object$pred)>1) {
    sinfo$pred.est <- sum(unlist(object$pred), na.rm=TRUE)
  } else {
    sinfo$pred.est <- object$pred[[1]]
  }
  
  # if we're just using the GAM variance, then we need to combine using
  # the delta method
  
  # setup everything to be multi-ddf compatible
  ddf <- object$dsm.object$ddf
  
  sinfo$detfct.cv <- c() #(phatInterval(detfunc)[2])[[1]]
  cvp.sq <- 0
  for(i in seq_along(ddf)) {
    
    this_ddf <- ddf[[i]]
    # if(all(class(this_ddf)!="fake_ddf")){
    #   ddf.summary <- summary(this_ddf)
    #
    # this_cvp.sq <- (ddf.summary$average.p.se/
    #                   ddf.summary$average.p)^2
    this_cvp.sq <- ((phatInterval(detfunc)[2])[[1]])^2
    cvp.sq <- cvp.sq + this_cvp.sq
    # }else{
    # this_cvp.sq <- NA
    # }
    # sinfo$detfct.cv <- c(sinfo$detfct.cv, sqrt(this_cvp.sq))
    sinfo$detfct.cv <- c(sqrt(this_cvp.sq))
    # }
    
    sinfo$gam.cv <- sinfo$se/sinfo$pred.est
    
    sinfo$cv <- sqrt(cvp.sq+sinfo$gam.cv^2)
    
    # total se
    sinfo$se <- sinfo$cv*sinfo$pred.est
    # }
    if(sinfo$varprop) {
      sinfo$model.check <- object$model.check
    }
    
  }
  
  class(sinfo) <- "summary.dsm.var"
  return(sinfo)
}