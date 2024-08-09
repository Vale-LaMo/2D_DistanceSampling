scalevalue<-function(key.scale, z){
  # scalevalue - computes scale of detection function for a given 
  #              set of scale covariates and parameters, using log link
  # Args:
  #  key.scale      scale parameters
  #  z              design matrix for scale covariates
  exp(z %*% key.scale)
}


keyfct.hz <- function(distance, key.scale, key.shape){
  return(1 - exp( - (distance/key.scale)^( - key.shape)))
}


integratepdf <- function(ddfobj, select, width, int.range, standardize=TRUE,
                         point=FALSE, left=0, doeachint=FALSE){
  # Make sure there is consistency between integration ranges and data
  # It is ok to have a single observation with multiple ranges or a single range
  # with multiple observations but otherwise the numbers must agree if both >1
  
  if(!is.matrix(int.range)){
    if(is.vector(int.range) && length(int.range)==2){
      int.range <- matrix(int.range, ncol=2, nrow=1)
    }else{
      stop("int.range is not a matrix and cannot be given the required matrix structure")
    }
  }
  
  ## select tells us which data we compute the integrals for,
  ## if it's null then we compute for all data in ddfobj$xmat
  if(is.null(select)){
    nobs <- nrow(ddfobj$xmat)
    index <- 1:nobs
  }else{
    nobs <- sum(select)
    index <- which(select)
  }
  
  # number of integration ranges must be 1 or match number of observations
  # OR only be 1 observation and many ranges
  if(nrow(int.range)>1 && nobs>1 && nrow(int.range)!=nobs){
    stop("\n Number of integration ranges (int.range) does not match number of observations\n")
  }
  
  ## Now compute the integrals
  
  # if there is only 1 integral to compute (no covariates/1 set of covariates
  # & only one set of integration ranges), that's easy
  if(nobs==1){
    return(gstdint(int.range[1, ], ddfobj=ddfobj, index=1, select=NULL,
                   width=width, standardize=standardize, point=point,
                   stdint=FALSE, left=left))
  }else{
    # if there are multiple covariates or multiple ranges
    
    # make int.range have nintegrals rows if we want it to
    # already checked above that this is okay
    # this allows us to simplify the code below
    if(nrow(int.range)==1){
      int.range <- t(replicate(nobs, int.range, simplify=TRUE))
    }
    
    ### find unique observations
    # need unique (model matrix)-(int.range) combinations
    # want them within the rows we selected to compute for.
    # Know from above that int.range has either nrow(data) rows or
    #  length(index) rows.
    if(ddfobj$type=="tpn"){
      ind <- 1:nrow(ddfobj$shape$dm[index, , drop=FALSE])
    }else{
      if(is.null(ddfobj$shape)){
        newdat <- cbind(ddfobj$scale$dm[index, , drop=FALSE], int.range)
      }else{
        if(ncol(ddfobj$shape$dm)>1){
          scale_dm <- ddfobj$shape$dm[index, , drop=FALSE]
          scale_dm[, "(Intercept)"] <- NULL
        }else{
          scale_dm <- NULL
        }
        newdat <- cbind(ddfobj$scale$dm[index, , drop=FALSE],
                        scale_dm, int.range)
      }
      u.rows <- mgcv::uniquecombs(newdat)
      uu.index <- sort(unique(attr(u.rows, "index")))
      u.index <- attr(u.rows, "index")
      
      # generate the indices that we want to calculate integrals for
      ind <- match(uu.index, u.index)
    }
    if(length(width)==1){
      width <- rep(width, nrow(int.range))
    }
    if(length(left)==1){
      left <- rep(left, nrow(int.range))
    }
    
    # calculate the integrals
    ints <- gstdint(int.range[ind, , drop=FALSE], ddfobj=ddfobj,
                    index=index[ind], select=NULL, width=width[ind],
                    standardize=standardize, point=point,
                    stdint=FALSE, left=left[ind], doeachint=doeachint)
    
    if(ddfobj$type=="tpn"){
      integrals <- ints
    }else{
      ## now rebuild the integrals and populate the return vector
      integrals <- ints[attr(u.rows, "index")]
    }
  }
  return(integrals)
}


gstdint <- function(x, ddfobj, index=NULL, select=NULL, width,
                    standardize=TRUE, point=FALSE, stdint=TRUE,
                    doeachint=FALSE, left=left){
  
  if(!is.matrix(x)){
    x <- matrix(x, ncol=2)
  }
  ## NB this calculates the integral of:
  ##       g(x)/w             for line transects
  ##       2*r*g(r)/width^2   for point transects
  
  # Set of observations for computation of detection function can
  # be specified with logical (select) and numeric (index) values.
  # Either or both can be specified although the latter is unlikely.
  if(is.null(select)){
    # use all
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm
    }else{
      # use only those with specific indices
      scale.dm <- ddfobj$scale$dm[index, , drop=FALSE]
    }
  }else{
    # Use those with select=TRUE
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm[select, , drop=FALSE]
    }else{
      # use the numeric index within those with select=TRUE
      scale.dm <- ddfobj$scale$dm[select, , drop=FALSE][index, , drop=FALSE]
    }
  }
  
  # when we have half-normal, key only use the exact analytic expression
  # for the integral using the error function/analytic expression
  if(ddfobj$type=="hn" & is.null(ddfobj$adjustment) & !doeachint){
    
    key.scale <- scalevalue(ddfobj$scale$parameters, scale.dm)
    
    if(point){
      # analytic expression for integral of 2*r*g(r)/width^2 when
      #  g(r) is half-normal
      int <- (2*(key.scale^2*exp(-x[ ,1]^2/(2*key.scale^2))-
                   key.scale^2*exp(-x[ ,2]^2/(2*key.scale^2))))/width^2
    }else{
      # define the error function in terms of pnorm
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # analytic expression for integral of g(x)/w when g(x) is half-normal
      
      # let's speed the above up a bit
      # when the left limit is zero (i.e. no left truncation) then just
      # need one call to erf
      if(all(x[,1]==0)){
        int <- (1/width)*(sqrt(pi/2)*key.scale*(erf(x[, 2]/
                                                      (key.scale*sqrt(2)))))
        
      }else{
        int <- (1/(width-left))*sqrt(pi/2)*key.scale*
          (-erf(x[, 1]/(key.scale*sqrt(2)))+
             erf(x[, 2]/(key.scale*sqrt(2))))
      }
    }
    return(int)
  }else if(ddfobj$type=="hr" & is.null(ddfobj$adjustment) &
           !doeachint & !point & all(x[, 1]==0)){
    
    # do the spline shortcut
    # note that this will only work for a given shape parameter
    # don't do this if you have left truncation, point transects or
    # adjustments
    
    # Documentation by Jeff Laake from a previous incarnation of this code:
    # uses the cgftab which is a spline fitted to a table of standardized
    # integrals and the value is interpolated from the spline for each
    # observation.
    # This used to speed up integration of the detection function with scale
    # covariates. The detection function is integrated at a series of points
    # from 0 to W and then a spline is fitted to the computed values which are
    # cumulative (integral from 0 to x < integral from 0 to x+dx fr dx>0). The
    # spline is then used to predict values of the integral which depend on the 
    # scale which can depend on observation specific covariates.
    
    # set up the integration grid (these are values of x/sigma)
    xx <- exp(0.05*(1:100))-1
    xx <- cbind(c(0, xx[1:99]), xx)
    # fix scale and shape matrices
    ddfobj$scale$dm <- matrix(1, nrow=nrow(xx))
    if(!is.null(ddfobj$shape)){
      ddfobj$shape$dm <- matrix(1, nrow=nrow(xx))
    }
    
    # Create cumulative sums of values of g(x) integrals from grid (xx)
    y <- cumsum(gstdint(xx, ddfobj=ddfobj, index=1:nrow(xx), width=width,
                        standardize=standardize, point=point, doeachint=TRUE,
                        left=left))
    # Return smoothed spline of cumulative integral values
    spp <- smooth.spline(c(0, xx[ ,2]), c(0, y))
    
    # get the scale parameter
    xscale <- scalevalue(ddfobj$scale$parameters, scale.dm)
    
    # do the integration via predict.smooth.spline
    integrals <- predict(spp, as.vector(x[, 2]/xscale))$y -
      predict(spp, as.vector(x[, 1]/xscale))$y
    
    # rescale for point or line
    if(!point){
      integrals <- xscale*integrals
    }else{
      integrals <- xscale^2*integrals
    }
    
    return(integrals)
    
  }else if(ddfobj$type=="tpn" & is.null(ddfobj$adjustment) & !doeachint){
    
    apex <- exp(ddfobj$shape$parameters)
    
    # left key
    ddfobj$xmat$.dummy_apex_side <- 0
    left.scale.dm <- setcov(ddfobj$xmat, ddfobj$scale$formula)
    left.scale <- scalevalue(ddfobj$scale$parameters, left.scale.dm)
    # right key
    ddfobj$xmat$.dummy_apex_side <- 1
    right.scale.dm <- setcov(ddfobj$xmat, ddfobj$scale$formula)
    right.scale <- scalevalue(ddfobj$scale$parameters, right.scale.dm)
    
    int1 <- int2 <- left.scale*0
    
    for(i in 1:length(int1)){
      int1[i] <- left.scale[i]*(pnorm(apex, mean=apex, sd=left.scale[i]) -
                                  pnorm(left[i], mean=apex, sd=left.scale[i]))
      int2[i] <- right.scale[i]*(pnorm(width[i], mean=apex, sd=right.scale[i]) -
                                   pnorm(apex, mean=apex, sd=right.scale[i]))
    }
    int <- int1+int2
    
    return(sqrt(pi*2)/(width-left)*int)
  }else{
    # loop over the integration ranges, calculating integrals
    res <- rep(NA, nrow(x))
    
    # duplicate width/left if necessary
    if(length(width) != nrow(x)){
      width <- rep(width, nrow(x))
    }
    if(length(left) != nrow(x)){
      left <- rep(left, nrow(x))
    }
    
    # wrapper around detection function to handle the case where g(x) < 0
    dpdf <- function(x, width, ddfobj, select, index, standardize, stdint,
                     point, left){
      v <- distpdf(x, width=width, ddfobj=ddfobj, select=select, index=index,
                   standardize=standardize, stdint=stdint, point=point,
                   left=left)
      v[v<1e-6] <- 0
      v
    }
    
    # now integrate for each observation
    xmatsave <- ddfobj$xmat
    for(i in 1:nrow(x)){
      ddfobj$xmat <- xmatsave[i, , drop=FALSE]
      res[i] <- integrate(dpdf, lower=x[i, 1], upper=x[i, 2], width=width[i],
                          ddfobj=ddfobj, select=select[i], index=index[i],
                          rel.tol=1e-7, standardize=standardize,
                          stdint=stdint, point=point, left=left[i])$value
    }
    return(res)
  }
}