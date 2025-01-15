model <- cds.hr$ddf
lower <- 0
dat <- model$data

width <- model$meta.data$width
left <- model$meta.data$left
ddfobj <- model$ds$aux$ddfobj
point <- model$ds$aux$point

if(is.null(model$ds$aux$int.range)){
    int.range <- c(0,width)
}else{
    int.range <- model$ds$aux$int.range
}

if(is.matrix(int.range)){
    max.range <- as.vector(int.range[1,])
    int.range <- int.range[2:dim(int.range)[1],]
    range.varies <- TRUE
}else{
    max.range <- int.range
    normalize <- FALSE
    range.varies <- FALSE
}

selected <- rep(TRUE, nrow(ddfobj$xmat))

if(is.matrix(int.range)){
    int.range <- int.range[selected, ]
}
  
xmat <- ddfobj$xmat[selected,]
if(!is.null(ddfobj$scale)){
    z <- ddfobj$scale$dm[, , drop=FALSE]
}else{
    z <- matrix(1, nrow=1, ncol=1)
}
  
if(length(model$fitted)==1){
    pdot <- rep(model$fitted, sum(as.numeric(selected)))
}else{
    pdot <- model$fitted[selected]
    Nhat <- sum(1/pdot)
}
  
zdim <- dim(z)[2]
n <- length(xmat$distance)
  
if(!is.null(breaks)){
    nc <- length(breaks)-1
}
  
if(is.null(nc)){
    nc <- round(sqrt(n), 0)
}
  
xgrid <- seq(0, width, length.out=101)
zgrid <- matrix(rep(z[1,], length(xgrid)), byrow=TRUE, ncol=sum(zdim))

  
breaks <- seq(0,width,width/nc)
# create intervals of distance (breaks) for the chosen number of classes (nc).
  if(is.null(breaks)){
    if(is.null(model$meta.data$binned)){
      binned <- FALSE
    }else{
      binned <- model$meta.data$binned
    }
    if(binned){
      breaks <- model$ds$aux$breaks
      nc <- length(breaks)-1
    }else{
      breaks <- c(max(0, (max.range[1])),
                  max.range[1]+((max.range[2]-max.range[1])/nc)*(1:nc))
      if(breaks[1]>left){
        breaks <- c(left, breaks)
        nc <- nc+1
      }
    }
  }
  
  
  # get the histogram object
  # hist.obj <- hist(dat[ ,vname][keep], breaks=breaks, plot=FALSE)
  hist.obj <- hist(dat$distance, breaks=seq(0,width,width/nc), plot=FALSE)
  # what's the top of the largest bar?
  ymax <- max(hist.obj$counts)
  
  # breaks=seq(0,79,79/5)
  # Rescaling for the histogram
  if(normalize & !point){
    bindata <- function(x, r, breaks){
      return(hist(r[r>=x[1] & r<=x[2]], breaks=breaks, plot=FALSE)$counts)
    }
    sumit<-function(x,n,wt){
      return(sum(x/(wt*n)))
    }
    expected.counts <- apply(int.range, 1, bindata,
                             r=(0:1000)*width/1001, breaks=breaks)
    expected.counts <- apply(expected.counts, 1, sumit, n=1001, wt=pdot)
  }else{
    if(!point){
      expected.counts <- (breaks[2:(nc+1)]-breaks[1:nc])*(Nhat/breaks[nc+1])
    }else{
      if(!pdf){
        expected.counts <- -apply(matrix(c(breaks[2:(nc+1)]^2, breaks[1:nc]^2),
                                         ncol=2, nrow=nc),
                                  1, diff)*(Nhat/breaks[nc+1]^2)
      }else{
        expected.counts <- sum(hist.obj$counts)
      }
    }
  }
  
  # rescale the histogram object by the expected counts
  # but only if we don't have point/pdf plots
  pdf=FALSE
  if(!(pdf & point)){
    hist.obj$density <- hist.obj$counts/expected.counts
    hist.obj$density[expected.counts==0] <- 0
  }
  hist.obj$equidist <- FALSE
  

  # histline(hist.obj$counts, breaks=breaks, lineonly=FALSE, ylim=c(0, ymax))#,
             # xlab=xlab, ylab="Frequency") #angle=angval1,
             # density=denval1, col=pl.col)
    # area under the histogram
    hist_area <- sum(hist.obj$density*diff(breaks))

    
    ## detection function
    
    point_vals <- detfct(xmat$distance, ddfobj, select=selected, width=width,
                           left=left)
    
    # set y labels, limits and tick marks (det.plot) depending on if we
    # are plotting PDF or df
    if(is.null(ylim)) ylim<-c(0, max(hist.obj$density, max(point_vals)))
    if(pdf){
      if(is.null(ylab)) ylab <- "Probability density"
      det.plot <- FALSE
    }else{
      if(is.null(ylab)) ylab <- "Detection probability"
      det.plot <- TRUE
    }
    
    # # plot the histogram
    # histline(hist.obj$density, breaks=breaks, lineonly=FALSE,
    #          # xlab=xlab, ylab=ylab,
    #          ylim=c(0, max(hist.obj$density, max(point_vals))))#,
    #          # angle=angval1, density=denval1, col=pl.col)
    hist_data <- data.frame(density=hist.obj$density, counts=hist.obj$counts, mids=hist.obj$mids)
    

      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- ddfobj$scale$dm[rep(1, length(xgrid)), ,drop=FALSE]
      }
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- ddfobj$shape$dm[rep(1, length(xgrid)), ,drop=FALSE]
      }
      
linevalues <- detfct(xgrid, ddfobj, width=width, left=left)

detection_data <- data.frame(observed_distances=xgrid,
                               fitted_values=linevalues)
cds_plot_data <- list(hist_data,detection_data)


