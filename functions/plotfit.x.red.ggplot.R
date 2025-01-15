plotfit.x.red.ggplot <- 
function(x,est,nclass=10,nint=100,
         plot=TRUE,title="",ymx=0.08,
         xaxislabel="Perpendicular distance (x)", showlegend=TRUE,
         species = "duiker", xminlogo, yminlogo,
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
    ## plot with base R
    # breaks=seq(0,w,length=(nclass+1))
    # hx=hist(x,breaks=breaks,plot=FALSE) # get hist bar heights
    # ymax=max(f.xfit,p.xfit.std,adbn,f.x,p.x.std,adbnTRUE,hx$density,ymx) 
    # ## Plot empirical distribution of sightings along x axis
    # hx=hist(x,breaks=breaks,freq=FALSE,ylim=c(0,ymax),
    #         main=title,xlab="Perpendicular distance (x)",ylab="pdf", ...)
    # ## Overlay predicted distribution of sightings along x axis, f(x)
    # lines(gridx,f.xfit,col="red",lwd=4)
    # # overlay fitted detection function p(x), scaled to have area=1
    # lines(gridx,p.xfit.std,lty=2,col="black",lwd=2)
    # # overlay fitted animal density function pi(x), scaled to have area=1
    # lines(gridx,adbn,col="black",lwd=2)
    # legend("topright",legend=c("f(x)","p(x)",expression(pi(x))),
    #        col=c("red","black","black"),lwd=c(2,2,2),lty=c(1,2,1))
  
    
    ## perplexity: plot with ggplot2
    library(ggplot2)
    library(ggh4x)
    library(grid)
    library(png)
    
    # Assuming you have the following variables defined:
    # x, gridx, f.xfit, p.xfit.std, adbn, w, nclass
    
    # Calculate breaks for the histogram
    breaks <- seq(0, w, length.out = nclass + 1)
    
    # Create a data frame for the histogram and the fitted lines
    data_hist <- data.frame(x = x)
    data_lines <- data.frame(gridx = gridx, f_xfit = f.xfit, p_xfit_std = p.xfit.std, adbn = adbn)
    
    # Calculate the maximum y value for scaling
    hx <- hist(x, breaks = breaks, plot = FALSE)
    ymax <- max(f.xfit, p.xfit.std, adbn, hx$density)
    
    # Add a column to identify each line for the legend
    data_lines_long <- data_lines %>%
      pivot_longer(cols = c(f_xfit, p_xfit_std, adbn),
                   names_to = "line_type",
                   values_to = "value")
    
    # Load the logo image
    logo_path <- paste("logos/", species, "_logo.png",sep = "")
    logo <- readPNG(logo_path)
    
    # Create the ggplot
    p <- ggplot(data_hist, aes(x = x)) +
      geom_histogram(aes(y = ..density..),  # Use density on the y-axis
                     breaks = breaks, 
                     color = "black", fill = "grey", boundary = 0, alpha = 0.5) +
      geom_line(data = data_lines_long, aes(x = gridx, y = value, color = line_type, linetype = line_type), size = 1) +
      # geom_line(data = data_lines, aes(x = gridx, y = f_xfit, colour = "f_xfit", linetype = "f_xfit")) +
      # geom_line(data = data_lines, aes(x = gridx, y = p_xfit_std, colour = "p_xfit_std", linetype = "p_xfit_std")) +
      # geom_line(data = data_lines, aes(x = gridx, y = adbn, colour = "adbn", linetype = "adbn")) +
      labs(title = "", 
           x = xaxislabel, 
           y = "Probability Density Function (pdf)") +
      scale_y_continuous(limits = c(0, ymx)) +  # Set limits for y-axis
      theme_classic(base_size = 10) +  # Using classic theme for a clean appearance
      scale_color_manual(name = "",
                         values = c("f_xfit" = "red", "p_xfit_std" = "black", "adbn" = "black"),
                         labels = c("f(x)", "p(x)", expression(pi(x)))) +
      scale_linetype_manual(name = "",
        values = c("f_xfit" = "solid", "p_xfit_std" = "dashed", "adbn" = "solid"),
                            labels = c("f(x)", "p(x)", expression(pi(x)))) +
      theme(legend.position=c(0.9, 0.8), legend.key.width=unit(1.5,"cm")) +
      # theme(
      #   # plot.title = element_text(hjust = 0.5),  # Center the title
      #   legend.position = "topright"  # Position of the legend
      # )
      guides(x = "axis_truncated", y = "axis_truncated") +
      annotation_custom(rasterGrob(logo, width = unit(2.5, "cm"), height = unit(2.5, "cm")), xmin=xminlogo, ymin=yminlogo)
    
    # Print the plot
    # if(showlegend==FALSE) print(p + theme(legend.position = "none")) else print(p)
    
  }
  invisible(list(gridx=gridx,p.xpifit=p.xpifit,mufit=mufit,
                 f.xfit=f.xfit,p.xfit=p.xfit,ptot=ptot,p.xfit.std=p.xfit.std,adbn=adbn
  ))
  return(if(showlegend==FALSE) print(p + theme(legend.position = "none")) else print(p))
}