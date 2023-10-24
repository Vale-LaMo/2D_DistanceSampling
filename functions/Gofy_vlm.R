GoFy_vlm <- 
function (fit, plot = FALSE, dotitle = FALSE) 
{
  ystart = fit$ystart
  w = fit$w
  hr = match.fun(fit$hr)
  analytic.F0 = TRUE
  if (fName == "h1") analytic.F0 = TRUE else analytic.F0 = FALSE
  ystart = fit$ystart
  pi.x = match.fun(fit$pi.x)
  b = fit$b
  logphi = fit$logphi
  x = fit$dat$x
  y = fit$dat$y
  n = length(x)
  if (analytic.F0) {
    zeros = (y == 0)
    Fy = rep(NA, length(x))
    if (length(zeros) > 0) {
      Fy[zeros] = HBhr(x[zeros], h1.to.HB(b))
      Fy[!zeros] = (1 - Sy(x = x[!zeros], y = y[!zeros], 
                           ymax = ystart, b = b, hfun = hr))
    }
    F0 = HBhr(x, h1.to.HB(b))
  }
  else {
    Fy = (1 - Sy(x = x, y = y, ymax = ystart, b = b, hfun = hr))
    F0 = (1 - Sy(x = x, y = rep(0, length(y)), ymax = ystart, 
                 b = b, hfun = hr))
  }
  Fy0 = Fy/F0
  Fy0.order = order(Fy0)
  yy = y[Fy0.order]
  xx = x[Fy0.order]
  cdf = Fy0[Fy0.order]
  e.cdf = order(cdf)/n
  dF = cdf - e.cdf
  worst = which(abs(dF) == max(abs(dF)))
  Dn = max(abs(dF)) * sqrt(n)
  p.cvm = goftest::cvm.test(Fy0)$p.value
  p.kolomogarov = 1 - p.kolomogarov(Dn)
  if (plot) {
    main = ""
    if (dotitle) 
      main = "Forward Dist. Q-Q Plot"
    plot(e.cdf, cdf, xlab = "Empirical Distribution Function", 
         ylab = "Cumulative Distribution Function", main = main, 
         xlim = c(1, 0), ylim = c(0, 1), pch = "+")
    lines(c(0, 1), c(0, 1))
    if (dotitle) 
      mtext(paste("p-values: Cramer-von Mises=", round(p.cvm, 
                                                       2), " ;\n                kolomogarov=", round(p.kolomogarov, 
                                                                                                     2)))
    points(e.cdf[worst], cdf[worst], col = "red")
  }
  pvals = c(p.cvm, p.kolomogarov)
  names(pvals) = c("Cramer-von Mises", "Kolmogarov-Smirnov")
  # return(data.frame(pvals = pvals, D.kolomogarov = Dn))
  return(list(pvals = pvals, D.kolomogarov = Dn))
}