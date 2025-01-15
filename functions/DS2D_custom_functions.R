# Generate parameters with constraints on the last two parameters
generate_params_with_constraint <- function(length.pars, starting_values, sd) {
  pars <- rnorm(length.pars, mean = starting_values, sd = sd)
  
  # Enforce the constraint: last two parameters must be non-positive
  while (pars[length(pars) - 1] > 0 || pars[length(pars)] > 0) {
    pars[length(pars) - 1] <- rnorm(1, mean = starting_values[length(pars) - 1], sd = sd)
    pars[length(pars)] <- rnorm(1, mean = starting_values[length(pars)], sd = sd)
  }
  
  return(pars)
}

#### Model fitting function
fit_model <- function(hr, pi.x,
                      starting_values = c(0.25,0.25,-4,-1), sd = 6) {
  y = data_trunc$forw_dist
  x = data_trunc$perp_dist
  
  # # these functions work: h.RE, h.IP, h.SS, h.okamura
  # # functions tested and working with h.RE: pi.sigmo, pi.CHN, pi.TN
  
  ystart = ceiling(max(y))
  w = ceiling(max(x))
  length.b <- if (identical(hr, h.RE)) 2 else 3  
  length.logphi <- if (identical(pi.x, pi.sigmo)) 2 else 3
  length.pars = length.b + length.logphi
  debug=FALSE
  
  FIT=list(); AICvalues=NULL
  for (m in 1:params$n_models) {
    # print(m)
    set.seed(m)
    pars <- generate_params_with_constraint(length.pars,
                                            starting_values, sd)
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
    
    if(!is.na(fit[1])) {
      if(any(is.nan(fit$corr))) {
        AICvalues=c(AICvalues, 1e12)
      } else {
        AICvalues=c(AICvalues, fit$AIC)
      }
    } else {
      AICvalues=c(AICvalues, 1e12)
    }
  }
  save(FIT, file = paste("output/FIT_", attr(hr, "fName"),"_", attr(pi.x, "fName"),"_", params$species_name, ".RData", sep = ""), compress = FALSE)
  return(list(FIT = FIT, AICvalues = AICvalues))
}



analyze_selected_models <- function(hr, pi.x, fit_model) {
  
  # Select models with modAIC < threshold
  good_models <- data.frame(m = 1:length(fit_model$FIT), modAIC = fit_model$AICvalues) %>%
    filter(modAIC < 1e12)
  
  # Extract the list of models to analyze further
  m_list <- good_models$m
  selected_models <- fit_model$FIT[m_list]
  
  # Initialize lists and vectors to store results
  gx <- list()
  gxp <- numeric(length(selected_models))
  CV.phat.values <- numeric(length(selected_models))
  
  # Analyze GoFx results
  for (i in seq_along(selected_models)) {
    # print(i)
    gx[[i]] <- tryCatch.W.E(GoFx(selected_models[[i]]))
    if (!"error" %in% class(gx[[i]]$value)) {
      gxp[i] <- GoFx(selected_models[[i]])$pvals[["Kolmogarov-Smirnov"]]
    } else {
      gxp[i] <- NA  # Handle errors by assigning NA
    }
    CV.phat.values[i] <- phatModels(list(selected_models[[i]]))$CV.phat
  }
  
  # Initialize lists and vectors for GoFy_mod analysis
  gy <- list()
  gyp <- numeric(length(selected_models))
  
  # Analyze GoFy_mod results
  for (i in seq_along(selected_models)) {
    # print(i)
    gy[[i]] <- tryCatch.W.E(GoFy_mod(selected_models[[i]]))
    if (!"error" %in% class(gy[[i]]$value)) {
      gyp[i] <- GoFy_mod(selected_models[[i]])$pvals[["Kolmogarov-Smirnov"]]
    } else {
      gyp[i] <- NA  # Handle errors by assigning NA
    }
  }
  
  # Initialize lists and vectors for BIC and AICc
  modBIC <- list()
  modAICc <- list()
  
  # Calculate BIC
  for (i in seq_along(selected_models)) {
    # print(i)
    modBIC[[i]] <- 2*selected_models[[i]]$value+(length(selected_models[[i]]$par))*log(dim(data_trunc)[1])
    modAICc[[i]] <- 2*selected_models[[i]]$value + length(selected_models[[i]]$par) +
      # correction term
      2*(2*(length(selected_models[[i]]$par))^2 + 
           2*(length(selected_models[[i]]$par)))/((dim(data_trunc)[1])-length(selected_models[[i]]$par)-1)
  }
  
  # Compile results into a data frame
  results <- data.frame(m = m_list,
                        modAIC = fit_model$AICvalues[m_list],
                        modAICc = unlist(modAICc),
                        modBIC = unlist(modBIC)) %>%
    mutate(GoFx_pvalue = gxp,
           GoFy_pvalue = gyp,
           CV.phat = CV.phat.values) %>%
    filter(round(GoFx_pvalue,1) >= 0.1 & round(GoFy_pvalue,1) >= 0.1) %>%
    arrange(modAIC) %>%
    mutate(h_func = attr(hr, "fName"),
           pi_func = attr(pi.x, "fName"))
  
  return(results)
}