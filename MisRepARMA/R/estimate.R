estimate <- function(data, tol, p_AR, q_MA, covars=NULL, misReport="U")
{
  if (!is.null(covars))
  {
    mod1  <- coef(lm(data~covars))
    data  <- data - covars %*% mod1[2:length(mod1)]
  }
  w0    <- 0.5
  q0    <- 0.5
  sd_x  <- sqrt(var(data)/(w0+(q0^2)*(1-w0)))
  invisible(capture.output(res1 <- try(normalmixEM(data, lambda=c(w0, (1-w0)), 
                                                   mu=c(mean(data), mean(data)/(w0+q0*(1-w0))),
                                                   sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE))) 
  if (class(res1)=="try-error") return(NA)
  if(misReport=="U")
  {
    q0   <- min(res1$mu)/max(res1$mu)
    w0   <- res1$lambda[which(res1$mu==min(res1$mu))]
    k1   <- which(res1$mu==min(res1$mu))
    k2   <- which(res1$mu==max(res1$mu))
  }
  if(misReport=="O")
  {
    q0   <- max(res1$mu)/min(res1$mu)
    w0   <- res1$lambda[which(res1$mu==max(res1$mu))]
    k1   <- which(res1$mu==max(res1$mu))
    k2   <- which(res1$mu==min(res1$mu))
  }
  q   <- q0
  w   <- w0
  eps <- tol + 1
  eps_ant <- tol + 2
  res <- res1
  i=1
  while (eps >= tol & eps < eps_ant)
  {
    eps_ant <- eps
    aux   <- ifelse(res$posterior[, k2]>0.5, 0, 1)
    x1    <- ifelse(aux==0, data, NA)   ### Misreported observations
    x2    <- ifelse(aux==1, data, NA)   ### Non-misreported observations
    q_ant <- q
    w_ant <- w
    invisible(capture.output(mod_infra   <- try(arima(x2, order=c(p_AR, 0, q_MA)), silent=TRUE)))
    invisible(capture.output(mod_noinfra <- try(arima(x1, order=c(p_AR, 0, q_MA)), silent=TRUE)))
    i=1
    while ((class(mod_infra)=="try-error" | class(mod_noinfra)=="try-error") & i < 10)
    {
      sd_x <- sqrt(var(data)/(w+(q^2)*(1-w)))
      invisible(capture.output(res1 <- try(normalmixEM(data, lambda=c(w, (1-w)), 
                                                       mu=c(mean(data), mean(data)/(w+q*(1-w))),
                                                       sigma=c(sd(data), sd_x), epsilon=1e-16, maxit=50000), silent=TRUE)))
      
      if(misReport=="U")
      {
        q <- min(res$mu)/max(res$mu)
        w <- res$lambda[which(res$mu==min(res$mu))]
        k1    <- which(res$mu==min(res$mu))
        k2    <- which(res$mu==max(res$mu))
      }
      if(misReport=="O")
      {
        q <- max(res$mu)/min(res$mu)
        w <- res$lambda[which(res$mu==max(res$mu))]
        k1    <- which(res$mu==max(res$mu))
        k2    <- which(res$mu==min(res$mu))
      }
      
      aux   <- ifelse(res$posterior[, k2]>0.5, 0, 1)
      x1    <- ifelse(aux==0, data, NA)   ### Misreported observations
      x2    <- ifelse(aux==1, data, NA)   ### Non-misreported observations
      invisible(capture.output(mod_infra   <- try(arima(x2, order=c(p_AR, 0, q_MA)), silent=TRUE)))
      invisible(capture.output(mod_noinfra <- try(arima(x1, order=c(p_AR, 0, q_MA)), silent=TRUE)))
      i<-i+1
    }
    if (class(mod_infra)!="try-error" & class(mod_noinfra)!="try-error")
    {
      q     <- as.numeric(coef(mod_infra)[names(coef(mod_infra))=="intercept"]/coef(mod_noinfra)[names(coef(mod_noinfra))=="intercept"])
      if (q < 0)  q <- sqrt(mod_infra$sigma2/mod_noinfra$sigma2)
      invisible(capture.output(res <- try(normalmixEM(data, lambda=c(w, (1-w)), 
                                                      mu=c(coef(mod_infra)[names(coef(mod_infra))=="intercept"],
                                                           coef(mod_noinfra)[names(coef(mod_noinfra))=="intercept"]),
                                                      sigma=c(sd(x2, na.rm=T), sd(x1, na.rm=T)), epsilon=1e-16, maxit=50000), silent=TRUE)))
      
      if(misReport=="U")
      {
        k1    <- which(res$mu==min(res$mu))
        k2    <- which(res$mu==max(res$mu))
      }
      if(misReport=="O")
      {
        k1    <- which(res$mu==max(res$mu))
        k2    <- which(res$mu==min(res$mu))
      }
      
      w     <- res$lambda[k1]
      eps   <- sqrt((q_ant-q)^2+(w_ant-w)^2)
    }else{
      return(NA)
    }
  }
  x_rec   <- ifelse(aux==0, data, data/q)
  mod_rec <- arima(x_rec, order=c(p_AR, 0, q_MA))
  sigma2  <- mod_rec$sigma2
  if (!is.null(covars))
  {
    estimates <- c(coef(mod_rec), mod1[2:length(mod1)], sigma2, q, w, AIC(mod_rec))
    names(estimates) <- c(names(coef(mod_rec)), names(mod1[2:length(mod1)]), "var", "q", "w", "AIC")
  }else{
    estimates <- c(coef(mod_rec), sigma2, q, w, AIC(mod_rec))
    names(estimates) <- c(names(coef(mod_rec)), "var", "q", "w", "AIC")
  }
  if (is.null(covars))
  {
    attr(estimates, "covars") <- NULL
  }else{
    attr(estimates, "covars") <- covars %*% mod1[2:length(mod1)]
  }
  attr(estimates, "z") <- aux
  attr(estimates, "q") <- q
  attr(estimates, "w") <- w
  attr(estimates, "var") <- sigma2
  attr(estimates, "AIC") <- AIC(mod_rec)
  return(estimates)
}