fitMisRepARMA <- function(y, tol, B, p_AR, q_MA, covars=NULL, misReport="U", method="freq", ...)
{
  Call <- match.call()
  params_post <- NA
  if (method=="freq")
  {
    class(params_post) <- "try-error"
    j <- 1
    orig.model <- estimate(as.numeric(y), tol, p_AR, q_MA, covars=covars, misReport)
    if (!is.null(covars))
    {
      while(class(params_post)=="try-error" & j < 10)
      {
        invisible(capture.output(params_post <- try(tsboot(y, statistic=estimate, R=B, sim="model", orig.t=TRUE, n.sim=length(y),
                                                           ran.gen=ran.genf, ran.args=list(ar=orig.model[grepl("ar", names(orig.model)) & !grepl("var", names(orig.model))],
                                                                                           ma=orig.model[grepl("ma", names(orig.model))],
                                                                                           intercept=orig.model[grepl("intercept", names(orig.model))],
                                                                                           var=attr(orig.model, "var"),
                                                                                           q=attr(orig.model, "q"),
                                                                                           w=attr(orig.model, "w"),
                                                                                           z=attr(orig.model, "z"),
                                                                                           misReport=misReport), 
                                                           ...,
                                                           tol=tol, p_AR=p_AR, q_MA=q_MA, covars=covars, misReport=misReport), silent=TRUE)))
        j <- j + 1
      }
    }else{
      while(class(params_post)=="try-error" & j < 10)
      {
        invisible(capture.output(params_post <- try(tsbootstrap(y, nb=B, statistic=estimate, 
                                                                type="stationary", tol=tol, p_AR=p_AR, q_MA=q_MA, covars=covars, misReport=misReport), silent=TRUE)))
        j <- j + 1
      }
      params_post$t  <- params_post$statistic
      params_post$t0 <- params_post$orig.statistic
    }
    colnames(params_post$t) <- names(orig.model)
    names(params_post$t0) <- names(orig.model)
    params_post <- list(data=y, t0=params_post$t0, t=params_post$t)
    class(params_post) <- c("fitMisRepARMA")
    attr(params_post, "covars") <- attr(orig.model, "covars")
    attr(params_post, "q") <- attr(orig.model, "q")
    attr(params_post, "w") <- attr(orig.model, "w")
    attr(params_post, "z") <- attr(orig.model, "z")
    attr(params_post, "Call") <- Call
    return(params_post)
  }else{
    if (method=="bayes")
    {
      stop("Bayesian method not yet implemented")
    }else{
      stop("Method not recognized. Choose either 'freq' or 'bayes'")
    }
  }
}
