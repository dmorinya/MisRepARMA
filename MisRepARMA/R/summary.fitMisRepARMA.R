summary.fitMisRepARMA <- function(object, ...)
{
  call <- attr(object, "Call")
  coef.p <- object$t0
  coef.p <- coef.p[1:(length(coef.p)-1)]
  sds <- apply(object$t, 2, sd, na.rm=T)
  sds <- sds[1:(length(sds)-1)]
  dn <- c("Estimate", "Std. Error")
  coef.table <- cbind(coef.p, sds)
  dimnames(coef.table) <- list(names(coef.p), dn)
  aic <- object$t0[length(object$t0)]
  ans <- list(coefficients = coef.table, aic=aic, call=call)
  class(ans) <- "summary.fitMisRepARMA"
  return(ans)
}