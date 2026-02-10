print.summary.fitMisRepARMA <- function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nCoefficients:\n")
    print(x$coefficients)
    
    cat("\nAIC: ", x$aic, "\n")
  }
